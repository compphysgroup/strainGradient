/*strainGradient.cc
 * Describe boundary value problem to solve.
 * --------------------------
 * Things may need modification:
 * 1.setting parameters' values   
 * 2.Reading NURBS geometry	  
 * 3.applying boundary condition  
 * 4.solver			  
 * 5.generating output data	
 *--------------------------
 *June 23, 2014  
 */

//#define SOLVER_MT
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <ctime>
#include <vector>
#include "hdf5Functions.h"
#include "igaClasses.h"
#include "vtkFunctions.h"
#include "sparseDataClass.h"
#include "mechanics.h"
#include "functionEvaluations.h"
#include "solutionClasses.h"
#include "base/table.h"
#include <Sacado.hpp>
#include <base/thread_management.h>
#include <base/multithread_info.h>
#include <tbb/task_scheduler_init.h>
#ifdef SOLVER_MT
#include "superLUMT_solver.h"
#else
#include "superLU_solver.h"
#endif
#include "parameters.h"
using namespace std;

#define DIMS 3 //dimention of problem
#define NUM_QUAD_POINTS 3 //NUM_QUAD_POINTS<=5 implemented
#define NUM_THREADS 1
#define FINITE_STRAIN_FLAG true

template <int dim>
class elasticity 
{
public:
  elasticity (NURBSMesh<dim>& _mesh, std::string& filePrefix, parametersClass* _params);
  ~elasticity(){}
  void run();
  parametersClass* params;
private:
  unsigned int numKnots;
  void solve();
  void apply_boundary_conditions();
  void condenseKandRHS();
  void apply_initial_values();
  void mark_boundaries();
  void mark_lines();
  void setup();
  void setupCellValues();
  void output (unsigned int _cycle);
  void assemble_system ();
  void assemble_system_interval (const typename std::vector<knotSpan<dim> >::iterator &begin, const typename std::vector<knotSpan<dim> >::iterator &end);
  NURBSMesh<dim>* mesh;
  std::vector<IGAValues<dim>*> cellValues;
  sparsityPattern sparsity_pattern;
  sparseMatrix system_matrix;
  denseVector  system_rhs, U, Un, dU;
  std::map<unsigned int, double> dirichletMap;
  std::map<double, unsigned int> outputMap;
  //Solver options
  unsigned int numIncrements, currentIncrement, currentIteration;
  //solution variables
  std::string filePrefix;
  solutionClass<dim> displacement, derivedValueP, derivedValueR;
  std::vector<solutionClass<dim>* > outputVariables, outputVariablesR;
  double dt;
  dealii::Threads::ThreadMutex     assembler_lock;
  double freeEnergy, interfaceEnergy, maxDisp, finalNorm, dirchletBC;
  std::ofstream energyFile;
};

/*elasticity:construction function.
 * Define input and output parameters/variables.
 */ 
template <int dim>
elasticity<dim>::elasticity (NURBSMesh<dim>& _mesh, std::string& _filePrefix, parametersClass* _params): mesh(&_mesh), filePrefix(_filePrefix), params(_params), \
													 displacement(_mesh, NODAL, VECTOR, std::string("u")), \
													 derivedValueP(_mesh, QUADRATURE, SCALAR, std::string("Pyy")), \
													 derivedValueR(_mesh, NODAL, VECTOR, std::string("R")) {
  numKnots=params->getInt("knots");
  //analysis variables
  //solution variables
  numIncrements=1; currentIncrement=0;//define load step
  dt=1.0/numIncrements;
  //output variables
  outputVariables.push_back(&displacement);
  outputVariables.push_back(&derivedValueP);
  outputVariables.push_back(&derivedValueR);
}

/*assemble_system
 * Initialize assemble process for "global Stiffiness Matrix" and "right hand term".  
 */
template <int dim>
void elasticity<dim>::assemble_system (){
  system_matrix=0.0; system_rhs=0.0; dU=0.0;
  freeEnergy=0.0; interfaceEnergy=0.0; maxDisp=0.0; dirchletBC=0.0;
  const unsigned int n_threads=dealii::multithread_info.n_default_threads;
  dealii::Threads::ThreadGroup<> threads;
  typedef typename std::vector<knotSpan<dim> >::iterator knotSpan_iterator;
  std::vector<std::pair<knotSpan_iterator,knotSpan_iterator> > thread_ranges = dealii::Threads::split_range<knotSpan_iterator> (mesh->knotSpanVector.begin(), mesh->knotSpanVector.end(), n_threads);
  printf("start assemble\n");
  for (unsigned int thread=0; thread<n_threads; ++thread){
    threads += dealii::Threads::new_thread (&elasticity<dim>::assemble_system_interval, *this, thread_ranges[thread].first, thread_ranges[thread].second);
  }
  threads.join_all ();
  printf("end assemble\n");
  printf("\nStrainEnergy:%12.4e, InterfaceEnergy:%12.6e\n", freeEnergy, interfaceEnergy);
  condenseKandRHS();//Adjust K and RHS to reflect BC's
}

/*assemble_system_interval
 * Assemble local matrix and local RHS by looping over every cell/element.
 * Using algorithmic (or automatic) differentiation (AD) to linearize Residual and compute the Jacobian matrix. 
 */
template <int dim>
void elasticity<dim>::assemble_system_interval (const typename std::vector<knotSpan<dim> >::iterator &begin, const typename std::vector<knotSpan<dim> >::iterator &end){
  //element loop
  IGAValues<dim> fe_values_base(mesh, dim, 2);//IGAValues(NURBSMesh,dofPerControlPoint,numberOfDerivatives)
  for (typename std::vector<knotSpan<dim> >::iterator cell=begin; cell<end; cell++){
    fe_values_base.reinit(*cell);
    IGAValues<dim>* fe_values=&fe_values_base;
    //IGAValues<dim>* fe_values=cellValues[cell->id];
    unsigned int n_q_points= fe_values->n_quadrature_points;
    unsigned int dofs_per_cell=fe_values->dofs_per_cell;
    denseMatrix local_matrix(dofs_per_cell, dofs_per_cell);
    denseVector local_rhs(dofs_per_cell);
    //AD variables
    dealii::Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); dealii::Table<1, double > ULocalConv(dofs_per_cell);
    dealii::Table<1, double> ULocalTemp(dofs_per_cell); 
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      ULocal[i]=U(cell->local_dof_indices[i]);
      ULocalTemp[i]=U(cell->local_dof_indices[i]);
      ULocal[i].diff (i, dofs_per_cell);
      ULocalConv[i]= Un(cell->local_dof_indices[i]);
    }
    //evaluate maximum displacement
    double maxDispLocal=0.0;
    dealii::Table<2,double> UTemp(n_q_points, dim);
    evaluateVectorFunction<double, dim>(*fe_values, 0, ULocalTemp, UTemp);
    for (unsigned int q=0; q<n_q_points; ++q){
      double uVal=0.0;
      for (unsigned int i=0; i<dim; ++i){
	uVal+=std::pow(UTemp[q][i],2);
      }
      uVal=std::pow(uVal, 0.5);      
      maxDispLocal=std::max(maxDispLocal, uVal);
    }    
    deformationMapwithGrad<Sacado::Fad::DFad<double>, dim> defMap(n_q_points);//structural variable defined in functionEvaluations.h from IGAbase
    getDeformationMapWithGradient<Sacado::Fad::DFad<double>, dim>(*fe_values, 0, ULocal, defMap, currentIteration);//function in mechanics.h
    dealii::Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); R.fill(0.0);//Residual for mechanics 
    double localfreeEnergy=0.0, localinterfaceEnergy=0.0, localdirchletBC=0.0;
    residualForMechanics(*fe_values, 0, FINITE_STRAIN_FLAG, ULocal, ULocalConv, R, currentIteration, localfreeEnergy, localinterfaceEnergy, localdirchletBC, *cell, numKnots, params, mesh, &derivedValueP);//function in mechanics.h
			
    //Residual(R) and Jacobian(R')
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      for (unsigned int j=0; j<dofs_per_cell; ++j){
	// R' by AD
	local_matrix(i,j)= R[i].fastAccessDx(j);
      }
      //R
      local_rhs(i) = -R[i].val();
    }

		
    //Global Assembly
    assembler_lock.acquire();
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      for (unsigned int j=0; j<dofs_per_cell; ++j){
	system_matrix(cell->local_dof_indices[i], cell->local_dof_indices[j])+=local_matrix(i,j);
      }
      system_rhs(cell->local_dof_indices[i]) += local_rhs(i);
    }
    freeEnergy+=localfreeEnergy; interfaceEnergy+=localinterfaceEnergy; 
    dirchletBC=std::max(dirchletBC, localdirchletBC);
    maxDisp=std::max(maxDisp, maxDispLocal);
    assembler_lock.release();
  }
}

/*solve
 * Using Newton method to solve nonlinear problem.
 * Define parameters for Solver.
 */
template <int dim>
void elasticity<dim>::solve (){
  double res=1, tol=1.0e-8, abs_tol=1.0e-11, initial_norm=0, current_norm=0;
  currentIteration=0;//interation mark in Newton method
  double oldNorm=1;
  while (true){
    if ((currentIteration>=20) or (res<1.0e-15)){printf ("Maximum number of iterations reached without convergence. \n"); break; exit (1);}
    double currentNorm=dU.l2_norm();
    //if (std::abs(currentNorm-oldNorm)<1.0e-14){printf ("Minimum change on solution norm: %12.e8\n", std::abs(currentNorm-oldNorm)); break; exit (1);}
    //oldNorm=currentNorm;
    if (current_norm>1/std::pow(tol,2)){printf ("\nNorm is too high. \n\n"); exit (1);}
    assemble_system();
    current_norm=system_rhs.l2_norm(); 
    initial_norm=std::max(initial_norm, current_norm);
    //if (initial_norm==0){printf ("Initial norm is zero.\n\n"); exit (1);}
    res=current_norm/initial_norm;
    printf ("Inc:%3u, Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e, UNorm: %12.8e \n", currentIncrement, currentIteration, current_norm, res, currentNorm); 
    finalNorm=res;
    if (res<tol || current_norm< abs_tol){printf ("Residual converged in %u iterations.\n", currentIteration); break;}
    //call solver in superLUMT_solver.h from IGAbase
    int status=luSolver(&system_matrix.nonZeroValues.at(0), &sparsity_pattern.columnIndices.at(0), &sparsity_pattern.rowIndex.at(0), &system_rhs.values.at(0), (int) system_rhs.size(), (int) system_rhs.size(), (int) sparsity_pattern.nnz, &dU.values.at(0), 0);
    if (status!=0) {printf("Solver exit status:%u\n", status); exit(-1);}
    U+=dU;
    ++currentIteration;
  }
  Un=U;
  //fill displacement Vector
  for (typename std::vector<controlPoint<dim> >::iterator i=mesh->controlPointVector.begin(); i<mesh->controlPointVector.end(); i++){
    for (unsigned int a=0; a<dim; a++){
      displacement(i->id, a)=Un(i->id*mesh->dofPerControlPoint + a);//the first dim nodal values are the displacement components by convention
    }
  }
}

/*apply_initial_values
 * Set initial "gussing" displacement value for Newton method.
 * Default setting is 0.
 */
template <int dim>
void elasticity<dim>::apply_initial_values(){
  Un=0.0; U=0.0;
  //fill results Vector
  for (typename std::vector<controlPoint<dim> >::iterator i=mesh->controlPointVector.begin(); i<mesh->controlPointVector.end(); i++){
    for (unsigned int a=0; a<dim; a++){
      displacement(i->id, a)=Un(i->id*mesh->dofPerControlPoint + a);
      derivedValueR(i->id, a)=0.0;//initialize output varialbe of Residual.
    }
  }
}

/*mark_boundaries
 * Mark six/four boundary surfaces.
 * Marking Rule:
 * X_i  boundaryFlags_Index  Default_Value  Activated_Value
 *  0		i*2+0		0		i*2+1
 *  1		i*2+1		0		i*2+2
 * -------------------------------------------------------
 * Boundary surfaces applied Neumann B.C must be activated. 
 */
template <int dim>
void elasticity<dim>::mark_boundaries(){
  for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
    for (unsigned int i=0;i<dim; i++){
      if (cell->endKnots[i][0]==0) cell->boundaryFlags[i*2+0]=i*2+1;
      if (cell->endKnots[i][1]==1) cell->boundaryFlags[i*2+1]=i*2+2;
    }
  }
}

/*mark_lines
 * Mark line that parallel to cell edge. 
 * Following code marks twelve/four boundary edges. 
 */
template <int dim>
void elasticity<dim>::mark_lines(){
  int lineMark[2]; 
  int Nmark=0;
  for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
    for (unsigned int i=0;i<dim; i++){
      Nmark=0;
      for(unsigned int j=0;j<dim;j++){
      	if(j!=i){
          lineMark[Nmark]=j;
          Nmark++;
        }
      }
      if(dim==3){
      	if ((cell->endKnots[lineMark[0]][0]==0)and(cell->endKnots[lineMark[1]][0]==0)) cell->lineFlags[i*4+0]=i*4+1;
      	if ((cell->endKnots[lineMark[0]][1]==1)and(cell->endKnots[lineMark[1]][0]==0)) cell->lineFlags[i*4+1]=i*4+2;
      	if ((cell->endKnots[lineMark[0]][0]==0)and(cell->endKnots[lineMark[1]][1]==1)) cell->lineFlags[i*4+2]=i*4+3;
      	if ((cell->endKnots[lineMark[0]][1]==1)and(cell->endKnots[lineMark[1]][1]==1)) cell->lineFlags[i*4+3]=i*4+4;
      }
     if(dim==2){
        if (cell->endKnots[lineMark[0]][0]==0) cell->lineFlags[i*2+0]=i*2+1;
        if (cell->endKnots[lineMark[0]][1]==1) cell->lineFlags[i*2+1]=i*2+2;
     }
   }
  }
}

/*apply_boundary_conditions
 * Following code sets dirichlet B.C of "crack" example.
 */
//Apply boundaries
template <int dim>
void elasticity<dim>::apply_boundary_conditions(){
  //Dirichlet map
  dirichletMap.clear(); outputMap.clear();
  unsigned int  controlPointDOF=-dim;
  for (typename std::vector<controlPoint<dim> >::iterator controlpoint=mesh->controlPointVector.begin(); controlpoint<mesh->controlPointVector.end(); controlpoint++){
    std::vector<double> coords(controlpoint->coords);
    controlPointDOF+=dim;
    double alpha=((double) currentIncrement)/numIncrements;
    //apply dirichlet B.C's
    //if (coords[dim-1]==0.0) {
      //if (coords[0]==0.0) dirichletMap[controlPointDOF+0]=0.0;
      //if (coords[1]==0.0) dirichletMap[controlPointDOF+1]=0.0;
      //if (coords[2]==0.0) dirichletMap[controlPointDOF+2]=0.0;
      //for (unsigned int i=0; i<dim; i++) dirichletMap[controlPointDOF+i]=0.0;  
      //}
    //if (coords[2]==10) {dirichletMap[controlPointDOF+1]=alpha*LOAD;}
    //crack BVP
    if (coords[0]==1.0) {
       dirichletMap[controlPointDOF+0]=0.0;
    }
    if ((coords[1]==0.0) and (coords[0]>=0.5)) {
       dirichletMap[controlPointDOF+1]=0.0;
    }
    if (coords[2]==0.0) {
       dirichletMap[controlPointDOF+2]=0.0;
    }
    //mark specific cell and the components(x,y,z) that have variable to output.
    if ((coords[1]==0.0) and (coords[2]==0.0)){
      outputMap[coords[0]]=controlPointDOF+1;
    }
    //Apply values to solution vector
    for (std::map<unsigned int, double>::iterator dof=dirichletMap.begin(); dof!=dirichletMap.end(); dof++){
      U(dof->first)=dof->second;
    }
  }
}

/*condenseKandRHS
 * Adjust K and RHS to reflect Dirichlet BC's.
 */
template <int dim>
void elasticity<dim>::condenseKandRHS(){
  //Storing Residual value for output before condensing RHS
  for (typename std::vector<controlPoint<dim> >::iterator i=mesh->controlPointVector.begin(); i<mesh->controlPointVector.end(); i++){
   for (unsigned int a=0; a<dim; a++){
     derivedValueR(i->id, a)=system_rhs(i->id*mesh->dofPerControlPoint + a);
   }
  }
  
  //Adjust K by applying dirichlet map which stores dof that has dirichlet BC's.   
  for (std::map<unsigned int, double>::iterator dof=dirichletMap.begin(); dof!=dirichletMap.end(); dof++){
    unsigned int i=dof->first;
    system_rhs(i)=0.0;
    for (unsigned int j=0; j<system_rhs.size(); j++){
      if ((i!=j) && (sparsity_pattern.nzMap.at(i).count(j)>0)) {system_matrix(i, j)=system_matrix(j, i)=0.0;}
    } 
  }
}

/*setupCellValues
 * Set up fe_value of each cell.
 * Normally no need here.
 */
template <int dim>
void elasticity<dim>::setupCellValues(){
  printf("setting up cell values\n");
  for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
    IGAValues<dim>* tempCellValues=new IGAValues<dim>(mesh, dim, 2);
    tempCellValues->reinit(*cell);
    cellValues.push_back(tempCellValues);
  }
}

/*setup
 * Initialize global data structures
 */
template <int dim>
void elasticity<dim>::setup (){	
  sparsity_pattern.init(mesh);
  system_matrix.reinit(sparsity_pattern);
  system_rhs.reinit(sparsity_pattern); U.reinit(sparsity_pattern); Un.reinit(sparsity_pattern); dU.reinit(sparsity_pattern); 
}

/*output
 * Generate vts file.
 */
template <int dim>
void elasticity<dim>::output (unsigned int _cycle){
  derivedValueP.projectQuadratureValues();
  
  //Generate output mesh
  char fileName[200];
  std::sprintf (fileName, "results/%s/%s/%u/lam-%4.2fmu-%4.2fmuSG-%4.2fl-%4.2fload-%4.2fG-%4.2fC-%4.2fWeak-%uFS-%s", \
		params->getString("bcType").c_str(), \
		params->getString("order").c_str(),  \
		params->getInt("knots"),  \
		params->getDouble("lambda"), \
		params->getDouble("mu"), \
		params->getDouble("muSG"), \
		params->getDouble("l"), \
		params->getDouble("load"), \
		params->getDouble("Gamma"), \
		params->getDouble("C"),			\
		params->getBool("enforceWeakBC"),	\
		(FINITE_STRAIN_FLAG)?"true":"false");
  writeMesh<dim>(fileName, _cycle, mesh, 2*numKnots+1, outputVariables);
}

/*run
 * 
 */
template <int dim>
void elasticity<dim>::run (){
  setup();
  mark_boundaries();
  //mark_lines();
  //setupCellValues();
  std::ofstream dataFileP, dataFileR;
  char fileNameP[200], fileNameR[200];

  std::sprintf (fileNameP, "results/%s/%s/%u/P.txt", \
		params->getString("bcType").c_str(), \
		params->getString("order").c_str(),  \
		params->getInt("knots"));
  std::sprintf (fileNameR, "results/%s/%s/%u/R.txt", \
		params->getString("bcType").c_str(), \
		params->getString("order").c_str(),  \
		params->getInt("knots"));
  dataFileP.open(fileNameP); dataFileR.open(fileNameR); 
  double values[]={0.0, 0.01, 0.1, 1.0, 10.0};//value for l.
  for (unsigned int j=0; j<1; ++j){
    params->setBool("enforceWeakBC", (j==1)? true:false);   
    for (unsigned int i=0; i<sizeof(values)/sizeof(double); ++i){
      params->setDouble("l",values[i]);
      apply_initial_values();
      //output(0);
      for (currentIncrement=1; currentIncrement<=numIncrements; ++currentIncrement){
	apply_boundary_conditions();
	solve(); 
	output(currentIncrement);
      }
      char dataP[100], dataR[100];
      if (i==0){
	std::sprintf(dataP, "l:%7.3e ", 0.0);
	std::sprintf(dataR, "l:%7.3e ", 0.0);
	dataFileP << dataP; dataFileR << dataR;
	for (std::map<double,unsigned int>::iterator it=outputMap.begin(); it!=outputMap.end(); ++it){
	  std::sprintf(dataP, "%12.5e ", it->first);
	  std::sprintf(dataR, "%12.5e ", it->first);
	  dataFileP << dataP; dataFileR << dataR;
	}
	dataFileP << "\n"; dataFileR << "\n";
      }
      //std::sprintf(energyValue, "l:%7.3e wBC:%s maxU:%12.6e SE:%12.6e GE:%12.6e Bnn:%12.6e Norm:%12.6e\n", params->getDouble("l"), (j==1)?"true":"fals", maxDisp, freeEnergy, interfaceEnergy, dirchletBC, finalNorm);
      std::sprintf(dataP, "l:%7.3e ", params->getDouble("l"));
      std::sprintf(dataR, "l:%7.3e ", params->getDouble("l"));
      dataFileP << dataP; dataFileR << dataR;
      for (std::map<double,unsigned int>::iterator it=outputMap.begin(); it!=outputMap.end(); ++it){
	std::sprintf(dataP, "%12.5e ", derivedValueP.projectedValues.at(it->second/dim));
	std::sprintf(dataR, "%12.5e ", derivedValueR(it->second/dim,1));
	dataFileP << dataP; dataFileR << dataR;
      }
      dataFileP << "\n"; dataFileR << "\n";
      dataFileP.flush(); dataFileR.flush();
    }
  }
  dataFileP.close(); dataFileR.close();
}

/*main
 * 
 */
int main(int argc, char *argv[]){
  std::clock_t start=std::clock();
  dealii::multithread_info.n_default_threads=NUM_THREADS;

  //set param values
  parametersClass params;
  params.setDouble("lambda", 1.0);
  params.setDouble("mu", 1.0);
  params.setDouble("muSG", 1.0);
  params.setDouble("l", 0.1);
  params.setDouble("load", 0.1);
  params.setDouble("Gamma", 0.0);
  params.setDouble("C", 5.0);
  params.setString("bcType", "tension");
  params.setString("order", "Quadratic");
  params.setBool("enforceWeakBC", true);
  params.setInt("knots", 50);
  //params.readInParameters("params.txt");
  
  //read environmental variables
  printf("reading environmental variables...\n");
  //char *bcType, *order, *knots, *enforceWeakBC;
  //bcType = getenv ("bcType"); params.setString("bcType", bcType); 
  //order = getenv ("order"); params.setString("order", order);
  //knots = getenv ("knots"); params.setInt("knots", std::atoi(knots));
  //enforceWeakBC = getenv ("enforceWeakBC"); params.setBool("enforceWeakBC", (std::atoi(enforceWeakBC)==1)?"true":"false");

  //NURBS file prefix
  char fileName[100];
  std::sprintf (fileName, "%uD%s%u", DIMS, params.getString("order").c_str(), params.getInt("knots"));
  std::string filePrefix(fileName);
		
  //Read NURBS geometry  
  nativeNURBSStructure<DIMS> geometry;
  char meshFile[100];
  std::sprintf (meshFile, "meshes/IGAMesh%s.h5", filePrefix.c_str());
  readHDF5<DIMS>(meshFile, geometry); 

  //Generate IGA mesh and data structures (control nodes, knot spans, basis functions, etc)
  NURBSMesh<DIMS> mesh(geometry, DIMS, NUM_QUAD_POINTS);
	
  //Problem
  elasticity<DIMS> problem(mesh, filePrefix, &params);
  problem.run();
	
  //Stats
  printf ("\nTime taken:%10.2e sec\n", (std::clock()-start)/((double)CLOCKS_PER_SEC));
}

