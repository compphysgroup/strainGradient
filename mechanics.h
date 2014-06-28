/*mechanics.h
 * Mechanics implementation on cell/element level.
 * Every function is on cell/element level.
 * Calculate F(deformation gradient tensor), P(Piola-Kirchhoff stress tensor), B(higher-order stress tensor) , freeEnergy.
 * Implement weak dirichlet condition.
 *--------------------------------------------
 * Things may need to take a look:
 * 1.material model
 * 2.Form of Newmann boundary condition
 */
#ifndef MECHANICS_H_
#define MECHANICS_H_
#include "igaClasses.h"
#include "base/table.h"
#include "functionEvaluations.h"
#include "supplementaryFunctions.h"
#include "solutionClasses.h"
#include "parameters.h"

/*getDeformationMapWithGradient
 * Calculate F, and it's inverse and derivatives.
 */
template <class T, int dim>
  void getDeformationMapWithGradient(IGAValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T>& ULocal, deformationMapwithGrad<T, dim>& defMap, const unsigned int iteration, int faceID=-1){
  //unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  //defaulting setting is faceID=-1 and this function is for normal quadrature point.
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;

  //evaluate dx/dX
  dealii::Table<3, T> GradU(n_q_points, dim, dim);
  dealii::Table<4, T> GradGradU(n_q_points, dim, dim, dim);
  evaluateVectorFunctionGradient<T, dim>(fe_values, DOF, ULocal, GradU, faceID);
  evaluateVectorFunctionSecondGradient<T, dim>(fe_values, DOF, ULocal, GradGradU, faceID);
	
  //Loop over quadrature points
  //F=GradU + I; derivatives of F=GradGradU
  for (unsigned int q=0; q<n_q_points; ++q){
    dealii::Table<2, T > Fq(dim, dim), invFq(dim, dim); T detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	defMap.F[q][i][j] = Fq[i][j] = (i==j) + GradU[q][i][j]; //F (as double value);
	for (unsigned int k=0; k<dim; ++k){
	  defMap.gradF[q][i][j][k] = GradGradU[q][i][j][k]; //gradF 
	}
      }
    }
		
    getInverse<T, dim>(Fq, invFq, detFq); //get inverse(F)
    defMap.detF[q] = detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
	defMap.invF[q][i][j] = invFq[i][j];
      }
    }
    //detF
    if (defMap.detF[q].val()<=1.0e-15 && iteration==0){
      printf("**************Non positive jacobian detected**************. Value: %12.4e\n", defMap.detF[q].val());
      for (unsigned int i=0; i<dim; ++i){
	for (unsigned int j=0; j<dim; ++j) printf("%12.6e  ", defMap.F[q][i][j].val());
	printf("\n"); exit(-1);
      }
      throw "Non positive jacobian detected";
    }
  }
}

/*evaluateStress
 * Calculate strain tensor E, P, B and free energy.
 */
template <class T, int dim>
  void evaluateStress(IGAValues<dim>& fe_values, const unsigned int DOF, bool finiteStrain, dealii::Table<3, T>& P, dealii::Table<4, T>& Beta, unsigned int currentIteration, double& freeEnergy, double& interfacEnergy, parametersClass* params, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, solutionClass<dim>* derivedValueP, unsigned int cellID,  int faceID=-1){
  //unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;

  deformationMapwithGrad<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
  getDeformationMapWithGradient<Sacado::Fad::DFad<double>, dim>(fe_values, 0, ULocal, defMap, currentIteration, faceID);
  
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){    
    //get F and it's derivatives explicitly
    T F[dim][dim], dF[dim][dim][dim];
    for (unsigned int i=0; i<dim; i++) {
      for (unsigned int j=0; j<dim; j++) {
	F[i][j]=defMap.F[q][i][j];	
	for (unsigned int k=0; k<dim; k++) {
	  dF[i][j][k]=defMap.gradF[q][i][j][k];
	}
      }
    }
    //Compute strain metric, E (E=0.5*(F^T*F-I)), it's derivatives (dE=2*F:dF) and some other values.
    T E[dim][dim], dE[dim][dim][dim], trE=0.0, trE2=0.0;
    for (unsigned int I=0; I<dim; I++){
      for (unsigned int J=0; J<dim; J++){
	E[I][J] = -0.5*(I==J);
	for (unsigned int k=0; k<dim; k++){
	  E[I][J] += 0.5*F[k][I]*F[k][J];
	}
	//compute derivatives of F
	for (unsigned int K=0; K<dim; K++){
	  dE[I][J][K]=0.0;
	  for (unsigned int k=0; k<dim; k++){
	    dE[I][J][K] += 0.5*(F[k][I]*dF[k][J][K]+F[k][J]*dF[k][I][K]);
	  }
	}
	trE2 += E[I][J]* E[I][J];//E:E
      }
      trE +=E[I][I];//E:E
    }

    //infinitesimal strain implementation( finiteStrain==false)
    //F=I (identity tensor); E=dx/dX
    if (!finiteStrain){
      trE=0.0; trE2=0.0;
      for (unsigned int i=0; i<dim; i++) {
	for (unsigned int j=0; j<dim; j++) {
	  F[i][j]=(i==j);	
	  E[i][j]=defMap.F[q][i][j] - (i==j);//defMap.F=I + GradU
	  for (unsigned int k=0; k<dim; k++) {
	    dE[i][j][k]=defMap.gradF[q][i][j][k];//GradGradU[q][i][j][k]
	  }
	  trE2 += E[i][j]*E[i][j];//E:E	
	}
	trE +=E[i][i]; //E:E  
      }
    }

    //compute P and Beta
    //standard St. Venant-Kirchhoff model for P
    //W (E, GradE) =1/2*(λ)*(Eaa)^2 + μ*(Eab * Eab) + 1/2*μ*l*(Eab,c * Eab,c)
    double lambda= params->getDouble("lambda");
    double mu=  params->getDouble("mu");
    double muSG= params->getDouble("muSG");
    double l= params->getDouble("l");
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int J=0; J<dim; ++J){
	P[q][i][J]=lambda*trE*F[i][J];
	for (unsigned int K=0; K<dim; ++K){
	  P[q][i][J]+= 2*mu*F[i][K]*E[K][J];
	  Beta[q][i][J][K]=0.0;
	  for (unsigned int A=0; A<dim; ++A){
	    P[q][i][J]+= mu*l*l*dF[i][K][A]*dE[J][K][A];
	    Beta[q][i][J][K]+=muSG*l*l*F[i][A]*dE[A][J][K];
	  }
	}
      }
    }
    
    //store certain values for output
    if (faceID==-1) {
      (*derivedValueP)(cellID,q)=P[q][dim-2][dim-2].val();
    }
    //store energies
    if (faceID==-1){
      double tempFreeEnergy=0.0, tempInterfaceEnergy=0.0;
      tempFreeEnergy= 0.5*lambda*std::pow(trE.val(),2)+ mu*trE2.val();
      for (unsigned int I=0; I<dim; ++I){
	for (unsigned int J=0; J<dim; ++J){
	  for (unsigned int K=0; K<dim; ++K){
	    tempInterfaceEnergy+=0.5*muSG*l*l*std::pow(dE[I][J][K].val(),2);
	  }	
	}
      }
      freeEnergy+=tempFreeEnergy*fe_values.JxW(q);
      interfacEnergy+=tempInterfaceEnergy*fe_values.JxW(q);
    }
  }
}

/*residualForMechanics
 * Mechanics residual R implementation
 * Applying Neumann boundary conditions 
 * Applying Weak dirichlet boundary condition
 */
template <int dim>
void residualForMechanics(IGAValues<dim>& fe_values, unsigned int DOF, bool finiteStrain, dealii::Table<1, Sacado::Fad::DFad<double> >& ULocal, dealii::Table<1, double>& ULocalConv, dealii::Table<1, Sacado::Fad::DFad<double> >& R, unsigned int currentIteration, double& freeEnergy, double& interfaceEnergy, double& dirchletBC, knotSpan<dim>& cell, unsigned int NumKnotInterval, parametersClass* params, NURBSMesh<dim>* mesh, solutionClass<dim>* derivedValueP){
  
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  //Temporary arrays
  dealii::Table<3,Sacado::Fad::DFad<double> > P (n_q_points, dim, dim);
  dealii::Table<4,Sacado::Fad::DFad<double> > Beta (n_q_points, dim, dim, dim);
  
  //evaluate mechanics parameters
  evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, DOF, finiteStrain, P, Beta, currentIteration, freeEnergy, interfaceEnergy, params, ULocal, derivedValueP, cell.id);
  
  //evaluate Residual
  for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
    R[dof] = 0;
    const unsigned int ck = fe_values.system_to_component_index(dof) - DOF;
    if (ck>=0 && ck<dim){
      // R = Grad(w)*P +GradGrad(w)*B
      for (unsigned int q=0; q<n_q_points; ++q){
	for (unsigned int J = 0; J < dim; J++){
	  R[dof] +=  fe_values.shape_grad(dof, q)[J]*P[q][ck][J]*fe_values.JxW(q);
	  for (unsigned int K = 0; K < dim; K++){
	    R[dof] +=  fe_values.shape_grad_grad(dof, q)[J][K]*Beta[q][ck][J][K]*fe_values.JxW(q);
	  }
	}
      }
    }
  }

  //Neumann  conditions
  double load= params->getDouble("load"); 
  const char* bcType= params->getString("bcType").c_str();
  
  //loop over boundary faces
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if (cell.boundaryFlags[faceID]>0){      
      for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
	//ck:specific dof of node. DOF make using flexible such at each node has many dof addition to three displacements.
	const unsigned int ck = fe_values.system_to_component_index(dof) - DOF;
	for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  if (std::strcmp(bcType,"bending")==0){
	    //bending along dim-1 direction
	    if ((cell.boundaryFlags[faceID]==2*dim) and (ck==dim-2)){
	      R[dof] += -fe_values.shape_value(dof, q, faceID)*load*fe_values.JxW(q, faceID);
	    }
	  }
	  else if (std::strcmp(bcType,"torsion")==0){
	    if (dim<3) throw "torsion B.C only valid for 3D problem";
	    //torque on z=1 face
	    if ((cell.boundaryFlags[faceID]==2*dim) and (ck<=dim-2)){
	      unsigned int fq=fe_values.quadratureMap[faceID][q]; //quadrature point index for face quadrature points
	      double x=fe_values.quadPointLocations[fq][0]-0.5, y=fe_values.quadPointLocations[fq][1]-0.5;
	      double r=std::sqrt(x*x+y*y);
	      double theta=std::atan2(y,x);
	      double loadX=load*r*std::sin(theta);
	      double loadY=-load*r*std::cos(theta);
	      if (ck==0) {
		R[dof] += -fe_values.shape_value(dof, q, faceID)*loadX*fe_values.JxW(q, faceID);
	      }
	      else if (ck==1) {
		R[dof] += -fe_values.shape_value(dof, q, faceID)*loadY*fe_values.JxW(q, faceID);
	      }
	    }
	  }
	  else if (std::strcmp(bcType,"tension")==0){
	    //tension along dim direction
	    if ((cell.boundaryFlags[faceID]==2*dim) and (ck==dim-1)){
	      R[dof] += -fe_values.shape_value(dof, q, faceID)*load*fe_values.JxW(q, faceID);
	    }
	  }
	  else if (std::strcmp(bcType,"line")==0){
	    if ((cell.boundaryFlags[1]>0) and (cell.boundaryFlags[3]>0)){ // line between faces (X=1,Y=1)
	      IGAValues<dim> fe_values_temp(mesh, dim, 0);
	      std::vector<std::vector<double> > quadPoints(2);
	      quadPoints[0].push_back(1.0); quadPoints[0].push_back(1.0); quadPoints[0].push_back(-1.0/std::pow(3, 0.5)); quadPoints[0].push_back(1.0);
	      quadPoints[1].push_back(1.0); quadPoints[1].push_back(1.0); quadPoints[1].push_back(1.0/std::pow(3, 0.5)); quadPoints[1].push_back(1.0);
	      fe_values_temp.reinit(cell, &quadPoints);
	      for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
		const unsigned int ck = fe_values.system_to_component_index(dof) - DOF;
		for (unsigned int q=0; q<2; ++q){
		  if (ck==dim-1){
		    R[dof] += -fe_values_temp.shape_value(dof, q)*load*(1.0/20);//fe_values_temp.JxW(q);
		  }
		}
	      }
	    }
	  }
	  else if  (std::strcmp(bcType,"crack")==0){
	    if ((cell.boundaryFlags[faceID]==2*(dim-1)) and (ck==dim-2)){
	      R[dof] += -fe_values.shape_value(dof, q, faceID)*load*fe_values.JxW(q, faceID);
	    }
	  }
	  else {
	    throw "unknown bcType";
	  }
	} 
      } 
    }
  }

  //Weak dirichlet condition grad(u).n=0
  //C: positive penalty parameter; he: characteristic mesh size parameter.
  if (params->getBool("enforceWeakBC")){ 
    double muSG= params->getDouble("muSG");
    double l= params->getDouble("l");
    double gamma= params->getDouble("Gamma");
    double C= params->getDouble("C");
    double he=1.0/NumKnotInterval;

    for (unsigned int faceID=0; faceID<2*dim; faceID++){
      if ((cell.boundaryFlags[faceID]==(dim-1)*2+1) or (cell.boundaryFlags[faceID]==(dim-1)*2+2)){
	//compute face normal (logic for computing n like this only works for cube geometries)
	std::vector<double> n(dim, 0);
	n[(cell.boundaryFlags[faceID]-1)/2]=std::pow(-1.0, (int)cell.boundaryFlags[faceID]%2);
	
	//Temporary arrays
	dealii::Table<3,Sacado::Fad::DFad<double> > PFace (fe_values.n_face_quadrature_points, dim, dim);
	dealii::Table<4,Sacado::Fad::DFad<double> > BetaFace (fe_values.n_face_quadrature_points, dim, dim, dim);
	
	//evaluate mechanics
	evaluateStress<Sacado::Fad::DFad<double>, dim>(fe_values, DOF, finiteStrain, PFace, BetaFace, currentIteration, freeEnergy, interfaceEnergy, params, ULocal, derivedValueP, cell.id, faceID);
 	
	//evaluate gradients on the faces
	dealii::Table<3,Sacado::Fad::DFad<double> > uij(fe_values.n_face_quadrature_points, dim, dim); uij.fill(0.0);
	dealii::Table<4,Sacado::Fad::DFad<double> > uijk(fe_values.n_face_quadrature_points, dim, dim, dim); uijk.fill(0.0);
	
	//Loop over quadrature points
	for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  for (unsigned int d=0; d<dofs_per_cell; ++d){
	    unsigned int i = fe_values.system_to_component_index(d) - DOF;
	    for (unsigned int j=0; j<dim; ++j){
	      uij[q][i][j]+=ULocal[d]*fe_values.shape_grad(d, q, faceID)[j];
	      for (unsigned int k=0; k<dim; ++k){
		uijk[q][i][j][k]+=ULocal[d]*fe_values.shape_grad_grad(d, q, faceID)[j][k];
	      }
	    }
	  }
	}

	//evaluate tensor multiplications 
	dealii::Table<2,Sacado::Fad::DFad<double> > uijn(fe_values.n_face_quadrature_points, dim); uijn.fill(0.0);
	dealii::Table<2,Sacado::Fad::DFad<double> > Betaijknn(fe_values.n_face_quadrature_points, dim); Betaijknn.fill(0.0);
	dealii::Table<2,double> wijn(fe_values.n_face_quadrature_points, dofs_per_cell); wijn.fill(0.0);
	dealii::Table<2,double> wijknn(fe_values.n_face_quadrature_points, dofs_per_cell); wijknn.fill(0.0);
	//evaluate uijn, Betaijknn
	for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  double tempdirchletBC=0.0;
	  for (unsigned int i=0; i<dim; ++i){
	    for (unsigned int j=0; j<dim; ++j){
	      uijn[q][i]+=uij[q][i][j]*n[j];
	      for (unsigned int k=0; k<dim; ++k){
		Betaijknn[q][i]+=BetaFace[q][i][j][k]*n[j]*n[k];
	      }
	    }
	    tempdirchletBC+=std::pow(uijn[q][i].val(),2.0);
	  }
	  dirchletBC=std::max(dirchletBC, std::pow(tempdirchletBC,0.5));
	  //gradun(fe_values.cell->id,q,i)= uijn[q][i].val();
	}
	//evaluate wijn, wijknn
	for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  for (unsigned int d=0; d<dofs_per_cell; ++d){
	    for (unsigned int j=0; j<dim; ++j){
	      wijn[q][d]+= fe_values.shape_grad(d, q, faceID)[j]*n[j];
	      for (unsigned int k=0; k<dim; ++k){
		wijknn[q][d]+= fe_values.shape_grad_grad(d, q, faceID)[j][k]*n[j]*n[k];
	      }
	    }
	  }	  
	}
	
	//Add the weak dirichlet terms
	for (unsigned int d=0; d<dofs_per_cell; ++d) {
	  unsigned int i = fe_values.system_to_component_index(d) - DOF;
	  if (i==dim-1){ //enforcing weak dirichlet only along dim direction
	    for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){	
	      //-mu*l*l*w_(ij)n_j*Beta_(ijk)n_jn_k
	      R[d] += -wijn[q][d]*Betaijknn[q][i]*fe_values.JxW(q, faceID);
	      //-mu*l*l*gamma*wijknn*uijn //WRONG form curretnly. Not used now as there is no point trying to make an unsymmetric euqation adjoint consistent. So gamma is always set to zero.
	      gamma=0.0;
	      R[d] += -gamma*wijknn[q][d]*uijn[q][i]*fe_values.JxW(q, faceID);
	      //mu*l*l*(C/he)*wijn*uijn
	      R[d] += (C/he)*wijn[q][d]*uijn[q][i]*fe_values.JxW(q, faceID);
	      //std::cout << uijn[q][i].val() << " ";
	    }
	  }
	} 
      }
    } 
  }
  
}

#endif /* MECHANICS_H_ */
