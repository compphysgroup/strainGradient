/*
 * solutionClasses.h
 *  Created on: Oct 28, 2011
 */
#ifndef SOLUTIONCLASSES_H_
#define SOLUTIONCLASSES_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include "nurbsClasses.h"
#include "sparseDataClass.h"
#ifdef SOLVER_MT
#include "superLUMT_solver.h"
#else
#include "superLU_solver.h"
#endif
using namespace std;

enum dataType {SCALAR, VECTOR, TENSOR};
enum dataLocation {NODAL, QUADRATURE};

template <int dim>
class solutionClass{
	private:
		NURBSMesh<dim>* mesh;
		
	public:
		solutionClass(NURBSMesh<dim>& _mesh, dataLocation location, dataType data, std::string  _variableName);
		~solutionClass();
		void createMassMatrix();
		void projectQuadratureValues();
		std::vector<double> values;
		std::vector<double> projectedValues;
		dataLocation datalocation;
		dataType datatype;
		//operator  overloading
		double& operator() (unsigned int _a) {
			if ((datalocation==NODAL) && (datatype==SCALAR)) return values.at(_a);
			else {printf("solutionClass: incompatible arguments\n"); exit(-1);}
		}
		double& operator() (unsigned int _a, unsigned int _b) {
			if (datalocation==NODAL){
				if (datatype==VECTOR) return values.at(_a*dim+_b);
				else if (datatype==TENSOR) return values.at(_a*dim*dim+_b);
				else {printf("solutionClass: incompatible datatype\n"); exit(-1);}
			}
			else if ((datalocation==QUADRATURE) && (datatype==SCALAR)) return values.at(_a*numQuadPoints +_b);
			else {printf("solutionClass: incompatible arguments\n"); exit(-1);}
		}
		double& operator() (unsigned int _a, unsigned int _b, unsigned int _c) {
			if (datalocation==QUADRATURE){
				if (datatype==VECTOR) return values.at(_a*numQuadPoints*dim+_b*dim+_c);
				else if (datatype==TENSOR) return values.at(_a*numQuadPoints*dim*dim+_b*dim*dim+_c);
				else {printf("solutionClass: incompatible datatype\n"); exit(-1);}
			}
			else {printf("solutionClass: incompatible arguments\n"); exit(-1);}
		}
		unsigned int numQuadPoints;
		unsigned int numVariablesPerPoint;
		std::string variableName;
		static sparseMatrix* massMatrix;
		static sparsityPattern* mass_sparsity_pattern;
};

template <int dim>
sparseMatrix* solutionClass<dim>::massMatrix=0;

template <int dim>
sparsityPattern* solutionClass<dim>::mass_sparsity_pattern=0;

template <int dim>
solutionClass<dim>::~solutionClass(){
	if (mass_sparsity_pattern) {delete mass_sparsity_pattern; mass_sparsity_pattern=0;}
	if (massMatrix) {delete massMatrix; massMatrix=0;}
}

template <int dim>
solutionClass<dim>::solutionClass(NURBSMesh<dim>& _mesh, dataLocation _datalocation, dataType _datatype, std::string _variableName): mesh(&_mesh), datalocation(_datalocation), datatype(_datatype), variableName(_variableName){
	numQuadPoints=std::pow(mesh->quadPtStencilSize,dim);
	if (datalocation==NODAL){
		if (datatype==SCALAR) {values.resize(mesh->controlPointVector.size(), 0.0); numVariablesPerPoint=1;}
		else if (datatype==VECTOR) {values.resize(mesh->controlPointVector.size()*dim, 0.0); numVariablesPerPoint=dim;}
		else if (datatype==TENSOR) {values.resize(mesh->controlPointVector.size()*dim*dim, 0.0); numVariablesPerPoint=dim*dim;}
		else {printf("unknown dataType\n"); exit(-1);}
	}
	else if (datalocation==QUADRATURE){
		if (datatype==SCALAR) {values.resize(mesh->knotSpanVector.size()*numQuadPoints, 0.0); numVariablesPerPoint=1;}
		else if (datatype==VECTOR) {values.resize(mesh->knotSpanVector.size()*dim*numQuadPoints, 0.0); numVariablesPerPoint=dim;}
		else if (datatype==TENSOR) {values.resize(mesh->knotSpanVector.size()*dim*dim*numQuadPoints, 0.0); numVariablesPerPoint=dim*dim;}
		else {printf("unknown dataType\n"); exit(-1);}
		if (mass_sparsity_pattern==0) {mass_sparsity_pattern=new sparsityPattern(); mass_sparsity_pattern->init(mesh, true);}
		if (massMatrix==0) createMassMatrix();
	}
	else {printf("unknown dataLocation\n"); exit(-1);}
}

template <int dim>
void solutionClass<dim>::createMassMatrix(){
	if (massMatrix) {printf("massMatrix already initialized\n"); exit(-1);}
	massMatrix=new sparseMatrix; 
	massMatrix->reinit(*mass_sparsity_pattern); *massMatrix=0;
	IGAValues<dim> fe_values(mesh, 1); 
		
	for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
		fe_values.reinit (*cell);
		unsigned int n_q_points= fe_values.n_quadrature_points;
		unsigned int dofs_per_cell=cell->local_mass_dof_indices.size();
		denseMatrix local_matrix(dofs_per_cell, dofs_per_cell);
			
		//Mass matrix
		for (unsigned int i=0; i<dofs_per_cell; ++i) {
			for (unsigned int j=0; j<dofs_per_cell; ++j){
				local_matrix(i,j)=0.0;
				for (unsigned int q=0; q<n_q_points; q++){
					local_matrix(i,j)+= fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q); // Mij= Ni*Nj
				}
				(*massMatrix)(cell->local_mass_dof_indices[i], cell->local_mass_dof_indices[j])+=local_matrix(i,j); //Global Assembly
			}
		}
	}
	//printf("mass matrix generated\n");
}

template <int dim>
void solutionClass<dim>::projectQuadratureValues(){
	if (datalocation!=QUADRATURE) {printf("datalocation!=QUADRATURE, so projection invalid\n"); exit(-1);}
	
	if (projectedValues.size()>0) projectedValues.clear();
	if (datatype==SCALAR) projectedValues.resize(mesh->controlPointVector.size(), 0.0);
	else if (datatype==VECTOR) projectedValues.resize(mesh->controlPointVector.size()*dim, 0.0);
	else if (datatype==TENSOR) projectedValues.resize(mesh->controlPointVector.size()*dim*dim, 0.0);
	
	//Create RHS
	denseVector system_rhs; system_rhs.reinit(*mass_sparsity_pattern);
	IGAValues<dim> fe_values(mesh, 1);
		
	//loop over all variables
	std::vector<double> tempValues(mesh->controlPointVector.size(), 0.0);
	for (unsigned int var=0; var<numVariablesPerPoint; var++){
		 system_rhs=0;
		 for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
			fe_values.reinit (*cell);
			unsigned int n_q_points= fe_values.n_quadrature_points;
			unsigned int dofs_per_cell=cell->local_mass_dof_indices.size();
			denseVector local_rhs(dofs_per_cell);
				
			//Mass matrix
			for (unsigned int i=0; i<dofs_per_cell; ++i) {
				local_rhs(i)=0.0;
				for (unsigned int q=0; q<n_q_points; q++){
					if (datatype==SCALAR) local_rhs(i)+= fe_values.shape_value(i,q)*((*this)(cell->id, q))*fe_values.JxW(q);
					else local_rhs(i)+= fe_values.shape_value(i,q)*((*this)(cell->id, q, var))*fe_values.JxW(q);
				}
				system_rhs(cell->local_mass_dof_indices[i]) += local_rhs(i); //Global Assembly
			}
		} 
		//project
		for(std::vector<double>::iterator it=tempValues.begin(); it<tempValues.end(); it++) *it=0.0; //initializing solution vector
		int status=luSolver(&massMatrix->nonZeroValues.at(0), &mass_sparsity_pattern->columnIndices.at(0), &mass_sparsity_pattern->rowIndex.at(0), &system_rhs.values.at(0), (int) system_rhs.size(), (int) system_rhs.size(), (int) mass_sparsity_pattern->nnz, &tempValues.at(0), 0);
		if (status!=0) {printf("Solver exit status:%u\n", status); exit(-1);}
		for (unsigned int a=0; a<mesh->controlPointVector.size(); a++){
			projectedValues.at(a*numVariablesPerPoint+var)=tempValues.at(a);
			if (a<10) std::cout << tempValues.at(a) << std::endl;
		}
	}
}


#endif /* SOLUTIONCLASSES_H_ */
