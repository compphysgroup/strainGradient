/*
 * nurbsClasses.h
 *  Created on: Oct 20, 2011
 */

#ifndef NURBSCLASSES_H_
#define NURBSCLASSES_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <set>
#include <lac/full_matrix.h>
#include <base/tensor.h>
#include "supplementaryFunctions.h"
class sparsityPattern;

//forward declaration
template <int dim>
class IGAValues;

template <int dim>
class NURBSMesh;

//Nurbs structure class
template <int dim>
class nativeNURBSStructure
{
	public:
		nativeNURBSStructure();
		std::vector<unsigned int> order;
		std::vector<unsigned int> numCtrlPts;
		std::vector<std::vector<double> > knots;
		std::vector<std::vector<double> > ctrlPtCoefs;
		std::vector<std::vector<unsigned int> > ctrlPtKnots;
};
template <int dim>
nativeNURBSStructure<dim>::nativeNURBSStructure(){order.resize(dim); numCtrlPts.resize(dim); knots.resize(dim); ctrlPtCoefs.resize(dim+1); ctrlPtKnots.resize(dim);}

//basis function class
template <int dim>
class basisFunction
{
	private:
		double evaluateN(unsigned int i, unsigned int localKnotIndex, unsigned int _p, double xi, unsigned int k=0);
	public:
		double weight;
		std::vector<unsigned int> p; //p_i: order  along each i^th direction
		std::vector<std::vector<double> > localKnots; //p_i+2 localKnots vector along each i^th direction
		basisFunction();
		bool pointInBasisSupport(std::vector<double>& coords);
		double value(std::vector<double>& coords, std::vector<unsigned int>& dervs);
};
template <int dim>
basisFunction<dim>::basisFunction(): weight(0.0) {p.resize(dim); localKnots.resize(dim);}

template <int dim>
bool basisFunction<dim>::pointInBasisSupport(std::vector<double>& coords){
	bool flag=true;
	for (unsigned int i=0; i<dim; i++){
		if  (!((coords[i]>=localKnots[i][0] && (coords[i]<=localKnots[i].back())))) flag=false;
	}
	return flag;
}

template <int dim>
double basisFunction<dim>::evaluateN(unsigned int i, unsigned int localKnotIndex, unsigned int _p, double xi, unsigned int k){
	double val=0.0;
	//Cox-DeBoor recursion formula
	if (_p==0){
		if ((xi>=localKnots[i][localKnotIndex]) && (xi<=localKnots[i][localKnotIndex+1]) && (k==0)) val=1;
		else val=0;
	}
	else{
		if (k>0){ 
			//recursion formula for derivative of basis function
			if (localKnots[i][localKnotIndex+_p]!=localKnots[i][localKnotIndex])     val+=   evaluateN(i, localKnotIndex, _p-1, xi, k-1)*(_p)/(localKnots[i][localKnotIndex+_p]-localKnots[i][localKnotIndex]);
			if (localKnots[i][localKnotIndex+_p+1]!=localKnots[i][localKnotIndex+1]) val+=-evaluateN(i, localKnotIndex+1, _p-1, xi, k-1)*(_p)/(localKnots[i][localKnotIndex+_p+1]-localKnots[i][localKnotIndex+1]);
		}
		else{
			//recursion formula for evaluation of basis function
			if (localKnots[i][localKnotIndex+_p]!=localKnots[i][localKnotIndex]) val+=evaluateN(i, localKnotIndex, _p-1, xi)*(xi-localKnots[i][localKnotIndex])/(localKnots[i][localKnotIndex+_p]-localKnots[i][localKnotIndex]);
			if (localKnots[i][localKnotIndex+_p+1]!=localKnots[i][localKnotIndex+1]) val+=evaluateN(i, localKnotIndex+1, _p-1, xi)*(localKnots[i][localKnotIndex+_p+1]-xi)/(localKnots[i][localKnotIndex+_p+1]-localKnots[i][localKnotIndex+1]);
		}
	}
	return val;
}

template <int dim>
double basisFunction<dim>::value(std::vector<double>& coords, std::vector<unsigned int>& dervs){ 
	double temp=1.0;
	for (unsigned int a=0; a<dim; a++){
		temp*=evaluateN(a, 0, p[a], coords[a], dervs[a]);
	}
	return temp*weight;
}

//control point class
template <int dim>
class controlPoint
{
	public:
		unsigned int id;
		controlPoint(unsigned int _id, nativeNURBSStructure<dim>& mesh, std::vector<unsigned int>& _dofs);
		std::vector<unsigned int> dofs;
		std::vector<double> coords;
		basisFunction<dim> basis; 
};
template <int dim>
controlPoint<dim>::controlPoint(unsigned int _id, nativeNURBSStructure<dim>& mesh, std::vector<unsigned int>& _dofs): id(_id), dofs(_dofs), basis(){
	coords.resize(dim); dofs.resize(dim);
	//fill basis variables
	basis.weight=mesh.ctrlPtCoefs[dim].at(id);
	for (unsigned int a=0; a<dim; a++) {
		//fill control point variables
		coords[a]=mesh.ctrlPtCoefs[a].at(id)/basis.weight; //the division by basis.weight is required as octave nurbs outputs coords multiplied by the weights. Not required if nurbs coefficients obtained from any other source.
		//fill basis variables
		basis.p[a]=mesh.order[a];
		for (unsigned int b=0; b<mesh.order[a]+2; b++) basis.localKnots[a].push_back(mesh.knots[a].at(mesh.ctrlPtKnots[a][id]+b));
	}
}

//knot span class
template <int dim>
class knotSpan
{
	public:
		unsigned int id;
		std::vector<unsigned int> boundaryFlags;
		std::vector<std::vector<double> > endKnots;
		std::vector<std::vector<double> > edgeCoords;
		std::vector<double> center;
		std::vector<controlPoint<dim>* > controlPoints;
		std::vector<unsigned int> local_dof_indices;
		std::vector<unsigned int> local_mass_dof_indices;
		knotSpan(unsigned int _id,  nativeNURBSStructure<dim>& mesh, NURBSMesh<dim>& nurbsmesh, std::vector<unsigned int>& tempIndices, IGAValues<dim>& fe_values); 
		bool pointInKnotSpan(std::vector<double>& coords);
};
template <int dim>
knotSpan<dim>::knotSpan(unsigned int _id,  nativeNURBSStructure<dim>& mesh, NURBSMesh<dim>& nurbsmesh, std::vector<unsigned int>& tempIndices, IGAValues<dim>& fe_values): \
	id(_id), boundaryFlags(2*dim, 0), endKnots(dim) {
	//fill endKnots
	for (unsigned int a=0; a<dim; a++) {endKnots[a].push_back(mesh.knots[a].at(tempIndices[a])); endKnots[a].push_back(mesh.knots[a].at(tempIndices[a]+1));}
	//fill controlPoints
	for (typename std::vector<controlPoint<dim> >::iterator a=nurbsmesh.controlPointVector.begin(); a<nurbsmesh.controlPointVector.end(); a++){ 
		std::vector<double> tempVec(dim);
		for (unsigned int b=0; b<dim; b++) tempVec[b]=0.5*(endKnots[b][0]+endKnots[b][1]);
		if (a->basis.pointInBasisSupport(tempVec)) {
			controlPoints.push_back(&(*a));
		}
	}

	//fill local_dof_indices
	for (typename std::vector<controlPoint<dim>*>::iterator a=controlPoints.begin(); a<controlPoints.end(); a++){ 
		//fill local_dof_indices for this knotSpan
		for (unsigned int c=0; c<nurbsmesh.dofPerControlPoint; c++){
			unsigned int aIndex=((*a)->id)*nurbsmesh.dofPerControlPoint+c;
			local_dof_indices.push_back(aIndex);
		}
		local_mass_dof_indices.push_back((*a)->id);
	}
	//evaluate edgeCoords and center
	std::vector<double> quadPtStencil; 
	quadPtStencil.push_back(-1.0); quadPtStencil.push_back(1.0); quadPtStencil.push_back(1.0); quadPtStencil.push_back(1.0);
	//fill QuadPoints
	std::vector<std::vector<double> > quadPoints(std::pow(2, dim));
	std::vector<unsigned int> indices;
	for (unsigned int quadPointIndex=0; quadPointIndex<quadPoints.size(); quadPointIndex++){
		resolve<dim>(quadPointIndex, 2, indices); double temp=1.0;
		for (unsigned int j=0; j<dim; j++){
			quadPoints[quadPointIndex].push_back(quadPtStencil.at(2*indices[j]));
			temp*=quadPtStencil.at(2*indices[j]+1); 
		}
		quadPoints[quadPointIndex].push_back(temp);
	}
	fe_values.reinit(*this, &quadPoints);
	//fill edgeCoords and center
	edgeCoords.resize(quadPoints.size());
	center.resize(dim, 0.0);
	for (unsigned int quadPointIndex=0; quadPointIndex<quadPoints.size(); quadPointIndex++){
		for (unsigned int i=0; i<dim; i++) {
			edgeCoords[quadPointIndex].push_back(fe_values.quadPointLocations[quadPointIndex][i]);
			center[i]+=fe_values.quadPointLocations[quadPointIndex][i]/quadPoints.size();
		}
	}
}

template <int dim>
bool knotSpan<dim>::pointInKnotSpan(std::vector<double>& coords){
	bool temp=true;
	for (unsigned int i=0; i<dim; i++){
		if  (!((coords[i]>=endKnots[i][0] && (coords[i]<=endKnots[i][1])))) temp=false;
	}
	return temp;
}

//NURBS mesh class
template <int dim>
class NURBSMesh
{
	public:
		unsigned int dofPerControlPoint, quadPtStencilSize;
		NURBSMesh(nativeNURBSStructure<dim>& mesh, unsigned int _dofPerControlPoint, unsigned int _quadPtStencilSize);
		std::vector<controlPoint<dim> > controlPointVector;
		std::vector<knotSpan<dim> > knotSpanVector;
};
template <int dim>
NURBSMesh<dim>::NURBSMesh(nativeNURBSStructure<dim>& mesh, unsigned int _dofPerControlPoint, unsigned int _quadPtStencilSize): \
	dofPerControlPoint(_dofPerControlPoint), quadPtStencilSize(_quadPtStencilSize){	
	//fill controlPointVector
	unsigned int dofID=0;
	for (unsigned int controlPointID=0; controlPointID<mesh.ctrlPtCoefs[0].size(); controlPointID++){
		std::vector<unsigned int> dofs(dim);
		for (unsigned int a=0; a<dofPerControlPoint; a++) dofs[a]=dofID++;
		//create controlPoint and push back into controlPointVector
		controlPointVector.push_back(*(new controlPoint<dim>(controlPointID, mesh, dofs)));
	}
	//fill knotSpanVector
	IGAValues<dim> fe_values(this); //used to compute knotSpan edgeCoords and center 
	unsigned int numKnots=1, knotSpanID=0;
	for (unsigned int a=0; a<dim; a++) numKnots*=mesh.knots[a].size();
	for (unsigned int tempID=0; tempID<numKnots; tempID++){
		std::vector<unsigned int> tempIndices(dim);
		if (dim==3){
			std::div_t divresult= std::div ((int) tempID, (int) mesh.knots[0].size()*mesh.knots[1].size()); tempIndices[2]=divresult.quot;
			divresult= std::div ((int) divresult.rem, (int) mesh.knots[0].size()); tempIndices[1]=divresult.quot; tempIndices[0]=divresult.rem;
		}
		else if (dim==2){
			std::div_t divresult= std::div ((int) tempID, (int) mesh.knots[0].size()); 
			tempIndices[1]=divresult.quot; tempIndices[0]=divresult.rem;
		}
		else tempIndices[0]=tempID;
		
		//check if span width is greater than zero in each direction 
		bool skipLoopLast=false, skipLoopRep=false;
		for (unsigned int a=0; a<dim; a++) {
			if (tempIndices[a]==(mesh.knots[a].size()-1))  skipLoopLast=true;
			else if (mesh.knots[a].at(tempIndices[a])==mesh.knots[a].at(tempIndices[a]+1))  skipLoopRep=true;
		}
		if (skipLoopLast || skipLoopRep) continue;
		
		//create knotSpan and push back into knotSpanVector
		knotSpanVector.push_back(*(new knotSpan<dim>(knotSpanID++, mesh, *this, tempIndices, fe_values)));
	}
	printf("IGA data structures generated\n");
}

#endif /* NURBSCLASSES_H_ */
