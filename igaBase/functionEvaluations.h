/*
 * functionEvaluations.h
 *
 *  Created on: May 1, 2011, Modified for IGA on Oct 25, 2011
 */

#ifndef FUNCTIONEVALUATIONS_H_
#define FUNCTIONEVALUATIONS_H_
#include "nurbsClasses.h"
#include "base/table.h"

template <class T, int dim>
  struct deformationMap{
  deformationMap(unsigned int n_q_points): F(n_q_points, dim, dim),  invF(n_q_points, dim, dim), detF(n_q_points){}
    dealii::Table<3, T> F, invF;
    dealii::Table<1, T> detF;
  };

template <class T, int dim>
  struct deformationMapwithGrad{
  deformationMapwithGrad(unsigned int n_q_points): F(n_q_points, dim, dim),  invF(n_q_points, dim, dim), gradF(n_q_points, dim, dim, dim), detF(n_q_points){}
    dealii::Table<3, T> F, invF;
    dealii::Table<4, T> gradF;
    dealii::Table<1, T> detF;
  };


template <class T, int dim>
  void evaluateScalarFunction(IGAValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T>& ULocal, dealii::Table<1, T>& U,  int faceID =-1){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    U[q]=0.0; //U
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      if (fe_values.system_to_component_index(k)==DOF){
	if (faceID>-1){
	  U[q]+=ULocal[k]*fe_values.shape_value(k, q, faceID); //U
	}
	else{
	  U[q]+=ULocal[k]*fe_values.shape_value(k, q); 
	}
      }
    }
  }
}

template <class T, int dim>
  void evaluateScalarFunctionGradient(IGAValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T>& ULocal, dealii::Table<2, T>& gradU, deformationMap<T, dim>& defMap, bool gradientInCurrentConfiguration, int faceID =-1){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  dealii::Table<1, T> refGradU(dim);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    refGradU.fill(0.0);
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      if (fe_values.system_to_component_index(k)==DOF){
	for (unsigned int i=0; i<dim; ++i){
	  if (faceID>-1){
	    refGradU[i]+=ULocal[k]*fe_values.shape_grad(k, q, faceID)[i]; //gradU
	  }
	  else{
	    refGradU[i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	  }
	}
      }
    }
    //Transform gradient to current configuration. gradW=(F^-T)*GradW
    for (unsigned int i=0; i<dim; ++i){
      if (gradientInCurrentConfiguration==false) gradU[q][i]=refGradU[i];
      else{
	gradU[q][i]=0.0;
	for (unsigned int j=0; j<dim; ++j){
	  gradU[q][i]+=defMap.invF[q][j][i]*refGradU[j];
	}
      }
    }
  }
}

template <class T, int dim>
  void evaluateVectorFunction(IGAValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T>& ULocal, dealii::Table<2, T>& U,  int faceID =-1){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  U.fill(0.0);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.system_to_component_index(k) - DOF;
      if (ck>=0 && ck<dim){
	if (faceID>-1){
	  U[q][ck]+=ULocal[k]*fe_values.shape_value(k, q, faceID); //U
	}
	else{
	  U[q][ck]+=ULocal[k]*fe_values.shape_value(k, q); //U
	}
      }
    }
  }
}

template <class T, int dim>
  void evaluateVectorFunctionGradient(IGAValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T>& ULocal, dealii::Table<3, T>& GradU, int faceID =-1){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  GradU.fill(0.0);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.system_to_component_index(k) - DOF;
      if (ck>=0 && ck<dim){
	for (unsigned int i=0; i<dim; ++i){
	  if (faceID>-1){
	    GradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad(k, q, faceID)[i]; //GradU
	  }
	  else{
	    GradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //GradU
	  }
	}
      }
    }
  }
}

template <class T, int dim>
  void evaluateVectorFunctionSecondGradient(IGAValues<dim>& fe_values, unsigned int DOF, dealii::Table<1, T>& ULocal, dealii::Table<4, T>& GradGradU, int faceID =-1){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  GradGradU.fill(0.0);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.system_to_component_index(k) - DOF;
      if (ck>=0 && ck<dim){
	for (unsigned int i=0; i<dim; ++i){
	  for (unsigned int j=0; j<dim; ++j){
	    if (faceID>-1){
	      GradGradU[q][ck][i][j]+=ULocal[k]*fe_values.shape_grad_grad(k, q, faceID)[i][j]; //GradGradU
	    }
	    else{
	      GradGradU[q][ck][i][j]+=ULocal[k]*fe_values.shape_grad_grad(k, q)[i][j]; //GradGradU
	    }
	  }
	}
      }
    }
  }
}

#endif /* FUNCTIONEVALUATIONS_H_ */

