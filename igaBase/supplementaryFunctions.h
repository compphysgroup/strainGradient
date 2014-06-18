/*
 * supplementaryFunctions.h
 *
 *  Created on: May 11, 2011, Modified for IGA on Oct 25, 2011
 */

#ifndef SUPPLEMENTARYFUNCTIONS_H_
#define SUPPLEMENTARYFUNCTIONS_H_
#include "base/table.h"
#include <Sacado.hpp>

template <int dim>
void resolve(unsigned int i, unsigned int div, std::vector<unsigned int>& indices){
	indices.clear();
	if (dim==1){
		indices.push_back(i);
	}
	else if (dim==2){
		indices.push_back(i%div);
		indices.push_back(i/div);
	}
	else if (dim==3){
		indices.push_back(i%div);
		indices.push_back((i%(div*div))/div);
		indices.push_back(i/(div*div));
	}
}

unsigned int factorial(unsigned int number) {
	int temp;
	if(number <= 1) return 1;
	temp = number * factorial(number - 1);
	return temp;
}

double Bin(unsigned int i, unsigned int j) {
	if (j>i) {printf("j>i in binomial coefficient computation\n"); exit(-1);}
	return ((double)factorial(i))/(factorial(j)*factorial(i-j));
}

template <class T, int dim>
T determinantOfMinor(unsigned int theRowHeightY, unsigned int theColumnWidthX, dealii::Table<2, T>& matrix){
  unsigned int x1 = theColumnWidthX == 0 ? 1 : 0;  /* always either 0 or 1 */
  unsigned int x2 = theColumnWidthX == 2 ? 1 : 2;  /* always either 1 or 2 */
  unsigned int y1 = theRowHeightY   == 0 ? 1 : 0;  /* always either 0 or 1 */
  unsigned int y2 = theRowHeightY   == 2 ? 1 : 2;  /* always either 1 or 2 */
  return matrix[y1][x1]*matrix[y2][x2] - matrix[y1][x2]*matrix[y2][x1];
}


template <class T, int dim>
void getInverse(dealii::Table<2, T>& matrix, dealii::Table<2, T>& invMatrix, T& det){
	if (dim==1){
		det=matrix[0][0];
		invMatrix[0][0]=1.0/det;
	}
	else if(dim==2){
		det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
		invMatrix[0][0]=matrix[1][1]/det;
		invMatrix[1][0]=-matrix[1][0]/det;
		invMatrix[0][1]=-matrix[0][1]/det;
		invMatrix[1][1]=matrix[0][0]/det;
	}
	else if(dim==3){
		det=  matrix[0][0]*determinantOfMinor<T, dim>(0, 0, matrix) - matrix[0][1]*determinantOfMinor<T, dim>(0, 1, matrix) +  matrix[0][2]*determinantOfMinor<T, dim>(0, 2, matrix);
		for (int y=0;  y< dim;  y++){
			for (int x=0; x< dim;  x++){
				invMatrix[y][x] = determinantOfMinor<T, dim>(x, y, matrix)/det;
				if( ((x + y) % 2)==1){invMatrix[y][x]*=-1;}
			}
		}
	}
	else throw "dim>3";
	if (std::abs(det)< 1.0e-15){
		printf("**************Near zero determinant in Matrix inversion***********************\n"); throw "Near zero determinant in Matrix inversion";
	}
}



#endif /* SUPPLEMENTARYFUNCTIONS_H_ */

