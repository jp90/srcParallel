/*
 * stencil.cpp
 *
 *  Created on: Nov 12, 2013
 *      Author: jan-philippwolf
 */

#include "gridfunction.h"
#include "stencil.h"
#include <iostream>
#include <fstream>


using namespace std;

/*!
 * Constructor
 * @param respresents the size of the stencile array (mostly 3)
 * @param stands for the grid distance
 */
Stencil::Stencil(int stencilwidth, const PointType& h) : h(h) {
	this->stencilwidth = stencilwidth;

	abs=false;

	stencil = new RealType* [stencilwidth];
	 for (int i=0;i<stencilwidth;i++){
		 stencil[i] = new RealType[stencilwidth];
	 }
	// Initialize grid values to zero
	 for (int i=0;i<stencilwidth;i++){
		 for (int j=0;j<stencilwidth;j++){
			 stencil[i][j]=0;
		 }
	 }
}

Stencil::~Stencil(){
	 for (int i=0;i<stencilwidth;i++){
		 delete[] stencil[i];
	 }
    delete[] stencil;
}
/*!
 * ApplyStencilOperator
 * applies the stencil to the whole grid
 * @param indices where to start reading
 * @param indices where to end reading
 * @param indices where to start writing
 * @param indices where to end writing
 * @param Source  of the values, which are combined with the values of the stencil
 * @param destination of the computet values
 */
void Stencil::ApplyStencilOperator(const MultiIndexType& gridreadbegin, const MultiIndexType& gridreadend,
		                           const MultiIndexType& gridwritebegin, const MultiIndexType& gridwriteend,
		                           GridFunction& sourcegridfunction, GridFunction& imagegridfunction){
	int a =int((stencilwidth-1)/2);
int count=0;
	for(int i=gridwritebegin[0];i<=gridwriteend[0];i++){
		for(int j=gridwritebegin[1];j<=gridwriteend[1];j++){
			RealType sum = 0.0;
			for(int k=0;k<stencilwidth;k++){
				for(int l=0;l<stencilwidth;l++){
							sum += sourcegridfunction.getGridFunction()[i+k-a][j+l-a]*stencil[k][l];
				}
			}
          if(abs) {
        	  if(sum<0.0)sum=-1.0*(sum);
          }
		imagegridfunction.getGridFunction()[i][j] = sum;
	}
	}
}
/*! Stencil for derivation of the U-values in x direction
	 */
void Stencil::setUxStencil(){
	abs=false;
	stencil[0][1]=-1.0/h[0];
	stencil[1][1]=1.0/h[0];
}
/*! Stencil for derivation of the U-values in y direction
 */
void Stencil::setUyStencil(){
	abs=false;
	stencil[1][1]=1.0/h[1];
	stencil[1][0]=-1.0/h[1];
}
/*! Stencil for derivation of P-values in x direction
 */
void Stencil::setPxStencil(){
	abs=false;
	stencil[2][1]=1.0/h[0];
	stencil[1][1]=-1.0/h[0];
}
/*! Stencil for derivation of P-values in y direction
 */
void Stencil::setPyStencil(){
	abs=false;
	stencil[1][2]=1.0/h[1];
	stencil[1][1]=-1.0/h[1];
}
/*! Stencil for the second derivation of P-values in x direction
 */
void Stencil::setPxxStencil(){
	abs=false;
	stencil[2][1]=1.0/(h[0]*h[0]);
	stencil[1][1]=-2.0/(h[0]*h[0]);
	stencil[0][1]=1.0/(h[0]*h[0]);
}
/*! Stencil for the second derivation of P-values in y direction
 */
void Stencil::setPyyStencil(){
	abs=false;
	stencil[1][2]=1.0/(h[1]*h[1]);
	stencil[1][1]=-2.0/(h[1]*h[1]);
	stencil[1][0]=1.0/(h[1]*h[1]);
}
/*! Stencil for the second derivation of U-values in x direction
 */
void Stencil::setUxxStencil(){
	abs=false;
	stencil[0][1]=1.0/(h[0]*h[0]);
	stencil[1][1]=-2.0/(h[0]*h[0]);
	stencil[2][1]=1.0/(h[0]*h[0]);
}
/*! Stencil for the second derivation of U-values in y direction
 */
void Stencil::setUyyStencil(){
	abs=false;
	stencil[1][0]=1.0/(h[1]*h[1]);
	stencil[1][1]=-2.0/(h[1]*h[1]);
	stencil[1][2]=1.0/(h[1]*h[1]);
}
/*! intermediate stencils for the derivation of squared U-values in x direction
 */
void Stencil::setUUx_1Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[2][1]=1.0/2.0;
}
void Stencil::setUUx_2Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[0][1]=1.0/2.0;
}
void Stencil::setUUx_3Stencil(){
	abs = true;
	stencil[1][1]=1.0/2.0;
	stencil[2][1]=1.0/2.0;
}
void Stencil::setUUx_4Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[2][1]=-1.0/2.0;
}
void Stencil::setUUx_5Stencil(){
	abs = true;
	stencil[1][1]=1.0/2.0;
	stencil[0][1]=1.0/2.0;
}
void Stencil::setUUx_6Stencil(){
	abs=false;
	stencil[1][1]=-1.0/2.0;
	stencil[0][1]=1.0/2.0;
}
/*! intermediate stencils for the derivation of squared V-values in y direction
 */
void Stencil::setVVy_1Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[1][2]=1.0/2.0;
}
void Stencil::setVVy_2Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[1][0]=1.0/2.0;
}
void Stencil::setVVy_3Stencil(){
	abs = true;
	stencil[1][1]=1.0/2.0;
	stencil[1][2]=1.0/2.0;
}
void Stencil::setVVy_4Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[1][2]=-1.0/2.0;
}
void Stencil::setVVy_5Stencil(){
	abs = true;
	stencil[1][1]=1.0/2.0;
	stencil[1][0]=1.0/2.0;
}
void Stencil::setVVy_6Stencil(){
	abs=false;
	stencil[1][1]=-1.0/2.0;
	stencil[1][0]=1.0/2.0;
}
/*! intermediate stencils for the derivation of "U multiplicated with V"-values in x direction
 */
void Stencil::setUVx_1Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[1][2]=1.0/2.0;
}
void Stencil::setUVx_2Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[2][1]=1.0/2.0;
}
void Stencil::setUVx_3Stencil(){
	abs=false;
	stencil[0][1]=1.0/2.0;
	stencil[0][2]=1.0/2.0;
}
void Stencil::setUVx_4Stencil(){
	abs=false;
	stencil[0][1]=1.0/2.0;
	stencil[1][1]=1.0/2.0;
}
void Stencil::setUVx_5Stencil(){
	abs=true;
	stencil[1][1]=1.0/2.0;
	stencil[1][2]=1.0/2.0;
}
void Stencil::setUVx_6Stencil(){
	abs=false;
	stencil[1][1]=1.0/2.0;
	stencil[2][1]=-1.0/2.0;
}
void Stencil::setUVx_7Stencil(){
	abs=true;
	stencil[0][1]=1.0/2.0;
	stencil[0][2]=1.0/2.0;
}
void Stencil::setUVx_8Stencil(){
	abs=false;
	stencil[0][1]=1.0/2.0;
	stencil[1][1]=-1.0/2.0;
}
/*! intermediate stencils for the derivation of "U multiplicated with V"-values in y direction
 */
void Stencil::setUVy_1Stencil(){
	abs=false;
	stencil[2][1]=1.0/2.0;
	stencil[1][1]=1.0/2.0;
}
void Stencil::setUVy_2Stencil(){
	abs=false;
	stencil[1][2]=1.0/2.0;
	stencil[1][1]=1.0/2.0;
}
void Stencil::setUVy_3Stencil(){
	abs=false;
	stencil[1][0]=1.0/2.0;
	stencil[2][0]=1.0/2.0;
}
void Stencil::setUVy_4Stencil(){
	abs=false;
	stencil[1][0]=1.0/2.0;
	stencil[1][1]=1.0/2.0;
}
void Stencil::setUVy_5Stencil(){
	abs=true;
	stencil[2][1]=1.0/2.0;
	stencil[1][1]=1.0/2.0;
}
void Stencil::setUVy_6Stencil(){
	abs=false;
	stencil[1][2]=-1.0/2.0;
	stencil[1][1]=1.0/2.0;
}
void Stencil::setUVy_7Stencil(){
	abs=true;
	stencil[1][0]=1.0/2.0;
	stencil[2][0]=1.0/2.0;
}
void Stencil::setUVy_8Stencil(){
	abs=false;
	stencil[1][0]=1.0/2.0;
	stencil[1][1]=-1.0/2.0;
}
/*! Stencil for derivation of the F-values in x direction
	 */
void Stencil::setFxStencil(){
	abs=false;
	stencil[1][1]=1.0/h[0];
	stencil[0][1]=-1.0/h[0];
}
/*! Stencil for derivation of the G-values in y direction
	 */
void Stencil::setGyStencil(){
	abs=false;
	stencil[1][1]=1.0/h[1];
	stencil[1][0]=-1.0/h[1];
}
