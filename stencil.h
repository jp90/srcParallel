/*
 * stencil.h
 *
 *  Created on: Nov 12, 2013
 *      Author: jan-philippwolf
 */

#ifndef STENCIL_H_
#define STENCIL_H_

#include "gridfunction.h"
#include "typedef.h"
/*! \brief Stencil class containing the stencils used for the computations of the derivatives
 *
 *	Stencil class containing the stencils used for the computations of the derivatives
 *  For every partial differention exists one stencile the capitalized letter indicats the source
 *  of the finit differential quotient, whereas the Number of x and y indicats the level and direction
 *  of differencation
 *
  */
class Stencil {
public:
	/*!
	 * Constructor
	 * @param respresents the size of the stencile array (mostly 3)
	 * @param stands for the grid distance
	 */
	Stencil(int stencilwidth, const PointType& h);

	~Stencil();
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
	void ApplyStencilOperator(const MultiIndexType& gridreadbegin,
			const MultiIndexType& gridreadend,
			const MultiIndexType& gridwritebegin,
			const MultiIndexType& gridwriteend,
			GridFunction& sourcegridFunction, GridFunction& imagegridFunction);
	/*! Stencil for derivation of the U-values in x direction
	 */
	void setUxStencil();
	/*! Stencil for derivation of the U-values in y direction
	 */
	void setUyStencil();
	/*! Stencil for derivation of P-values in x direction
	 */
	void setPxStencil();
	/*! Stencil for derivation of P-values in y direction
	 */
	void setPyStencil();
	/*! Stencil for the second derivation of P-values in x direction
	 */
	void setPxxStencil();
	/*! Stencil for the second derivation of P-values in y direction
	 */
	void setPyyStencil();
	/*! Stencil for the second derivation of U-values in x direction
	 */
	void setUxxStencil();
	/*! Stencil for the second derivation of U-values in y direction
	 */
	void setUyyStencil();
	/*! intermediate stencils for the derivation of squared U-values in x direction
	 */
	void setUUx_1Stencil();
	void setUUx_2Stencil();
	void setUUx_3Stencil();
	void setUUx_4Stencil();
	void setUUx_5Stencil();
	void setUUx_6Stencil();
	/*! intermediate stencils for the derivation of squared V-values in y direction
	 */
	void setVVy_1Stencil();
	void setVVy_2Stencil();
	void setVVy_3Stencil();
	void setVVy_4Stencil();
	void setVVy_5Stencil();
	void setVVy_6Stencil();
	/*! intermediate stencils for the derivation of "U multiplicated with V"-values in x direction
	 */
	void setUVx_1Stencil();
	void setUVx_2Stencil();
	void setUVx_3Stencil();
	void setUVx_4Stencil();
	void setUVx_5Stencil();
	void setUVx_6Stencil();
	void setUVx_7Stencil();
	void setUVx_8Stencil();
	/*! intermediate stencils for the derivation of "U multiplicated with V"-values in y direction
	 */
	void setUVy_1Stencil();
	void setUVy_2Stencil();
	void setUVy_3Stencil();
	void setUVy_4Stencil();
	void setUVy_5Stencil();
	void setUVy_6Stencil();
	void setUVy_7Stencil();
	void setUVy_8Stencil();
	/*! Stencil for derivation of the F-values in x direction
	 */
	void setFxStencil();
	/*! Stencil for derivation of the G-values in y direction
	 */
	void setGyStencil();

	StencilType stencil;
	int stencilwidth;
	const PointType& h;
	/*! true: the absolute value of the computed stencil-value
	 * false: the signed values is used
		 */
	bool abs;

};
#endif /* STENCIL_H_ */
