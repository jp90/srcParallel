/*
 * derivatives.h
 *
 *  Created on: Nov 14, 2013
 *      Author: jan-philippwolf
 */

#ifndef DERIVATIVES_H_
#define DERIVATIVES_H_

#include "gridfunction.h"
#include "stencil.h"

	//! Computes the first derivative of u in x-direction.
	/*!
	 * \param output Stores the result
	 * \param u The gridfunction u
	 * \param h Grid distance in x- and y-direction
	 */
void Ux(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the first derivative of u in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Uy(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the second derivative of u in x-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Uxx(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the second derivative of u in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Uyy(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the first derivative of p in x-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Px(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the first derivative of p in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Py(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the second derivative of p in x-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Pxx(GridFunction& output, GridFunction& p, const PointType& h);

//! Computes the second derivative of p in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Pyy(GridFunction& output, GridFunction& p, const PointType& h);

//! Computes the second derivative of v in x-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Vxx(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the second derivative of v in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Vyy(GridFunction& output, GridFunction& u, const PointType& h);

//! Computes the first derivative of u*u in x-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void UUx(GridFunction& output, GridFunction& u, const RealType alpha,
		const PointType& h);

//! Computes the first derivative of v*v in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void VVy(GridFunction& output, GridFunction& u, const RealType alpha,
		const PointType& h);

//! Computes the derivative of u*v in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void UVy(GridFunction& output, GridFunction& u, GridFunction& v,
		const RealType alpha, const PointType& h);

//! Computes the derivative of u*v in x-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void UVx(GridFunction& output, GridFunction& u, GridFunction& v,
		const RealType alpha, const PointType& h);

//! Computes the derivative of F in x-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Fx(GridFunction& output, GridFunction& f, const PointType& h);

//! Computes the derivative of G in y-direction.
/*!
 * \param output Stores the result
 * \param u The gridfunction u
 * \param h Grid distance in x- and y-direction
 */
void Gy(GridFunction& output, GridFunction& g, const PointType& h);

void UTX(GridFunction& output, GridFunction& u, GridFunction& T,
		const RealType gamma, const PointType& h);

void VTY(GridFunction& output, GridFunction& v, GridFunction& T,
		const RealType gamma, const PointType& h);
void TXX(GridFunction& output, GridFunction& T, const PointType& h);
void TYY(GridFunction& output, GridFunction& T, const PointType& h);

#endif /* DERIVATIVES_H_ */
