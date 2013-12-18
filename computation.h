/*
 * computation.h
 *
 *  Created on: 12.11.2013
 *      Author: David
 */

#include "gridfunction.h"
#include "derivatives.h"
#include "IO.hpp"

#ifndef COMPUTATION_H_
#define COMPUTATION_H_

class Computation {
public:

	Computation(IO& SimIO);
	/*! \brief Computation of time step size
	 *
	 * Computes the time step size deltaT using the maximum velocities uMax and vMax.
	 * @param uMax maximum velocity stored in u
	 * @param vMax maximum velocity stored in v
	 * @return deltaT time step size
	 */
	RealType computeTimesstep(RealType uMax, RealType vMax);
	/*! \brief Computes the new velocites u and v
	 *
	 * Computes the new velocities u and v using the pressure p, time step size deltaT and f,g.
	 * @param u velocity in x-direction
	 * @param v velocity in y-direction
	 * @param f
	 * @param g
	 * @param p pressure
	 * @param deltaT time step size
	 */
	void computeNewVelocities(GridFunction& u, GridFunction& v, GridFunction& f,
			GridFunction& g, GridFunction& p, RealType deltaT);
	/*! \brief Computes the momentums f and g
	 *
	 * Computes the momentums f and g.
	 * @param u velocity in x-direction
	 * @param v velocity in y-direction
	 * @param f
	 * @param g
	 * @param gx external force in x-direction
	 * @param gy external force in y-direction
	 * @param deltaT
	 */
	void computeMomentumEquations(GridFunction& f, GridFunction& g,
			GridFunction& u, GridFunction& v, GridFunction& gx,
			GridFunction& gy, RealType& deltaT);
	/*! \brief Sets the boundary values for u.
	 *
	 * Sets the boundary values for u.
	 * @param u velocity in x-direction
	 */
	void setBoundaryU(GridFunction& u);
	/*! \brief Sets the boundary values for v.
	 *
	 * Sets the boundary values for v.
	 * @param v velocity in y-direction
	 */

	void setBoundaryV(GridFunction& v);

	/*! \brief Sets the boundary values for p.
	 *
	 * Sets the boundary values for p.
	 * @param p pressure
	 */
void setBoundaryP(GridFunction& p);
/*! \brief Sets the boundary values for f.
 *
 * Sets the boundary values for f.
 * @param f
 */

void setBoundaryF(GridFunction& f, GridFunction& u);
/*! \brief Sets the boundary values for g.
 *
 * Sets the boundary values for g.
 * @param g
 */

void setBoundaryG(GridFunction& g, GridFunction& v);
/*! \brief Computes the right hands side (rhs) of the pressure equation.
 *
 * Computes the right hands side (rhs) of the pressure equation.
 * @param rhs value of right hand side
 * @param f
 * @param g
 * @param deltaT time step size
 */

	void computeRighthandSide(GridFunction& rhs, GridFunction& f, GridFunction& g,RealType deltaT);

	IO SimIO;

};

#endif /* COMPUTATION_H_ */
