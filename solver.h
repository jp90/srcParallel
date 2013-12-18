/*
 * solver.h
 *
 *  Created on: 12.11.2013
 *      Author: David
 */
#include "gridfunction.h"
#include "IO.hpp"
#include "derivatives.h"

#ifndef SOLVER_H_
#define SOLVER_H_
/*! \brief Solver class containing the Solver of the pressure
 *
 *	Solver class containing the Solverfunction of the pressure with the chess board algorithm.
 *  The related functions are SORCycle_Black and SORCycle_White. Also the function computeResidual
 *  is given and computs the resuduum over a subspace
 */
class Solver {
public:
	/*!
	 * Constructor of the Solver
	 * @param SimIO structur containing inputvalues
	 * */
	Solver(IO& SimIO);
	/*!
	 * Constructor of the Solver
	 * @param Gridfunction of the pressure (call by reference)
	 * @param Gridfunction of the right hand side (call by reference)
	 * @return returns the Residuum (RealType)
	 * */
	RealType computeResidual(GridFunction& p,
			GridFunction& rhs);
	/*!
	 * SOR Cycle over first half half of the subspace (begins with the field (1,1)
	 * @param Gridfunction of the pressure (call by reference)
	 * @param Gridfunction of the right hand side (call by reference)
	 * */
	void SORCycle_Black(GridFunction& p, GridFunction& rhs);
	/*!
	 * SOR Cycle over the other half of the subspace (begins with the field (1,2)
	 * @param Gridfunction of the pressure (call by reference)
	 * @param Gridfunction of the right hand side (call by reference)
	 * */
	void SORCycle_White(GridFunction& p, GridFunction& rhs);

	IO SimIO;
};

#endif /* SOLVER_H_ */
