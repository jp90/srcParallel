/*
 * solver.cpp
 *
 *  Created on: 12.11.2013
 *      Author: David
 */

#include "solver.h"
#include <iostream>
#include "math.h"
using namespace std;

//Using the default copy constructor of the SimIO class
Solver::Solver(IO& SimIO) :
		SimIO(SimIO) {
}

//computes residuum over a subspace
RealType Solver::computeResidual(GridFunction& p, GridFunction& rhs) {
	PointType delta;
	delta[0] = SimIO.para.deltaX;
	delta[1] = SimIO.para.deltaY;

	MultiIndexType begin, end;
	begin[0] = 1;
	end[0] = p.griddimension[0] - 2;
	begin[1] = 1;
	end[1] = p.griddimension[1] - 2;

	//temporary gridfunctions
	GridFunction branch_1(p.griddimension);
	Pxx(branch_1, p, delta);
	GridFunction branch_2(p.griddimension);
	Pyy(branch_2, p, delta);

	branch_1.AddToGridFunction(begin, end, 1.0, branch_2);
	branch_1.AddToGridFunction(begin, end, -1.0, rhs);
	branch_1.MultiplyGridFunctions(begin, end, branch_1);
	branch_1.ScaleGridFunction(begin, end,
			1.0 / (SimIO.para.nfc));

	//equation for the residuum
	RealType res = 0.0;
	for (int i = 1; i <= SimIO.para.iMax; i++) {
		for (int j = 1; j <= SimIO.para.jMax; j++) {
			if (SimIO.geometry_field[i][j] > 15) {
				res += branch_1.getGridFunction()[i][j];
			}
		}
	}
	return sqrt(res);

}

//one cycle of the solver, starts with (1,2)
void Solver::SORCycle_Black(GridFunction& p, GridFunction& rhs) {

	for (int i = 1; i <= SimIO.para.iMax; i++) {
		for (int j = 2 - i % 2; j <= SimIO.para.jMax; j += 2) {
			if (SimIO.geometry_field[i][j] > 15) {
				p.getGridFunction()[i][j] =
						(1 - SimIO.para.omg) * p.getGridFunction()[i][j]
								+ (SimIO.para.omg)
										/ (2.0
												* (1.0
														/ (SimIO.para.deltaX
																* SimIO.para.deltaX)
														+ 1.0
																/ (SimIO.para.deltaY
																		* SimIO.para.deltaY)))
										* ((p.getGridFunction()[i + 1][j]
												+ p.getGridFunction()[i - 1][j])
												/ (SimIO.para.deltaX
														* SimIO.para.deltaX)
												+ (p.getGridFunction()[i][j + 1]
														+ p.getGridFunction()[i][j
																- 1])
														/ (SimIO.para.deltaY
																* SimIO.para.deltaY)
												- rhs.getGridFunction()[i][j]);
			}
		}
	}

}
//one cycle of the solver, starts with (1,1)
void Solver::SORCycle_White(GridFunction& p, GridFunction& rhs) {

	for (int i = 1; i <= SimIO.para.iMax; i++) {
		for (int j = 1 + i % 2; j <= SimIO.para.jMax; j += 2) {
			if (SimIO.geometry_field[i][j] > 15) {
				p.getGridFunction()[i][j] =
						(1 - SimIO.para.omg) * p.getGridFunction()[i][j]
								+ (SimIO.para.omg)
										/ (2.0
												* (1.0
														/ (SimIO.para.deltaX
																* SimIO.para.deltaX)
														+ 1.0
																/ (SimIO.para.deltaY
																		* SimIO.para.deltaY)))
										* ((p.getGridFunction()[i + 1][j]
												+ p.getGridFunction()[i - 1][j])
												/ (SimIO.para.deltaX
														* SimIO.para.deltaX)
												+ (p.getGridFunction()[i][j + 1]
														+ p.getGridFunction()[i][j
																- 1])
														/ (SimIO.para.deltaY
																* SimIO.para.deltaY)
												- rhs.getGridFunction()[i][j]);
			}
		}
	}
}
