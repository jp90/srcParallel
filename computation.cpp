/*
 * computation.cpp
 *
 *  Created on: 12.11.2013
 *      Author: David
 */

#include "computation.h"
#include <iostream>

using namespace std;

Computation::Computation(IO& SimIO) :
		SimIO(SimIO) {
}

RealType Computation::computeTimesstep(RealType uMax, RealType vMax) {
	RealType deltaT, min, c, d;
	c = SimIO.para.re
			/ (2.0
					* (1.0 / (SimIO.para.deltaX * SimIO.para.deltaX)
							+ 1.0 / (SimIO.para.deltaY * SimIO.para.deltaY)));
	d = SimIO.para.re * SimIO.para.Pr
			/ (2.0
					* (1.0 / (SimIO.para.deltaX * SimIO.para.deltaX)
							+ 1.0 / (SimIO.para.deltaY * SimIO.para.deltaY)));

	if ((c < SimIO.para.deltaX / uMax) && (c < SimIO.para.deltaY / vMax)
			&& (c < d))
		min = c;

	else if ((d < SimIO.para.deltaX / uMax) && (d < SimIO.para.deltaY / vMax)
			&& (d < c))
		min = d;

	else if ((SimIO.para.deltaX / uMax < c)
			&& (SimIO.para.deltaX / uMax < SimIO.para.deltaY / vMax)
			&& (SimIO.para.deltaX / uMax < d))
		min = SimIO.para.deltaX / uMax;

	else
		min = SimIO.para.deltaY / vMax;

	deltaT = SimIO.para.tau * min;
	return deltaT;
}

void Computation::computeNewVelocities(GridFunction& u, GridFunction& v,
		GridFunction& f, GridFunction& g, GridFunction& p, RealType deltaT) {
	GridFunction branch_1(p.griddimension);
	MultiIndexType begin, end;
	begin[0] = 1;
	end[0] = u.griddimension[0] - 3;
	begin[1] = 1;
	end[1] = u.griddimension[1] - 2;
	// u Update
	PointType delta;
	delta[0] = SimIO.para.deltaX;
	delta[1] = SimIO.para.deltaY;
	Px(branch_1, p, delta);
	f.AddToGridFunction(begin, end, -1.0 * deltaT, branch_1);
	u.SetGridFunctionFluidU(begin, end, 1.0, f,SimIO.geometry_field);
	// v Update
	begin[0] = 1;
	end[0] = v.griddimension[0] - 2;
	begin[1] = 1;
	end[1] = v.griddimension[1] - 3;
	GridFunction branch_2(p.griddimension);
	Py(branch_2, p, delta);
	g.AddToGridFunction(begin, end, -1.0 * deltaT, branch_2);
	v.SetGridFunctionFluidV(begin, end, 1.0, g, SimIO.geometry_field);

}
;

void Computation::computeMomentumEquations(GridFunction& f, GridFunction& g,
		GridFunction& u, GridFunction& v, GridFunction& t, RealType& deltaT) {
	PointType h;
	h[0] = SimIO.para.deltaX;
	h[1] = SimIO.para.deltaY;
	RealType alpha = SimIO.para.alpha;

	MultiIndexType begin, end;
	begin[0] = 1;
	end[0] = u.griddimension[0] - 3;
	begin[1] = 1;
	end[1] = u.griddimension[1] - 2;

	// Term F

	Uxx(f, u, h);

	GridFunction branch_2(u.griddimension);
	Uyy(branch_2, u, h);
	f.AddToGridFunction(begin, end, 1.0, branch_2);
	f.ScaleGridFunction(begin, end, 1.0 / SimIO.para.re);
	GridFunction branch_3(u.griddimension);
	UUx(branch_3, u, alpha, h);
	GridFunction branch_4(u.griddimension);
	UVy(branch_4, u, v, alpha, h);
	branch_3.AddToGridFunction(begin, end, 1.0, branch_4);

	f.AddToGridFunction(begin, end, -1.0, branch_3);
	// KILL branch 2-4

	f.ScaleGridFunction(begin, end, deltaT);
	f.AddToGridFunction(begin, end, 1.0, u);

	// - beta deltaT (tx stencil) gx
	GridFunction branch_5(u.griddimension);

	Stencil stencil_1(3, h);
	stencil_1.setTxStencil();
	stencil_1.ApplyStencilOperator(begin, end, begin, end, t, branch_5);

	branch_5.ScaleGridFunction(begin, end,
			deltaT * SimIO.para.beta * SimIO.para.gx);

	f.AddToGridFunction(begin, end, -1.0, branch_5);

	//Term G
	begin[0] = 1;
	end[0] = v.griddimension[0] - 2;
	begin[1] = 1;
	end[1] = v.griddimension[1] - 3;

	Vxx(g, v, h);
	GridFunction branch_6(v.griddimension);
	Vyy(branch_6, v, h);
	g.AddToGridFunction(begin, end, 1.0, branch_6);
	g.ScaleGridFunction(begin, end, 1.0 / SimIO.para.re);

	GridFunction branch_7(v.griddimension);
	UVx(branch_7, u, v, alpha, h);
	GridFunction branch_8(v.griddimension);
	VVy(branch_8, v, alpha, h);
	branch_7.AddToGridFunction(begin, end, 1.0, branch_8);

	g.AddToGridFunction(begin, end, -1.0, branch_7);

	g.ScaleGridFunction(begin, end, deltaT);
	g.AddToGridFunction(begin, end, 1.0, v);

	// - beta deltaT (ty stencil) gy
	GridFunction branch_9(u.griddimension);

	Stencil stencil_2(3, h);
	stencil_2.setTyStencil();
	stencil_2.ApplyStencilOperator(begin, end, begin, end, t, branch_9);

	branch_9.ScaleGridFunction(begin, end,
			deltaT * SimIO.para.beta * SimIO.para.gy);

	g.AddToGridFunction(begin, end, -1.0, branch_9);
}

void Computation::setBoundaryU(GridFunction& u) {

	for (int i = 1; i < u.griddimension[0] - 2; i++) {
		for (int j = 1; j < u.griddimension[1] - 2; j++) {
			switch (SimIO.geometry_field[i][j]) {
			//hindernis, Ÿberall hindernis
			case 0:
				u.getGridFunction()[i][j] = 0.0;
				break;
				//hindernis, N Fluid
			case 1:
				u.getGridFunction()[i - 1][j] = -u.getGridFunction()[i - 1][j
						+ 1];
				u.getGridFunction()[i][j] = -u.getGridFunction()[i][j + 1];
				break;
				//hindernis, S Fluid
			case 2:
				u.getGridFunction()[i - 1][j] = -u.getGridFunction()[i - 1][j
						- 1];
				u.getGridFunction()[i][j] = -u.getGridFunction()[i][j - 1];
				break;
			case 3:
				cout << "Forbidden geometry, sucker!";
				break;
				//hindernis, W Fluid
			case 4:
				u.getGridFunction()[i - 1][j] = 0.0;
				break;
				//hindernis, NW Fluid
			case 5:
				u.getGridFunction()[i][j] = -u.getGridFunction()[i][j + 1];
				u.getGridFunction()[i - 1][j] = 0.0;
				break;
				//hindernis, SW Fluid
			case 6:
				u.getGridFunction()[i][j] = -u.getGridFunction()[i][j - 1];
				u.getGridFunction()[i - 1][j] = 0.0;
				break;
			case 7:
				cout << "Forbidden geometry, sucker!";
				break;
				//hindernis, O Fluid
			case 8:
				u.getGridFunction()[i][j] = 0.0;
				break;
				//hindernis, NO Fluid
			case 9:
				u.getGridFunction()[i][j] = 0.0;
				u.getGridFunction()[i - 1][j] = -u.getGridFunction()[i - 1][j
						+ 1];
				break;
				//hindernis, SO Fluid
			case 10:
				u.getGridFunction()[i][j] = 0.0;
				u.getGridFunction()[i - 1][j] = -u.getGridFunction()[i - 1][j
						- 1];
				break;

			case 11:
				cout << "Forbidden geometry, sucker!";
				break;
			case 12:
				cout << "Forbidden geometry, sucker!";
				break;
			case 13:
				cout << "Forbidden geometry, sucker!";
				break;
			case 14:
				cout << "Forbidden geometry, sucker!";
				break;
			case 15:
				cout << "Forbidden geometry, sucker!";
				break;

			}
		}
	}
	MultiIndexType begin, end;
	if (SimIO.para.world_rank == 0) {
		// u_0,j = 0
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 1;
		end[1] = u.griddimension[1] - 2;

		if (SimIO.para.WL == 1) {
			u.SetGridFunction(begin, end, 0.0);
		} else if (SimIO.para.WL == 2) {
			u.SetGridFunction(begin, end, 0.0);
		} else if (SimIO.para.WL == 3) {
			MultiIndexType Offset;
			Offset[0] = 1;
			Offset[1] = 0;
			u.SetGridFunction(begin, end, 1.0, Offset);
		} else if (SimIO.para.WL == 4) {
			MultiIndexType Offset;
			Offset[0] = 1;
			Offset[1] = 0;
			u.SetGridFunction(begin, end, -1.0, u, Offset, 2.0);
			//u.SetGridFunction(begin,end,1.0);
		}

	}
	if (SimIO.para.world_rank == 1) {
		//u_iMax,j = 0
		begin[0] = u.griddimension[0] - 2;
		end[0] = u.griddimension[0] - 2;
		begin[1] = 1;
		end[1] = u.griddimension[1] - 2;
		if (SimIO.para.WR == 1) {
			u.SetGridFunction(begin, end, 0.0);
		} else if (SimIO.para.WR == 2) {
			u.SetGridFunction(begin, end, 0.0);
		} else if (SimIO.para.WR == 3) {
			MultiIndexType Offset;
			Offset[0] = -1;
			Offset[1] = 0;
			u.SetGridFunction(begin, end, 1.0, Offset);
		}
	}
	// u_i,0
	begin[0] = 1;
	end[0] = u.griddimension[0] - 2;
	begin[1] = 0;
	end[1] = 0;
	if (SimIO.para.WU == 1) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = 1;
		u.SetGridFunction(begin, end, -1.0, Offset);
	} else if (SimIO.para.WU == 2) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = 1;
		u.SetGridFunction(begin, end, 1.0, Offset);
	} else if (SimIO.para.WU == 3) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = 1;
		u.SetGridFunction(begin, end, 1.0, Offset);
	}
	// u_i,jMax+1
	begin[0] = 1;
	end[0] = u.griddimension[0] - 2;
	begin[1] = u.griddimension[1] - 1;
	end[1] = u.griddimension[1] - 1;

	if (SimIO.para.WO == 1) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = -1;
		u.SetGridFunction(begin, end, -1.0, Offset);
		//u.SetGridFunction(begin, end, -1.0, u, Offset, 2.0);
	} else if (SimIO.para.WO == 2) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = -1;
		u.SetGridFunction(begin, end, 1.0, Offset);
	} else if (SimIO.para.WO == 3) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = -1;
		u.SetGridFunction(begin, end, 1.0, Offset);
	}
}
void Computation::setBoundaryV(GridFunction& v) {
	MultiIndexType begin, end;
	for (int i = 1; i < v.griddimension[0] - 2; i++) {
		for (int j = 1; j < v.griddimension[1] - 2; j++) {
			switch (SimIO.geometry_field[i][j]) {
			//hindernis, Ÿberall hindernis
			case 0:
				v.getGridFunction()[i][j] = 0.0;
				break;
				//hindernis, N Fluid
			case 1:
				v.getGridFunction()[i][j] = 0.0;
				break;
				//hindernis, S Fluid
			case 2:
				v.getGridFunction()[i][j - 1] = 0.0;
				break;
			case 3:
				cout << "Forbidden geometry, sucker!";
				break;
				//hindernis, W Fluid
			case 4:
				v.getGridFunction()[i][j - 1] = -v.getGridFunction()[i - 1][j
						- 1];
				v.getGridFunction()[i][j] = -v.getGridFunction()[i - 1][j];
				break;
				//hindernis, NW Fluid
			case 5:
				v.getGridFunction()[i][j] = 0.0;
				v.getGridFunction()[i][j - 1] = -v.getGridFunction()[i - 1][j
						- 1];

				break;
				//hindernis, SW Fluid
			case 6:
				v.getGridFunction()[i][j] = -v.getGridFunction()[i - 1][j];
				v.getGridFunction()[i][j - 1] = 0.0;

				break;
			case 7:
				cout << "Forbidden geometry, sucker!";
				break;
				//hindernis, O Fluid
			case 8:
				v.getGridFunction()[i][j - 1] = -v.getGridFunction()[i + 1][j
						- 1];
				v.getGridFunction()[i][j] = -v.getGridFunction()[i + 1][j];

				break;
				//hindernis, NO Fluid
			case 9:
				v.getGridFunction()[i][j] = 0;
				v.getGridFunction()[i][j - 1] = -v.getGridFunction()[i + 1][j
						- 1];

				break;
				//hindernis, SO Fluid
			case 10:
				v.getGridFunction()[i][j] = -v.getGridFunction()[i + 1][j];
				v.getGridFunction()[i][j - 1] = 0.0;

				break;

			case 11:
				cout << "Forbidden geometry, sucker!";
				break;
			case 12:
				cout << "Forbidden geometry, sucker!";
				break;
			case 13:
				cout << "Forbidden geometry, sucker!";
				break;
			case 14:
				cout << "Forbidden geometry, sucker!";
				break;
			case 15:
				cout << "Forbidden geometry, sucker!";
				break;

			}
		}
	}
	// v_i,0 = 0
	begin[0] = 1;
	end[0] = v.griddimension[0] - 2;
	begin[1] = 0;
	end[1] = 0;
	if (SimIO.para.WU == 1) {
		v.SetGridFunction(begin, end, 0.0);
	} else if (SimIO.para.WU == 2) {
		v.SetGridFunction(begin, end, 0.0);
	} else if (SimIO.para.WU == 3) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = 1;
		v.SetGridFunction(begin, end, 1.0, Offset);
	}

	// v_i,jMax =0
	begin[0] = 1;
	end[0] = v.griddimension[0] - 2;
	begin[1] = v.griddimension[1] - 2;
	end[1] = v.griddimension[1] - 2;
	if (SimIO.para.WO == 1) {
		v.SetGridFunction(begin, end, 0.0);
	} else if (SimIO.para.WO == 2) {
		v.SetGridFunction(begin, end, 0.0);
	} else if (SimIO.para.WO == 3) {
		MultiIndexType Offset;
		Offset[0] = 0;
		Offset[1] = -1;
		v.SetGridFunction(begin, end, 1.0, Offset);
	}

	if (SimIO.para.world_rank == 0) {
		// v_0,j = -v_1,j
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 1;
		end[1] = v.griddimension[1] - 2;
		if (SimIO.para.WL == 1) {
			MultiIndexType Offset;
			Offset[0] = 1;
			Offset[1] = 0;
			v.SetGridFunction(begin, end, -1.0, Offset);
		} else if (SimIO.para.WL == 2) {
			MultiIndexType Offset;
			Offset[0] = 1;
			Offset[1] = 0;
			v.SetGridFunction(begin, end, 1.0, Offset);
		} else if (SimIO.para.WL == 3) {
			MultiIndexType Offset;
			Offset[0] = 1;
			Offset[1] = 0;
			v.SetGridFunction(begin, end, 1.0, Offset);
		} else if (SimIO.para.WL == 4) {
			MultiIndexType Offset;
			Offset[0] = 1;
			Offset[1] = 0;
			v.SetGridFunction(begin, end, 0.0);
		}
	}
	if (SimIO.para.world_rank == 1) {
		// v_iMax+1,j = -v_iMax,j
		begin[0] = v.griddimension[0] - 1;
		end[0] = v.griddimension[0] - 1;
		begin[1] = 1;
		end[1] = v.griddimension[1] - 2;
		if (SimIO.para.WR == 1) {
			MultiIndexType Offset;
			Offset[0] = -1;
			Offset[1] = 0;
			v.SetGridFunction(begin, end, -1.0, Offset);
		} else if (SimIO.para.WR == 2) {
			MultiIndexType Offset;
			Offset[0] = -1;
			Offset[1] = 0;
			v.SetGridFunction(begin, end, 1.0, Offset);
		} else if (SimIO.para.WR == 3) {
			MultiIndexType Offset;
			Offset[0] = -1;
			Offset[1] = 0;
			v.SetGridFunction(begin, end, 1.0, Offset);
		}
	}
}
void Computation::setBoundaryP(GridFunction& p) {
	MultiIndexType begin, end;
	for (int i = 1; i < p.griddimension[0] - 2; i++) {
		for (int j = 1; j < p.griddimension[1] - 2; j++) {
			switch (SimIO.geometry_field[i][j]) {
			//hindernis, Ÿberall hindernis
			case 0:
				p.getGridFunction()[i][j] = 0.0;
				break;
				//hindernis, N Fluid
			case 1:
				p.getGridFunction()[i][j] = p.getGridFunction()[i][j + 1];
				break;
				//hindernis, S Fluid
			case 2:
				p.getGridFunction()[i][j] = p.getGridFunction()[i][j - 1];
				break;
			case 3:
				cout << "Forbidden geometry, sucker!";
				break;
				//hindernis, W Fluid
			case 4:
				p.getGridFunction()[i][j] = p.getGridFunction()[i - 1][j];
				break;
				//hindernis, NW Fluid
			case 5:
				p.getGridFunction()[i][j] = 1.0
						/ (SimIO.para.deltaX * SimIO.para.deltaX
								+ SimIO.para.deltaY * SimIO.para.deltaY)
						* (SimIO.para.deltaX * SimIO.para.deltaX
								* p.getGridFunction()[i][j + 1]
								+ SimIO.para.deltaY * SimIO.para.deltaY
										* p.getGridFunction()[i + 1][j]);
				break;
				//hindernis, SW Fluid
			case 6:
				p.getGridFunction()[i][j] = 1.0
						/ (SimIO.para.deltaX * SimIO.para.deltaX
								+ SimIO.para.deltaY * SimIO.para.deltaY)
						* (SimIO.para.deltaX * SimIO.para.deltaX
								* p.getGridFunction()[i][j - 1]
								+ SimIO.para.deltaY * SimIO.para.deltaY
										* p.getGridFunction()[i + 1][j]);
				break;
			case 7:
				cout << "Forbidden geometry, sucker!";
				break;
				//hindernis, O Fluid
			case 8:
				p.getGridFunction()[i][j] = p.getGridFunction()[i + 1][j];
				break;
				//hindernis, NO Fluid
			case 9:
				p.getGridFunction()[i][j] = 1.0
						/ (SimIO.para.deltaX * SimIO.para.deltaX
								+ SimIO.para.deltaY * SimIO.para.deltaY)
						* (SimIO.para.deltaX * SimIO.para.deltaX
								* p.getGridFunction()[i][j + 1]
								+ SimIO.para.deltaY * SimIO.para.deltaY
										* p.getGridFunction()[i - 1][j]);
				break;
				//hindernis, SO Fluid
			case 10:
				p.getGridFunction()[i][j] = 1.0
						/ (SimIO.para.deltaX * SimIO.para.deltaX
								+ SimIO.para.deltaY * SimIO.para.deltaY)
						* (SimIO.para.deltaX * SimIO.para.deltaX
								* p.getGridFunction()[i][j - 1]
								+ SimIO.para.deltaY * SimIO.para.deltaY
										* p.getGridFunction()[i - 1][j]);
				break;

			case 11:
				cout << "Forbidden geometry, sucker!";
				break;
			case 12:
				cout << "Forbidden geometry, sucker!";
				break;
			case 13:
				cout << "Forbidden geometry, sucker!";
				break;
			case 14:
				cout << "Forbidden geometry, sucker!";
				break;
			case 15:
				cout << "Forbidden geometry, sucker!";
				break;

			}
		}
	}
	if (SimIO.para.world_rank == 0) {
		// p_0,j = p_1,j
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 1;
		end[1] = p.griddimension[1] - 2;
		MultiIndexType Offset;
		Offset[0] = 1;
		Offset[1] = 0;
		p.SetGridFunction(begin, end, 1.0, Offset);
	}

	if (SimIO.para.world_rank == 1) {
		// p_iMax+1,j = p_iMax,j
		begin[0] = p.griddimension[0] - 1;
		end[0] = p.griddimension[0] - 1;
		begin[1] = 1;
		end[1] = p.griddimension[1] - 2;
		MultiIndexType Offset;
		Offset[0] = -1;
		Offset[1] = 0;
		p.SetGridFunction(begin, end, 1.0, Offset);
	}
	// p_i,0 = p_i,1
	begin[0] = 1;
	end[0] = p.griddimension[0] - 2;
	begin[1] = 0;
	end[1] = 0;
	MultiIndexType Offset;
	Offset[0] = 0;
	Offset[1] = 1;
	p.SetGridFunction(begin, end, 1.0, Offset);

	// p_i,jMax+1 = p_i,jMax
	begin[0] = 1;
	end[0] = p.griddimension[0] - 2;
	begin[1] = p.griddimension[1] - 1;
	end[1] = p.griddimension[1] - 1;
	Offset[0] = 0;
	Offset[1] = -1;
	p.SetGridFunction(begin, end, 1.0, Offset);

}
void Computation::setBoundaryF(GridFunction& f, GridFunction& u) {
	for (int i = 1; i < f.griddimension[0] - 2; i++) {
		for (int j = 1; j < f.griddimension[1] - 2; j++) {
			switch (SimIO.geometry_field[i][j]) {
			//hindernis, Ÿberall hindernis
			case 0:
				f.getGridFunction()[i][j] = 0.0;
				break;
				//hindernis, N Fluid

			case 4:
				f.getGridFunction()[i - 1][j] = u.getGridFunction()[i - 1][j];
				break;
				//hindernis, NW Fluid
			case 5:
				f.getGridFunction()[i - 1][j] = u.getGridFunction()[i - 1][j];

				break;
				//hindernis, SW Fluid
			case 6:
				f.getGridFunction()[i - 1][j] = u.getGridFunction()[i - 1][j];

				break;
				//hindernis, O Fluid
			case 8:
				f.getGridFunction()[i][j] = u.getGridFunction()[i][j];

				break;
				//hindernis, NO Fluid
			case 9:

				f.getGridFunction()[i][j] = u.getGridFunction()[i][j];

				break;
				//hindernis, SO Fluid
			case 10:
				f.getGridFunction()[i][j] = u.getGridFunction()[i][j];
				break;

			case 11:
				cout << "Forbidden geometry, sucker!";
				break;
			case 12:
				cout << "Forbidden geometry, sucker!";
				break;
			case 13:
				cout << "Forbidden geometry, sucker!";
				break;
			case 14:
				cout << "Forbidden geometry, sucker!";
				break;
			case 15:
				cout << "Forbidden geometry, sucker!";
				break;

			}
		}
	}
	MultiIndexType begin, end;
	begin[0] = 0;
	end[0] = 0;
	begin[1] = 1;
	end[1] = f.griddimension[1] - 2;
	f.SetGridFunction(begin, end, 1.0, u);

// F_iMax,j=u_iMax+1,j
	begin[0] = f.griddimension[0] - 2;
	end[0] = f.griddimension[0] - 2;
	begin[1] = 1;
	end[1] = f.griddimension[1] - 2;
	f.SetGridFunction(begin, end, 1.0, u);
}
void Computation::setBoundaryG(GridFunction& g, GridFunction& v) {
	for (int i = 1; i < g.griddimension[0] - 2; i++) {
		for (int j = 1; j < g.griddimension[1] - 2; j++) {
			switch (SimIO.geometry_field[i][j]) {
			//hindernis, Ÿberall hindernis
			case 0:
				g.getGridFunction()[i][j] = 0.0;
				break;
				//hindernis, N Fluid

			case 1:
				g.getGridFunction()[i][j] = v.getGridFunction()[i][j];
				break;
				//hindernis, NW Fluid
			case 2:
				g.getGridFunction()[i][j - 1] = v.getGridFunction()[i][j - 1];

				break;
				//hindernis, SW Fluid
			case 5:
				g.getGridFunction()[i][j] = v.getGridFunction()[i][j];

				break;
				//hindernis, O Fluid
			case 6:
				g.getGridFunction()[i][j - 1] = v.getGridFunction()[i][j - 1];

				break;
				//hindernis, NO Fluid
			case 9:

				g.getGridFunction()[i][j] = v.getGridFunction()[i][j];

				break;
				//hindernis, SO Fluid
			case 10:
				g.getGridFunction()[i][j - 1] = v.getGridFunction()[i][j - 1];
				break;

			case 11:
				cout << "Forbidden geometry, sucker!";
				break;
			case 12:
				cout << "Forbidden geometry, sucker!";
				break;
			case 13:
				cout << "Forbidden geometry, sucker!";
				break;
			case 14:
				cout << "Forbidden geometry, sucker!";
				break;
			case 15:
				cout << "Forbidden geometry, sucker!";
				break;

			}
		}
	}
	MultiIndexType begin, end;
	begin[0] = 1;
	end[0] = g.griddimension[0] - 2;
	begin[1] = 0;
	end[1] = 0;
	g.SetGridFunction(begin, end, 1.0, v);

	begin[0] = 1;
	end[0] = g.griddimension[0] - 2;
	begin[1] = g.griddimension[1] - 2;
	end[1] = g.griddimension[1] - 2;
	g.SetGridFunction(begin, end, 1.0, v);
}

void Computation::computeRighthandSide(GridFunction& rhs, GridFunction& f,
		GridFunction& g, RealType deltaT) {

	GridFunction branch_1(g.griddimension);

	MultiIndexType begin, end;

	begin[0] = 1;
	end[0] = f.griddimension[0] - 2;
	begin[1] = 1;
	end[1] = f.griddimension[1] - 2;
	PointType delta;
	delta[0] = SimIO.para.deltaX;
	delta[1] = SimIO.para.deltaY;
	Fx(rhs, f, delta);
	Gy(branch_1, g, delta);

	rhs.AddToGridFunction(begin, end, 1.0, branch_1);
	rhs.ScaleGridFunction(begin, end, 1.0 / deltaT);

}

void Computation::ComputeTemperature(GridFunction& T, GridFunction& u,
		GridFunction& v, RealType deltaT) {
	GridFunction branch_1(T.griddimension);
	GridFunction branch_2(T.griddimension);
	MultiIndexType begin, end;
//missing wï¿½rmequelle
	begin[0] = 1;
	end[0] = T.griddimension[0] - 2;
	begin[1] = 1;
	end[1] = T.griddimension[1] - 2;
	PointType delta;
	delta[0] = SimIO.para.deltaX;
	delta[1] = SimIO.para.deltaY;

	TXX(branch_1, T, delta);
	TYY(branch_2, T, delta);
	branch_1.AddToGridFunction(begin, end, 1.0, branch_2);
	branch_1.ScaleGridFunction(begin, end,
			1.0 / (SimIO.para.Pr * SimIO.para.re));

	GridFunction branch_3(T.griddimension);
	GridFunction branch_4(T.griddimension);
	UTX(branch_3, u, T, SimIO.para.gamma, delta);
	VTY(branch_4, v, T, SimIO.para.gamma, delta);
	branch_3.AddToGridFunction(begin, end, 1.0, branch_4);

	branch_1.AddToGridFunction(begin, end, -1.0, branch_3);
	branch_1.ScaleGridFunction(begin, end, deltaT);
	T.AddToGridFunction(begin, end, 1.0, branch_1);
}

void Computation::ComputeHeatfunction(GridFunction& h, GridFunction& t,
		GridFunction& u, RealType deltaT) {
	for (int i = 0; i <= h.griddimension[0] - 2; i++) {
		for (int j = 1; j <= h.griddimension[1] - 2; j++) {
			h.getGridFunction()[i][j] = h.getGridFunction()[i][j - 1]
					+ deltaT
							* (SimIO.para.re * SimIO.para.Pr
									* u.getGridFunction()[i][j]
									* (t.getGridFunction()[i + 1][j]
											+ t.getGridFunction()[i][j]) / 2.0
									- (t.getGridFunction()[i + 1][j]
											- t.getGridFunction()[i][j])
											/ SimIO.para.deltaX);
		}
	}

}

void Computation::setBoundaryTD(GridFunction& T, RealType (*TO)(RealType),
		RealType (*TU)(RealType), RealType (*TL)(RealType),
		RealType (*TR)(RealType)) {
	MultiIndexType begin, end;
	if (SimIO.para.world_rank == 0) {
		// p_0,j = p_1,j
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 1;
		end[1] = T.griddimension[1] - 2;

		T.SetGridFunction(begin, end, TL, true, SimIO.para.deltaY);
		T.ScaleGridFunction(begin, end, 2.0);
		MultiIndexType Offset;
		Offset[0] = 1;
		Offset[1] = 0;
		T.AddToGridFunction(begin, end, -1.0, T, Offset);

	}

	if (SimIO.para.world_rank == 1) {
		// T_imax+1,j
		begin[0] = T.griddimension[0] - 1;
		end[0] = T.griddimension[0] - 1;
		begin[1] = 1;
		end[1] = T.griddimension[1] - 2;
		MultiIndexType Offset;
		Offset[0] = -1;
		Offset[1] = 0;
		T.SetGridFunction(begin, end, TR, true, SimIO.para.deltaY);
		T.ScaleGridFunction(begin, end, 2.0);
		T.AddToGridFunction(begin, end, -1.0, T, Offset);
	}
	/*
	 // T_i,0

	 begin[0] = 1;
	 end[0] = T.griddimension[0] - 2;
	 begin[1] = 0;
	 end[1] = 0;
	 MultiIndexType Offset;
	 Offset[0] = 0;
	 Offset[1] = 1;
	 T.SetGridFunction(begin, end, TU, false, SimIO.para.deltaX);
	 T.ScaleGridFunction(begin, end, 2.0);
	 T.AddToGridFunction(begin, end, -1.0, T, Offset);

	 // T_i,jmax+1

	 begin[0] = 1;
	 end[0] = T.griddimension[0] - 2;
	 begin[1] = T.griddimension[1] - 1;
	 end[1] = T.griddimension[1] - 1;
	 MultiIndexType Offset;
	 Offset[0] = 0;
	 Offset[1] = -1;
	 T.SetGridFunction(begin, end, TO, false, SimIO.para.deltaX);
	 T.ScaleGridFunction(begin, end, 2.0);
	 T.AddToGridFunction(begin, end, -1.0, T, Offset);*/

}
void Computation::setBoundaryTN(GridFunction& T, RealType (*TO)(RealType),
		RealType (*TU)(RealType), RealType (*TL)(RealType),
		RealType (*TR)(RealType)) {
	MultiIndexType begin, end;
	// T_i,0

	begin[0] = 1;
	end[0] = T.griddimension[0] - 2;
	begin[1] = 0;
	end[1] = 0;
	MultiIndexType Offset;
	Offset[0] = 0;
	Offset[1] = 1;
	T.SetGridFunction(begin, end, TU, false, SimIO.para.deltaX);
	T.ScaleGridFunction(begin, end, SimIO.para.deltaY);
	T.AddToGridFunction(begin, end, 1.0, T, Offset);

	// T_i,jmax+1

	begin[0] = 1;
	end[0] = T.griddimension[0] - 2;
	begin[1] = T.griddimension[1] - 1;
	end[1] = T.griddimension[1] - 1;
	//MultiIndexType Offset;
	Offset[0] = 0;
	Offset[1] = -1;
	T.SetGridFunction(begin, end, TO, false, SimIO.para.deltaX);
	T.ScaleGridFunction(begin, end, SimIO.para.deltaY);
	T.AddToGridFunction(begin, end, 1.0, T, Offset);
}
