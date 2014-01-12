//============================================================================
// Name        : NumSimFluid.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <iostream>
#include "IO.hpp"
#include "gridfunction.h"
#include "stencil.h"
#include "derivatives.h"
#include "computation.h"
#include "solver.h"
#include <mpi.h>

#include <unistd.h>
#include "communication.h"
using namespace std;

RealType const_zero(RealType c) {
	return 0.0;
}
RealType const_one(RealType c) {
	return 1.0;
}
int main() {

	RealType (*TO)(RealType);
	RealType (*TU)(RealType);
	RealType (*TL)(RealType);
	RealType (*TR)(RealType);

	TO = &const_zero;
	TU = &const_zero;
	TL = &const_one;
	TR = &const_zero;

	int s = 1;

	char input[] = "./srcParallel/inputvals.txt";
	char output[] = "./srcParallel/inputvals.txt";
	IO SimIO(input, output);
	char geometry_path[] ="./srcParallel/wirbel_geometrie.csv";
	SimIO.readGeometry(geometry_path);


	IndexType n = 0;
	RealType time = 0.0;

	MPI_Init(NULL, NULL);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	MultiIndexType global_grid;
	global_grid[0] = SimIO.para.iMax + 2;
	global_grid[1] = SimIO.para.jMax + 2;

	//PARA:
	SimIO.para.iMax = SimIO.para.iMax / 2;
	SimIO.para.world_rank = world_rank;

	//initialize u,v,p
	MultiIndexType begin, end;
	//Initialize u-velocity
	begin[0] = 1;
	end[0] = SimIO.para.iMax - 1;
	begin[1] = 1;
	end[1] = SimIO.para.jMax;
	GridFunction u(SimIO.para.iMax + 2, SimIO.para.jMax + 2);
	u.SetGridFunction(begin, end, SimIO.para.ui);
	//Initialize v-velocity
	GridFunction v(SimIO.para.iMax + 2, SimIO.para.jMax + 2);
	begin[0] = 1;
	end[0] = SimIO.para.iMax;
	begin[1] = 1;
	end[1] = SimIO.para.jMax - 1;
	v.SetGridFunction(begin, end, SimIO.para.vi);
	//Initialize pressure
	GridFunction p(SimIO.para.iMax + 2, SimIO.para.jMax + 2);
	begin[0] = 1;
	end[0] = SimIO.para.iMax;
	begin[1] = 1;
	end[1] = SimIO.para.jMax;
	p.SetGridFunction(begin, end, SimIO.para.pi);
	//Initialize temperature
	GridFunction t(SimIO.para.iMax + 2, SimIO.para.jMax + 2);
	begin[0] = 1;
	end[0] = SimIO.para.iMax;
	begin[1] = 1;
	end[1] = SimIO.para.jMax;
	t.SetGridFunction(begin, end, SimIO.para.TI);
	//Initialize heat
	GridFunction h(SimIO.para.iMax + 2, SimIO.para.jMax + 2);

	GridFunction f(u.griddimension);
	GridFunction g(v.griddimension);

	GridFunction rhs(p.griddimension);

	Computation computer(SimIO);
	Solver solve(SimIO);

	PointType delta;

	//u.Grid_Print();
	//return 9;
	delta[0] = SimIO.para.deltaX;
	delta[1] = SimIO.para.deltaY;
	//Start Main Loop

	computer.setBoundaryTD(t, TO, TU, TL, TR);
	computer.setBoundaryTN(t, TO, TU, TL, TR);
	computer.setBoundaryU(u);
	computer.setBoundaryV(v);
//	if(world_rank==0){sleep(s);}
	cout << world_rank << ": nach boundary" << endl;
	p.Grid_Print();
//	if(world_rank==1){sleep(s);}
	cout << world_rank << ": nach boundary" << endl;
	p.Grid_Print();

	Communication communication(world_rank);

	int count = 0;

	while ((time < SimIO.para.tEnd)) {
		SimIO.writeVTKMasterfile(global_grid, u.getGridFunction(),
				v.getGridFunction(), p.getGridFunction(), delta, n);
		SimIO.writeVTKSlavefile(global_grid, u.getGridFunction(),
				v.getGridFunction(), p.getGridFunction(), t.getGridFunction(),
				h.getGridFunction(), delta, world_rank, n, n);

		// compute timestep size deltaT
		RealType uMax = u.MaxValueGridFunction(begin, end);
		RealType vMax = v.MaxValueGridFunction(begin, end);
		RealType deltaT = computer.computeTimesstep(uMax, vMax);

		double deltaTmin;

		MPI_Reduce(&deltaT, &deltaTmin, 1, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
		MPI_Bcast(&deltaTmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//	if(world_rank==1){sleep(s);}
		//		cout << world_rank <<": bevor boundary" << endl;
		//			p.Grid_Print();
		// set boundary values
		computer.setBoundaryU(u);
		computer.setBoundaryV(v);
		//	if(world_rank==1){sleep(s);}
		//				cout << world_rank <<": nach boundary" << endl;
		//					p.Grid_Print();

		computer.ComputeTemperature(t, u, v, deltaT);
	    computer.setBoundaryTD(t, TO, TU, TL, TR);
		computer.setBoundaryTN(t, TO, TU, TL, TR);
		communication.ExchangeTValues(t);
		computer.ComputeHeatfunction(h, t, u, deltaT);
		/*
		 if(world_rank==1){sleep(s);}
		 cout << world_rank << ": U" << endl;
		 u.Grid_Print();

		 if(world_rank==1){sleep(s);}
		 cout << world_rank << ": V" << endl;
		 v.Grid_Print();

		 if(world_rank==1){sleep(s);}
		 cout << world_rank << ": P" << endl;
		 p.Grid_Print();

		 if(world_rank==1){sleep(s);}
		 cout << world_rank << ": T" << endl;
		 t.Grid_Print();
		 */

		// Compute F and G
		computer.computeMomentumEquations(f, g, u, v, t, deltaTmin);

		computer.setBoundaryF(f, u);
		computer.setBoundaryG(g, v);

		/*	if(world_rank==1)sleep(s);
		 cout<<"u"<<endl;
		 p.Grid_Print();

		 cout<<"v"<<endl;
		 p.Grid_Print(); */
		//if(world_rank==1){(s);}
		cout << world_rank << ": f" << endl;

		//f.Grid_Print();
		//if(world_rank==1){sleep(s);}
		cout << world_rank << ": t" << endl;

		//t.Grid_Print();

		// set right hand side of p equation
		computer.computeRighthandSide(rhs, f, g, deltaTmin);
		//	rhs.Grid_Print();
		//SOR-SOLVER
		int it = 0;
		RealType Residuum = SimIO.para.eps + 1.0;
		RealType Residuum_local = 0.0;
		while ((it < SimIO.para.iterMax) && (Residuum > SimIO.para.eps)) {
			//  cout << "Computing pressure" <<endl;
			//Set boundary
			computer.setBoundaryP(p);
			// SOR Cycle
			solve.SORCycle_Black(p, rhs);
			// p.Grid_Print();
			// 		if(world_rank==1){sleep(s);}
			// 		cout << world_rank <<": after black" << endl;
			//		p.Grid_Print();

			communication.ExchangePValues(p);

			solve.SORCycle_White(p, rhs);

			communication.ExchangePValues(p);
			//	if(world_rank==1){sleep(s);}
			//           		cout << world_rank <<": after white" << endl;
			//      		p.Grid_Print();
			// if(world_rank==0) sleep(s);
			//p.Grid_Print();
			//return 0;
			Residuum_local = solve.computeResidual(p, rhs);

			MPI_Reduce(&Residuum_local, &Residuum, 1, MPI_DOUBLE, MPI_SUM, 0,
			MPI_COMM_WORLD);
			MPI_Bcast(&Residuum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if (world_rank == 0) {
				cout << "Current Residuum: ";
				cout << Residuum << endl;
				cout << "it= " << it << "; n=" << n << endl;
			}
			it++;
		}
		//if (Residuum > 10.0)
		//	return 0;

		// Update velocites u and v
		computer.computeNewVelocities(u, v, f, g, p, deltaTmin);
		//	if(world_rank==1){sleep(s);}
		//	cout << world_rank <<": bevor austausch" << endl;
		//	p.Grid_Print();
		communication.ExchangeUVValues(u, v);

//		if(world_rank==1){sleep(s);}
//		cout << world_rank <<": nach austausch" << endl;
		//	p.Grid_Print();
		cout << "Prozessor " << world_rank << " finished" << endl;

		n++;
		//	cout << "t= " << t<< endl;
		time += deltaTmin;
		count++;
	}
	cout << "laeuft!";
	MPI_Finalize();
}
