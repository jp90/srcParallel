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
//#include <unistd.h>
#include "communication.h"
using namespace std;

int main() {

	int s = 0;

	char input[] = "./srcParallel/inputvals.txt";
	char output[] = "./srcParallel/inputvals.txt";
	IO SimIO(input, output);

	IndexType n = 0;
	RealType t = 0.0;

	MPI_Init(NULL,NULL);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	MultiIndexType global_grid;
		global_grid[0] = SimIO.para.iMax + 2;
		global_grid[1] = SimIO.para.jMax + 2;

	//PARA:
	SimIO.para.iMax = SimIO.para.iMax/2;
	SimIO.para.world_rank = world_rank;

	//initialize u,v,p
	MultiIndexType begin, end;
	//Initialize u-velocity
	begin[0] = 1;
	end[0] = SimIO.para.iMax - 1;
	begin[1] = 1;
	end[1] = SimIO.para.jMax;
	GridFunction u(SimIO.para.iMax+2, SimIO.para.jMax+2);
	u.SetGridFunction(begin, end, SimIO.para.ui);
    //Initialize v-velocity
	GridFunction v(SimIO.para.iMax+2, SimIO.para.jMax+2);
	begin[0] = 1;
	end[0] = SimIO.para.iMax;
	begin[1] = 1;
	end[1] = SimIO.para.jMax -1;
	v.SetGridFunction(begin, end, SimIO.para.vi);
	//Initialize pressure
	GridFunction p(SimIO.para.iMax+2, SimIO.para.jMax+2);
	begin[0] = 1;
	end[0] = SimIO.para.iMax;
	begin[1] = 1;
	end[1] = SimIO.para.jMax;
	p.SetGridFunction(begin, end, SimIO.para.pi);

	GridFunction gx(u.griddimension);
	GridFunction gy(v.griddimension);

	GridFunction f(u.griddimension);
	GridFunction g(v.griddimension);

	GridFunction rhs(p.griddimension);

	Computation computer(SimIO);
	Solver solve(SimIO);

	PointType delta;
    //computer.setBoundaryU(u);
    //p.Grid_Print();
	delta[0]=SimIO.para.deltaX;
	delta[1]=SimIO.para.deltaY;
	//Start Main Loop

	computer.setBoundaryV(v);
//	if(world_rank==0){sleep(s);}
					cout << world_rank <<": nach boundary" << endl;
						p.Grid_Print();
//	if(world_rank==1){sleep(s);}
					cout << world_rank <<": nach boundary" << endl;
						p.Grid_Print();

	Communication communication(world_rank);


int count=0;
	while ((t < SimIO.para.tEnd)){
		SimIO.writeVTKMasterfile(global_grid,u.getGridFunction(),v.getGridFunction(),p.getGridFunction(),delta,n);
		SimIO.writeVTKSlavefile(global_grid,u.getGridFunction(),v.getGridFunction(),p.getGridFunction(),delta,world_rank,n,n);

		// compute timestep size deltaT
		RealType uMax = u.MaxValueGridFunction(begin, end);
		RealType vMax = v.MaxValueGridFunction(begin, end);
		RealType deltaT = computer.computeTimesstep(uMax, vMax);

		double deltaTmin;

		//cout << "was soll des?" << endl;
		MPI_Reduce(&deltaT,&deltaTmin,1,MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Bcast(&deltaTmin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//	if(world_rank==1){sleep(s);}
	//		cout << world_rank <<": bevor boundary" << endl;
	//			p.Grid_Print();
		// set boundary values
		computer.setBoundaryU(u);
		computer.setBoundaryV(v);
	//	if(world_rank==1){sleep(s);}
	//				cout << world_rank <<": nach boundary" << endl;
	//					p.Grid_Print();

		// Compute F and G
		computer.computeMomentumEquations(f, g, u, v, gx, gy, deltaTmin);

		computer.setBoundaryF(f,u);
		computer.setBoundaryG(g,v);

	/*	if(world_rank==1)sleep(s);
		cout<<"u"<<endl;
		p.Grid_Print();

		cout<<"v"<<endl;
		p.Grid_Print(); */
		//if(world_rank==1){(s);}
				cout << world_rank <<": f" << endl;

		f.Grid_Print();
		//if(world_rank==1){sleep(s);}
		cout << world_rank <<": g" << endl;

		g.Grid_Print();


		// set right hand side of p equation
		computer.computeRighthandSide(rhs,f,g,deltaTmin);
	//	rhs.Grid_Print();
		//SOR-SOLVER
		int it =0;
		RealType Residuum = SimIO.para.eps+1.0;
		RealType Residuum_local = 0.0;
		while ((it < SimIO.para.iterMax) && (Residuum > SimIO.para.eps)) {
          //  cout << "Computing pressure" <<endl;
            //Set boundary
            computer.setBoundaryP(p);
            // SOR Cycle
            solve.SORCycle_Black(p,rhs);
           // p.Grid_Print();
          // 		if(world_rank==1){sleep(s);}
           // 		cout << world_rank <<": after black" << endl;
            //		p.Grid_Print();

            communication.ExchangePValues(p);


            solve.SORCycle_White(p,rhs);

            communication.ExchangePValues(p);
        //	if(world_rank==1){sleep(s);}
         //           		cout << world_rank <<": after white" << endl;
              //      		p.Grid_Print();
           // if(world_rank==0) sleep(s);
            //p.Grid_Print();
            //return 0;
			Residuum_local = solve.computeResidual(p,rhs);

			MPI_Reduce(&Residuum_local,&Residuum,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&Residuum,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			if(world_rank==0){
			cout << "Current Residuum: ";
			cout << Residuum<<endl;
			cout << "it= "<<it << "n="<<n<< endl;}
			it++;
		}
		if (Residuum>10.0)
			return 0;

		// Update velocites u and v
		computer.computeNewVelocities(u, v, f, g, p, deltaTmin);
	//	if(world_rank==1){sleep(s);}
	//	cout << world_rank <<": bevor austausch" << endl;
	//	p.Grid_Print();
		communication.ExchangeUVValues(u,v);

//		if(world_rank==1){sleep(s);}
//		cout << world_rank <<": nach austausch" << endl;
	//	p.Grid_Print();
		cout <<"Prozessor " << world_rank <<" finished" <<endl;

		n++;
	//	cout << "t= " << t<< endl;
		t += deltaTmin;
		count++;
	}
	cout << "laeuft!";
	MPI_Finalize();
}
