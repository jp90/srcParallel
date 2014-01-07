/*
 * communciation.cpp
 *
 *  Created on: 2013-11-27
 *      Author: sina
 */

#include "communication.h"
//#include <unistd.h>
Communication::Communication(int world_rank) {
	this->world_rank = world_rank;
}
;

void Communication::ExchangePValues(GridFunction& p) {


	if (world_rank == 0) {
		GridFunction package(1, p.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = p.griddimension[1] - 2;
		MultiIndexType Offset;
		Offset[0] = p.griddimension[0] - 2;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, p, Offset);
	//	package.Grid_Print();

		int test = MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_P,MPI_COMM_WORLD);
		if(test!=MPI_SUCCESS){
		cout <<"0: MPI: Sending failed!"<<endl;}
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_P,MPI_COMM_WORLD,&status);
	//	package.Grid_Print();

				begin[0] = p.griddimension[0] - 1;
				end[0] = p.griddimension[0] - 1;
				begin[1] = 0;
				end[1] = p.griddimension[1] - 1;
				Offset[0] = -1*(p.griddimension[0] - 1);
				Offset[1] = 0;
				p.SetGridFunction(begin,end,1.0,package,Offset);
	}
	if (world_rank == 1) {
		GridFunction package(1, p.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = p.griddimension[1] - 1;
		MultiIndexType Offset;
		Offset[0] = 1;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, p, Offset);
//		package.Grid_Print();
		MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_P,MPI_COMM_WORLD);
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_P,MPI_COMM_WORLD,&status);
	//	package.Grid_Print();

		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = p.griddimension[1] - 1;
		Offset[0] = 0;
		Offset[1] = 0;
		p.SetGridFunction(begin,end,1.0,package);

		}

}

void Communication::ExchangeTValues(GridFunction& t) {


	if (world_rank == 0) {
		GridFunction package(1, t.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = t.griddimension[1] - 2;
		MultiIndexType Offset;
		Offset[0] = t.griddimension[0] - 2;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, t, Offset);
	//	package.Grid_Print();

		int test = MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_P,MPI_COMM_WORLD);
		if(test!=MPI_SUCCESS){
		cout <<"0: MPI: Sending failed!"<<endl;}
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_P,MPI_COMM_WORLD,&status);
	//	package.Grid_Print();

				begin[0] = t.griddimension[0] - 1;
				end[0] = t.griddimension[0] - 1;
				begin[1] = 0;
				end[1] = t.griddimension[1] - 1;
				Offset[0] = -1*(t.griddimension[0] - 1);
				Offset[1] = 0;
				t.SetGridFunction(begin,end,1.0,package,Offset);
	}
	if (world_rank == 1) {
		GridFunction package(1, t.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = t.griddimension[1] - 1;
		MultiIndexType Offset;
		Offset[0] = 1;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, t, Offset);
//		package.Grid_Print();
		MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_P,MPI_COMM_WORLD);
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_P,MPI_COMM_WORLD,&status);
	//	package.Grid_Print();

		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = t.griddimension[1] - 1;
		Offset[0] = 0;
		Offset[1] = 0;
		t.SetGridFunction(begin,end,1.0,package);

		}

}


void Communication::ExchangeUVValues(GridFunction& u,GridFunction& v) {

	//Exchange U

	if (world_rank == 0) {
		GridFunction package(1, u.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = u.griddimension[1] - 1;
		MultiIndexType Offset;
		Offset[0] = u.griddimension[0] - 3;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, u, Offset);
		//package.Grid_Print();

		int test = MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_U,MPI_COMM_WORLD);
		if(test!=MPI_SUCCESS){
		cout <<"0: MPI: Sending failed!"<<endl;}
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_U,MPI_COMM_WORLD,&status);

				begin[0] = u.griddimension[0] - 2;
				end[0] = u.griddimension[0] - 2;
				begin[1] = 0;
				end[1] = u.griddimension[1] - 1;
				Offset[0] = -1*(u.griddimension[0] - 2);
				Offset[1] = 0;
				u.SetGridFunction(begin,end,1.0,package,Offset);
			//	u.Grid_Print();
	}
	if (world_rank == 1) {
		GridFunction package(1, u.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = u.griddimension[1] - 1;
		MultiIndexType Offset;
		Offset[0] = 1;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, u, Offset);
//		package.Grid_Print();
		MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_U,MPI_COMM_WORLD);
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_U,MPI_COMM_WORLD,&status);
	//	package.Grid_Print();

		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = u.griddimension[1] - 1;
		Offset[0] = 0;
		Offset[1] = 0;
		u.SetGridFunction(begin,end,1.0,package);

		}



		// Exchange V

	if (world_rank == 0) {
		GridFunction package(1, v.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = v.griddimension[1] - 1;
		MultiIndexType Offset;
		Offset[0] = v.griddimension[0] - 2;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, v, Offset);
	//	package.Grid_Print();

		int test = MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_V,MPI_COMM_WORLD);
		if(test!=MPI_SUCCESS){
		cout <<"0: MPI: Sending failed!"<<endl;}
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,1,COM_V,MPI_COMM_WORLD,&status);
	//	package.Grid_Print();

				begin[0] = v.griddimension[0] - 1;
				end[0] = v.griddimension[0] - 1;
				begin[1] = 0;
				end[1] = v.griddimension[1] - 1;
				Offset[0] = -1*(v.griddimension[0] - 1);
				Offset[1] = 0;
				v.SetGridFunction(begin,end,1.0,package,Offset);
	}
	if (world_rank == 1) {
		GridFunction package(1, v.griddimension[1]);
		MultiIndexType begin, end;
		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = v.griddimension[1] - 1;
		MultiIndexType Offset;
		Offset[0] = 1;
		Offset[1] = 0;
		package.SetGridFunction(begin, end, 1.0, v, Offset);
//		package.Grid_Print();
		MPI_Send(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_V,MPI_COMM_WORLD);
		MPI_Status  status;
		MPI_Recv(&package.getGridFunction()[0][0],package.griddimension[1],MPI_DOUBLE,0,COM_V,MPI_COMM_WORLD,&status);
	//	package.Grid_Print();

		begin[0] = 0;
		end[0] = 0;
		begin[1] = 0;
		end[1] = v.griddimension[1] - 1;
		Offset[0] = 0;
		Offset[1] = 0;
		v.SetGridFunction(begin,end,1.0,package);

		}



}
