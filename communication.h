/*
 * communication.h
 *
 *  Created on: 2013-11-27
 *      Author: sina
 */

#include "gridfunction.h"
#include "mpi.h"


#ifndef COMMUNICATION_H_
#define COMMUNICATION_H_

const int COM_P = 1;
const int COM_U = 2;
const int COM_V = 4;
using namespace std;

class Communication {
public:



	void ExchangePValues(GridFunction& p);

	void ExchangeUVValues(GridFunction& u,GridFunction& v);
	Communication(int world_rank);
//	MPI_Comm comm_cart;
//	int dims[size];
//	int periods[size];
//	int reorder;
int world_rank;

};


#endif /* COMMUNICATION_H_ */
