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


	/*! \brief Exchanges p values
	 *
	 * Exchanges the p values between the domain's borders
	 * @param p pressure
	 */
	void ExchangePValues(GridFunction& p);
	void ExchangeTValues(GridFunction& t);
	/*! \brief Exchanges p values
	 *
	 * Exchanges the velocity values u and v between the domain's borders
	 * @param u velocity in x-direction
	 * @param v velocity in y-direction
	 */
	void ExchangeUVValues(GridFunction& u,GridFunction& v);
	Communication(int world_rank);

int world_rank;

};


#endif /* COMMUNICATION_H_ */
