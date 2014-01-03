//! The class implements the IO
/*!
 * @author sina,jp,dave
 * @date 2013
 */

#include "typedef.h"
#include "simparam.h"

#ifndef GRIDFUNCTION_HPP_
#define GRIDFUNCTION_HPP_

#include <iostream>
#include <fstream>

class GridFunction {
public:
	/*! \brief Allocates a grid
	 *
	 *  Allocates a grid [0, dimX - 1][0, dimY -1] with the data type 'RealType'.
	 *  @param dimX grid dimension in x-direction
	 *  @param dimY grid dimension in y-direction
	 */
	GridFunction(int dimX, int dimY);

	/*! \brief Allocates a grid
	 *
	 *  Allocates a grid [0, griddimension[0] - 1][0, griddimension[1] -1] with the data type 'RealType'.
	 *  @param griddimension vector holding the grid dimensions in x- and y-direction
	 */
	GridFunction(const MultiIndexType griddimension);

	/*! \brief Deallocates a grid
	 *
	 *   Deallocates memory which was allocated before with the method GridFunction.
	 */
	~GridFunction();

	/*! \brief Points on the gridfunction
	 *
	 *   Returns a pointer on the gridfunction.
	 */
	GridFunctionType getGridFunction();

	/*! \brief Sets entries to a certain value
	 *
	 *  Sets all entries in a rectangular domain to a given value.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param value value which is set in the given domain
	 */
	void SetGridFunction(const MultiIndexType& begin, const MultiIndexType& end,
			RealType value);

	/*! \brief Sets entries to a certain value
	 *
	 *  Takes the entries in a certain domain of the Gridfunction and multiplies them with a given factor.
	 *  The entries within another given domain are set to these values.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param factor factor with which the entries are multiplied
	 *  @param offset determines the domain which is multiplied with the given factor
	 */
	void SetGridFunction(const MultiIndexType& begin, const MultiIndexType& end,
			RealType factor, MultiIndexType& offset);

	/*! \brief Scales entries by a certain value
	 *
	 *  Scales all entries in a rectangular domain by a certain value.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param factor value by which the entries in the given domain are multiplied
	 */
	void ScaleGridFunction(const MultiIndexType& begin,
			const MultiIndexType& end, RealType factor);

	/*! \brief Sets entries to a certain value
	 *
	 *  Takes entries of a certain domain of a given Gridfunction, multiplies them  with a given value
	 *  and writes the results into the original Gridfunction.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param factor value by which the entries in the given domain are multiplied
	 *  @param sourcegridfunction source of the values multiplied by a given factor
	 */
	void SetGridFunction(const MultiIndexType& begin, const MultiIndexType& end,
			RealType factor, GridFunction& sourcegridFunction);

	/*! \brief Sets entries to a certain value
	 *
	 *  Takes entries of a certain domain of a given Gridfunction, multiplies them  with a given value
	 *  and writes the results into another given domain in the original Gridfunction.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param factor value by which the entries in the given domain are multiplied
	 *  @param sourcegridfunction source of the values multiplied by a given factor
	 *  @param offset determines the domain which is multiplied with the given factor
	 */
	void SetGridFunction(const MultiIndexType& begin, const MultiIndexType& end,
			RealType factor, GridFunction& sourcegridFunction,
			MultiIndexType& offset);

	/*! \brief Sets entries to a certain value
	 *
	 *  Takes entries of a certain domain of a given Gridfunction, multiplies them  with a given value,
	 *  adds a constant value and writes the results into another given domain in the original Gridfunction.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param factor value by which the entries in the given domain are multiplied
	 *  @param sourcegridfunction source of the values multiplied by a given factor
	 *  @param offset determines the domain which is multiplied with the given factor
	 *  @param constant value added to te result of the multiplication
	 */
	void SetGridFunction(const MultiIndexType& begin, const MultiIndexType& end,
			RealType factor, GridFunction& sourcegridFunction,
			MultiIndexType& offset, RealType constant);

	/*! \brief Sets entries to a certain value
	 *
	 *  Takes entries of a certain domain of a given Gridfunction, multiplies them  with a given value
	 *  and writes the results into another given domain in the original Gridfunction.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param factor value by which the entries in the given domain are multiplied
	 */
	void AddToGridFunction(const MultiIndexType& begin,
			const MultiIndexType& end, RealType factor,
			GridFunction& sourcegridFunction);

	/*! \brief Returns the maximum value
	 *
	 *  Returns the maximum value within the given gridfunction.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @return Maximum value
	 */
	RealType MaxValueGridFunction(const MultiIndexType& begin,
			const MultiIndexType& end);

	/*! \brief Sets entries to a certain value
	 *
	 *  Multiplies the values of the given gridfunction by the values in the sourcegridfunction.
	 *  @param begin grid coordinates of the top left corner of the domain where the grid function is set
	 *  @param end grid coordinates of the bottom right corner of the domain where the grid function is set
	 *  @param sourcegridfunction source of the values multiplied by a given factor
	 */
	void MultiplyGridFunctions(const MultiIndexType& begin,
			const MultiIndexType& end, GridFunction& sourcegridFunction);

	/*! \brief Prints a gridfunction
	 *
	 *  Prints the given gridfunction as a matrix.
	 */
	void Grid_Print();

	GridFunctionType gridfunction;
	MultiIndexType griddimension;

};

#endif
