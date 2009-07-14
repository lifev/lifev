/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-07-09

  Copyright (C) 2009 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file BCInterfaceFunctionFile.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-09
 */

#ifndef __BCInterfaceFunctionFile_H
#define __BCInterfaceFunctionFile_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

#include <sstream>

#include <lifemc/lifefem/BCInterfaceFunction.hpp>





// ===================================================
//! Namespaces & Enums
// ===================================================
namespace LifeV {





/*!
 * \class BCInterfaceFunctionFile
 * \brief LifeV bcFunction wrapper for BCInterface
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  This class is an interface between BCInterface and SpiritParser. It allows to construct LifeV
 *  functions type for boundary conditions, using a GetPot file containing a function string and a
 *  table of discrete data (for example a Flux or a Pressure depending on time).
 *
 *
 *
 *  <b>DETAILS:</b>
 *
 *  The constructor of the class takes a string contains the GetPot file name.
 *  The GetPot file has the following structure:
 *
 *  function:  contains the expression of the function to use (as described in the BCInterfaceFunction class).
 *
 *  variables: contains the list of variables and coefficients present in the function.
 *             The first one is the variable and should be sorted in a growing order,
 *             while all the others are coefficients.
 *
 *  data:      contains the discrete list of variable and related coefficients. It must
 *             have n columns, where n is the number of variables (and coefficients).
 *
 *	scale:     Contains n optional coefficients (Default = 1), which multiply each column
 *	           in the data table.
 *
 *	NOTE:
 *	During the execution, if the value of the variable (usually the time) is not present in the 'data' table,
 *	the class linearly interpolates the value between the two closest one. Moreover, if the value of the variable is higher
 *	than anyone present in the 'data' table, the class linearly extrapolates the value using the last two values in the table.
 *
 *   <b>EXAMPLE OF DATAFILE</b>
 *
 *  function	= '(0,0,q)'
 *  variables	=   't				q'
 *  scale		= 	'1				1'
 *  data		=  '0.000000000		1.00
 *  				0.333333333		2.00
 *  				0.666666666		3.00
 *  				1.000000000		4.00'
 *
 */
class BCInterfaceFunctionFile : public BCInterfaceFunction
//     :
//     public LifeV::Application
{
public:

	// ===================================================
	//! Public functions
	// ===================================================

	/** @name Constructors & Destructor
     */
    //@{

    //! Constructor
	/*!
	 * \param baseString	- function string
	 * \param comV 			- vector of components
	 */
	BCInterfaceFunctionFile( const std::string& baseString, const BCComV& comV );

	//! Copy constructor
	/*!
	 * \param function		- BCInterfaceFunctionFile
	 */
	BCInterfaceFunctionFile( const BCInterfaceFunctionFile& function );

    //! Destructor
    ~BCInterfaceFunctionFile() {}

    //@}



    /** @name Methods
     */
    //@{

    //! Operator =
    /*!
     * \param function - BCInterfaceFunctionFile
     */
    BCInterfaceFunctionFile& operator=( const BCInterfaceFunctionFile& function );

    //! Compare function
    /*!
     * \param fileName
     * \param comV
     */
    bool compare( const std::string& fileName, const BCComV& comV );

    //! loadData
    void loadData( void );

    //@}

private:

	// ===================================================
	//! Private variables
	// ===================================================

	std::string											M_fileName;
	std::vector<std::string>							M_variables;
	std::vector<Real>									M_scale;
	std::map< std::string, std::vector<Real> >			M_data;
	std::vector<Real>::iterator							M_dataIterator;



	// ===================================================
	//! Private functions
	// ===================================================

	/** @name Private functions
	 */
	//@{

	//! Linear interpolation (extrapolation) between two values of the data.
	inline void dataInterpolation( void );

    //@}
};

} // Namespace LifeV

#endif /* __BCInterfaceFunctionFile_H */
