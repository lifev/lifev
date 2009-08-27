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
 *  table of discrete data (for example a discrete Flux or Pressure depending on time).
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
template <typename Operator>
class BCInterfaceFunctionFile : public virtual BCInterfaceFunction<Operator>
//     :
//     public LifeV::Application
{
public:

	typedef BCInterfaceFunction<Operator>				super;

	/** @name Constructors & Destructor
     */
    //@{

	//! Empty Constructor
	BCInterfaceFunctionFile( void );

    //! Constructor
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
	BCInterfaceFunctionFile( const BCInterfaceData<Operator>& data );

	//! Copy constructor
	/*!
	 * \param function			- BCInterfaceFunctionFile
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
     * \param function			- BCInterfaceFunctionFile
     */
    virtual BCInterfaceFunctionFile& operator=( const BCInterfaceFunctionFile& function );

    //! Set data
    /*!
	 * \param data				- BC data loaded from GetPot file
	 */
    virtual void setData( const BCInterfaceData<Operator>& data );

	//! Compare function
	/*!
	 * \param data				- BC data loaded from GetPot file
	 */
    virtual bool compare( const BCInterfaceData<Operator>& data );

    //@}

protected:

	std::string												M_fileName;

private:

	std::vector<std::string>								M_variables;
	std::vector<Real>										M_scale;
	std::map< std::string, std::vector<Real> >				M_data;
	std::vector<Real>::iterator								M_dataIterator;

	/** @name Private functions
	 */
	//@{

    //! loadData
    inline void loadData( const GetPot& dataFile );

	//! Linear interpolation (extrapolation) between two values of the data.
	inline void dataInterpolation( void );

    //@}
};

//! Factory create function
template <typename Operator>
inline BCInterfaceFunction<Operator>* createFunctionFile()
{
	return new BCInterfaceFunctionFile<Operator>();
}



// ===================================================
//! Constructors
// ===================================================
template <typename Operator>
BCInterfaceFunctionFile<Operator>::BCInterfaceFunctionFile( ) :
	BCInterfaceFunction<Operator>	( ),
	M_fileName						( ),
	M_variables						( ),
	M_scale							( ),
	M_data							( ),
	M_dataIterator					( )
{

#ifdef DEBUG
	Debug( 5022 ) << "BCInterfaceFunctionFile::BCInterfaceFunctionFile( void )" << "\n";
#endif

}



template <typename Operator>
BCInterfaceFunctionFile<Operator>::BCInterfaceFunctionFile( const BCInterfaceData<Operator>& data ) :
	BCInterfaceFunction<Operator>	( ),
	M_fileName						( ),
	M_variables						( ),
	M_scale							( ),
	M_data							( ),
	M_dataIterator					( )
{

#ifdef DEBUG
	Debug( 5022 ) << "BCInterfaceFunctionFile::BCInterfaceFunctionFile( data )" << "\n";
#endif

	this->setData( data );
}



template <typename Operator>
BCInterfaceFunctionFile<Operator>::BCInterfaceFunctionFile( const BCInterfaceFunctionFile& function ) :
	BCInterfaceFunction<Operator>	( function ),
	M_fileName						( function.M_fileName ),
	M_variables						( function.M_variables ),
	M_scale							( function.M_scale ),
	M_data							( function.M_data ),
	M_dataIterator					( function.M_dataIterator )
{
}





// ===================================================
//! Methods
// ===================================================
template <typename Operator>
BCInterfaceFunctionFile<Operator>&
BCInterfaceFunctionFile<Operator>::operator=( const BCInterfaceFunctionFile& function )
{
	if ( this != &function )
	{
		super::operator=( function );
		M_fileName		= function.M_fileName;
		M_variables		= function.M_variables;
		M_scale			= function.M_scale;
		M_data			= function.M_data;
		M_dataIterator	= function.M_dataIterator;
	}

	return *this;
}



template <typename Operator>
void
BCInterfaceFunctionFile<Operator>::setData( const BCInterfaceData<Operator>& data )
{
	M_fileName = data.get_baseString();	//The base string contains the file name

	//Load data from file
	GetPot dataFile( M_fileName );
	loadData( dataFile );

#ifdef DEBUG
	Debug( 5022 ) << "BCInterfaceFunctionFile::setData             fileName: " << M_fileName  << "\n";
	Debug( 5022 ) << "                                             function: " << data.get_baseString()  << "\n";
#endif

	//Create a new data container with the correct base string
	BCInterfaceData<Operator> newData = data;
	newData.set_baseString( dataFile( "function", "Undefined" ) ); // Now it contains the real base string
	super::setData( newData );
}



template <typename Operator>
bool
BCInterfaceFunctionFile<Operator>::compare( const BCInterfaceData<Operator>& data )
{
	return M_fileName.compare( data.get_baseString() ) == 0 && super::M_comV == data.get_comV();
}





// ===================================================
//! Private functions
// ===================================================
template <typename Operator>
inline void
BCInterfaceFunctionFile<Operator>::loadData( const GetPot& dataFile )
{

	//Set variables
	UInt variablesNumber = dataFile.vector_variable_size( "variables" );

	M_variables.clear();
	M_variables.reserve( variablesNumber );

	M_scale.clear();
	M_scale.reserve( variablesNumber );

	for ( UInt j(0) ; j < variablesNumber ; ++j )
	{
		M_variables.push_back( dataFile( "variables", "unknown", j ) );
		M_scale.push_back( dataFile( "scale", 1.0, j ) );
	}

#ifdef DEBUG
	std::stringstream output;
	output << "BCInterfaceFunctionFile::loadData           variables: ";
	for ( UInt j(0) ; j < variablesNumber ; ++j )
		output << M_variables[j] << "  ";

	output << "\n                                                           scale: ";
	for ( UInt j(0) ; j < variablesNumber ; ++j )
		output << M_scale[j] << "  ";

	Debug( 5022 ) << output.str() << "\n";
#endif



	//Load data
	UInt dataLines = dataFile.vector_variable_size( "data" ) / variablesNumber;

	M_data.clear();
	for ( UInt j(0) ; j < variablesNumber ; ++j )
		M_data[M_variables[j]].reserve( dataLines );

	for ( UInt i(0) ; i < dataLines ; ++i )
		for ( UInt j(0) ; j < variablesNumber ; ++j )
			M_data[ M_variables[j] ].push_back( M_scale[j] * dataFile( "data", 0.0, i*variablesNumber + j ) );

#ifdef DEBUG
	output.str("");
	output << "                                                 data:";
	for ( UInt i(0) ; i < dataLines ; ++i )
	{
		if (i > 0)
			output << "                                                                 ";

		for ( UInt j(0) ; j < variablesNumber ; ++j )
			output << " " << M_data[ M_variables[j] ][i];
		output << "\n";
	}
	Debug( 5022 ) << output.str();
#endif



	//Initialize iterator
	M_dataIterator = M_data[M_variables[0]].begin();
}




template <typename Operator>
inline void
BCInterfaceFunctionFile<Operator>::dataInterpolation( void )
{
	//Get variable
	Real X = super::M_parser->getVariable( M_variables[0] );

#ifdef DEBUG
	Debug( 5022 ) << "                                                    variable: " << X  << "\n";
#endif



	//Move Iterator
	for ( ; ; )
	{

#ifdef DEBUG
		Debug( 5022 ) << "                                       iterator  position   : " << static_cast<Real> ( M_dataIterator - M_data[ M_variables[0] ].begin() )  << "\n";
		Debug( 5022 ) << "                                       variable (position)  : " << *M_dataIterator << "\n";
		Debug( 5022 ) << "                                       variable (position+1): " << *(M_dataIterator+1) << "\n";
#endif

		if (X >= *M_dataIterator && X <= *(M_dataIterator+1) )
			break;

		if ( X > *M_dataIterator )
		{
			if ( M_dataIterator+1 == M_data[ M_variables[0] ].end() )
				break;
			else
				++M_dataIterator;
		}
		else
		{
			if ( M_dataIterator == M_data[ M_variables[0] ].begin() )
				break;
			else
				--M_dataIterator;
		}
	}



	//Linear interpolation (extrapolation if X > xB)
	Real xA, xB, A, B;
	for ( UInt j(1), position = static_cast<UInt> (M_dataIterator - M_data[ M_variables[0] ].begin()) ; j < static_cast<UInt> (M_variables.size()) ; ++j )
	{
		xA	= M_data[ M_variables[0] ][position];
		xB	= M_data[ M_variables[0] ][position+1];
		A 	= M_data[ M_variables[j] ][position];
		B 	= M_data[ M_variables[j] ][position+1];

		super::M_parser->setVariable( M_variables[j], A+(B-A)/(xB-xA)*(X-xA) );

#ifdef DEBUG
		Debug( 5022 ) << "                                                          " << M_variables[j] << " = " << A+(B-A)/(xB-xA)*(X-xA) << "\n";
#endif
	}
}

} // Namespace LifeV

#endif /* __BCInterfaceFunctionFile_H */
