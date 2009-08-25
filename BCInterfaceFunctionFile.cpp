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
   \file BCInterfaceFunctionFile.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-09
 */

#include <lifemc/lifefem/BCInterfaceFunctionFile.hpp>

namespace LifeV {

// ===================================================
//! Constructor & Destructor
// ===================================================
BCInterfaceFunctionFile::BCInterfaceFunctionFile( ) :
	BCInterfaceFunction			( ),
	M_fileName					( ),
	M_variables					( ),
	M_scale						( ),
	M_data						( ),
	M_dataIterator				( )
{

#ifdef DEBUG
	Debug( 5022 ) << "BCInterfaceFunctionFile::BCInterfaceFunctionFile( void )" << "\n";
#endif

}

BCInterfaceFunctionFile::BCInterfaceFunctionFile( const BCInterfaceData& data ) :
	BCInterfaceFunction			( ),
	M_fileName					( ),
	M_variables					( ),
	M_scale						( ),
	M_data						( ),
	M_dataIterator				( )
{

#ifdef DEBUG
	Debug( 5022 ) << "BCInterfaceFunctionFile::BCInterfaceFunctionFile( data )" << "\n";
#endif

	this->setData( data );
}



BCInterfaceFunctionFile::BCInterfaceFunctionFile( const BCInterfaceFunctionFile& function ) :
	BCInterfaceFunction	( function ),
	M_fileName			( function.M_fileName ),
	M_variables			( function.M_variables ),
	M_scale				( function.M_scale ),
	M_data				( function.M_data ),
	M_dataIterator		( function.M_dataIterator )
{
}





// ===================================================
//! Methods
// ===================================================
BCInterfaceFunctionFile&
BCInterfaceFunctionFile::operator=( const BCInterfaceFunctionFile& function )
{
	if ( this != &function )
	{
		BCInterfaceFunction::operator=( function );
		M_fileName		= function.M_fileName;
		M_variables		= function.M_variables;
		M_scale			= function.M_scale;
		M_data			= function.M_data;
		M_dataIterator	= function.M_dataIterator;
	}

	return *this;
}



void
BCInterfaceFunctionFile::setData( const BCInterfaceData& data )
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
	BCInterfaceData newData = data;
	newData.set_baseString( dataFile( "function", "Undefined" ) ); // Now it contains the real base string
	BCInterfaceFunction::setData( newData );
}



bool
BCInterfaceFunctionFile::compare( const BCInterfaceData& data )
{
	return M_fileName.compare( data.get_baseString() ) == 0 && M_comV == data.get_comV();
}





// ===================================================
//! Private functions
// ===================================================
void
BCInterfaceFunctionFile::loadData( const GetPot& dataFile )
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





inline void
BCInterfaceFunctionFile::dataInterpolation( void )
{
	//Get variable
	Real X = M_parser->getVariable( M_variables[0] );

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

		M_parser->setVariable( M_variables[j], A+(B-A)/(xB-xA)*(X-xA) );

#ifdef DEBUG
		Debug( 5022 ) << "                                                          " << M_variables[j] << " = " << A+(B-A)/(xB-xA)*(X-xA) << "\n";
#endif
	}
}

} // Namespace LifeV
