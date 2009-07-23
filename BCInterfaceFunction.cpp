/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-04-06

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
   \file BCInterfaceFunction.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-06
 */

#include <lifemc/lifefem/BCInterfaceFunction.hpp>

namespace LifeV {

// ===================================================
//! Constructor & Destructor
// ===================================================
BCInterfaceFunction::BCInterfaceFunction( ) :
	M_baseString				( ),
	M_comV						( ),
	M_base						( ),
	M_parser					( ),
	M_mapID						( )
{

#ifdef DEBUG
	Debug( 5021 ) << "BCInterfaceFunction::BCInterfaceFunction( void )" << "\n";
#endif

}



BCInterfaceFunction::BCInterfaceFunction( const BCInterfaceData& data ) :
	M_baseString				( ),
	M_comV						( ),
	M_base						( ),
	M_parser					( ),
	M_mapID						( )
{

#ifdef DEBUG
	Debug( 5021 ) << "BCInterfaceFunction::BCInterfaceFunction" << "\n";
#endif

	this->setData( data );
}



BCInterfaceFunction::BCInterfaceFunction( const BCInterfaceFunction& function ) :
	M_baseString	( function.M_baseString ),
	M_comV			( function.M_comV ),
	M_base			( function.M_base ),
	M_parser		( function.M_parser ),
	M_mapID			( function.M_mapID )
{
}





// ===================================================
//! Methods
// ===================================================
BCInterfaceFunction&
BCInterfaceFunction::operator=( const BCInterfaceFunction& function )
{
    if ( this != &function )
    {
    	M_baseString	= function.M_baseString;
    	M_comV			= function.M_comV;
    	M_base			= function.M_base;
    	M_parser 		= function.M_parser;
    	M_mapID			= function.M_mapID;
    }

	return *this;
}



void
BCInterfaceFunction::setData( const BCInterfaceData& data )
{

#ifdef DEBUG
	Debug( 5022 ) << "BCInterfaceFunction::setData" << "\n";
#endif

	M_comV			= data.get_comV();
	M_baseString	= data.get_baseString();

	//boost::shared_ptr<SpiritParser> emptyParser( );
	//if ( M_parser == emptyParser )
	if ( M_parser )
		M_parser->setString( M_baseString );
	else
		M_parser.reset( new SpiritParser( M_baseString ) ); // INVERTITI

	setFunction();
}



bool
BCInterfaceFunction::compare( const BCInterfaceData& data )
{
	return M_baseString.compare( data.get_baseString() ) == 0 && M_comV == data.get_comV();
}





// ===================================================
//! Private functions
// ===================================================
void
BCInterfaceFunction::setFunction( void )
{
	/*
	 * MODE          COMPONENT     FUNCTION      |      COMV.SIZE     ARGUMENTS     INTERFACEFUNCTION
	 * ------------------------------------------|---------------------------------------------------
	 *                                           |
	 * COMPONENT     2             x*y*z         |      1             1             Function
	 * FULL          3             x*y*z         |      1             1             Function
	 * FULL          1             x*y*z         |      1             1             Function
	 * FULL          3             (y*z,x*z,x*y) |      1             3             FunctionID
	 * FULL          2             (x,y)         |      1             2             FunctionID
	 * COMPONENT     '1 3'         (x,y)         |      2             2             FunctionID
	 */

	UInt arguments = M_parser->countSubstring( "," ) + 1;

#ifdef DEBUG
	Debug( 5021 ) << "BCInterfaceFunction::setFunction            arguments: " << arguments  << "\n";
#endif

	if ( arguments == 1 )
		M_base.setFunction( boost::bind(&BCInterfaceFunction::Function, this, _1, _2, _3, _4, _5) );
	else
	{
		//Create the ID map
		if ( M_comV.size() > 1 )	// Component
			for ( ID i(0) ; i < static_cast<ID>( M_comV.size() ) ; ++i )
				M_mapID[ M_comV[i] ] = i+1;
		else					// if ( M_comV.front() == arguments )  Full
			for ( ID i(1) ; i <= M_comV.front() ; ++i )
				M_mapID[ i ] = i;

		M_base.setFunction( boost::bind(&BCInterfaceFunction::FunctionID, this, _1, _2, _3, _4, _5) );
	}
}



Real
BCInterfaceFunction::Function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/ )
{
	M_parser->setVariable( "t", t );
	M_parser->setVariable( "x", x );
	M_parser->setVariable( "y", y );
	M_parser->setVariable( "z", z );

	this->dataInterpolation();

	this->addFSIVariables( t );

#ifdef DEBUG
	Debug( 5021 ) << "BCInterfaceFunction::Function: " << "\n";
	Debug( 5021 ) << "                                                           x: " << x  << "\n";
	Debug( 5021 ) << "                                                           y: " << y  << "\n";
	Debug( 5021 ) << "                                                           z: " << z  << "\n";
	Debug( 5021 ) << "                                                           t: " << t  << "\n";
#endif

	return M_parser->evaluate( 1 );
}



Real
BCInterfaceFunction::FunctionID( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id )
{
	M_parser->setVariable( "t", t );
	M_parser->setVariable( "x", x );
	M_parser->setVariable( "y", y );
	M_parser->setVariable( "z", z );

	this->dataInterpolation();

	this->addFSIVariables( t );

#ifdef DEBUG
	Debug( 5021 ) << "BCInterfaceFunction::FunctionID: " << "\n";
	Debug( 5021 ) << "                                                           x: " << x  << "\n";
	Debug( 5021 ) << "                                                           y: " << y  << "\n";
	Debug( 5021 ) << "                                                           z: " << z  << "\n";
	Debug( 5021 ) << "                                                           t: " << t  << "\n";
	Debug( 5021 ) << "                                                          id: " << id  << "\n";
	Debug( 5021 ) << "                                                evaluate(" << M_mapID[id] << ") : " << M_parser->evaluate( M_mapID[id] )  << "\n";
#endif

	return M_parser->evaluate( M_mapID[id] );
}

} // Namespace LifeV
