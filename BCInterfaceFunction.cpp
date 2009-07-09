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
BCInterfaceFunction::BCInterfaceFunction( const std::string& baseString, const BCComV& comV ) :
	M_baseString				( baseString ),
	M_comV						( comV ),
	M_base						( ),
	M_parser					( new SpiritParser( baseString ) )
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

#ifdef DEBUG
	Debug( 5021 ) << "BCInterfaceFunction::BCInterfaceFunction: " << "\n";
#endif

	UInt arguments = M_parser->countSubstring( "," ) + 1;

#ifdef DEBUG
	Debug( 5021 ) << "                                                   arguments: " << arguments  << "\n";
#endif

	if ( arguments == 1 )
		M_base.setFunction( boost::bind(&BCInterfaceFunction::Function, this, _1, _2, _3, _4, _5) );
	else
	{
		//Create the ID map
		if ( comV.size() > 1 )	// Component
			for ( ID i(0) ; i < static_cast<ID>( comV.size() ) ; ++i )
				M_mapID[ comV[i] ] = i+1;
		else					// if ( comV.front() == arguments )  Full
			for ( ID i(1) ; i <= comV.front() ; ++i )
				M_mapID[ i ] = i;

		M_base.setFunction( boost::bind(&BCInterfaceFunction::FunctionID, this, _1, _2, _3, _4, _5) );
	}
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



bool
BCInterfaceFunction::compare( const std::string& baseString, const BCComV& comV )
{
	return M_baseString.compare( baseString ) == 0 && M_comV == comV;
}





// ===================================================
//! Private functions
// ===================================================
Real
BCInterfaceFunction::Function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*id*/ )
{
	M_parser->setVariable( "t", t );
	M_parser->setVariable( "x", x );
	M_parser->setVariable( "y", y );
	M_parser->setVariable( "z", z );

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

#ifdef DEBUG
	Debug( 5021 ) << "BCInterfaceFunction::FunctionID: " << "\n";
	Debug( 5021 ) << "                                                           x: " << x  << "\n";
	Debug( 5021 ) << "                                                           y: " << y  << "\n";
	Debug( 5021 ) << "                                                           z: " << z  << "\n";
	Debug( 5021 ) << "                                                           t: " << t  << "\n";
	Debug( 5021 ) << "                                                evaluate(" << id<< ") : " << M_parser->evaluate( id )  << "\n";
#endif

	return M_parser->evaluate( M_mapID[id] );
}

} // Namespace LifeV
