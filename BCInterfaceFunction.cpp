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
BCInterfaceFunction::BCInterfaceFunction( const std::string& baseString ) :
	M_base						( ),
	M_parser					( new SpiritParser( baseString ) )
{
    Debug( 5021 ) << "BCInterfaceFunction::BCInterfaceFunction: " << "\n";

	buildFunctionBase();
}



BCInterfaceFunction::BCInterfaceFunction( const BCInterfaceFunction& function ) :
	M_base		( function.M_base ),
	M_parser	( function.M_parser )
{
}



BCInterfaceFunction&
BCInterfaceFunction::operator=( const BCInterfaceFunction& function )
{
    if ( this != &function )
    {
    	M_base		= function.M_base;
    	M_parser 	= function.M_parser;
    }

	return *this;
}


// ===================================================
//! Private functions
// ===================================================
void
BCInterfaceFunction::buildFunctionBase( void )
{
	M_base.setFunction( getFunction() );
}



BCInterfaceFunction::function_type
BCInterfaceFunction::getFunction( void )
{
	return boost::bind(&BCInterfaceFunction::Function, this, _1, _2, _3, _4, _5);
}



Real
BCInterfaceFunction::Function( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id )
{
	M_parser->setVariable( "t", t );
	M_parser->setVariable( "x", x );
	M_parser->setVariable( "y", y );
	M_parser->setVariable( "z", z );

	Debug( 5021 ) << "BCInterfaceFunction::Function: " << "\n";
	Debug( 5021 ) << "                                                           x: " << x  << "\n";
	Debug( 5021 ) << "                                                           y: " << y  << "\n";
	Debug( 5021 ) << "                                                           z: " << z  << "\n";
	Debug( 5021 ) << "                                                           t: " << t  << "\n";
	Debug( 5021 ) << "                                                evaluate(" << id<< ") : " << M_parser->evaluate( id )  << "\n";

	return M_parser->evaluate( id );
}

} // Namespace LifeV
