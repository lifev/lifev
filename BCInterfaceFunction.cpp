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





// ===================================================
//! Constructor & Destructor
// ===================================================
BCInterfaceFunction::BCInterfaceFunction( const std::string& baseString, const std::string& stringSeparator ) :
	M_mapID						( ),
	M_base						( )
{
	//Set default mapID
	M_mapID[1] = 0;
	M_mapID[2] = 0;
	M_mapID[3] = 0;

	std::vector<std::string> baseStringVector;
	boost::split( baseStringVector, baseString, boost::is_any_of(stringSeparator) );

    Debug( 5021 ) << "BCInterfaceFunction::BCInterfaceFunction: " << "\n";
	for ( UInt i = 0 ; i < baseStringVector.size() ; ++i )
	{
		boost::shared_ptr<SpiritParser> parser( new SpiritParser( baseStringVector[i] ) );
		M_parserVector.push_back( parser );

		M_mapID[i+1] = i;

	    Debug( 5021 ) << "                                          baseStringVector["<< i << "]: " << baseStringVector[i]  << "\n";
	    Debug( 5021 ) << "                                                   M_mapID["<< i << "]: " << i 					<< "\n";
	}

	buildFunctionBase();
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
	UInt i = M_mapID[id];

	M_parserVector[i]->setVariable( "t", t );
	M_parserVector[i]->setVariable( "x", x );
	M_parserVector[i]->setVariable( "y", y );
	M_parserVector[i]->setVariable( "z", z );

	Debug( 5021 ) << "BCInterfaceFunction::Function: " << "\n";
	Debug( 5021 ) << "                                                  M_mapID["<< id << "]: " << i  << "\n";
	Debug( 5021 ) << "                                                           x: " << x  << "\n";
	Debug( 5021 ) << "                                                           y: " << y  << "\n";
	Debug( 5021 ) << "                                                           z: " << z  << "\n";
	Debug( 5021 ) << "                                                           t: " << t  << "\n";

	return M_parserVector[i]->evaluate();
}
