/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-07-17

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
   \file BCInterfaceData.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-17
 */

#include <lifemc/lifefem/BCInterfaceData.hpp>

namespace LifeV {

BCInterfaceData::BCInterfaceData( ) :
	M_name			( ),
	M_flag			( ),
	M_type			( ),
	M_mode			( ),
	M_comV			( ),
	M_baseString	( )
{}



BCInterfaceData::BCInterfaceData( const BCInterfaceData& data ) :
		M_name					( data.M_name ),
		M_flag					( data.M_flag ),
		M_type					( data.M_type ),
		M_mode					( data.M_mode ),
		M_comV					( data.M_comV ),
		M_baseString			( data.M_baseString )
{}


BCInterfaceData&
BCInterfaceData::operator=( const BCInterfaceData& data )
{
	if ( this != &data )
	{
		M_name			= data.M_name;
		M_flag			= data.M_flag;
		M_type			= data.M_type;
		M_mode			= data.M_mode;
		M_comV			= data.M_comV;
		M_baseString	= data.M_baseString;
	}

	return *this;
}





// ===================================================
//! Methods
// ===================================================
void
BCInterfaceData::reset_comV( const UInt& components )
{
	M_comV.clear();
	M_comV.reserve( components );
}

} // Namespace LifeV
