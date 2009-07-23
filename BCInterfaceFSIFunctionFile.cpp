/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-07-22

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
   \file BCInterfaceFSIFunctionFile.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-07-22
 */

#include <lifemc/lifefem/BCInterfaceFSIFunctionFile.hpp>

namespace LifeV {

// ===================================================
//! Constructor & Destructor
// ===================================================
BCInterfaceFSIFunctionFile::BCInterfaceFSIFunctionFile( const boost::shared_ptr<FSIOperator>& oper ) :
	BCInterfaceFunctionFile		( ),
	BCInterfaceFSIFunction		( oper )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunctionFile::BCInterfaceFSIFunctionFile( void )" << "\n";
#endif

}



BCInterfaceFSIFunctionFile::BCInterfaceFSIFunctionFile( const BCInterfaceData& data,
														const boost::shared_ptr<FSIOperator>& oper ) :
	BCInterfaceFunctionFile		( ),
	BCInterfaceFSIFunction		( oper )
{

#ifdef DEBUG
	Debug( 5025 ) << "BCInterfaceFSIFunctionFile::BCInterfaceFSIFunctionFile" << "\n";
#endif

	this->setData( data );
}



BCInterfaceFSIFunctionFile::BCInterfaceFSIFunctionFile( const BCInterfaceFSIFunctionFile& function ) :
	BCInterfaceFunction			( function ), //To avoid warning!
	BCInterfaceFunctionFile		( function ),
	BCInterfaceFSIFunction		( function )
{
}





// ===================================================
//! Methods
// ===================================================
BCInterfaceFSIFunctionFile&
BCInterfaceFSIFunctionFile::operator=( const BCInterfaceFSIFunctionFile& function )
{
	if ( this != &function )
	{
		BCInterfaceFunctionFile::operator=( function );
		BCInterfaceFSIFunction::operator=( function );
	}

	return *this;
}



void
BCInterfaceFSIFunctionFile::setData( const BCInterfaceData& data )
{

#ifdef DEBUG
	Debug( 5024 ) << "BCInterfaceFSIFunctionFile::setData" << "\n";
#endif

	BCInterfaceFunctionFile::setData( data );
	BCInterfaceFSIFunction::setData( data );
}



bool
BCInterfaceFSIFunctionFile::compare( const BCInterfaceData& data )
{
	return M_fileName.compare( data.get_baseString() ) == 0 && M_comV == data.get_comV() && M_flag == data.get_flag();
}

} // Namespace LifeV
