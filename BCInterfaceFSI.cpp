/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-04-23

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
   \file BCInterfaceFSI.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-23
 */

#include <lifemc/lifefem/BCInterfaceFSI.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
BCInterfaceFSI<FSIOperator>::BCInterfaceFSI( ) :
	M_operator							( ),
	M_baseString						( ),
	M_base								( )
{

#ifdef DEBUG
	Debug( 5029 ) << "BCInterfaceFSI::BCInterfaceFSI( void )" << "\n";
#endif

}



BCInterfaceFSI<FSIOperator>::BCInterfaceFSI( const BCInterfaceData<FSIOperator>& data ) :
	M_operator							( ),
	M_baseString						( ),
	M_base								( )
{

#ifdef DEBUG
	Debug( 5029 ) << "BCInterfaceFSI::BCInterfaceFSI( data )" << "\n";
#endif

	this->setData( data );
}



BCInterfaceFSI<FSIOperator>::BCInterfaceFSI( const BCInterfaceFSI& fsi ) :
	M_operator							( fsi.M_operator ),
	M_baseString						( fsi.M_baseString ),
	M_base								( fsi.M_base ),
	M_mapMethod							( fsi.M_mapMethod ),
	M_mapFunction						( fsi.M_mapFunction )
{
}





// ===================================================
//! Methods
// ===================================================
BCInterfaceFSI<FSIOperator>&
BCInterfaceFSI<FSIOperator>::operator=( const BCInterfaceFSI& fsi )
{
    if ( this != &fsi )
    {
    	M_operator							= fsi.M_operator;
    	M_baseString						= fsi.M_baseString;
    	M_base								= fsi.M_base;
    	M_mapMethod							= fsi.M_mapMethod;
    	M_mapFunction						= fsi.M_mapFunction;
    }

	return *this;
}



void
BCInterfaceFSI<FSIOperator>::setData( const BCInterfaceData<FSIOperator>& data )
{

#ifdef DEBUG
	Debug( 5029 ) << "BCInterfaceFSIFunctionFile::setData" << "\n";
#endif

	M_operator		= data.get_operator();
	M_baseString	= data.get_baseString();

	this->checkMethod();
}



bool
BCInterfaceFSI<FSIOperator>::compare( const BCInterfaceData<FSIOperator>& data )
{
	return	M_baseString.compare( data.get_baseString() ) == 0; //&& add compare for Operator!
}





// ===================================================
//! Private functions
// ===================================================
inline void
BCInterfaceFSI<FSIOperator>::checkMethod( void )
{
	//Set mapMethod
	M_mapMethod["exactJacobian"]	= EXACTJACOBIAN;
	M_mapMethod["fixedPoint"]		= FIXEDPOINT;
	M_mapMethod["monolithic"]		= MONOLITHIC;
	M_mapMethod["steklovPoincare"]	= STEKLOVPOINCARE;

	switch ( M_mapMethod[M_operator->method()] )
	{
			case EXACTJACOBIAN :

#ifdef DEBUG
				Debug( 5029 ) << "BCInterfaceFSI::checkMethod                            exactJacobian" << "\n";
#endif

				checkFunction<exactJacobian>();

				break;

			case FIXEDPOINT :

#ifdef DEBUG
				Debug( 5029 ) << "BCInterfaceFSI::checkMethod                            fixedPoint" << "\n";
#endif

				checkFunction<fixedPoint>();

				break;

			case MONOLITHIC :

#ifdef DEBUG
				Debug( 5029 ) << "BCInterfaceFSI::checkMethod                            monolithic" << "\n";
#endif

				//checkFunction<monolithic>();

				break;

			case STEKLOVPOINCARE :

#ifdef DEBUG
				Debug( 5029 ) << "BCInterfaceFSI::checkMethod                            steklovPoincare" << "\n";
#endif

				//checkFunction<steklovPoincare>();

				break;
	}
}

} // Namespace LifeV
