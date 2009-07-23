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
//! Constructor & Destructor
// ===================================================
BCInterfaceFSI::BCInterfaceFSI( const BCInterfaceData& data,
								const boost::shared_ptr<FSIOperator>& oper ) :
	M_baseString			( data.get_baseString() ),
	M_FSIOperator			( oper ),
	M_base					( )
{
#ifdef DEBUG
	Debug( 5023 ) << "BCInterfaceFSI::BCInterfaceFSI" << "\n";
#endif

	//Set mapMethod
	M_mapMethod["exactJacobian"]						= EXACTJACOBIAN;
	M_mapMethod["fixedPoint"]							= FIXEDPOINT;
	M_mapMethod["monolithic"]							= MONOLITHIC;
	M_mapMethod["steklovPoincare"]						= STEKLOVPOINCARE;

	//Set mapFunction
	M_mapFunction["DerFluidLoadToFluid"]				= DerFluidLoadToFluid;
	M_mapFunction["DerFluidLoadToStructure"]			= DerFluidLoadToStructure;
	M_mapFunction["DerHarmonicExtensionVelToFluid"]		= DerHarmonicExtensionVelToFluid;
	M_mapFunction["DerStructureDispToSolid"]			= DerStructureDispToSolid;
	M_mapFunction["FluidInterfaceDisp"]					= FluidInterfaceDisp;
	M_mapFunction["FluidLoadToStructure"]				= FluidLoadToStructure;
	M_mapFunction["HarmonicExtensionVelToFluid"] 		= HarmonicExtensionVelToFluid;
	M_mapFunction["SolidLoadToStructure"]				= SolidLoadToStructure;
	M_mapFunction["StructureDispToHarmonicExtension"]	= StructureDispToHarmonicExtension;
	M_mapFunction["StructureDispToSolid"]				= StructureDispToSolid;
	M_mapFunction["StructureToFluid"]					= StructureToFluid;

	checkMethod( );
}



BCInterfaceFSI::BCInterfaceFSI( const BCInterfaceFSI& fsiOperator ) :
	M_baseString	( fsiOperator.M_baseString ),
	M_FSIOperator	( fsiOperator.M_FSIOperator ),
	M_base			( fsiOperator.M_base ),
	M_mapMethod		( fsiOperator.M_mapMethod ),
	M_mapFunction	( fsiOperator.M_mapFunction )
{
}



BCInterfaceFSI&
BCInterfaceFSI::operator=( const BCInterfaceFSI& fsiOperator )
{
    if ( this != &fsiOperator )
    {
    	M_baseString	= fsiOperator.M_baseString;
    	M_FSIOperator	= fsiOperator.M_FSIOperator;
    	M_base			= fsiOperator.M_base;
    	M_mapMethod		= fsiOperator.M_mapMethod;
    	M_mapFunction	= fsiOperator.M_mapFunction;
    }

	return *this;
}





// ===================================================
//! Private functions
// ===================================================
inline void
BCInterfaceFSI::checkMethod( void )
{
	switch ( M_mapMethod[M_FSIOperator->method()] )
	{
			case EXACTJACOBIAN :

#ifdef DEBUG
				Debug( 5023 ) << "BCInterfaceFSI::checkMethod                            exactJacobian" << "\n";
#endif

				checkFunction<exactJacobian>();

				break;

			case FIXEDPOINT :

#ifdef DEBUG
				Debug( 5023 ) << "BCInterfaceFSI::checkMethod                            fixedPoint" << "\n";
#endif

				checkFunction<fixedPoint>();

				break;

			case MONOLITHIC :

#ifdef DEBUG
				Debug( 5023 ) << "BCInterfaceFSI::checkMethod                            monolithic" << "\n";
#endif

				//checkFunction<monolithic>();

				break;

			case STEKLOVPOINCARE :

#ifdef DEBUG
				Debug( 5023 ) << "BCInterfaceFSI::checkMethod                            steklovPoincare" << "\n";
#endif

				//checkFunction<steklovPoincare>();

				break;
	}
}

} // Namespace LifeV
