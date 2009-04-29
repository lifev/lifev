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
   \file BCInterfaceFSIOperator.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-23
 */

#include <lifemc/lifefem/BCInterfaceFSIOperator.hpp>





// ===================================================
//! Constructor & Destructor
// ===================================================
BCInterfaceFSIOperator::BCInterfaceFSIOperator( const std::string& baseString,
												const boost::shared_ptr<FSIOperator>& oper ) :
	M_baseString			( baseString ),
	M_FSIOperator			( oper ),
	M_base					( )
{
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





// ===================================================
//! Private functions
// ===================================================
inline void
BCInterfaceFSIOperator::checkMethod( void )
{
	switch ( M_mapMethod[M_FSIOperator->method()] )
	{
			case EXACTJACOBIAN :

				Debug( 5022 ) << "BCInterfaceFSIOperator::checkMethod   -> exactJacobian" << "\n";

				checkFunction<exactJacobian>();

				break;

			case FIXEDPOINT :

				Debug( 5022 ) << "BCInterfaceFSIOperator::checkMethod   -> fixedPoint" << "\n";

				checkFunction<fixedPoint>();

				break;

			case MONOLITHIC :

				Debug( 5022 ) << "BCInterfaceFSIOperator::checkMethod   -> monolithic" << "\n";

				//checkFunction<monolithic>();

				break;

			case STEKLOVPOINCARE :

				Debug( 5022 ) << "BCInterfaceFSIOperator::checkMethod   -> steklovPoincare" << "\n";

				//checkFunction<steklovPoincare>();

				break;
	}
}
