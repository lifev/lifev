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
BCInterfaceFSIOperator::BCInterfaceFSIOperator( const std::string& baseString, const boost::shared_ptr<FSIOperator>& oper ) :
	M_baseString			( baseString ),
	M_FSIOperator			( oper ),
	M_base					( )
{
	//Set mapFSI
	M_mapFSI["DerFluidLoadToFluid"]					= DerFluidLoadToFluid;
	M_mapFSI["DerFluidLoadToStructure"]				= DerFluidLoadToStructure;
	M_mapFSI["DerHarmonicExtensionVelToFluid"]		= DerHarmonicExtensionVelToFluid;
	M_mapFSI["DerStructureDispToSolid"]				= DerStructureDispToSolid;
	M_mapFSI["FluidInterfaceDisp"]					= FluidInterfaceDisp;
	M_mapFSI["FluidLoadToStructure"]				= FluidLoadToStructure;
	M_mapFSI["HarmonicExtensionVelToFluid"] 		= HarmonicExtensionVelToFluid;
	M_mapFSI["SolidLoadToStructure"]				= SolidLoadToStructure;
	M_mapFSI["StructureDispToHarmonicExtension"]	= StructureDispToHarmonicExtension;
	M_mapFSI["StructureDispToSolid"]				= StructureDispToSolid;
	M_mapFSI["StructureToFluid"]					= StructureToFluid;

	checkFSIList( );
}





// ===================================================
//! Private functions
// ===================================================
inline void
BCInterfaceFSIOperator::checkFSIList( void )
{
	switch ( M_mapFSI[M_baseString] )
	{
		case DerFluidLoadToFluid :

			std::cout << "BUILD DerFluidLoadToFluid" << std::endl;

			break;

		case DerFluidLoadToStructure :

			std::cout << "BUILD DerFluidLoadToStructure" << std::endl;

			break;

		case DerHarmonicExtensionVelToFluid :

			std::cout << "BUILD DerHarmonicExtensionVelToFluid" << std::endl;

			break;

		case DerStructureDispToSolid :

			std::cout << "BUILD DerStructureDispToSolid" << std::endl;

			break;

		case FluidInterfaceDisp :

			std::cout << "BUILD FluidInterfaceDisp" << std::endl;

			break;

		case FluidLoadToStructure :

			std::cout << "BUILD FluidLoadToStructure" << std::endl;

			break;

		case HarmonicExtensionVelToFluid :

			std::cout << "BUILD HarmonicExtensionVelToFluid" << std::endl;

			break;

		case SolidLoadToStructure :

			std::cout << "BUILD SolidLoadToStructure" << std::endl;
			M_FSIOperator->setSolidLoadToStructure( M_FSIOperator->minusSigmaFluidRepeated() );

			M_base = M_FSIOperator->bcvSolidLoadToStructure();

			break;

		case StructureDispToHarmonicExtension :

			std::cout << "BUILD StructureDispToHarmonicExtension" << std::endl;

			break;

		case StructureDispToSolid :

			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFSIList        StructureDispToSolid" << "\n";

			break;

		case StructureToFluid :

			Debug( 5022 ) << "BCInterfaceFSIOperator::checkFSIList            StructureToFluid" << "\n";

			M_FSIOperator->setStructureToFluid( M_FSIOperator->veloFluidMesh() );

			M_base = M_FSIOperator->bcvStructureToFluid();

			break;
	}
}
