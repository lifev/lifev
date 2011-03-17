//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the BCInterface3DFSI class
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/BCInterface3DFSI.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
BCInterface3DFSI< FSIOperator >::BCInterface3DFSI() :
        M_FSIFunction   (),
        M_name          (),
        M_flag          (),
        M_type          (),
        M_mode          (),
        M_comV          ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFunctionFSI::BCInterface3DFunctionFSI()" << "\n";
#endif

}

BCInterface3DFSI< FSIOperator >::BCInterface3DFSI( const data_Type& data ) :
        M_FSIFunction   (),
        M_name          (),
        M_flag          (),
        M_type          (),
        M_mode          (),
        M_comV          ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFunctionFSI::BCInterface3DFunctionFSI( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Methods
// ===================================================
void
BCInterface3DFSI< FSIOperator >::exportData( data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFunctionFSIFunctionFile::exportData" << "\n";
#endif

    data.setName( M_name );
    data.setFlag( M_flag );
    data.setType( M_type );
    data.setMode( M_mode );
    data.setComV( M_comV );
}

void
BCInterface3DFSI< FSIOperator >::assignFunction( const boost::shared_ptr< FSIOperator >& physicalSolver, BCVectorInterface& base )
{
    //Set mapMethod
    std::map< std::string, FSIMethod > mapMethod;

    mapMethod["exactJacobian"]   = EXACTJACOBIAN;
    mapMethod["fixedPoint"]      = FIXEDPOINT;
    mapMethod["monolithicGE"]    = MONOLITHIC_GE;
    mapMethod["monolithicGI"]    = MONOLITHIC_GI;

    switch ( mapMethod[physicalSolver->data().method()] )
    {
    case EXACTJACOBIAN:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkMethod                            exactJacobian" << "\n";
#endif

        checkFunction< FSIExactJacobian > ( physicalSolver, base );

        break;

    case FIXEDPOINT:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkMethod                            fixedPoint" << "\n";
#endif

        checkFunction< FSIFixedPoint > ( physicalSolver, base );

        break;

    case MONOLITHIC_GE:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkMethod                            monolithicGE" << "\n";
#endif

        checkFunction< FSIMonolithicGE >( physicalSolver, base );

        break;

    case MONOLITHIC_GI:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkMethod                            monolithicGI" << "\n";
#endif

        checkFunction< FSIMonolithicGI >( physicalSolver, base );

        break;

    default:

        std::cout << " !!! Warning:" << mapMethod[physicalSolver->data().method()] << " not assigned !!!" << std::endl;

    }
}

// ===================================================
// Set Methods
// ===================================================
void
BCInterface3DFSI< FSIOperator >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface3DFunctionFSI::setData" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, FSIFunction > mapFunction;
    mapFunction["DerFluidLoadToFluid"]              = DerFluidLoadToFluid;
    mapFunction["DerFluidLoadToStructure"]          = DerFluidLoadToStructure;
    mapFunction["DerHarmonicExtensionVelToFluid"]   = DerHarmonicExtensionVelToFluid;
    mapFunction["DerStructureDispToSolid"]          = DerStructureDispToSolid;
    mapFunction["FluidInterfaceDisp"]               = FluidInterfaceDisp;
    mapFunction["FluidLoadToStructure"]             = FluidLoadToStructure;
    mapFunction["HarmonicExtensionVelToFluid"]      = HarmonicExtensionVelToFluid;
    mapFunction["SolidLoadToStructure"]             = SolidLoadToStructure;
    mapFunction["StructureDispToHarmonicExtension"] = StructureDispToHarmonicExtension;
    mapFunction["StructureDispToSolid"]             = StructureDispToSolid;
    mapFunction["StructureToFluid"]                 = StructureToFluid;

    // Retrieving the strings
    M_FSIFunction = mapFunction[ data.baseString() ];

    M_name = data.name();
    M_flag = data.flag();
    M_type = data.type();
    M_mode = data.mode();
    M_comV = data.comV();
}

} // Namespace LifeV
