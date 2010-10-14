//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief BCInterface_FSI
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 23-04-2009
 */

#include <lifemc/lifesolver/BCInterface_FSI.hpp>

namespace LifeV {

// ===================================================
// Constructors
// ===================================================
BCInterface_FSI< FSIOperator >::BCInterface_FSI() :
    M_FSIFunction   (),
    M_name          (),
    M_flag          (),
    M_type          (),
    M_mode          (),
    M_comV          (),
    M_base          ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface_FSI::BCInterface_FSI()" << "\n";
#endif

}

BCInterface_FSI< FSIOperator >::BCInterface_FSI( const Data_Type& data ) :
    M_FSIFunction   (),
    M_name          (),
    M_flag          (),
    M_type          (),
    M_mode          (),
    M_comV          (),
    M_base          ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface_FSI::BCInterface_FSI( data )" << "\n";
#endif

    this->SetData( data );
}

BCInterface_FSI< FSIOperator >::BCInterface_FSI( const BCInterface_FSI& fsi ) :
    M_FSIFunction   ( fsi.M_FSIFunction ),
    M_name          ( fsi.M_name ),
    M_flag          ( fsi.M_flag ),
    M_type          ( fsi.M_type ),
    M_mode          ( fsi.M_mode ),
    M_comV          ( fsi.M_comV ),
    M_base          ( fsi.M_base )
{
}

// ===================================================
// Methods
// ===================================================
BCInterface_FSI< FSIOperator >&
BCInterface_FSI< FSIOperator >::operator=( const BCInterface_FSI& fsi )
{
    if ( this != &fsi )
    {
        M_FSIFunction   = fsi.M_FSIFunction;
        M_name          = fsi.M_name;
        M_flag          = fsi.M_flag;
        M_type          = fsi.M_type;
        M_mode          = fsi.M_mode;
        M_comV          = fsi.M_comV;
        M_base          = fsi.M_base;
    }

    return *this;
}

void
BCInterface_FSI< FSIOperator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface_FSIFunctionFile::setData" << "\n";
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
    M_FSIFunction = mapFunction[ data.GetBaseString() ];

    M_name = data.GetName();
    M_flag = data.GetFlag();
    M_type = data.GetType();
    M_mode = data.GetMode();
    M_comV = data.GetComV();
}

void
BCInterface_FSI< FSIOperator >::ExportData( Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface_FSIFunctionFile::ExportData" << "\n";
#endif

    data.SetName( M_name );
    data.SetFlag( M_flag );
    data.SetType( M_type );
    data.SetMode( M_mode );
    data.SetComV( M_comV );
}

void
BCInterface_FSI< FSIOperator >::CheckMethod( const boost::shared_ptr< FSIOperator >& Oper )
{
    //Set mapMethod
    std::map< std::string, FSIMethod > mapMethod;

    mapMethod["exactJacobian"]   = EXACTJACOBIAN;
    mapMethod["fixedPoint"]      = FIXEDPOINT;
    mapMethod["monolithicGE"]    = MONOLITHIC_GE;
    mapMethod["monolithicGI"]    = MONOLITHIC_GI;
    mapMethod["steklovPoincare"] = STEKLOVPOINCARE;

    switch ( mapMethod[Oper->data().method()] )
    {
        case EXACTJACOBIAN:

    #ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            exactJacobian" << "\n";
    #endif

            CheckFunction< exactJacobian > ( Oper );

            break;

        case FIXEDPOINT:

    #ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            fixedPoint" << "\n";
    #endif

            CheckFunction< fixedPoint > ( Oper );

            break;

        case MONOLITHIC_GE:

    #ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            monolithicGE" << "\n";
    #endif

            CheckFunction< MonolithicGE >( Oper );

            break;

        case MONOLITHIC_GI:

    #ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            monolithicGI" << "\n";
    #endif

            CheckFunction< MonolithicGI >( Oper );

            break;

        case STEKLOVPOINCARE:

    #ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            steklovPoincare" << "\n";
    #endif

            //CheckFunction< steklovPoincare >( Oper );

            break;
    }
}

// ===================================================
// Get Methods
// ===================================================
BCInterface_FSI< FSIOperator >::BCFunction_Type&
BCInterface_FSI< FSIOperator >::GetBase()
{
    return *M_base;
}

} // Namespace LifeV
