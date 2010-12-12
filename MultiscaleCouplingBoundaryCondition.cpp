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
 *  @brief File containing the MultiScale Coupling BoundaryCondition
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleCouplingBoundaryCondition.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingBoundaryCondition::MultiscaleCouplingBoundaryCondition() :
        MS_Coupling_Type       (),
        M_fileName             (),
        M_list                 (),
        M_listSize             ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8210 ) << "MultiscaleCouplingBoundaryCondition::MultiscaleCouplingBoundaryCondition() \n";
#endif

    M_type = BoundaryCondition;
}

// ===================================================
// MultiScale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingBoundaryCondition::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8210 ) << "MultiscaleCouplingBoundaryCondition::SetupData() \n";
#endif

    MS_Coupling_Type::setupData( fileName );

    M_fileName = fileName;
    GetPot dataFile( fileName );

    //Load the list of boundary conditions
    M_listSize = dataFile.vector_variable_size( "boundary_conditions/list" );

    M_list.reserve( M_listSize );
    for ( UInt i( 0 ); i < M_listSize; ++i )
        M_list.push_back( dataFile( "boundary_conditions/list", " ", i ) );
}

void
MultiscaleCouplingBoundaryCondition::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8210 ) << "MultiscaleCouplingBoundaryCondition::SetupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first  = 0;

    //Create local vectors
    createLocalVectors();

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:

            applyBoundaryConditions3D< MultiscaleModelFluid3D > ( i );

            break;

        case FSI3D:

            applyBoundaryConditions3D< MultiscaleModelFSI3D > ( i );

            break;

        case OneDimensional:

            applyBoundaryConditions1D< MultiscaleModel1D > ( i );

            break;

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }
}

void
MultiscaleCouplingBoundaryCondition::showMe()
{
    if ( M_displayer->isLeader() )
    {
        MS_Coupling_Type::showMe();

        std::cout << "List size           = " << M_listSize << std::endl;
        std::cout << "List                = ";
        for ( UInt i( 0 ); i < M_listSize; ++i )
            std::cout << M_list[i] << " ";
        std::cout << std::endl << std::endl << std::endl << std::endl;
    }
}

// ===================================================
// Private MultiScale PhysicalCoupling Implementation
// ===================================================
MS_ModelsVector_Type
MultiscaleCouplingBoundaryCondition::listOfPerturbedModels( const UInt& /*localCouplingVariableID*/ )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8210 ) << "MultiscaleCouplingBoundaryCondition::GetListOfPerturbedModels() \n";
#endif

    MS_ModelsVector_Type emptyList;

    return emptyList;
}

void
MultiscaleCouplingBoundaryCondition::displayCouplingValues( std::ostream& output )
{
    Real FlowRate(0), Pressure(0), DynamicPressure(0);
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        switch ( M_models[i]->type() )
        {
        case Fluid3D:
        {
            FlowRate        = MS_DynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            Pressure        = MS_DynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryPressure( M_flags[i] );
            DynamicPressure = MS_DynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        case FSI3D:
        {
            FlowRate        = MS_DynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            Pressure        = MS_DynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryPressure( M_flags[i] );
            DynamicPressure = MS_DynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        case OneDimensional:
        {
            FlowRate        = MS_DynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            Pressure        = MS_DynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryPressure( M_flags[i] );
            DynamicPressure = MS_DynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_globalData->dataTime()->getTime() << "    " << M_models[i]->ID()
            << "    " << M_flags[i]
            << "    " << FlowRate
            << "    " << "NaN                   "
            << "    " << Pressure
            << "    " << DynamicPressure << std::endl;
    }
}

} // Namespace LifeV
