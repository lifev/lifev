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
 *  @brief MultiScale Coupling BoundaryCondition
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 02-09-2009
 */

#include <lifemc/lifesolver/MS_Coupling_BoundaryCondition.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition() :
    super       (),
    M_FileName  (),
    M_list      (),
    M_listSize  ()
{

#ifdef DEBUG
    Debug( 8210 ) << "MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition() \n";
#endif

    M_type = BoundaryCondition;
}

MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition( const MS_Coupling_BoundaryCondition& BoundaryCondition ) :
    super       ( BoundaryCondition ),
    M_FileName  ( BoundaryCondition.M_FileName ),
    M_list      ( BoundaryCondition.M_list ),
    M_listSize  ( BoundaryCondition.M_listSize )
{

#ifdef DEBUG
    Debug( 8210 ) << "MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition( BoundaryCondition ) \n";
#endif

}

// ===================================================
// Operators
// ===================================================
MS_Coupling_BoundaryCondition&
MS_Coupling_BoundaryCondition::operator=( const MS_Coupling_BoundaryCondition& boundaryCondition )
{
    if ( this != &boundaryCondition )
    {
        super::operator=( boundaryCondition );
        M_FileName     = boundaryCondition.M_FileName;
        M_list         = boundaryCondition.M_list;
        M_listSize     = boundaryCondition.M_listSize;
    }
    return *this;
}

// ===================================================
// MultiScale PhysicalCoupling Implementation
// ===================================================
void
MS_Coupling_BoundaryCondition::SetupData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8210 ) << "MS_Coupling_BoundaryCondition::SetupData() \n";
#endif

    super::SetupData( FileName );

    M_FileName = FileName;
    GetPot DataFile( FileName );

    //Load the list of boundary conditions
    M_listSize = DataFile.vector_variable_size( "boundary_conditions/list" );

    M_list.reserve( M_listSize );
    for ( UInt i( 0 ); i < M_listSize; ++i )
        M_list.push_back( DataFile( "boundary_conditions/list", " ", i ) );
}

void
MS_Coupling_BoundaryCondition::SetupCoupling()
{

#ifdef DEBUG
    Debug( 8210 ) << "MS_Coupling_BoundaryCondition::SetupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first  = 0;

    //Create local vectors
    CreateLocalVectors();

    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:

                ApplyBoundaryConditions3D< MS_Model_Fluid3D > ( i );
                ApplyDeltaBoundaryConditions3D< MS_Model_Fluid3D > ( i );

                break;

            case FSI3D:

                ApplyBoundaryConditions3D< MS_Model_FSI3D > ( i );
                //ApplyDeltaBoundaryConditions3D< MS_Model_FSI3D > ( i );

                break;

            case OneDimensional:

                ApplyBoundaryConditions1D< MS_Model_1D > ( i );
                //ApplyDeltaBoundaryConditions1D< MS_Model_1D > ( i );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
}

ModelsVector_Type
MS_Coupling_BoundaryCondition::GetListOfPerturbedModels( const UInt& /*LocalCouplingVariableID*/ )
{

#ifdef DEBUG
    Debug( 8210 ) << "MS_Coupling_BoundaryCondition::GetListOfPerturbedModels() \n";
#endif

    ModelsVector_Type emptyList;

    return emptyList;
}

void
MS_Coupling_BoundaryCondition::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "List size           = " << M_listSize << std::endl;
        std::cout << "List                = ";
        for ( UInt i( 0 ); i < M_listSize; ++i )
            std::cout << M_list[i] << " ";
        std::cout << std::endl << std::endl << std::endl << std::endl;
    }
}

void
MS_Coupling_BoundaryCondition::DisplayCouplingValues( std::ostream& output )
{
    Real FlowRate(0), Pressure(0), DynamicPressure(0);
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
    {
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                FlowRate        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );
                Pressure        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryPressure( M_flags[i] );
                DynamicPressure = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryDynamicPressure( M_flags[i] );

                break;
            }

            case FSI3D:
            {
                FlowRate        = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );
                Pressure        = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryPressure( M_flags[i] );
                DynamicPressure = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryDynamicPressure( M_flags[i] );

                break;
            }

            case OneDimensional:
            {
                FlowRate        = MS_DynamicCast< MS_Model_1D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );
                Pressure        = MS_DynamicCast< MS_Model_1D >( M_models[i] )->GetBoundaryPressure( M_flags[i] );
                DynamicPressure = MS_DynamicCast< MS_Model_1D >( M_models[i] )->GetBoundaryDynamicPressure( M_flags[i] );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_PhysicalData->GetDataTime()->getTime() << "    " << M_models[i]->GetID()
                                                                       << "    " << M_flags[i]
                                                                       << "    " << FlowRate
                                                                       << "    " << "NaN                   "
                                                                       << "    " << Pressure
                                                                       << "    " << DynamicPressure << std::endl;
    }
}

} // Namespace LifeV
