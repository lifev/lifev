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
 *  @brief MultiScale Coupling Stress
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-10-2009
 */

#include <lifemc/lifesolver/MS_Coupling_Stress.hpp>

namespace LifeV {

std::map< std::string, stressTypes > stressMap;

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Coupling_Stress::MS_Coupling_Stress() :
    super              (),
    M_baseStress       (),
    M_baseDeltaStress  (),
    M_stressType       ()
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::MS_Coupling_Stress() \n";
#endif

    M_type = Stress;

    //Set type of stress coupling: StaticPressure, TotalPressure
    stressMap["StaticPressure"]     = StaticPressure;
    stressMap["TotalPressure"]      = TotalPressure;
    stressMap["LagrangeMultiplier"] = LagrangeMultiplier;
}

MS_Coupling_Stress::MS_Coupling_Stress( const MS_Coupling_Stress& Stress ) :
    super             ( Stress ),
    M_baseStress      ( Stress.M_baseStress ),
    M_baseDeltaStress ( Stress.M_baseDeltaStress ),
    M_stressType      ( Stress.M_stressType )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::MS_Coupling_Stress( Stress ) \n";
#endif

}

// ===================================================
// Operators
// ===================================================
MS_Coupling_Stress&
MS_Coupling_Stress::operator=( const MS_Coupling_Stress& Stress )
{
    if ( this != &Stress )
    {
        super::operator=( Stress );
        M_baseStress      = Stress.M_baseStress;
        M_baseDeltaStress = Stress.M_baseDeltaStress;
        M_stressType      = Stress.M_stressType;
    }
    return *this;
}

// ===================================================
// MultiScale PhysicalCoupling Implementation
// ===================================================
void
MS_Coupling_Stress::SetupData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::SetupData() \n";
#endif

    super::SetupData( FileName );

    GetPot DataFile( FileName );

    //Set type of stress coupling
    M_stressType = stressMap[DataFile( "MultiScale/stressType", "StaticPressure" )];
}

void
MS_Coupling_Stress::SetupCoupling()
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::SetupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first = GetModelsNumber();

    //Create local vectors
    CreateLocalVectors();

    //Set Functions
    M_baseStress.setFunction( boost::bind( &MS_Coupling_Stress::FunctionStress, this, _1, _2, _3, _4, _5 ) );

    M_baseDeltaStress.setFunction( boost::bind( &MS_Coupling_Stress::FunctionDeltaStress, this, _1, _2, _3, _4, _5 ) );

    // Impose stress
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:

                ImposeStress< MS_Model_Fluid3D > ( i );
                ImposeDeltaStress< MS_Model_Fluid3D > ( i );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
}

void
MS_Coupling_Stress::InitializeCouplingVariables()
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::InitializeCouplingVariables() \n";
#endif

    *M_LocalCouplingVariables      = 0.;
    *M_LocalDeltaCouplingVariables = 0.;

    // Compute the Stress coupling variable as an average of all the stresses
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                ( *M_LocalCouplingVariables )[0] += MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryStress( M_flags[i], M_stressType );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
    ( *M_LocalCouplingVariables )[0] /= GetModelsNumber();

    // Compute Flux coupling variables
    for ( UInt i( 1 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {

                ( *M_LocalCouplingVariables )[i] = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
}

void
MS_Coupling_Stress::ExportCouplingResiduals( VectorType& CouplingResiduals )
{
#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::ExportCouplingResiduals()  \n";
#endif

    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:

                ( *M_LocalCouplingResiduals )[i] = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    for ( UInt i( 1 ); i < GetModelsNumber(); ++i )
    {
        ( *M_LocalCouplingResiduals )[0]  += ( *M_LocalCouplingVariables )[i];
        ( *M_LocalCouplingResiduals )[i]  -= ( *M_LocalCouplingVariables )[i];
    }

    ExportCouplingVector( *M_LocalCouplingResiduals, CouplingResiduals );
}

ModelsVector_Type
MS_Coupling_Stress::GetListOfPerturbedModels( const UInt& LocalCouplingVariableID )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::GetListOfPerturbedModels( LocalCouplingVariableID ) \n";
#endif

    ModelsVector_Type perturbedModelsList(0);

    if ( LocalCouplingVariableID == 0 )
    {
        perturbedModelsList.resize( GetModelsNumber() );

        for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
            perturbedModelsList[i] = M_models[i];
    }

    return perturbedModelsList;
}

void
MS_Coupling_Stress::InsertJacobianConstantCoefficients( MatrixType& Jacobian )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::InsertJacobianConstantCoefficients( Jacobian )  \n";
#endif

    UInt Row    = M_couplingIndex.second;
    UInt Column = M_couplingIndex.second;

    if ( M_comm->MyPID() == 0 )
        for ( UInt i( 1 ); i < GetModelsNumber(); ++i )
        {
            Jacobian.set_mat_inc( Row,     Column + i,  1 );
            Jacobian.set_mat_inc( Row + i, Column + i, -1 );
        }
}

void
MS_Coupling_Stress::InsertJacobianDeltaCoefficients( MatrixType& Jacobian, const UInt& Column, const UInt& ID, bool& SolveLinearSystem )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::InsertJacobianDeltaCoefficients( Jacobian, Column, ID, LinearSystemSolved )  \n";
#endif

    // Definitions
    Real Coefficient  = 0;
    UInt Row          = 0;
    UInt ModelLocalID = GetModelLocalID( ID );

    switch ( M_models[ModelLocalID]->GetType() )
    {
        case Fluid3D:

            Coefficient = MS_DynamicCast< MS_Model_Fluid3D >( M_models[ModelLocalID] )->GetBoundaryDeltaFlux( M_flags[ModelLocalID], SolveLinearSystem );

            break;

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[ModelLocalID] );
    }

    // Add coefficient to the matrix
    Row = M_couplingIndex.second + ModelLocalID;
    if ( M_comm->MyPID() == 0 )
        Jacobian.set_mat_inc( Row, Column, Coefficient );

#ifdef DEBUG
    Debug( 8220 ) << "J(" << Row << "," << Column << ") = " << Coefficient  << "\n";
#endif

}

void
MS_Coupling_Stress::DisplayCouplingValues( std::ostream& output )
{
    Real FlowRate(0), Stress(0), Pressure(0), DynamicPressure(0);
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
    {
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                FlowRate        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );
                Stress          = ( *M_LocalCouplingVariables )[0];
                Pressure        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryPressure( M_flags[i] );
                DynamicPressure = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryDynamicPressure( M_flags[i] );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_dataPhysics->GetDataTime()->getTime() << "    " << M_models[i]->GetID()
                                                                      << "    " << M_flags[i]
                                                                      << "    " << FlowRate
                                                                      << "    " << Stress
                                                                      << "    " << Pressure
                                                                      << "    " << DynamicPressure << std::endl;
    }
}

void
MS_Coupling_Stress::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Stress Type         = " << Enum2String( M_stressType, stressMap ) << std::endl;
        std::cout << "Coupling Stress     = " << ( *M_LocalCouplingVariables )[0] << std::endl << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================
Real
MS_Coupling_Stress::FunctionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ )
{
    return ( *M_LocalCouplingVariables )[0];
}

Real
MS_Coupling_Stress::FunctionDeltaStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalDeltaCouplingVariables )[0];
}

} // Namespace LifeV
