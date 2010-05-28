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
 *  @brief MultiScale Coupling FluxStress
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 24-08-2009
 */

#include <lifemc/lifesolver/MS_Coupling_FluxStress.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Coupling_FluxStress::MS_Coupling_FluxStress() :
    super                (),
    M_baseFlux3D         (),
    M_baseStress3D       (),
    M_baseDeltaFlux3D    (),
    M_baseDeltaStress3D  (),
    M_baseFlux1D         (),
    M_baseStress1D       (),
    M_baseDeltaFlux1D    (),
    M_baseDeltaStress1D  (),
    M_stressType         ()
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::MS_Coupling_FluxStress() \n";
#endif

    M_type = FluxStress;

    //Set type of stress coupling: StaticPressure, TotalPressure, LagrangeMultiplier
    stressMap["StaticPressure"]     = StaticPressure;
    stressMap["TotalPressure"]      = TotalPressure;
    stressMap["LagrangeMultiplier"] = LagrangeMultiplier;
}

MS_Coupling_FluxStress::MS_Coupling_FluxStress( const MS_Coupling_FluxStress& FluxStress ) :
    super                ( FluxStress ),
    M_baseFlux3D         ( FluxStress.M_baseFlux3D ),
    M_baseStress3D       ( FluxStress.M_baseStress3D ),
    M_baseDeltaFlux3D    ( FluxStress.M_baseDeltaFlux3D ),
    M_baseDeltaStress3D  ( FluxStress.M_baseDeltaStress3D ),
    M_baseFlux1D         ( FluxStress.M_baseFlux1D ),
    M_baseStress1D       ( FluxStress.M_baseStress1D ),
    M_baseDeltaFlux1D    ( FluxStress.M_baseDeltaFlux1D ),
    M_baseDeltaStress1D  ( FluxStress.M_baseDeltaStress1D ),
    M_stressType         ( FluxStress.M_stressType )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::MS_Coupling_FluxStress( FluxStress ) \n";
#endif

}

// ===================================================
// Operators
// ===================================================
MS_Coupling_FluxStress&
MS_Coupling_FluxStress::operator=( const MS_Coupling_FluxStress& FluxStress )
{
    if ( this != &FluxStress )
    {
        super::operator=( FluxStress );
        M_baseFlux3D         = FluxStress.M_baseFlux3D;
        M_baseStress3D       = FluxStress.M_baseStress3D;
        M_baseDeltaFlux3D    = FluxStress.M_baseDeltaFlux3D;
        M_baseDeltaStress3D  = FluxStress.M_baseDeltaStress3D;
        M_baseFlux1D         = FluxStress.M_baseFlux1D;
        M_baseStress1D       = FluxStress.M_baseStress1D;
        M_baseDeltaFlux1D    = FluxStress.M_baseDeltaFlux1D;
        M_baseDeltaStress1D  = FluxStress.M_baseDeltaStress1D;
        M_stressType         = FluxStress.M_stressType;
    }
    return *this;
}

// ===================================================
// MultiScale PhysicalCoupling Implementation
// ===================================================
void
MS_Coupling_FluxStress::SetupData( const std::string& FileName )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::SetupData() \n";
#endif

    super::SetupData( FileName );

    GetPot DataFile( FileName );

    //Set type of stress coupling
    M_stressType = stressMap[DataFile( "MultiScale/stressType", "StaticPressure" )];
}

void
MS_Coupling_FluxStress::SetupCoupling()
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::SetupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first  = 2;

    //Create local vectors
    CreateLocalVectors();

    //Set Functions
    M_baseFlux3D.setFunction  ( boost::bind( &MS_Coupling_FluxStress::FunctionFlux3D,   this, _1, _2, _3, _4, _5 ) );
    M_baseStress3D.setFunction( boost::bind( &MS_Coupling_FluxStress::FunctionStress3D, this, _1, _2, _3, _4, _5 ) );

    M_baseDeltaFlux3D.setFunction  ( boost::bind( &MS_Coupling_FluxStress::FunctionDeltaFlux3D,   this, _1, _2, _3, _4, _5 ) );
    M_baseDeltaStress3D.setFunction( boost::bind( &MS_Coupling_FluxStress::FunctionDeltaStress3D, this, _1, _2, _3, _4, _5 ) );

    M_baseFlux1D.setFunction  ( boost::bind( &MS_Coupling_FluxStress::FunctionFlux1D,   this, _1 ) );
    M_baseStress1D.setFunction( boost::bind( &MS_Coupling_FluxStress::FunctionStress1D, this, _1 ) );

    M_baseDeltaFlux1D.setFunction  ( boost::bind( &MS_Coupling_FluxStress::FunctionDeltaFlux1D,   this, _1 ) );
    M_baseDeltaStress1D.setFunction( boost::bind( &MS_Coupling_FluxStress::FunctionDeltaStress1D, this, _1 ) );

    // Impose FlowRate
    switch ( M_models[0]->GetType() )
    {
        case Fluid3D:

            ImposeFlux3D< MS_Model_Fluid3D > ( 0 );
            ImposeDeltaFlux3D< MS_Model_Fluid3D > ( 0 );

            break;

        case FSI3D:

            ImposeFlux3D< MS_Model_FSI3D > ( 0 );
            //ImposeDeltaFlux3D< MS_Model_FSI3D > ( 0 );

            break;

        case OneDimensional:

            ImposeFlux1D< MS_Model_1D > ( 0 );
            //ImposeDeltaFlux1D< MS_Model_1D > ( 0 );

            break;

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[0] );
    }

    // Impose stress
    for ( UInt i( 1 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:

                ImposeStress3D< MS_Model_Fluid3D > ( i );
                ImposeDeltaStress3D< MS_Model_Fluid3D > ( i );

                break;

            case FSI3D:

                ImposeStress3D< MS_Model_FSI3D > ( i );
                //ImposeDeltaStress3D< MS_Model_FSI3D > ( i );

                break;

            case OneDimensional:

                ImposeStress1D< MS_Model_1D > ( i );
                //ImposeDeltaStress1D< MS_Model_1D > ( i );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
}

void
MS_Coupling_FluxStress::InitializeCouplingVariables()
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::InitializeCouplingVariables() \n";
#endif

    *M_LocalCouplingVariables      = 0.;
    *M_LocalDeltaCouplingVariables = 0.;

    // Compute the FlowRate
    for ( UInt i( 1 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                ( *M_LocalCouplingVariables )[0] -= MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;
            }

            case FSI3D:
            {
                ( *M_LocalCouplingVariables )[0] -= MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;
            }

            case OneDimensional:
            {
                ( *M_LocalCouplingVariables )[0] -= MS_DynamicCast< MS_Model_1D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    // Compute the Stress
    switch ( M_models[0]->GetType() )
    {
        case Fluid3D:
        {
            ( *M_LocalCouplingVariables )[1] = MS_DynamicCast< MS_Model_Fluid3D >( M_models[0] )->GetBoundaryStress( M_flags[0], M_stressType );

            break;
        }

        case FSI3D:
        {
            ( *M_LocalCouplingVariables )[1] = MS_DynamicCast< MS_Model_FSI3D >( M_models[0] )->GetBoundaryStress( M_flags[0], M_stressType );

            break;
        }

        case OneDimensional:
        {
            ( *M_LocalCouplingVariables )[1] = MS_DynamicCast< MS_Model_1D >( M_models[0] )->GetBoundaryStress( M_flags[0], M_stressType );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[0] );
    }
}

void
MS_Coupling_FluxStress::ExportCouplingResiduals( VectorType& CouplingResiduals )
{
#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::ExportCouplingResiduals() \n";
#endif

    *M_LocalCouplingResiduals = 0.;

    // Compute the FlowRate
    for ( UInt i( 1 ); i < GetModelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                ( *M_LocalCouplingResiduals )[0] -= MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;
            }

            case FSI3D:
            {
                ( *M_LocalCouplingResiduals )[0] -= MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;
            }

            case OneDimensional:
            {
                ( *M_LocalCouplingResiduals )[0] -= MS_DynamicCast< MS_Model_1D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    // Compute the Stress
    switch ( M_models[0]->GetType() )
    {
        case Fluid3D:
        {
            ( *M_LocalCouplingResiduals )[1] = MS_DynamicCast< MS_Model_Fluid3D >( M_models[0] )->GetBoundaryStress( M_flags[0], M_stressType );

            break;
        }

        case FSI3D:
        {
            ( *M_LocalCouplingResiduals )[1] = MS_DynamicCast< MS_Model_FSI3D >( M_models[0] )->GetBoundaryStress( M_flags[0], M_stressType );

            break;
        }

        case OneDimensional:
        {
            ( *M_LocalCouplingResiduals )[1] = MS_DynamicCast< MS_Model_1D >( M_models[0] )->GetBoundaryStress( M_flags[0], M_stressType );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[0] );
    }

    *M_LocalCouplingResiduals -= *M_LocalCouplingVariables;

    ExportCouplingVector( *M_LocalCouplingResiduals, CouplingResiduals );
}

ModelsVector_Type
MS_Coupling_FluxStress::GetListOfPerturbedModels( const UInt& LocalCouplingVariableID )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::GetListOfPerturbedModels( LocalCouplingVariableID ) \n";
#endif

    ModelsVector_Type perturbedModelsList(1);

    if ( LocalCouplingVariableID == 0 )
        perturbedModelsList[0] = M_models[0];
    else
    {
        perturbedModelsList.resize( GetModelsNumber() - 1 );

        for ( UInt i( 1 ); i < GetModelsNumber(); ++i )
            perturbedModelsList[i-1] = M_models[i];
    }

    return perturbedModelsList;
}

void
MS_Coupling_FluxStress::InsertJacobianConstantCoefficients( MatrixType& Jacobian )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::InsertJacobianConstantCoefficients( Jacobian )  \n";
#endif

    UInt Row = M_couplingIndex.second;

    if ( M_comm->MyPID() == 0 )
    {
        Jacobian.set_mat_inc( Row,     Row,     -1 );
        Jacobian.set_mat_inc( Row + 1, Row + 1, -1 );
    }
}

void
MS_Coupling_FluxStress::InsertJacobianDeltaCoefficients( MatrixType& Jacobian, const UInt& Column, const UInt& ID, bool& SolveLinearSystem )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::InsertJacobianDeltaCoefficients( Jacobian, Column, ID, LinearSystemSolved )  \n";
#endif

    // Definitions
    Real Coefficient  = 0;
    UInt Row          = 0;
    UInt ModelLocalID = GetModelLocalID( ID );

    // Compute the coefficient
    switch ( M_models[ModelLocalID]->GetType() )
    {
        case Fluid3D:
        {

            if ( ModelLocalID == 0 ) // DeltaSigma coefficient
                Coefficient =  MS_DynamicCast< MS_Model_Fluid3D >( M_models[ModelLocalID] )->GetBoundaryDeltaStress( M_flags[ModelLocalID], SolveLinearSystem, M_stressType );
            else                     // DeltaFlux coefficient
                Coefficient = -MS_DynamicCast< MS_Model_Fluid3D >( M_models[ModelLocalID] )->GetBoundaryDeltaFlux(  M_flags[ModelLocalID], SolveLinearSystem );

            break;
        }

        case FSI3D:
        {

            if ( ModelLocalID == 0 ) // DeltaSigma coefficient
                Coefficient =  MS_DynamicCast< MS_Model_FSI3D >( M_models[ModelLocalID] )->GetBoundaryDeltaStress( M_flags[ModelLocalID], SolveLinearSystem, M_stressType );
            else                     // DeltaFlux coefficient
                Coefficient = -MS_DynamicCast< MS_Model_FSI3D >( M_models[ModelLocalID] )->GetBoundaryDeltaFlux(  M_flags[ModelLocalID], SolveLinearSystem );

            break;
        }

        case OneDimensional:
        {

            if ( ModelLocalID == 0 ) // DeltaSigma coefficient
                Coefficient =  MS_DynamicCast< MS_Model_1D >( M_models[ModelLocalID] )->GetBoundaryDeltaStress( M_flags[ModelLocalID], SolveLinearSystem, M_stressType );
            else                     // DeltaFlux coefficient
                Coefficient = -MS_DynamicCast< MS_Model_1D >( M_models[ModelLocalID] )->GetBoundaryDeltaFlux(  M_flags[ModelLocalID], SolveLinearSystem );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[ModelLocalID] );
    }

    // Compute the row
    if ( ModelLocalID == 0 )
        Row = M_couplingIndex.second + 1;
    else
        Row = M_couplingIndex.second;

    // Add coefficient to the matrix
    if ( M_comm->MyPID() == 0 )
        Jacobian.set_mat_inc( Row, Column, Coefficient );

#ifdef DEBUG
    Debug( 8230 ) << "J(" << Row << "," << Column << ") = " << Coefficient  << "\n";
#endif
}

void
MS_Coupling_FluxStress::DisplayCouplingValues( std::ostream& output )
{
    Real FlowRate(0), Stress(0), Pressure(0), DynamicPressure(0);
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
    {
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                FlowRate        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );
                Stress          = ( *M_LocalCouplingVariables )[1];
                Pressure        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryPressure( M_flags[i] );
                DynamicPressure = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->GetBoundaryDynamicPressure( M_flags[i] );

                break;
            }

            case FSI3D:
            {
                FlowRate        = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );
                Stress          = ( *M_LocalCouplingVariables )[1];
                Pressure        = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryPressure( M_flags[i] );
                DynamicPressure = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->GetBoundaryDynamicPressure( M_flags[i] );

                break;
            }

            case OneDimensional:
            {
                FlowRate        = MS_DynamicCast< MS_Model_1D >( M_models[i] )->GetBoundaryFlowRate( M_flags[i] );
                Stress          = ( *M_LocalCouplingVariables )[1];
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
                                                                       << "    " << Stress
                                                                       << "    " << Pressure
                                                                       << "    " << DynamicPressure << std::endl;
    }
}

void
MS_Coupling_FluxStress::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Stress Type         = " << Enum2String( M_stressType, stressMap ) << std::endl;
        std::cout << "Coupling FlowRate   = " << ( *M_LocalCouplingVariables )[0] << std::endl
                  << "Coupling Stress     = " << ( *M_LocalCouplingVariables )[1] << std::endl << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================
Real
MS_Coupling_FluxStress::FunctionFlux3D( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalCouplingVariables )[0];
}

Real
MS_Coupling_FluxStress::FunctionStress3D( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalCouplingVariables )[1];
}

Real
MS_Coupling_FluxStress::FunctionDeltaFlux3D( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalDeltaCouplingVariables )[0];
}

Real
MS_Coupling_FluxStress::FunctionDeltaStress3D( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalDeltaCouplingVariables )[1];
}

Real
MS_Coupling_FluxStress::FunctionFlux1D( const Real& /*t*/ )
{
    return ( *M_LocalCouplingVariables )[0];
}

Real
MS_Coupling_FluxStress::FunctionStress1D( const Real& /*t*/ )
{
    return ( *M_LocalCouplingVariables )[1];
}

Real
MS_Coupling_FluxStress::FunctionDeltaFlux1D( const Real& /*t*/ )
{
    return ( *M_LocalDeltaCouplingVariables )[0];
}

Real
MS_Coupling_FluxStress::FunctionDeltaStress1D( const Real& /*t*/ )
{
    return ( *M_LocalDeltaCouplingVariables )[1];
}

} // Namespace LifeV
