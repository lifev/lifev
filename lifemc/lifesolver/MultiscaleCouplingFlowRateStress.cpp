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
 *  @brief File containing the MultiScale Coupling FlowRateStress
 *
 *  @date 24-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleCouplingFlowRateStress.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingFlowRateStress::MultiscaleCouplingFlowRateStress() :
        multiscaleCoupling_Type     (),
        M_baseFlowRate3D            (),
        M_baseStress3D              (),
        M_baseFlowRate1D            (),
        M_baseStress1D              (),
        M_stressType                ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::MultiscaleCouplingFlowRateStress() \n";
#endif

    M_type = FlowRateStress;
}

// ===================================================
// MultiScale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingFlowRateStress::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::SetupData() \n";
#endif

    multiscaleCoupling_Type::setupData( fileName );

    GetPot dataFile( fileName );

    //Set type of stress coupling
    M_stressType = multiscaleStressesMap[dataFile( "MultiScale/stressType", "StaticPressure" )];
}

void
MultiscaleCouplingFlowRateStress::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::setupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first  = 2;

    //Create local vectors
    createLocalVectors();

    //Set Functions
    M_baseFlowRate3D.setFunction  ( boost::bind( &MultiscaleCouplingFlowRateStress::functionFlowRate,   this, _1, _2, _3, _4, _5 ) );
    M_baseStress3D.setFunction( boost::bind( &MultiscaleCouplingFlowRateStress::functionStress, this, _1, _2, _3, _4, _5 ) );

    M_baseFlowRate1D.setFunction  ( boost::bind( &MultiscaleCouplingFlowRateStress::functionFlowRate,   this, _1, _1, _1, _1, _1 ) );
    M_baseStress1D.setFunction( boost::bind( &MultiscaleCouplingFlowRateStress::functionStress, this, _1, _1, _1, _1, _1 ) );

    // Impose flowRate
    switch ( M_models[0]->type() )
    {
    case Fluid3D:

        imposeFlowRate3D< MultiscaleModelFluid3D > ( 0 );

        break;

    case FSI3D:

        imposeFlowRate3D< MultiscaleModelFSI3D > ( 0 );

        break;

    case OneDimensional:

        imposeFlowRate1D< MultiscaleModel1D > ( 0 );

        break;

    default:

        if ( M_displayer->isLeader() )
            switchErrorMessage( M_models[0] );
    }

    // Impose stress
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:

            imposeStress3D< MultiscaleModelFluid3D > ( i );

            break;

        case FSI3D:

            imposeStress3D< MultiscaleModelFSI3D > ( i );

            break;

        case OneDimensional:

            imposeStress1D< MultiscaleModel1D > ( i );

            break;

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }
}

void
MultiscaleCouplingFlowRateStress::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::InitializeCouplingVariables() \n";
#endif

    *M_localCouplingVariables[0] = 0.;

    // Compute the FlowRate
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:
        {
            ( *M_localCouplingVariables[0] )[0] -= multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        case FSI3D:
        {
            ( *M_localCouplingVariables[0] )[0] -= multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        case OneDimensional:
        {
            ( *M_localCouplingVariables[0] )[0] -= multiscaleDynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }

    // Compute the Stress
    switch ( M_models[0]->type() )
    {
    case Fluid3D:
    {
        ( *M_localCouplingVariables[0] )[1] = multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[0] )->boundaryStress( M_flags[0], M_stressType );

        break;
    }

    case FSI3D:
    {
        ( *M_localCouplingVariables[0] )[1] = multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[0] )->boundaryStress( M_flags[0], M_stressType );

        break;
    }

    case OneDimensional:
    {
        ( *M_localCouplingVariables[0] )[1] = multiscaleDynamicCast< MultiscaleModel1D >( M_models[0] )->boundaryStress( M_flags[0], M_stressType );

        break;
    }

    default:

        if ( M_displayer->isLeader() )
            switchErrorMessage( M_models[0] );
    }

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8230 ) << "C(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingVariables[0] )[i]  << "\n";
#endif

}

void
MultiscaleCouplingFlowRateStress::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::ExportCouplingResiduals() \n";
#endif

    *M_localCouplingResiduals = 0.;

    // Compute the FlowRate
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:
        {
            ( *M_localCouplingResiduals )[0] -= multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        case FSI3D:
        {
            ( *M_localCouplingResiduals )[0] -= multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        case OneDimensional:
        {
            ( *M_localCouplingResiduals )[0] -= multiscaleDynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }

    // Compute the Stress
    switch ( M_models[0]->type() )
    {
    case Fluid3D:
    {
        ( *M_localCouplingResiduals )[1] = multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[0] )->boundaryStress( M_flags[0], M_stressType );

        break;
    }

    case FSI3D:
    {
        ( *M_localCouplingResiduals )[1] = multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[0] )->boundaryStress( M_flags[0], M_stressType );

        break;
    }

    case OneDimensional:
    {
        ( *M_localCouplingResiduals )[1] = multiscaleDynamicCast< MultiscaleModel1D >( M_models[0] )->boundaryStress( M_flags[0], M_stressType );

        break;
    }

    default:

        if ( M_displayer->isLeader() )
            switchErrorMessage( M_models[0] );
    }

    *M_localCouplingResiduals -= *M_localCouplingVariables[0];

    exportCouplingVector( *M_localCouplingResiduals, couplingResiduals );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8230 ) << "R(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif
}

void
MultiscaleCouplingFlowRateStress::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleCoupling_Type::showMe();

        std::cout << "Stress Type         = " << Enum2String( M_stressType, multiscaleStressesMap ) << std::endl;
        std::cout << "Coupling FlowRate   = " << ( *M_localCouplingVariables[0] )[0] << std::endl
                  << "Coupling Stress     = " << ( *M_localCouplingVariables[0] )[1] << std::endl << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private MultiScale PhysicalCoupling Implementation
// ===================================================
multiscaleModelsVector_Type
MultiscaleCouplingFlowRateStress::listOfPerturbedModels( const UInt& localCouplingVariableID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::GetListOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    multiscaleModelsVector_Type perturbedModelsList(1);

    if ( localCouplingVariableID == 0 )
        perturbedModelsList[0] = M_models[0];
    else
    {
        perturbedModelsList.resize( modelsNumber() - 1 );

        for ( UInt i( 1 ); i < modelsNumber(); ++i )
            perturbedModelsList[i-1] = M_models[i];
    }

    return perturbedModelsList;
}

void
MultiscaleCouplingFlowRateStress::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::InsertJacobianConstantCoefficients( jacobian )  \n";
#endif

    UInt row = M_couplingIndex.second;

    if ( M_comm->MyPID() == 0 )
    {
        jacobian.addToCoefficient( row,     row,     -1 );
        jacobian.addToCoefficient( row + 1, row + 1, -1 );
    }
}

void
MultiscaleCouplingFlowRateStress::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::InsertJacobianDeltaCoefficients( jacobian, column, ID, LinearSystemSolved )  \n";
#endif

    // Definitions
    Real coefficient  = 0;
    UInt row          = 0;
    UInt modelLocalID = modelGlobalToLocalID( ID );

    // Compute the coefficient
    switch ( M_models[modelLocalID]->type() )
    {
    case Fluid3D:
    {

        if ( modelLocalID == 0 ) // DeltaSigma coefficient
            coefficient =  multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[modelLocalID] )->boundaryDeltaStress( M_flags[modelLocalID], solveLinearSystem, M_stressType );
        else                     // DeltaFlowRate coefficient
            coefficient = -multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

        break;
    }

    case FSI3D:
    {

        if ( modelLocalID == 0 ) // DeltaSigma coefficient
            coefficient =  multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[modelLocalID] )->boundaryDeltaStress( M_flags[modelLocalID], solveLinearSystem, M_stressType );
        else                     // DeltaFlowRate coefficient
            coefficient = -multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

        break;
    }

    case OneDimensional:
    {

        if ( modelLocalID == 0 ) // DeltaSigma coefficient
            coefficient =  multiscaleDynamicCast< MultiscaleModel1D >( M_models[modelLocalID] )->boundaryDeltaStress( M_flags[modelLocalID], solveLinearSystem, M_stressType );
        else                     // DeltaFlowRate coefficient
            coefficient = -multiscaleDynamicCast< MultiscaleModel1D >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

        break;
    }

    default:

        if ( M_displayer->isLeader() )
            switchErrorMessage( M_models[modelLocalID] );
    }

    // Compute the row
    if ( modelLocalID == 0 )
        row = M_couplingIndex.second + 1;
    else
        row = M_couplingIndex.second;

    // Add coefficient to the matrix
    if ( M_comm->MyPID() == 0 )
        jacobian.addToCoefficient( row, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "J(" << row << "," << column << ") = " << coefficient  << "\n";
#endif

}

void
MultiscaleCouplingFlowRateStress::displayCouplingValues( std::ostream& output )
{
    Real flowRate(0), stress(0), pressure(0), dynamicPressure(0);
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        switch ( M_models[i]->type() )
        {
        case Fluid3D:
        {
            flowRate        = multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            stress          = ( *M_localCouplingVariables[0] )[1];
            pressure        = multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryPressure( M_flags[i] );
            dynamicPressure = multiscaleDynamicCast< MultiscaleModelFluid3D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        case FSI3D:
        {
            flowRate        = multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            stress          = ( *M_localCouplingVariables[0] )[1];
            pressure        = multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryPressure( M_flags[i] );
            dynamicPressure = multiscaleDynamicCast< MultiscaleModelFSI3D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        case OneDimensional:
        {
            flowRate        = multiscaleDynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            stress          = ( *M_localCouplingVariables[0] )[1];
            pressure        = multiscaleDynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryPressure( M_flags[i] );
            dynamicPressure = multiscaleDynamicCast< MultiscaleModel1D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_globalData->dataTime()->time() << "    " << M_models[i]->ID()
            << "    " << M_flags[i]
            << "    " << flowRate
            << "    " << stress
            << "    " << pressure
            << "    " << dynamicPressure << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================
Real
MultiscaleCouplingFlowRateStress::functionFlowRate( const Real& t, const Real&, const Real& , const Real&, const UInt& )
{
    multiscaleVector_Type interpolatedCouplingVariables( *M_localCouplingVariables[0] );

    timeContainer_Type timeContainer( M_timeInterpolationOrder + 1 );
    for ( UInt i(0) ; i <= M_timeInterpolationOrder ; ++i )
        timeContainer[i] = M_globalData->dataTime()->time() - i * M_globalData->dataTime()->timeStep();

    interpolateCouplingVariables( timeContainer, t, interpolatedCouplingVariables );

    return interpolatedCouplingVariables[0];
}

Real
MultiscaleCouplingFlowRateStress::functionStress( const Real& t, const Real&, const Real&, const Real&, const UInt& )
{
    multiscaleVector_Type interpolatedCouplingVariables( *M_localCouplingVariables[0] );

    timeContainer_Type timeContainer( M_timeInterpolationOrder + 1 );
    for ( UInt i(0) ; i <= M_timeInterpolationOrder ; ++i )
        timeContainer[i] = M_globalData->dataTime()->time() - i * M_globalData->dataTime()->timeStep();

    interpolateCouplingVariables( timeContainer, t, interpolatedCouplingVariables );

    return interpolatedCouplingVariables[1];
}

} // Namespace Multiscale
} // Namespace LifeV
