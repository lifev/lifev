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
 *  @brief File containing the Multiscale Coupling FlowRateStress
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
        multiscaleCoupling_Type     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::MultiscaleCouplingFlowRateStress() \n";
#endif

    M_type = FlowRateStress;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingFlowRateStress::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::setupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first = 2;

    //Create local vectors
    createLocalVectors();

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

    case Windkessel0D:

        imposeFlowRate0D< MultiscaleModelWindkessel0D > ( 0 );

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

        case Windkessel0D:

            imposeStress0D< MultiscaleModelWindkessel0D > ( i );

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
        ( *M_localCouplingVariables[0] )[0] -= multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );

    // Compute the Stress
    ( *M_localCouplingVariables[0] )[1] = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] );


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
        ( *M_localCouplingResiduals )[0] -= multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );

    // Compute the Stress
    ( *M_localCouplingResiduals )[1] = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] );

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

        std::cout << "Coupling FlowRate   = " << ( *M_localCouplingVariables[0] )[0] << std::endl
                  << "Coupling Stress     = " << ( *M_localCouplingVariables[0] )[1] << std::endl << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private Multiscale PhysicalCoupling Implementation
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
    if ( modelLocalID == 0 ) // DeltaSigma coefficient
        coefficient =  multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaStress( M_flags[modelLocalID], solveLinearSystem );
    else                     // DeltaFlowRate coefficient
        coefficient = -multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

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
    Real flowRate(0), stress(0);
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        flowRate = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );
        stress   = ( *M_localCouplingVariables[0] )[1];

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_globalData->dataTime()->time() << "    " << M_models[i]->ID()
            << "    " << M_flags[i]
            << "    " << flowRate
            << "    " << stress << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================
Real
MultiscaleCouplingFlowRateStress::functionFlowRate( const Real& t, const Real&, const Real& , const Real&, const UInt& )
{
    multiscaleVector_Type interpolatedCouplingVariables( *M_localCouplingVariables[0] );
    interpolateCouplingVariables( t, interpolatedCouplingVariables );

    return interpolatedCouplingVariables[0];
}

Real
MultiscaleCouplingFlowRateStress::functionStress( const Real& t, const Real&, const Real&, const Real&, const UInt& )
{
    multiscaleVector_Type interpolatedCouplingVariables( *M_localCouplingVariables[0] );
    interpolateCouplingVariables( t, interpolatedCouplingVariables );

    return interpolatedCouplingVariables[1];
}

} // Namespace Multiscale
} // Namespace LifeV
