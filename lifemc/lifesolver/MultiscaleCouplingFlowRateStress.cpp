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

    if ( myModelsNumber() > 0 )
    {
        // Set the number of coupling variables
        M_couplingIndex.first = 2;

        // Impose flow rate boundary condition on the first model
        if ( myModel( 0 ) )
        {
            M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, 0 ) );
            multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->imposeBoundaryFlowRate( M_flags[0], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
        }

        // Impose stress boundary condition on all the other models
        for ( UInt i( 1 ); i < modelsNumber(); ++i )
            if ( myModel( i ) )
            {
                M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, 1 ) );
                multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->imposeBoundaryStress( M_flags[i], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
            }
    }

    // Create local vectors
    createLocalVectors();
}

void
MultiscaleCouplingFlowRateStress::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::initializeCouplingVariables() \n";
#endif

    // Compute the flow rate coupling variable summing the flow rate of all the models but the first one
    Real localSum( 0 );
    Real globalSum( 0 );

    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        if ( myModel( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            if ( isModelLeaderProcess( i ) )
                localSum -= myValue;
        }
    
    // Sum the flow rate on all the models (but the first one)
    M_comm->SumAll( &localSum, &globalSum, 1 );
    if ( myModelsNumber() > 0 )
        localCouplingVariables( 0 )[0] = globalSum;

    // Compute the stress on the first model, then broadcast it with the others
    localSum  = 0;
    globalSum = 0;
    
    if ( myModel( 0 ) )
    {
        Real myValue = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] );
        if ( isModelLeaderProcess( 0 ) )
            localSum = myValue;
    }

    M_comm->SumAll( &localSum, &globalSum, 1 );
    if ( myModelsNumber() > 0 )
        localCouplingVariables( 0 )[1] = globalSum;

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8230 ) << "C(" << M_couplingIndex.second + i << ") = " << localCouplingVariables( 0 )[i]  << "\n";
#endif

}

void
MultiscaleCouplingFlowRateStress::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::exportCouplingResiduals() \n";
#endif

    // Reset coupling residual
    *M_localCouplingResiduals = 0.;

    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        if ( myModel( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            if ( isModelLeaderProcess( i ) )
                ( *M_localCouplingResiduals )[0] -= myValue;
        }

    if ( myModelsNumber() > 0 )
        if ( myModel( 0 ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] );
            if ( isModelLeaderProcess( 0 ) )
            {
                ( *M_localCouplingResiduals )[1] = myValue;
                *M_localCouplingResiduals -= localCouplingVariables( 0 );
            }
        }

    exportCouplingVector( couplingResiduals, *M_localCouplingResiduals );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8230 ) << "R(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif
}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
multiscaleModelsContainer_Type
MultiscaleCouplingFlowRateStress::listOfPerturbedModels( const UInt& localCouplingVariableID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::listOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    multiscaleModelsContainer_Type perturbedModelsList(0);

    if ( localCouplingVariableID == 0 )
    {
        if ( myModel(localCouplingVariableID) )
        {
            perturbedModelsList.reserve( 1 );
            perturbedModelsList.push_back( M_models[0] );
        }
    }
    else
    {
        perturbedModelsList.reserve( myModelsNumber() );

        for ( UInt i( 1 ); i < modelsNumber(); ++i )
            if ( myModel(i) )
                perturbedModelsList.push_back( M_models[i] );
    }

    return perturbedModelsList;
}

void
MultiscaleCouplingFlowRateStress::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    // The constant coefficients are added by the leader process of model 0.
    if ( myModel( 0 ) )
        if ( isModelLeaderProcess( 0 ) )
        {
            UInt row = M_couplingIndex.second;

            jacobian.addToCoefficient( row,     row,     -1 );
            jacobian.addToCoefficient( row + 1, row + 1, -1 );
        }
}

void
MultiscaleCouplingFlowRateStress::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    // Model global to local conversion
    UInt modelLocalID = modelGlobalToLocalID( ID );
    if ( myModel( modelLocalID ) )
    {
        Real coefficient = 0;
        UInt row         = M_couplingIndex.second;

        // Compute the coefficient
        if ( modelLocalID == 0 )
        {
            row += 1;
            coefficient =  multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaStress( M_flags[modelLocalID], solveLinearSystem );
        }
        else
            coefficient = -multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

        // Add the coefficient to the matrix
        if ( isModelLeaderProcess( modelLocalID ) )
        {
            jacobian.addToCoefficient( row, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
            Debug( 8230 ) << "J(" << row << "," << column << ") = " << coefficient  << "\n";
#endif
        }
    }

}

} // Namespace Multiscale
} // Namespace LifeV
