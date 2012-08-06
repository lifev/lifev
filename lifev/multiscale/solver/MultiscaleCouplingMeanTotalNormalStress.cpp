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
 *  @brief File containing the Multiscale Coupling Stress
 *
 *  @date 20-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleCouplingMeanTotalNormalStress.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingMeanTotalNormalStress::MultiscaleCouplingMeanTotalNormalStress() :
        multiscaleCoupling_Type     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingMeanTotalNormalStress::MultiscaleCouplingMeanTotalNormalStress() \n";
#endif

    M_type = MeanTotalNormalStress;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanTotalNormalStress::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingMeanTotalNormalStress::setupCoupling() \n";
#endif

    if ( myModelsNumber() > 0 )
    {
        // Set the number of coupling variables
        M_couplingVariablesNumber = modelsNumber()+1;

        // Impose stress boundary condition on all the models
        for ( UInt i( 0 ); i < modelsNumber(); ++i )
            if ( myModel( i ) )
            {
                M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, i ) );
                multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->imposeBoundaryStress( M_flags[i], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
            }
    }

    // Create local vectors
    createLocalVectors();
}

void
MultiscaleCouplingMeanTotalNormalStress::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingMeanTotalNormalStress::initializeCouplingVariables() \n";
#endif

    // Compute the mean total normal stress coupling variable as an average of all the stresses
    Real localSum( 0 );
    Real globalSum( 0 );

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        if ( myModel( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryTotalStress( M_flags[i] );
            if ( isModelLeaderProcess( i ) )
                localSum += myValue;
        }

    M_comm->SumAll( &localSum, &globalSum, 1 );
    if ( myModelsNumber() > 0 )
        localCouplingVariables( 0 )[M_couplingVariablesNumber-1] = globalSum / modelsNumber();

    // Compute the mean normal stress coupling variables on the other models
    localSum  = 0;
    globalSum = 0;

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        if ( myModel( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryStress( M_flags[i] );
            if ( isModelLeaderProcess( i ) )
                localSum = myValue;
        }

        // We use the SumAll() instead of the Broadcast() because this way we don't need the id of the leader process.
        M_comm->SumAll( &localSum, &globalSum, 1 );
        if ( myModelsNumber() > 0 )
            localCouplingVariables( 0 )[i] = globalSum;

        localSum  = 0;
        globalSum = 0;
    }

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingVariablesNumber; ++i )
        Debug( 8220 ) << "C(" << M_couplingVariablesOffset + i << ") = " << localCouplingVariables( 0 )[i]  << "\n";
#endif

}

void
MultiscaleCouplingMeanTotalNormalStress::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingMeanTotalNormalStress::exportCouplingResiduals()  \n";
#endif

    // Reset coupling residual
    *M_localCouplingResiduals = 0.;

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        if ( myModel( i ) )
        {
            Real myValueFlowRate    = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            Real myValueTotalStress = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryTotalStress( M_flags[i] );
            if ( isModelLeaderProcess( i ) )
            {
                ( *M_localCouplingResiduals )[0] += myValueFlowRate;
                ( *M_localCouplingResiduals )[i+1] = myValueTotalStress - localCouplingVariables( 0 )[M_couplingVariablesNumber-1];
            }
        }

    exportCouplingVector( couplingResiduals, *M_localCouplingResiduals );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingVariablesNumber; ++i )
        Debug( 8220 ) << "R(" << M_couplingVariablesOffset + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif

}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
multiscaleModelsContainer_Type
MultiscaleCouplingMeanTotalNormalStress::listOfPerturbedModels( const UInt& localCouplingVariableID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingMeanTotalNormalStress::listOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    multiscaleModelsContainer_Type perturbedModelsList(0);

    if ( localCouplingVariableID < M_couplingVariablesNumber - 1 )
        if ( myModel(localCouplingVariableID) )
        {
            perturbedModelsList.reserve( 1 );
            perturbedModelsList.push_back( M_models[localCouplingVariableID] );
        }
    return perturbedModelsList;
}

void
MultiscaleCouplingMeanTotalNormalStress::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingMeanTotalNormalStress::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    // The constant coefficients are added by the leader process of model 0.
    if ( myModel( 0 ) )
        if ( isModelLeaderProcess( 0 ) )
        {
            UInt row    = M_couplingVariablesOffset + 1;
            UInt column = M_couplingVariablesOffset + M_couplingVariablesNumber - 1;

            for ( UInt i( 0 ); i < modelsNumber(); ++i )
            {
                jacobian.addToCoefficient( row + i, column, -1 );
            }
        }
}

void
MultiscaleCouplingMeanTotalNormalStress::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingMeanTotalNormalStress::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    // Model global to local conversion
    UInt modelLocalID = modelGlobalToLocalID( ID );
    if ( myModel( modelLocalID ) )
    {
        // Compute the coefficient
        Real coefficientFlowRate              = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );
        Real coefficientMeanTotalNormalStress = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaTotalStress( M_flags[modelLocalID], solveLinearSystem );

        // Add the coefficient to the matrix
        if ( isModelLeaderProcess( modelLocalID ) )
        {
            UInt row = M_couplingVariablesOffset + 1 + modelLocalID;
            jacobian.addToCoefficient( M_couplingVariablesOffset, column, coefficientFlowRate );
            jacobian.addToCoefficient( row, column, coefficientMeanTotalNormalStress );

#ifdef HAVE_LIFEV_DEBUG
            Debug( 8220 ) << "J(" << M_couplingVariablesOffset << "," << column << ") = " << coefficientFlowRate  << "\n";
            Debug( 8220 ) << "J(" << row << "," << column << ") = " << coefficientMeanTotalNormalStress  << "\n";
#endif
        }
    }

}

} // Namespace Multiscale
} // Namespace LifeV
