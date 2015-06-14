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
 *  @brief File containing the multiscale mean normal stress coupling class
 *
 *  @date 07-08-2012
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/couplings/MultiscaleCouplingMeanNormalStress.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingMeanNormalStress::MultiscaleCouplingMeanNormalStress() :
    multiscaleCoupling_Type     ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::MultiscaleCouplingMeanNormalStress() \n";
#endif

    M_type = MeanNormalStress;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanNormalStress::setupCouplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::setupCouplingVariablesNumber() \n";
#endif

    M_couplingVariablesNumber = M_flowRateInterfaces + 1;
}

void
MultiscaleCouplingMeanNormalStress::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::setupCoupling() \n";
#endif

    if ( myModelsNumber() > 0 )
    {
        // Impose flow rate boundary conditions
        for ( UInt i ( 0 ); i < static_cast<UInt>(M_flowRateInterfaces); ++i )
            if ( myModel ( i ) )
            {
                M_localCouplingFunctions.push_back ( MultiscaleCouplingFunction ( this, i ) );
                multiscaleDynamicCast< MultiscaleInterface > ( M_models[i] )->imposeBoundaryFlowRate ( M_boundaryIDs[i], boost::bind ( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
            }

        // Impose stress boundary conditions
        for ( UInt i ( M_flowRateInterfaces ); i < modelsNumber(); ++i )
            if ( myModel ( i ) )
            {
                M_localCouplingFunctions.push_back ( MultiscaleCouplingFunction ( this, M_flowRateInterfaces ) );
                multiscaleDynamicCast< MultiscaleInterface > ( M_models[i] )->imposeBoundaryMeanNormalStress ( M_boundaryIDs[i], boost::bind ( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
            }
    }
}

void
MultiscaleCouplingMeanNormalStress::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::initializeCouplingVariables() \n";
#endif

    // Compute the flow rate coupling variables on the first M_flowRateInterfaces models
    Real localSum ( 0 );
    Real globalSum ( 0 );

    for ( UInt i ( 0 ); i < static_cast<UInt>(M_flowRateInterfaces); ++i )
    {
        if ( myModel ( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterface > ( M_models[i] )->boundaryFlowRate ( M_boundaryIDs[i] );
            if ( isModelLeaderProcess ( i ) )
            {
                localSum = myValue;
            }
        }

        // We use the SumAll() instead of the Broadcast() because this way we don't need the id of the leader process.
        M_comm->SumAll ( &localSum, &globalSum, 1 );
        if ( myModelsNumber() > 0 )
        {
            localCouplingVariables ( 0 ) [i] = globalSum;
        }

        localSum  = 0;
        globalSum = 0;
    }

    // Compute the mean normal stress coupling variable as an average of all the stresses
    for ( UInt i ( 0 ); i < modelsNumber(); ++i )
        if ( myModel ( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterface > ( M_models[i] )->boundaryMeanNormalStress ( M_boundaryIDs[i] );
            if ( isModelLeaderProcess ( i ) )
            {
                localSum += myValue;
            }
        }

    M_comm->SumAll ( &localSum, &globalSum, 1 );
    if ( myModelsNumber() > 0 )
    {
        localCouplingVariables ( 0 ) [M_flowRateInterfaces] = globalSum / modelsNumber();
    }
}

void
MultiscaleCouplingMeanNormalStress::computeCouplingResiduals()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::computeCouplingResiduals()  \n";
#endif

    // Reset coupling residual
    *M_localCouplingResiduals = 0.;

    if ( myModelsNumber() > 0 )
    {
        for ( UInt i ( 0 ); i < static_cast<UInt>(M_flowRateInterfaces); ++i )
            if ( myModel ( i ) )
            {
                Real myValueStress = multiscaleDynamicCast< MultiscaleInterface > ( M_models[i] )->boundaryMeanNormalStress ( M_boundaryIDs[i] );
                if ( isModelLeaderProcess ( i ) )
                {
                    ( *M_localCouplingResiduals ) [0]  += localCouplingVariables ( 0 ) [i];
                    ( *M_localCouplingResiduals ) [i + 1] = myValueStress - localCouplingVariables ( 0 ) [M_flowRateInterfaces];
                }
            }

        for ( UInt i ( M_flowRateInterfaces ); i < modelsNumber(); ++i )
            if ( myModel ( i ) )
            {
                Real myValueFlowRate = multiscaleDynamicCast< MultiscaleInterface > ( M_models[i] )->boundaryFlowRate ( M_boundaryIDs[i] );
                if ( isModelLeaderProcess ( i ) )
                {
                    ( *M_localCouplingResiduals ) [0]  += myValueFlowRate;
                }
            }
    }
}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanNormalStress::exportListOfPerturbedModels ( const UInt& localCouplingVariableID, multiscaleModelsContainer_Type& perturbedModelsList )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::exportListOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    if ( static_cast<Int>(localCouplingVariableID) < M_flowRateInterfaces )
    {
        if ( myModel (localCouplingVariableID) )
        {
            perturbedModelsList.reserve ( 1 );
            perturbedModelsList.push_back ( M_models[localCouplingVariableID] );
        }
    }
    else
    {
        perturbedModelsList.reserve ( modelsNumber() - M_flowRateInterfaces );
        for ( UInt i ( M_flowRateInterfaces ); i < modelsNumber(); ++i )
            if ( myModel ( i ) )
            {
                perturbedModelsList.push_back ( M_models[i] );
            }
    }
}

void
MultiscaleCouplingMeanNormalStress::insertJacobianConstantCoefficients ( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    // The constant coefficients are added by the leader process of model 0.
    if ( myModel ( 0 ) )
        if ( isModelLeaderProcess ( 0 ) )
            for ( UInt i ( 0 ); i < static_cast<UInt>(M_flowRateInterfaces); ++i )
            {
                jacobian.addToCoefficient ( M_couplingVariablesOffset,     M_couplingVariablesOffset + i,                     1 );
                jacobian.addToCoefficient ( M_couplingVariablesOffset + 1 + i, M_couplingVariablesOffset + M_flowRateInterfaces, -1 );
            }
}

void
MultiscaleCouplingMeanNormalStress::insertJacobianDeltaCoefficients ( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8220 ) << "MultiscaleCouplingMeanNormalStress::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    // Model global to local conversion
    UInt modelLocalID = modelGlobalToLocalID ( ID );
    if ( myModel ( modelLocalID ) )
    {
        Real row ( 0 );
        Real coefficient ( 0 );

        if ( static_cast<Int>(modelLocalID) >= M_flowRateInterfaces )
        {
            row = M_couplingVariablesOffset;
            coefficient = multiscaleDynamicCast< MultiscaleInterface > ( M_models[modelLocalID] )->boundaryDeltaFlowRate ( M_boundaryIDs[modelLocalID], solveLinearSystem );
        }
        else
        {
            row = M_couplingVariablesOffset + 1 + modelLocalID;
            coefficient = multiscaleDynamicCast< MultiscaleInterface > ( M_models[modelLocalID] )->boundaryDeltaMeanNormalStress ( M_boundaryIDs[modelLocalID], solveLinearSystem );
        }

        // Add the coefficient to the matrix
        if ( isModelLeaderProcess ( modelLocalID ) )
        {
            jacobian.addToCoefficient ( row, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
            debugStream ( 8220 ) << "J(" << row << "," << column << ") = " << coefficient << "\n";
#endif
        }
    }

}

} // Namespace Multiscale
} // Namespace LifeV
