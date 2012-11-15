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
 *  @date 11-10-2012
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleCouplingMeanNormalStressArea.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingMeanNormalStressArea::MultiscaleCouplingMeanNormalStressArea() :
        multiscaleCoupling_Type     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::MultiscaleCouplingMeanNormalStressArea() \n";
#endif

    M_type = MeanNormalStressArea;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanNormalStressArea::setupCouplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::setupCouplingVariablesNumber() \n";
#endif

    M_couplingVariablesNumber = M_flowRateInterfaces + 1 + 1;
}

void
MultiscaleCouplingMeanNormalStressArea::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::setupCoupling() \n";
#endif

    super_Type::setupCoupling();

    // Preliminary checks
    if ( myModelsNumber() > 0 )
    {
        if ( modelsNumber() > 2 )
            std::cerr << "!!! ERROR: MultiscaleCouplingMeanNormalStressArea does not work with more than two models !!!" << std::endl;

        //TODO: add a check on the type of the two models: they must be a FSI3D and a 1D model.
    }

    if ( myModelsNumber() > 0 )
    {
        // Impose area boundary condition on the FSI3D model
        for ( UInt i( 0 ); i < 2; ++i )
            if ( myModel( i ) )
                if ( M_models[i]->type() == FSI3D )
                {
                    M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, M_flowRateInterfaces + 1 ) );
                    multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->imposeBoundaryArea( M_boundaryIDs[i], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );

                    break;
                }
    }
}

void
MultiscaleCouplingMeanNormalStressArea::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::initializeCouplingVariables() \n";
#endif

    super_Type::initializeCouplingVariables();

    // Local variables initialization
    Real localSum( 0 );
    Real globalSum( 0 );

    // Compute the area coupling variable as an average of the two areas
    for ( UInt i( 0 ); i < 2; ++i )
        if ( myModel( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryArea( M_boundaryIDs[i] );
            if ( isModelLeaderProcess( i ) )
                localSum += myValue;
        }

    M_comm->SumAll( &localSum, &globalSum, 1 );
    if ( myModelsNumber() > 0 )
        localCouplingVariables( 0 )[M_flowRateInterfaces + 1] = globalSum / 2;

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingVariablesNumber; ++i )
        Debug( 8230 ) << "C(" << M_couplingVariablesOffset + i << ") = " << localCouplingVariables( 0 )[i]  << "\n";
#endif

}

void
MultiscaleCouplingMeanNormalStressArea::computeCouplingResiduals()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::computeCouplingResiduals()  \n";
#endif

    super_Type::computeCouplingResiduals();

    if ( myModelsNumber() > 0 )
    {
        // Impose area boundary condition on the FSI3D model
        for ( UInt i( 0 ); i < 2; ++i )
            if ( myModel( i ) )
                if ( M_models[i]->type() == OneDimensional )
                {
                    Real myValueArea = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryArea( M_boundaryIDs[i] );
                    if ( isModelLeaderProcess( i ) )
                    {
                        ( *M_localCouplingResiduals )[0]  = myValueArea - localCouplingVariables( 0 )[M_flowRateInterfaces+1];
                    }
                }
    }

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingVariablesNumber; ++i )
        Debug( 8230 ) << "R(" << M_couplingVariablesOffset + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif

}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
multiscaleModelsContainer_Type
MultiscaleCouplingMeanNormalStressArea::listOfPerturbedModels( const UInt& localCouplingVariableID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::listOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    multiscaleModelsContainer_Type perturbedModelsList(0);

    if ( localCouplingVariableID < M_couplingVariablesNumber - 1 )
    {
        if ( myModel(localCouplingVariableID) )
        {
            perturbedModelsList.reserve( 1 );
            perturbedModelsList.push_back( M_models[localCouplingVariableID] );
        }
    }
    else
    {
        perturbedModelsList.reserve( modelsNumber() - M_flowRateInterfaces );
        for ( UInt i( M_flowRateInterfaces ); i < modelsNumber(); ++i )
            if ( myModel( i ) )
                perturbedModelsList.push_back( M_models[i] );
    }

    return perturbedModelsList;
}

void
MultiscaleCouplingMeanNormalStressArea::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    // The constant coefficients are added by the leader process of model 0.
    if ( myModel( 0 ) )
        if ( isModelLeaderProcess( 0 ) )
            for ( UInt i( 0 ); i < M_flowRateInterfaces; ++i )
            {
                jacobian.addToCoefficient( M_couplingVariablesOffset,     M_couplingVariablesOffset + i,                     1 );
                jacobian.addToCoefficient( M_couplingVariablesOffset+1+i, M_couplingVariablesOffset + M_flowRateInterfaces, -1 );
            }
}

void
MultiscaleCouplingMeanNormalStressArea::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingMeanNormalStressArea::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    // Model global to local conversion
    UInt modelLocalID = modelGlobalToLocalID( ID );
    if ( myModel( modelLocalID ) )
    {
        Real row( 0 );
        Real coefficient( 0 );

        if ( modelLocalID >= M_flowRateInterfaces )
        {
            row = M_couplingVariablesOffset;
            coefficient = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_boundaryIDs[modelLocalID], solveLinearSystem );
        }
        else
        {
            row = M_couplingVariablesOffset + 1 + modelLocalID;
            coefficient = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaMeanNormalStress( M_boundaryIDs[modelLocalID], solveLinearSystem );
        }

        // Add the coefficient to the matrix
        if ( isModelLeaderProcess( modelLocalID ) )
        {
            jacobian.addToCoefficient( row, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
            Debug( 8230 ) << "J(" << row << "," << column << ") = " << coefficient << "\n";
#endif
        }
    }

}

} // Namespace Multiscale
} // Namespace LifeV
