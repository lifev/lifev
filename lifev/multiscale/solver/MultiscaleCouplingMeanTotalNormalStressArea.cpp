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

#include <lifev/multiscale/solver/MultiscaleCouplingMeanTotalNormalStressArea.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingMeanTotalNormalStressArea::MultiscaleCouplingMeanTotalNormalStressArea() :
        multiscaleCoupling_Type     ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::MultiscaleCouplingMeanTotalNormalStressArea() \n";
#endif

    M_type = MeanTotalNormalStressArea;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanTotalNormalStressArea::setupCouplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::setupCouplingVariablesNumber() \n";
#endif

    M_couplingVariablesNumber = 4;
}

void
MultiscaleCouplingMeanTotalNormalStressArea::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::setupCoupling() \n";
#endif

    super_Type::setupCoupling();

    // Preliminary checks
    if ( myModelsNumber() > 0 )
    {
        if ( modelsNumber() > 2 )
            std::cerr << "!!! ERROR: MultiscaleCouplingMeanTotalNormalStressArea does not work with more than two models !!!" << std::endl;

        //TODO: add a check on the type of the two models: they must be a FSI3D and a 1D model.
    }

    if ( myModelsNumber() > 0 )
    {
        // Impose area boundary condition on the FSI3D model
        for ( UInt i( 0 ); i < 2; ++i )
            if ( myModel( i ) )
                if ( M_models[i]->type() == FSI3D )
                {
                    M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, 3 ) );
                    multiscaleDynamicCast< MultiscaleInterface >( M_models[i] )->imposeBoundaryArea( M_boundaryIDs[i], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );

                    break;
                }
    }
}

void
MultiscaleCouplingMeanTotalNormalStressArea::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::initializeCouplingVariables() \n";
#endif

    super_Type::initializeCouplingVariables();

    // Local variables initialization
    Real localSum( 0 );
    Real globalSum( 0 );

    // Compute the area coupling variable as an average of the two areas
    for ( UInt i( 0 ); i < 2; ++i )
        if ( myModel( i ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterface >( M_models[i] )->boundaryArea( M_boundaryIDs[i] );
            if ( isModelLeaderProcess( i ) )
                localSum += myValue;
        }

    M_comm->SumAll( &localSum, &globalSum, 1 );
    if ( myModelsNumber() > 0 )
        localCouplingVariables( 0 )[3] = globalSum / 2;
}

void
MultiscaleCouplingMeanTotalNormalStressArea::computeCouplingResiduals()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::computeCouplingResiduals()  \n";
#endif

    super_Type::computeCouplingResiduals();

    if ( myModelsNumber() > 0 )
    {
        // Impose area boundary condition on the FSI3D model
        for ( UInt i( 0 ); i < 2; ++i )
            if ( myModel( i ) )
                if ( M_models[i]->type() == FSI1D )
                {
                    Real myValueArea = multiscaleDynamicCast< MultiscaleInterface >( M_models[i] )->boundaryArea( M_boundaryIDs[i] );
                    if ( isModelLeaderProcess( i ) )
                    {
                        ( *M_localCouplingResiduals )[3]  = myValueArea - localCouplingVariables( 0 )[3];
                    }
                }
    }
}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanTotalNormalStressArea::exportListOfPerturbedModels( const UInt& localCouplingVariableID, multiscaleModelsContainer_Type& perturbedModelsList )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::exportListOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    if ( localCouplingVariableID == 3 )
    {
        for ( UInt i( 0 ); i < 2; ++i )
            if ( myModel( i ) )
                if ( M_models[i]->type() == FSI3D )
                    perturbedModelsList.push_back( M_models[i] );
    }
    else
        super_Type::exportListOfPerturbedModels( localCouplingVariableID, perturbedModelsList );
}

void
MultiscaleCouplingMeanTotalNormalStressArea::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    super_Type::insertJacobianConstantCoefficients( jacobian );

    // The constant coefficients are added by the leader process of model 0.
    if ( myModel( 0 ) )
        if ( isModelLeaderProcess( 0 ) )
            jacobian.addToCoefficient( M_couplingVariablesOffset + 3, M_couplingVariablesOffset + 3, -1 );
}

void
MultiscaleCouplingMeanTotalNormalStressArea::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8260 ) << "MultiscaleCouplingMeanTotalNormalStressArea::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    super_Type::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem );

    // Model global to local conversion
    UInt modelLocalID = modelGlobalToLocalID( ID );
    if ( myModel( modelLocalID ) )
        if ( M_models[modelLocalID]->type() == FSI1D )
        {
            Real row( 0 );
            Real coefficient( 0 );

            row = M_couplingVariablesOffset + 3;
            coefficient = multiscaleDynamicCast< MultiscaleInterface >( M_models[modelLocalID] )->boundaryDeltaArea( M_boundaryIDs[modelLocalID], solveLinearSystem );

            // Add the coefficient to the matrix
            if ( isModelLeaderProcess( modelLocalID ) )
            {
                jacobian.addToCoefficient( row, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
                debugStream( 8260 ) << "J(" << row << "," << column << ") = " << coefficient << "\n";
#endif
            }
        }

}

} // Namespace Multiscale
} // Namespace LifeV
