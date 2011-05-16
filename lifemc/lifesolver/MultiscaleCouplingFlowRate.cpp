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
 *  @brief File containing the Multiscale Coupling FlowRate
 *
 *  @date 04-04-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleCouplingFlowRate.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingFlowRate::MultiscaleCouplingFlowRate() :
        multiscaleCoupling_Type     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRate::MultiscaleCouplingFlowRate() \n";
#endif

    M_type = FlowRate;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingFlowRate::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRate::setupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first = modelsNumber();

    //Create local vectors
    createLocalVectors();

    // Impose FlowRate boundary condition on all the models
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, i ) );
        multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->imposeBoundaryFlowRate( M_flags[i], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
    }
}

void
MultiscaleCouplingFlowRate::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRate::initializeCouplingVariables() \n";
#endif

    localCouplingVariables( 0 ) = 0.;

    // Compute Flow rate coupling variables
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        localCouplingVariables( 0 )[i] = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8240 ) << "C(" << M_couplingIndex.second + i << ") = " << localCouplingVariables( 0 )[i]  << "\n";
#endif

}

void
MultiscaleCouplingFlowRate::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRate::exportCouplingResiduals()  \n";
#endif

    *M_localCouplingResiduals = 0.;

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        ( *M_localCouplingResiduals )[0] += localCouplingVariables( 0 )[i];

    for ( UInt i( 1 ); i < modelsNumber(); ++i )
    {
        ( *M_localCouplingResiduals )[i]  = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] )
                                          - multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryStress( M_flags[i] );
    }

    exportCouplingVector( *M_localCouplingResiduals, couplingResiduals );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8240 ) << "R(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif

}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
multiscaleModelsContainer_Type
MultiscaleCouplingFlowRate::listOfPerturbedModels( const UInt& localCouplingVariableID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRate::listOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    multiscaleModelsContainer_Type perturbedModelsList(1);

    perturbedModelsList[0] = M_models[localCouplingVariableID];

    return perturbedModelsList;
}

void
MultiscaleCouplingFlowRate::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRate::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    UInt row    = M_couplingIndex.second;
    UInt column = M_couplingIndex.second;

    if ( M_comm->MyPID() == 0 )
        for ( UInt i( 0 ); i < modelsNumber(); ++i )
            jacobian.addToCoefficient( row, column + i, 1 );
}

void
MultiscaleCouplingFlowRate::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRate::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    // Definitions
    UInt modelLocalID = modelGlobalToLocalID( ID );
    Real coefficient  = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaStress( M_flags[modelLocalID], solveLinearSystem );
    UInt row          = M_couplingIndex.second + modelLocalID;

    // Add coefficient to the matrix
    if ( modelLocalID == 0 )
    {
        for ( UInt i( 1 ); i < modelsNumber(); ++i )
        {
            if ( M_comm->MyPID() == 0 )
                jacobian.addToCoefficient( row + i, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
            Debug( 8240 ) << "J(" << row + i << "," << column << ") = " << coefficient  << "\n";
#endif

        }
    }
    else
    {
        if ( M_comm->MyPID() == 0 )
            jacobian.addToCoefficient( row, column, -coefficient );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 8240 ) << "J(" << row << "," << column << ") = " << -coefficient  << "\n";
#endif

    }
}

} // Namespace Multiscale
} // Namespace LifeV
