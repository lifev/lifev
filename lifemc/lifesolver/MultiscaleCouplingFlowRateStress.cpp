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

    // Impose stress boundary condition on the first model
    M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, 0 ) );
    multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->imposeBoundaryFlowRate( M_flags[0], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );

    // Impose stress boundary condition on all the other models
    M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, 1 ) );
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->imposeBoundaryStress( M_flags[i], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
}

void
MultiscaleCouplingFlowRateStress::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::initializeCouplingVariables() \n";
#endif

    localCouplingVariables( 0 ) = 0.;

    // Compute the FlowRate
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        localCouplingVariables( 0 )[0] -= multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );

    // Compute the Stress
    localCouplingVariables( 0 )[1] = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] );


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

    *M_localCouplingResiduals = 0.;

    // Compute the FlowRate
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        ( *M_localCouplingResiduals )[0] -= multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );

    // Compute the Stress
    ( *M_localCouplingResiduals )[1] = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] );

    *M_localCouplingResiduals -= localCouplingVariables( 0 );

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

        std::cout << "Coupling FlowRate   = " << ( localCouplingVariables( 0 ) )[0] << std::endl
                  << "Coupling Stress     = " << ( localCouplingVariables( 0 ) )[1] << std::endl << std::endl;
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
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::listOfPerturbedModels( localCouplingVariableID ) \n";
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
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::insertJacobianConstantCoefficients( jacobian )  \n";
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
    Debug( 8230 ) << "MultiscaleCouplingFlowRateStress::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    // Definitions
    UInt modelLocalID = modelGlobalToLocalID( ID );
    Real coefficient  = 0;
    UInt row          = M_couplingIndex.second;

    // Compute the coefficient
    if ( modelLocalID == 0 )
    {
        row += 1;
        coefficient =  multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaStress( M_flags[modelLocalID], solveLinearSystem );
    }
    else
        coefficient = -multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

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
        stress   = ( localCouplingVariables( 0 ) )[1];

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_globalData->dataTime()->time() << "    " << M_models[i]->ID()
            << "    " << M_flags[i]
            << "    " << flowRate
            << "    " << stress << std::endl;
    }
}

} // Namespace Multiscale
} // Namespace LifeV
