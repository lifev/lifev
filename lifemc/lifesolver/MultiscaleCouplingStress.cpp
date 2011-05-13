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

#include <lifemc/lifesolver/MultiscaleCouplingStress.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingStress::MultiscaleCouplingStress() :
        multiscaleCoupling_Type     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingStress::MultiscaleCouplingStress() \n";
#endif

    M_type = Stress;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingStress::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingStress::setupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first = modelsNumber();

    //Create local vectors
    createLocalVectors();

    // Impose stress boundary condition on all the models
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        M_localCouplingFunctions.push_back( MultiscaleCouplingFunction( this, 0 ) );
        multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->imposeBoundaryStress( M_flags[i], boost::bind( &MultiscaleCouplingFunction::function, M_localCouplingFunctions.back(), _1, _2, _3, _4, _5 ) );
    }
}

void
MultiscaleCouplingStress::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingStress::initializeCouplingVariables() \n";
#endif

    localCouplingVariables( 0 ) = 0.;

    // Compute the Stress coupling variable as an average of all the stresses
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        localCouplingVariables( 0 )[0] += multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryStress( M_flags[i] );

    localCouplingVariables( 0 )[0] /= modelsNumber();

    // Compute Flow rate coupling variables
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        localCouplingVariables( 0 )[i] = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8220 ) << "C(" << M_couplingIndex.second + i << ") = " << localCouplingVariables( 0 )[i]  << "\n";
#endif

}

void
MultiscaleCouplingStress::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingStress::exportCouplingResiduals()  \n";
#endif

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        ( *M_localCouplingResiduals )[i] = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );

    for ( UInt i( 1 ); i < modelsNumber(); ++i )
    {
        ( *M_localCouplingResiduals )[0]  += localCouplingVariables( 0 )[i];
        ( *M_localCouplingResiduals )[i]  -= localCouplingVariables( 0 )[i];
    }

    exportCouplingVector( *M_localCouplingResiduals, couplingResiduals );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8220 ) << "R(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif

}

void
MultiscaleCouplingStress::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleCoupling_Type::showMe();

        std::cout << "Coupling Stress     = " << localCouplingVariables( 0 )[0] << std::endl;
        for ( UInt i( 1 ); i < modelsNumber(); ++i )
            std::cout << "Coupling FlowRate(" << i << ")= " << localCouplingVariables( 0 )[i] << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
multiscaleModelsContainer_Type
MultiscaleCouplingStress::listOfPerturbedModels( const UInt& localCouplingVariableID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingStress::listOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    multiscaleModelsContainer_Type perturbedModelsList(0);

    if ( localCouplingVariableID == 0 )
    {
        perturbedModelsList.resize( modelsNumber() );

        for ( UInt i( 0 ); i < modelsNumber(); ++i )
            perturbedModelsList[i] = M_models[i];
    }

    return perturbedModelsList;
}

void
MultiscaleCouplingStress::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingStress::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    UInt row    = M_couplingIndex.second;
    UInt column = M_couplingIndex.second;

    if ( M_comm->MyPID() == 0 )
        for ( UInt i( 1 ); i < modelsNumber(); ++i )
        {
            jacobian.addToCoefficient( row,     column + i,  1 );
            jacobian.addToCoefficient( row + i, column + i, -1 );
        }
}

void
MultiscaleCouplingStress::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MultiscaleCouplingStress::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    // Definitions
    UInt modelLocalID = modelGlobalToLocalID( ID );
    Real coefficient  = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );;
    UInt row          = M_couplingIndex.second + modelLocalID;

    // Add coefficient to the matrix
    if ( M_comm->MyPID() == 0 )
        jacobian.addToCoefficient( row, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "J(" << row << "," << column << ") = " << coefficient  << "\n";
#endif

}

void
MultiscaleCouplingStress::displayCouplingValues( std::ostream& output )
{
    Real flowRate(0), stress(0);
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        flowRate = multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[i] )->boundaryFlowRate( M_flags[i] );
        stress   = localCouplingVariables( 0 )[0];

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_globalData->dataTime()->time() << "    " << M_models[i]->ID()
            << "    " << M_flags[i]
            << "    " << flowRate
            << "    " << stress << std::endl;
    }
}

} // Namespace Multiscale
} // Namespace LifeV
