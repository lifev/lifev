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
 *  @brief File containing the Multiscale Physical Coupling
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleCoupling.hpp>

namespace LifeV
{
namespace Multiscale
{

std::map< std::string, couplings_Type > multiscaleCouplingsMap;

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCoupling::MultiscaleCoupling() :
        M_ID                          (),
        M_type                        (),
        M_models                      (),
        M_couplingName                (),
        M_flags                       (),
        M_globalData                  (),
        M_couplingIndex               (),
        M_localCouplingFunctions      (),
        M_localCouplingVariables      (),
        M_localCouplingResiduals      (),
        M_timeInterpolationOrder      ( 1 ),
        M_perturbedCoupling           ( false ),
        M_comm                        ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::MultiscaleCoupling() \n";
#endif

}

// ===================================================
// Multiscale PhysicalCoupling Virtual Methods
// ===================================================
void
MultiscaleCoupling::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::setupData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    // Read Multiscale parameters
    M_couplingName = dataFile( "Multiscale/couplingName", "couplingName" );
    M_timeInterpolationOrder = dataFile( "Multiscale/timeInterpolationOrder", 0 );

    // Set the size of the local coupling variables
    M_localCouplingVariables.reserve( M_timeInterpolationOrder + 1 );
}

void
MultiscaleCoupling::showMe()
{
    std::cout << "Coupling id         = " << M_ID << std::endl
              << "Coupling name       = " << M_couplingName << std::endl
              << "Coupling type       = " << enum2String( M_type, multiscaleCouplingsMap ) << std::endl << std::endl;

    std::cout << "Models number       = " << modelsNumber() << std::endl;
    std::cout << "Models ID(s)        = ";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        std::cout << M_models[i]->ID() << " ";
    std::cout << std::endl;
    std::cout << "Models type(s)      = ";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        std::cout << enum2String( M_models[i]->type(), multiscaleModelsMap ) << " ";
    std::cout << std::endl;
    std::cout << "Flags list          = ";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        std::cout << M_flags[i] << " ";
    std::cout << std::endl << std::endl;
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleCoupling::createCouplingMap( MapEpetra& couplingMap )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::createCouplingMap( couplingMap ) \n";
#endif

    M_couplingIndex.second = couplingMap.map( Unique )->NumGlobalElements();

    couplingMap += localCouplingVariables( 0 ).map();
}

void
MultiscaleCoupling::extrapolateCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::extrapolateCouplingVariables() \n";
#endif

    // Extrapolate the coupling variables at the next time
    multiscaleVector_Type extrapolatedCouplingVariables( localCouplingVariables( 0 ) );
    interpolateCouplingVariables( M_globalData->dataTime()->nextTime(), extrapolatedCouplingVariables );

    // If we have not yet enough samples for interpolation, we add a new one
    UInt couplingVariablesSize( M_localCouplingVariables.size() );
    if ( couplingVariablesSize <= M_timeInterpolationOrder )
    {
        ++couplingVariablesSize;
        M_localCouplingVariables.push_back( multiscaleVectorPtr_Type ( new VectorEpetra( localCouplingVariables( 0 ) ) ) );
    }

    // Updating database
    for ( UInt i(1) ; i < couplingVariablesSize ; ++i )
        localCouplingVariables( couplingVariablesSize-i ) = localCouplingVariables( couplingVariablesSize-i-1 );

    localCouplingVariables( 0 ) = extrapolatedCouplingVariables;

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8200 ) << "C(" << M_couplingIndex.second + i << ") = " << ( localCouplingVariables( 0 ) )[i]  << "\n";
#endif

}

void
MultiscaleCoupling::interpolateCouplingVariables( const Real& t, multiscaleVector_Type& interpolatedCouplingVariables ) const
{
    // Coupling variables size
    UInt couplingVariablesSize( M_localCouplingVariables.size() );

    // Time container for interpolation
    timeContainer_Type timeContainer( couplingVariablesSize, 0 );
    for ( UInt i(0) ; i < couplingVariablesSize ; ++i )
        timeContainer[i] = M_globalData->dataTime()->time() - i * M_globalData->dataTime()->timeStep();

    // Lagrange interpolation
    interpolatedCouplingVariables *= 0;
    Real base(1);

    for ( UInt i(0) ; i < M_localCouplingVariables.size() ; ++i )
    {
        base = 1;
        for ( UInt j(0) ; j < M_localCouplingVariables.size() ; ++j )
            if ( j != i )
                base *= (t - timeContainer[j]) / (timeContainer[i] - timeContainer[j]);

        interpolatedCouplingVariables += localCouplingVariables( i ) * base;
    }
}

void
MultiscaleCoupling::exportJacobian( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::exportJacobian( jacobian ) \n";
#endif

    // Definitions
    bool solveLinearSystem;                          // Flag to avoid multiple solution of the same linear system
    multiscaleModelsContainer_Type perturbedModelsList; // List of perturbed model

    // Insert constant values in the jacobian (due to this coupling condition)
    insertJacobianConstantCoefficients( jacobian );

    // Set as perturbed
    M_perturbedCoupling = 0;

    // Loop on all the local coupling variables that should be perturbed
    for ( UInt column(M_couplingIndex.second) ; M_perturbedCoupling < static_cast< Int > ( M_couplingIndex.first ); ++M_perturbedCoupling, ++column )
    {
        // Build the list of models affected by the perturbation of the variable associated with this column
        perturbedModelsList = listOfPerturbedModels( M_perturbedCoupling );

        // Loop on all the models, that are influenced by the perturbation of the coupling variable
        for ( multiscaleModelsContainerIterator_Type j = perturbedModelsList.begin() ; j < perturbedModelsList.end() ; ++j )
        {
            solveLinearSystem = true;

            // Loop on all the couplings (boundary flags) that connect the j-model
            for ( UInt k(0) ; k < ( *j )->couplingsNumber() ; ++k )
                ( *j )->coupling( k )->insertJacobianDeltaCoefficients( jacobian, column, ( *j )->ID(), solveLinearSystem );
        }
    }

    // Set as unperturbed
    M_perturbedCoupling = -1;
}

void
MultiscaleCoupling::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::saveSolution() \n";
#endif

    std::ofstream output;
    output << std::scientific << std::setprecision( 15 );

    if ( M_comm->MyPID() == 0 )
    {
        std::string filename = multiscaleProblemFolder + "Step_" + number2string( multiscaleProblemStep ) + "_Coupling_" + number2string( M_ID ) + ".mfile";

        if ( M_globalData->dataTime()->isFirstTimeStep() )
        {
            output.open( filename.c_str(), std::ios::trunc );
            output << "% Coupling Type: " << enum2String( M_type, multiscaleCouplingsMap ) << std::endl << std::endl;
            output << "% TIME                     ID   FLAG FLOW RATE                STRESS" << std::endl;
        }
        else
            output.open( filename.c_str(), std::ios::app );
    }

    displayCouplingValues( output );

    if ( M_comm->MyPID() == 0 )
        output.close();
}

// ===================================================
// Get Methods
// ===================================================
void
MultiscaleCoupling::addFlagID( const UInt& flagID )
{
    M_flags.push_back( M_models.back()->flag( flagID ) );
}

UInt
MultiscaleCoupling::modelGlobalToLocalID( const UInt& ID ) const
{
    for ( UInt localID( 0 ); localID < modelsNumber(); ++localID )
        if ( M_models[localID]->ID() == ID )
            return localID;

    return 0;
}

// ===================================================
// Protected Methods
// ===================================================
void
MultiscaleCoupling::createLocalVectors()
{
    // Build a repeated list of GlobalElements
    std::vector<Int> myGlobalElements( M_couplingIndex.first );
    for ( UInt i = 0 ; i < static_cast< UInt > ( myGlobalElements.size() ) ; ++i )
        myGlobalElements[i] = i;

    // Build a repeated map for the couplings
    MapEpetra map( -1, static_cast< Int > ( myGlobalElements.size() ), &myGlobalElements[0], M_comm );

    // Create local repeated vectors
    M_localCouplingVariables.push_back( multiscaleVectorPtr_Type ( new VectorEpetra( map, Repeated ) ) );
    M_localCouplingResiduals.reset( new VectorEpetra( map, Repeated ) );
}

void
MultiscaleCoupling::resetCouplingHistory()
{
    // Reset coupling variable history
    for ( UInt i( 0 ) ; i < M_localCouplingVariables.size() ; ++i )
        localCouplingVariables( i ) = 0;
}

void
MultiscaleCoupling::importCouplingVector( const multiscaleVector_Type& globalVector, multiscaleVector_Type& localVector )
{
    Real value(0);
    for ( UInt i(0) ; i < M_couplingIndex.first ; ++i )
    {
        if ( M_comm->MyPID() == 0 )
            value = globalVector[ M_couplingIndex.second + i ];

        M_comm->Broadcast( &value, 1, 0 );

        localVector[i] = value;
    }
}

void
MultiscaleCoupling::exportCouplingVector( const multiscaleVector_Type& localVector, multiscaleVector_Type& globalVector )
{
    for ( UInt i(0) ; i < M_couplingIndex.first ; ++i )
        if ( M_comm->MyPID() == 0 )
            globalVector[ M_couplingIndex.second + i ] = localVector[i];
}

void
MultiscaleCoupling::switchErrorMessage( const multiscaleModelPtr_Type& model )
{
    multiscaleErrorCheck( ModelType, "Invalid model type ["  + enum2String( model->type(), multiscaleModelsMap ) +
                        "] for coupling type [" + enum2String( M_type, multiscaleCouplingsMap ) +"]\n" );
}

} // Namespace Multiscale
} // Namespace LifeV
