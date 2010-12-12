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
 *  @brief File containing the MultiScale Physical Coupling
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleCoupling.hpp>

namespace LifeV
{

std::map< std::string, couplings_Type > MS_couplingsMap;

UInt MultiscaleCoupling::M_couplingsNumber = 0;

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
        M_localCouplingVariables      (),
        M_localCouplingResiduals      (),
        M_timeInterpolationOrder      ( 1 ),
        M_perturbedCoupling           ( false ),
        M_comm                        (),
        M_displayer                   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::MultiscaleCoupling() \n";
#endif

    M_ID = M_couplingsNumber++;
}

// ===================================================
// MultiScale PhysicalCoupling Virtual Methods
// ===================================================
void
MultiscaleCoupling::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::SetupData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    // Read multiscale parameters
    M_couplingName = dataFile( "MultiScale/couplingName", "couplingName" );
    M_timeInterpolationOrder = dataFile( "MultiScale/timeInterpolationOrder", 0 );

    // Set the size of the local coupling variables
    M_localCouplingVariables.reserve( M_timeInterpolationOrder + 1 );
}

void
MultiscaleCoupling::showMe()
{
    std::cout << "Coupling id         = " << M_ID << std::endl
              << "Coupling name       = " << M_couplingName << std::endl
              << "Coupling type       = " << Enum2String( M_type, MS_couplingsMap ) << std::endl << std::endl;

    std::cout << "Models number       = " << modelsNumber() << std::endl;
    std::cout << "Models ID(s)        = ";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        std::cout << M_models[i]->ID() << " ";
    std::cout << std::endl;
    std::cout << "Models type(s)      = ";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        std::cout << Enum2String( M_models[i]->type(), MS_modelsMap ) << " ";
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
MultiscaleCoupling::createCouplingMap( EpetraMap& couplingMap )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::CreateCouplingMap( couplingMap ) \n";
#endif

    M_couplingIndex.second = couplingMap.getMap( Unique )->NumGlobalElements();

    couplingMap += M_localCouplingVariables[0]->getMap();
}

void
MultiscaleCoupling::importCouplingVariables( const MS_Vector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::ImportCouplingVariables( couplingVariables ) \n";
#endif

    importCouplingVector( couplingVariables, *M_localCouplingVariables[0] );
}

void
MultiscaleCoupling::exportCouplingVariables( MS_Vector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::ExportCouplingVariables( couplingVariables ) \n";
#endif

    exportCouplingVector( *M_localCouplingVariables[0], couplingVariables );
}

void
MultiscaleCoupling::extrapolateCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::ExtrapolateCouplingVariables() \n";
#endif

    MS_Vector_Type extrapolatedCouplingVariables( *M_localCouplingVariables[0] );
    UInt couplingVariablesSize( M_localCouplingVariables.size() );

    // Time container for interpolation
    timeContainer_Type timeContainer( couplingVariablesSize, 0 );
    for ( UInt i(0) ; i < couplingVariablesSize ; ++i )
        timeContainer[i] = M_globalData->dataTime()->getTime() - i * M_globalData->dataTime()->getTimeStep();

    // Interpolate the coupling variables at the next time
    interpolateCouplingVariables( timeContainer, M_globalData->dataTime()->getNextTime(), extrapolatedCouplingVariables );

    // If we have not yet enough samples for interpolation, we add a new one
    if ( couplingVariablesSize <= M_timeInterpolationOrder )
    {
        ++couplingVariablesSize;
        M_localCouplingVariables.push_back( MS_Vector_PtrType ( new EpetraVector( *M_localCouplingVariables[0] ) ) );
    }

    // Updating database
    for ( UInt i(1) ; i < couplingVariablesSize ; ++i )
        *M_localCouplingVariables[couplingVariablesSize-i] = *M_localCouplingVariables[couplingVariablesSize-i-1];

    *M_localCouplingVariables[0] = extrapolatedCouplingVariables;

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8200 ) << "C(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingVariables[0] )[i]  << "\n";
#endif

}

bool
MultiscaleCoupling::isPerturbed() const
{
    return M_perturbedCoupling == -1 ? false : true;
}

void
MultiscaleCoupling::exportJacobian( MS_Matrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::exportJacobian( jacobian ) \n";
#endif

    // Definitions
    bool solveLinearSystem;                   // Flag to avoid multiple solution of the same linear system
    MS_ModelsVector_Type perturbedModelsList; // List of perturbed model

    // Insert constant values in the jacobian (due to this coupling condition)
    insertJacobianConstantCoefficients( jacobian );

    // Set as perturbed
    M_perturbedCoupling = 0;

    // Loop on all the local coupling variables that should be perturbed
    for ( UInt column(M_couplingIndex.second) ; M_perturbedCoupling < static_cast< Int > ( M_couplingIndex.first ); ++M_perturbedCoupling, ++column )
    {
        // Build the list of models affected by the perturbation
        perturbedModelsList = listOfPerturbedModels( M_perturbedCoupling );

        // Loop on all the models, that could be influenced by the perturbation of the coupling variable
        for ( MS_ModelsVector_Iterator j = perturbedModelsList.begin() ; j < perturbedModelsList.end() ; ++j )
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
    Debug( 8200 ) << "MultiscaleCoupling::SaveSolution() \n";
#endif

    std::ofstream output;
    output << std::scientific << std::setprecision( 15 );

    if ( M_comm->MyPID() == 0 )
    {
        std::string filename = MS_ProblemFolder + "Step_" + number2string( MS_ProblemStep ) + "_Coupling_" + number2string( M_ID ) + ".mfile";

        if ( M_globalData->dataTime()->isFirstTimeStep() )
        {
            output.open( filename.c_str(), std::ios::trunc );
            output << "% Coupling Type: " << Enum2String( M_type, MS_couplingsMap ) << std::endl << std::endl;
            output << "% TIME                     ID   FLAG FLUX                     STRESS                    S. PRESSURE              D. PRESSURE" << std::endl;
        }
        else
            output.open( filename.c_str(), std::ios::app );
    }

    displayCouplingValues( output );

    if ( M_comm->MyPID() == 0 )
        output.close();
}

void
MultiscaleCoupling::clearModelsList()
{
    M_models.clear();
}

// ===================================================
// Set Methods
// ===================================================
void
MultiscaleCoupling::setID( const UInt& id )
{
    M_ID = id;
}

void
MultiscaleCoupling::addModel( const MS_Model_PtrType& model )
{
    M_models.push_back( model );
}

void
MultiscaleCoupling::addFlag( const BCFlag& flag )
{
    M_flags.push_back( flag );
}

void
MultiscaleCoupling::addFlagID( const UInt& flagID )
{
    M_flags.push_back( M_models.back()->flag( flagID ) );
}

void
MultiscaleCoupling::setGlobalData( const MS_GlobalDataContainer_PtrType& globalData )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::SetupData( globalData ) \n";
#endif

    M_globalData = globalData;
}

void
MultiscaleCoupling::setCommunicator( const MS_Comm_PtrType& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MultiscaleCoupling::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

// ===================================================
// Get Methods
// ===================================================
const UInt&
MultiscaleCoupling::ID() const
{
    return M_ID;
}

const couplings_Type&
MultiscaleCoupling::type() const
{
    return M_type;
}

const std::string&
MultiscaleCoupling::couplingName() const
{
    return M_couplingName;
}

UInt
MultiscaleCoupling::modelsNumber() const
{
    return static_cast< UInt > ( M_models.size() );
}

UInt
MultiscaleCoupling::modelGlobalToLocalID( const UInt& ID ) const
{
    for ( UInt localID( 0 ); localID < modelsNumber(); ++localID )
        if ( M_models[localID]->ID() == ID )
            return localID;

    return -1;
}

MS_Model_PtrType
MultiscaleCoupling::model( const UInt& localID ) const
{
    return M_models[localID];
}

const BCFlag&
MultiscaleCoupling::flag( const UInt& localID ) const
{
    return M_flags[localID];
}

const UInt&
MultiscaleCoupling::couplingVariablesNumber() const
{
    return M_couplingIndex.first;
}

const Int&
MultiscaleCoupling::perturbedCoupling() const
{
    return M_perturbedCoupling;
}

const MS_Vector_Type&
MultiscaleCoupling::residual() const
{
    return *M_localCouplingResiduals;
}

const UInt&
MultiscaleCoupling::timeInterpolationOrder() const
{
    return M_timeInterpolationOrder;
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
    EpetraMap map( -1, static_cast< Int > ( myGlobalElements.size() ), &myGlobalElements[0], 0, M_comm );

    // Create local repeated vectors
    M_localCouplingVariables.push_back( MS_Vector_PtrType ( new EpetraVector( map, Repeated ) ) );
    M_localCouplingResiduals.reset( new EpetraVector( map, Repeated ) );
}

void
MultiscaleCoupling::importCouplingVector( const MS_Vector_Type& globalVector, MS_Vector_Type& localVector )
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
MultiscaleCoupling::exportCouplingVector( const MS_Vector_Type& localVector, MS_Vector_Type& globalVector )
{
    for ( UInt i(0) ; i < M_couplingIndex.first ; ++i )
        if ( M_comm->MyPID() == 0 )
            globalVector[ M_couplingIndex.second + i ] = localVector[i];
}

void
MultiscaleCoupling::interpolateCouplingVariables( const timeContainer_Type& timeContainer, const Real& t,
                                                   MS_Vector_Type& interpolatedCouplingVariables )
{
    // Lagrange interpolation
    interpolatedCouplingVariables *= 0;
    Real base(1);

    for ( UInt i(0) ; i < M_localCouplingVariables.size() ; ++i )
    {
        base = 1;
        for ( UInt j(0) ; j < M_localCouplingVariables.size() ; ++j )
            if ( j != i )
                base *= (t - timeContainer[j]) / (timeContainer[i] - timeContainer[j]);

        interpolatedCouplingVariables += *M_localCouplingVariables[i] * base;
    }
}

void
MultiscaleCoupling::switchErrorMessage( const MS_Model_PtrType& model )
{
    MS_ErrorCheck( MS_ModelType, "Invalid model type ["  + Enum2String( model->type(), MS_modelsMap ) +
                   "] for coupling type [" + Enum2String( M_type, MS_couplingsMap ) +"]\n" );
}

} // Namespace LifeV
