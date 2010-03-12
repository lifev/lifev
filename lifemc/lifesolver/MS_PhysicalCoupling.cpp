//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Physical Coupling
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 02-09-2009
 */

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>

namespace LifeV {

std::map< std::string, couplingsTypes > couplingsMap;

UInt MS_PhysicalCoupling::M_couplingsNumber = 0;

// ===================================================
// Constructors & Destructor
// ===================================================
MS_PhysicalCoupling::MS_PhysicalCoupling() :
    M_ID                          (),
    M_type                        (),
    M_dataFile                    (),
    M_models                      (),
    M_couplingName                (),
    M_flags                       (),
    M_couplingIndex               (),
    M_LocalCouplingVariables      (),
    M_LocalCouplingResiduals      (),
    M_LocalDeltaCouplingVariables (),
    M_dataPhysics                 (),
    M_dataTime                    (),
    M_comm                        (),
    M_displayer                   ()
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::MS_PhysicalCoupling() \n";
#endif

    M_ID = M_couplingsNumber++;
}

MS_PhysicalCoupling::MS_PhysicalCoupling( const MS_PhysicalCoupling& coupling ) :
    M_ID                          ( coupling.M_ID ),
    M_type                        ( coupling.M_type ),
    M_dataFile                    ( coupling.M_dataFile ),
    M_models                      ( coupling.M_models ),
    M_couplingName                ( coupling.M_couplingName ),
    M_flags                       ( coupling.M_flags ),
    M_couplingIndex               ( coupling.M_couplingIndex ),
    M_LocalCouplingVariables      ( coupling.M_LocalCouplingVariables ),
    M_LocalCouplingResiduals      ( coupling.M_LocalCouplingResiduals ),
    M_LocalDeltaCouplingVariables ( coupling.M_LocalDeltaCouplingVariables ),
    M_dataPhysics                 ( coupling.M_dataPhysics ),
    M_dataTime                    ( coupling.M_dataTime ),
    M_comm                        ( coupling.M_comm ),
    M_displayer                   ( coupling.M_displayer )

{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::MS_PhysicalCoupling( coupling ) \n";
#endif

    M_ID = M_couplingsNumber++;
}

// ===================================================
// Operators
// ===================================================
MS_PhysicalCoupling&
MS_PhysicalCoupling::operator=( const MS_PhysicalCoupling& coupling )
{
    if ( this != &coupling )
    {
        M_ID                          = coupling.M_ID;
        M_type                        = coupling.M_type;
        M_dataFile                    = coupling.M_dataFile;
        M_models                      = coupling.M_models;
        M_couplingName                = coupling.M_couplingName;
        M_flags                       = coupling.M_flags;
        M_couplingIndex               = coupling.M_couplingIndex;
        M_LocalCouplingVariables      = coupling.M_LocalCouplingVariables;
        M_LocalCouplingResiduals      = coupling.M_LocalCouplingResiduals;
        M_LocalDeltaCouplingVariables = coupling.M_LocalDeltaCouplingVariables;
        M_dataPhysics                 = coupling.M_dataPhysics;
        M_dataTime                    = coupling.M_dataTime;
        M_comm                        = coupling.M_comm;
        M_displayer                   = coupling.M_displayer;
    }
    return *this;
}

// ===================================================
// MultiScale PhysicalCoupling Virtual Methods
// ===================================================
void
MS_PhysicalCoupling::ShowMe()
{
    std::cout << "Coupling id         = " << M_ID << std::endl
              << "Coupling name       = " << M_couplingName << std::endl
              << "Coupling type       = " << Enum2String( M_type, couplingsMap ) << std::endl << std::endl;

    std::cout << "Models number       = " << GetModelsNumber() << std::endl;
    std::cout << "Models type(s)      = ";
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        std::cout << Enum2String( M_models[i]->GetType(), modelsMap ) << " ";
    std::cout << std::endl;
    std::cout << "Flags list          = ";
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        std::cout << M_flags[i] << " ";
    std::cout << std::endl << std::endl;
}

// ===================================================
// Methods
// ===================================================
void
MS_PhysicalCoupling::CreateCouplingMap( EpetraMap& couplingMap )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::CreateCouplingMap( couplingMap ) \n";
#endif

    M_couplingIndex.second = couplingMap.getMap( Unique )->NumGlobalElements();

    couplingMap += M_LocalCouplingVariables->getMap();
}

void
MS_PhysicalCoupling::ImportCouplingVariables( const VectorType& CouplingVariables )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ImportCouplingVariables( CouplingVariables ) \n";
#endif

    ImportCouplingVector( CouplingVariables, *M_LocalCouplingVariables );
}

void
MS_PhysicalCoupling::ExportCouplingVariables( VectorType& CouplingVariables )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ExportCouplingVariables( CouplingVariables ) \n";
#endif

    ExportCouplingVector( *M_LocalCouplingVariables, CouplingVariables );
}

void
MS_PhysicalCoupling::ExportJacobian( MatrixType& Jacobian )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ExportJacobian( Jacobian ) \n";
#endif

    // Definitions
    bool SolveLinearSystem;                  // Flag to avoid multiple solution of the same linear system
    ModelsVector_Type ListOfPerturbedModels; // List of perturbed model

    // Insert constant values in the Jacobian (due to this coupling condition)
    InsertJacobianConstantCoefficients( Jacobian );

    // Loop on all the local coupling variables that should be perturbed
    for ( UInt i(0), column(M_couplingIndex.second) ; i < M_couplingIndex.first ; ++i, ++column )
    {
        ( *M_LocalDeltaCouplingVariables )[i] = 1.;

        // Build the list of models affected by the perturbation
        ListOfPerturbedModels = GetListOfPerturbedModels( i );

        // Loop on all the models, that could be influenced by the perturbation of the coupling variable
        for ( ModelsVector_Iterator j = ListOfPerturbedModels.begin() ; j < ListOfPerturbedModels.end() ; ++j )
        {
            SolveLinearSystem = true;

            // Loop on all the couplings (boundary flags) that connect the j-model
            for ( UInt k(0) ; k < ( *j )->GetCouplingsNumber() ; ++k )
                ( *j )->GetCoupling( k )->InsertJacobianDeltaCoefficients( Jacobian, column, ( *j )->GetID(), SolveLinearSystem );
        }

        ( *M_LocalDeltaCouplingVariables )[i] = 0.;
    }
}

void
MS_PhysicalCoupling::SaveSolution()
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SaveSolution() \n";
#endif

    std::ofstream output;
    output << std::scientific << std::setprecision( 6 );

    if ( M_comm->MyPID() == 0 )
    {
        std::string filename = MS_ProblemName + "Coupling_" + number2string( M_ID ) + ".mfile";

        if ( M_dataTime->isFirstTimeStep() )
        {
            output.open( filename.c_str(), std::ios::trunc );
            output << "% Coupling Type: " << Enum2String( M_type, couplingsMap ) << std::endl << std::endl;
            output << "% TIME            ID   FLAG FLUX            STRESS           S. PRESSURE     D. PRESSURE" << std::endl;
        }
        else
            output.open( filename.c_str(), std::ios::app );
    }

    DisplayCouplingValues( output );

    if ( M_comm->MyPID() == 0 )
        output.close();
}

void
MS_PhysicalCoupling::ClearModelsList()
{
    M_models.clear();
}

// ===================================================
// Set Methods
// ===================================================
void
MS_PhysicalCoupling::SetID( const UInt& id )
{
    M_ID = id;
}

void
MS_PhysicalCoupling::SetDataFile( const std::string& dataFile )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SetDataFile( dataFile ) \n";
#endif

    M_dataFile = GetPot( dataFile );

    // Read multiscale parameters
    M_couplingName = M_dataFile( "MultiScale/couplingName", "couplingName" );

    UInt componentSize = M_dataFile.vector_variable_size( "MultiScale/couplingFlags" );

    M_flags.reserve( componentSize );
    for ( UInt j( 0 ); j < componentSize; ++j )
        M_flags.push_back( M_dataFile( "MultiScale/couplingFlags", 0, j ) ); // flags
}

void
MS_PhysicalCoupling::AddModel( const Model_ptrType& model )
{
    M_models.push_back( model );
}

void
MS_PhysicalCoupling::SetData( const boost::shared_ptr< MS_PhysicalData >& dataPhysics,
                              const boost::shared_ptr< DataTime >& dataTime )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SetData( dataPhysics, dataTime ) \n";
#endif

    M_dataPhysics = dataPhysics;
    M_dataTime    = dataTime;
}

void
MS_PhysicalCoupling::AddFlag( const BCFlag& flag )
{
    M_flags.push_back( flag );
}

void
MS_PhysicalCoupling::AddFlagID( const UInt& flagID )
{
    M_flags.push_back( M_models.back()->GetFlag( flagID ) );
}

void
MS_PhysicalCoupling::SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm )
{

#ifdef DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm.get() ) );
}

// ===================================================
// Get Methods
// ===================================================
const UInt&
MS_PhysicalCoupling::GetID() const
{
    return M_ID;
}

const couplingsTypes&
MS_PhysicalCoupling::GetType() const
{
    return M_type;
}

const std::string&
MS_PhysicalCoupling::GetCouplingName() const
{
    return M_couplingName;
}

UInt
MS_PhysicalCoupling::GetModelsNumber() const
{
    return static_cast< UInt > ( M_models.size() );
}

UInt
MS_PhysicalCoupling::GetModelLocalID( const UInt& ID ) const
{
    for ( UInt localID( 0 ); localID < GetModelsNumber(); ++localID )
        if ( M_models[localID]->GetID() == ID )
            return localID;

    return -1;
}

Model_ptrType
MS_PhysicalCoupling::GetModel( const UInt& LocalID ) const
{
    return M_models[LocalID];
}

const UInt&
MS_PhysicalCoupling::GetCouplingVariablesNumber() const
{
    return M_couplingIndex.first;
}

// ===================================================
// Protected Methods
// ===================================================
void
MS_PhysicalCoupling::CreateLocalVectors()
{
    // Build a repeated list of GlobalElements
    std::vector<int> MyGlobalElements( M_couplingIndex.first );
    for ( UInt i = 0 ; i < static_cast< UInt > ( MyGlobalElements.size() ) ; ++i )
        MyGlobalElements[i] = i;

    // Build a repeated map for the couplings
    EpetraMap map( -1, static_cast< int > ( MyGlobalElements.size() ), &MyGlobalElements[0], 0, *M_comm );

    // Create local repeated vectors
    M_LocalCouplingVariables.reset     ( new EpetraVector( map, Repeated ) );
    M_LocalCouplingResiduals.reset     ( new EpetraVector( map, Repeated ) );
    M_LocalDeltaCouplingVariables.reset( new EpetraVector( map, Repeated ) );
}

void
MS_PhysicalCoupling::ImportCouplingVector( const VectorType& globalVector, VectorType& localVector )
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
MS_PhysicalCoupling::ExportCouplingVector( const VectorType& localVector, VectorType& globalVector )
{
    for ( UInt i(0) ; i < M_couplingIndex.first ; ++i )
        if ( M_comm->MyPID() == 0 )
            globalVector[ M_couplingIndex.second + i ] = localVector[i];
}

void
MS_PhysicalCoupling::switchErrorMessage( const Model_ptrType& model )
{
    MS_ErrorCheck( MS_ModelType, "Invalid model type ["  + Enum2String( model->GetType(), modelsMap ) +
                                 "] for coupling type [" + Enum2String( M_type, couplingsMap ) +"]\n" );
}

} // Namespace LifeV
