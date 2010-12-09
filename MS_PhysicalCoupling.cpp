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

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>

namespace LifeV
{

std::map< std::string, couplingsTypes > MS_couplingsMap;

UInt MS_PhysicalCoupling::M_couplingsNumber = 0;

// ===================================================
// Constructors & Destructor
// ===================================================
MS_PhysicalCoupling::MS_PhysicalCoupling() :
        M_ID                          (),
        M_type                        (),
        M_models                      (),
        M_couplingName                (),
        M_flags                       (),
        M_globalData                  (),
        M_couplingIndex               (),
        M_LocalCouplingVariables      (),
        M_LocalCouplingResiduals      (),
        M_timeInterpolationOrder      ( 1 ),
        M_perturbedCoupling           ( false ),
        M_comm                        (),
        M_displayer                   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::MS_PhysicalCoupling() \n";
#endif

    M_ID = M_couplingsNumber++;
}

// ===================================================
// MultiScale PhysicalCoupling Virtual Methods
// ===================================================
void
MS_PhysicalCoupling::SetupData( const std::string& FileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SetupData( FileName ) \n";
#endif

    GetPot DataFile( FileName );

    // Read multiscale parameters
    M_couplingName = DataFile( "MultiScale/couplingName", "couplingName" );
    M_timeInterpolationOrder = DataFile( "MultiScale/timeInterpolationOrder", 0 );

    // Set the size of the local coupling variables
    M_LocalCouplingVariables.reserve( M_timeInterpolationOrder + 1 );
}

void
MS_PhysicalCoupling::ShowMe()
{
    std::cout << "Coupling id         = " << M_ID << std::endl
              << "Coupling name       = " << M_couplingName << std::endl
              << "Coupling type       = " << Enum2String( M_type, MS_couplingsMap ) << std::endl << std::endl;

    std::cout << "Models number       = " << GetModelsNumber() << std::endl;
    std::cout << "Models ID(s)        = ";
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        std::cout << M_models[i]->GetID() << " ";
    std::cout << std::endl;
    std::cout << "Models type(s)      = ";
    for ( UInt i( 0 ); i < GetModelsNumber(); ++i )
        std::cout << Enum2String( M_models[i]->GetType(), MS_modelsMap ) << " ";
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

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::CreateCouplingMap( couplingMap ) \n";
#endif

    M_couplingIndex.second = couplingMap.getMap( Unique )->NumGlobalElements();

    couplingMap += M_LocalCouplingVariables[0]->getMap();
}

void
MS_PhysicalCoupling::ImportCouplingVariables( const MS_Vector_Type& CouplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ImportCouplingVariables( CouplingVariables ) \n";
#endif

    ImportCouplingVector( CouplingVariables, *M_LocalCouplingVariables[0] );
}

void
MS_PhysicalCoupling::ExportCouplingVariables( MS_Vector_Type& CouplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ExportCouplingVariables( CouplingVariables ) \n";
#endif

    ExportCouplingVector( *M_LocalCouplingVariables[0], CouplingVariables );
}

void
MS_PhysicalCoupling::ExtrapolateCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ExtrapolateCouplingVariables() \n";
#endif

    MS_Vector_Type ExtrapolatedCouplingVariables( *M_LocalCouplingVariables[0] );
    UInt couplingVariablesSize( M_LocalCouplingVariables.size() );

    // Time container for interpolation
    TimeContainer_Type timeContainer( couplingVariablesSize, 0 );
    for ( UInt i(0) ; i < couplingVariablesSize ; ++i )
        timeContainer[i] = M_globalData->GetDataTime()->getTime() - i * M_globalData->GetDataTime()->getTimeStep();

    // Interpolate the coupling variables at the next time
    InterpolateCouplingVariables( timeContainer, M_globalData->GetDataTime()->getNextTime(), ExtrapolatedCouplingVariables );

    // If we have not yet enough samples for interpolation, we add a new one
    if ( couplingVariablesSize <= M_timeInterpolationOrder )
    {
        ++couplingVariablesSize;
        M_LocalCouplingVariables.push_back( MS_Vector_PtrType ( new EpetraVector( *M_LocalCouplingVariables[0] ) ) );
    }

    // Updating database
    for ( UInt i(1) ; i < couplingVariablesSize ; ++i )
        *M_LocalCouplingVariables[couplingVariablesSize-i] = *M_LocalCouplingVariables[couplingVariablesSize-i-1];

    *M_LocalCouplingVariables[0] = ExtrapolatedCouplingVariables;

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8200 ) << "C(" << M_couplingIndex.second + i << ") = " << ( *M_LocalCouplingVariables[0] )[i]  << "\n";
#endif

}

bool
MS_PhysicalCoupling::IsPerturbed() const
{
    return M_perturbedCoupling == -1 ? false : true;
}

void
MS_PhysicalCoupling::ExportJacobian( MS_Matrix_Type& Jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::ExportJacobian( Jacobian ) \n";
#endif

    // Definitions
    bool SolveLinearSystem;                     // Flag to avoid multiple solution of the same linear system
    MS_ModelsVector_Type ListOfPerturbedModels; // List of perturbed model

    // Insert constant values in the Jacobian (due to this coupling condition)
    InsertJacobianConstantCoefficients( Jacobian );

    // Set as perturbed
    M_perturbedCoupling = 0;

    // Loop on all the local coupling variables that should be perturbed
    for ( UInt column(M_couplingIndex.second) ; M_perturbedCoupling < static_cast< Int > ( M_couplingIndex.first ); ++M_perturbedCoupling, ++column )
    {
        // Build the list of models affected by the perturbation
        ListOfPerturbedModels = GetListOfPerturbedModels( M_perturbedCoupling );

        // Loop on all the models, that could be influenced by the perturbation of the coupling variable
        for ( MS_ModelsVector_Iterator j = ListOfPerturbedModels.begin() ; j < ListOfPerturbedModels.end() ; ++j )
        {
            SolveLinearSystem = true;

            // Loop on all the couplings (boundary flags) that connect the j-model
            for ( UInt k(0) ; k < ( *j )->GetCouplingsNumber() ; ++k )
                ( *j )->GetCoupling( k )->InsertJacobianDeltaCoefficients( Jacobian, column, ( *j )->GetID(), SolveLinearSystem );
        }
    }

    // Set as unperturbed
    M_perturbedCoupling = -1;
}

void
MS_PhysicalCoupling::SaveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SaveSolution() \n";
#endif

    std::ofstream output;
    output << std::scientific << std::setprecision( 15 );

    if ( M_comm->MyPID() == 0 )
    {
        std::string filename = MS_ProblemFolder + "Step_" + number2string( MS_ProblemStep ) + "_Coupling_" + number2string( M_ID ) + ".mfile";

        if ( M_globalData->GetDataTime()->isFirstTimeStep() )
        {
            output.open( filename.c_str(), std::ios::trunc );
            output << "% Coupling Type: " << Enum2String( M_type, MS_couplingsMap ) << std::endl << std::endl;
            output << "% TIME                     ID   FLAG FLUX                     STRESS                    S. PRESSURE              D. PRESSURE" << std::endl;
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
MS_PhysicalCoupling::AddModel( const MS_Model_PtrType& model )
{
    M_models.push_back( model );
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
MS_PhysicalCoupling::SetGlobalData( const MS_GlobalDataContainer_PtrType& globalData )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SetupData( globalData ) \n";
#endif

    M_globalData = globalData;
}

void
MS_PhysicalCoupling::SetCommunicator( const MS_Comm_PtrType& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8200 ) << "MS_PhysicalCoupling::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
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

MS_Model_PtrType
MS_PhysicalCoupling::GetModel( const UInt& LocalID ) const
{
    return M_models[LocalID];
}

const BCFlag&
MS_PhysicalCoupling::GetFlag( const UInt& LocalID ) const
{
    return M_flags[LocalID];
}

const UInt&
MS_PhysicalCoupling::GetCouplingVariablesNumber() const
{
    return M_couplingIndex.first;
}

const Int&
MS_PhysicalCoupling::GetPerturbedCoupling() const
{
    return M_perturbedCoupling;
}

const MS_Vector_Type&
MS_PhysicalCoupling::GetResidual() const
{
    return *M_LocalCouplingResiduals;
}

const UInt&
MS_PhysicalCoupling::GetTimeInterpolationOrder() const
{
    return M_timeInterpolationOrder;
}

// ===================================================
// Protected Methods
// ===================================================
void
MS_PhysicalCoupling::CreateLocalVectors()
{
    // Build a repeated list of GlobalElements
    std::vector<Int> MyGlobalElements( M_couplingIndex.first );
    for ( UInt i = 0 ; i < static_cast< UInt > ( MyGlobalElements.size() ) ; ++i )
        MyGlobalElements[i] = i;

    // Build a repeated map for the couplings
    EpetraMap map( -1, static_cast< Int > ( MyGlobalElements.size() ), &MyGlobalElements[0], 0, M_comm );

    // Create local repeated vectors
    M_LocalCouplingVariables.push_back( MS_Vector_PtrType ( new EpetraVector( map, Repeated ) ) );
    M_LocalCouplingResiduals.reset( new EpetraVector( map, Repeated ) );
}

void
MS_PhysicalCoupling::ImportCouplingVector( const MS_Vector_Type& globalVector, MS_Vector_Type& localVector )
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
MS_PhysicalCoupling::ExportCouplingVector( const MS_Vector_Type& localVector, MS_Vector_Type& globalVector )
{
    for ( UInt i(0) ; i < M_couplingIndex.first ; ++i )
        if ( M_comm->MyPID() == 0 )
            globalVector[ M_couplingIndex.second + i ] = localVector[i];
}

void
MS_PhysicalCoupling::InterpolateCouplingVariables( const TimeContainer_Type& timeContainer,
                                                   const Real& t,
                                                   MS_Vector_Type& interpolatedCouplingVariables )
{
    // Lagrange interpolation
    interpolatedCouplingVariables *= 0;
    Real base(1);

    for ( UInt i(0) ; i < M_LocalCouplingVariables.size() ; ++i )
    {
        base = 1;
        for ( UInt j(0) ; j < M_LocalCouplingVariables.size() ; ++j )
            if ( j != i )
                base *= (t - timeContainer[j]) / (timeContainer[i] - timeContainer[j]);

        interpolatedCouplingVariables += *M_LocalCouplingVariables[i] * base;
    }
}

void
MS_PhysicalCoupling::switchErrorMessage( const MS_Model_PtrType& model )
{
    MS_ErrorCheck( MS_ModelType, "Invalid model type ["  + Enum2String( model->GetType(), MS_modelsMap ) +
                   "] for coupling type [" + Enum2String( M_type, MS_couplingsMap ) +"]\n" );
}

} // Namespace LifeV
