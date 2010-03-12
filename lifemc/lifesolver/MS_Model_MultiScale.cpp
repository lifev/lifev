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
 *  @brief MultiScale Model MultiScale
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-03-2009
 */


#include <lifemc/lifesolver/MS_Model_MultiScale.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_MultiScale::MS_Model_MultiScale() :
    M_modelsList        (),
    M_couplingsList     ()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::MS_Model_MultiScale() \n";
#endif

    M_type = MultiScale;
}

MS_Model_MultiScale::MS_Model_MultiScale( const MS_Model_MultiScale& multiscale ) :
    super               ( multiscale ),
    M_modelsList        ( multiscale.M_modelsList ),
    M_couplingsList     ( multiscale.M_couplingsList )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::MS_Model_MultiScale( multiscale ) \n";
#endif

}

MS_Model_MultiScale::~MS_Model_MultiScale()
{

    #ifdef DEBUG
        Debug( 8110 ) << "MS_Model_MultiScale::~MS_Model_MultiScale( ) \n";
    #endif

        // Disconnect models and couplings to allow their destruction
        for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
            ( *i )->ClearCouplingsList();

        for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
            ( *i )->ClearModelsList();
}

// ===================================================
// Operators
// ===================================================
MS_Model_MultiScale&
MS_Model_MultiScale::operator=( const MS_Model_MultiScale& multiscale )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::operator=( multiscale ) \n";
#endif

    if ( this != &multiscale )
    {
        super::operator=( multiscale );
        M_modelsList      = multiscale.M_modelsList;
        M_couplingsList   = multiscale.M_couplingsList;
    }
    return *this;
}

// ===================================================
// MultiScale Physical Model
// ===================================================
void
MS_Model_MultiScale::SetupData()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SetupData() \n";
#endif

    // Load
    loadModels();
    loadCouplings();
    loadGeometry();

    // SetData & SetCommunicator
    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->SetData( M_dataPhysics, M_dataTime );

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->SetData( M_dataPhysics, M_dataTime );

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->SetCommunicator( M_comm );

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->SetCommunicator( M_comm );

    // SetupData
    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->SetupData();

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->SetupData();
}

void
MS_Model_MultiScale::SetupModel()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SetupModel() \n";
#endif

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->SetupCoupling();

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->SetupModel();
}

void
MS_Model_MultiScale::BuildSystem()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::BuildSystem() \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->BuildSystem();

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->InitializeCouplingVariables();
}

void
MS_Model_MultiScale::UpdateSystem()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::UpdateSystem() \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->UpdateSystem();
}

void
MS_Model_MultiScale::SolveSystem()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SolveSystem() \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->SolveSystem();
}

void
MS_Model_MultiScale::SaveSolution()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SaveSolution() \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->SaveSolution();

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->SaveSolution();
}

void
MS_Model_MultiScale::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();
        std::cout << "Models number       = " << M_modelsList.size() << std::endl
                  << "Couplings number    = " << M_couplingsList.size() << std::endl << std::endl;

        std::cout << "==================== Models Information =====================" << std::endl << std::endl;
    }

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        ( *i )->ShowMe();

    if ( M_displayer->isLeader() )
        std::cout << "=================== Couplings Information ===================" << std::endl << std::endl;

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->ShowMe();

    //MPI Barrier
    M_comm->Barrier();
}

// ===================================================
// Methods
// ===================================================
void
MS_Model_MultiScale::CreateCouplingMap( EpetraMap& couplingMap )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::CreateCouplingMap( couplingMap ) \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->CreateCouplingMap( couplingMap );

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->CreateCouplingMap( couplingMap );
}

void
MS_Model_MultiScale::InitializeCouplingVariables()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::InitializeCouplingVariables() \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->InitializeCouplingVariables();

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->InitializeCouplingVariables();
}

void
MS_Model_MultiScale::ImportCouplingVariables( const VectorType& CouplingVariables )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ImportCouplingVariables( CouplingVariables ) \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ImportCouplingVariables( CouplingVariables );

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->ImportCouplingVariables( CouplingVariables );
}

void
MS_Model_MultiScale::ExportCouplingVariables( VectorType& CouplingVariables )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ExportCouplingVariables( CouplingVariables ) \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ExportCouplingVariables( CouplingVariables);

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->ExportCouplingVariables( CouplingVariables );
}

void
MS_Model_MultiScale::ExportCouplingResiduals( VectorType& CouplingResiduals )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ExportCouplingResiduals( CouplingResiduals ) \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ExportCouplingResiduals( CouplingResiduals );

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->ExportCouplingResiduals( CouplingResiduals );
}

void
MS_Model_MultiScale::ExportJacobian( MatrixType& Jacobian )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ExportJacobian() \n";
#endif

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ExportJacobian( Jacobian );

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        ( *i )->ExportJacobian( Jacobian );
}

// ===================================================
// Get Methods
// ===================================================
UInt
MS_Model_MultiScale::GetCouplingVariablesNumber()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::GetCouplingVariablesNumber() \n";
#endif

    UInt CouplingVariablesNumber = 0;

    for ( ModelsVector_ConstIterator i = M_modelsList.begin(); i < M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            CouplingVariablesNumber += ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->GetCouplingVariablesNumber();

    for ( CouplingsVector_ConstIterator i = M_couplingsList.begin(); i < M_couplingsList.end(); ++i )
        CouplingVariablesNumber += ( *i )->GetCouplingVariablesNumber();

    return CouplingVariablesNumber;
}

// ===================================================
// Private Methods
// ===================================================
inline void
MS_Model_MultiScale::loadModels()
{
    UInt id;
    modelsTypes model;

    UInt columnNumber = 3.0;
    UInt linesNumber = M_dataFile.vector_variable_size( "Problem/models" ) / columnNumber;

    std::string path = M_dataFile( "Problem/modelsPath", "./" );
    M_modelsList.resize( linesNumber );
    for ( UInt i( 0 ); i < linesNumber; ++i )
    {
        id    =           M_dataFile( "Problem/models", 0, i * columnNumber );
        model = modelsMap[M_dataFile( "Problem/models", "undefined", i * columnNumber + 1 )];

        M_modelsList[id] = Model_ptrType( FactoryModels::instance().createObject( model ) );
        M_modelsList[id]->SetDataFile( path + M_dataFile( "Problem/models", "undefined", i * columnNumber + 2 ) + ".dat" );
    }
}

inline void
MS_Model_MultiScale::loadCouplings()
{
    UInt id;
    couplingsTypes coupling;

    std::vector< UInt > modelsIDVector;
    std::vector< UInt > flagsIDVector;

    UInt columnNumber = 5.0;
    UInt linesNumber = M_dataFile.vector_variable_size( "Problem/couplings" ) / columnNumber;

    std::string path = M_dataFile( "Problem/couplingsPath", "./" );
    M_couplingsList.resize( linesNumber );
    for ( UInt i( 0 ); i < linesNumber; ++i )
    {
        id       =              M_dataFile( "Problem/couplings", 0, i * columnNumber );
        coupling = couplingsMap[M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 1 )];

        M_couplingsList[id] = Coupling_ptrType( FactoryCouplings::instance().createObject( coupling ) );
        M_couplingsList[id]->SetDataFile( path + M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 2 ) + ".dat" );

        modelsIDVector = string2numVect< UInt > ( M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 3 ) );
        flagsIDVector  = string2numVect< UInt > ( M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 4 ) );
        for ( UInt j( 0 ); j < static_cast< UInt > ( modelsIDVector.size() ); ++j )
        {
            M_couplingsList[id]->AddModel( M_modelsList[modelsIDVector[j]] );
            M_couplingsList[id]->AddFlagID( flagsIDVector[j] );
            M_modelsList[modelsIDVector[j]]->AddCoupling( M_couplingsList[id] );
        }
    }
}

inline void
MS_Model_MultiScale::loadGeometry()
{
    UInt id;

    boost::array< Real, NDIM > geometryScale, geometryRotate, geometryTranslate;

    UInt columnNumber = 10.0;
    UInt linesNumber = M_dataFile.vector_variable_size( "Problem/offset" ) / columnNumber;

    for ( UInt i( 0 ); i < linesNumber; ++i )
    {
        id = M_dataFile( "Problem/offset", 0, i * columnNumber );

        geometryScale[0]     = M_geometryScale[0]     * M_dataFile( "Problem/offset", 1., i * columnNumber + 1 );
        geometryScale[1]     = M_geometryScale[1]     * M_dataFile( "Problem/offset", 1., i * columnNumber + 2 );
        geometryScale[2]     = M_geometryScale[2]     * M_dataFile( "Problem/offset", 1., i * columnNumber + 3 );

        geometryRotate[0]    = M_geometryRotate[0]    + M_dataFile( "Problem/offset", 0., i * columnNumber + 4 ) * Pi / 180;
        geometryRotate[1]    = M_geometryRotate[1]    + M_dataFile( "Problem/offset", 0., i * columnNumber + 5 ) * Pi / 180;
        geometryRotate[2]    = M_geometryRotate[2]    + M_dataFile( "Problem/offset", 0., i * columnNumber + 6 ) * Pi / 180;

        geometryTranslate[0] = M_geometryTranslate[0] + M_dataFile( "Problem/offset", 0., i * columnNumber + 7 );
        geometryTranslate[1] = M_geometryTranslate[1] + M_dataFile( "Problem/offset", 0., i * columnNumber + 8 );
        geometryTranslate[2] = M_geometryTranslate[2] + M_dataFile( "Problem/offset", 0., i * columnNumber + 9 );

        M_modelsList[id]->SetGeometry( geometryScale, geometryRotate, geometryTranslate );
    }
}

template< typename number >
inline std::vector< number >
MS_Model_MultiScale::string2numVect( const std::string& string )
{
    //Split the string
    std::vector< std::string > stringVector;
    boost::split( stringVector, string, boost::is_any_of( "," ) );

    //Convert to the right type
    std::vector< number > numberVector;
    for ( UInt i( 0 ); i < static_cast< UInt > ( stringVector.size() ); ++i )
        numberVector.push_back( static_cast< number > ( std::atoi( stringVector[i].c_str() ) ) );

    return numberVector;
}

} // Namespace LifeV
