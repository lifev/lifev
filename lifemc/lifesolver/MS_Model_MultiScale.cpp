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

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Model_MultiScale::MS_Model_MultiScale() :
        M_modelsList        (),
        M_couplingsList     ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::MS_Model_MultiScale() \n";
#endif

    M_type = MultiScale;
}

MS_Model_MultiScale::MS_Model_MultiScale( const MS_Model_MultiScale& multiscale ) :
        super               ( multiscale ),
        M_modelsList        ( multiscale.M_modelsList ),
        M_couplingsList     ( multiscale.M_couplingsList )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::MS_Model_MultiScale( multiscale ) \n";
#endif

}

MS_Model_MultiScale::~MS_Model_MultiScale()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::~MS_Model_MultiScale( ) \n";
#endif

    // Disconnect models and couplings to allow their destruction
    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->ClearCouplingsList();

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->ClearModelsList();
}

// ===================================================
// Operators
// ===================================================
MS_Model_MultiScale&
MS_Model_MultiScale::operator=( const MS_Model_MultiScale& multiscale )
{

#ifdef HAVE_LIFEV_DEBUG
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
MS_Model_MultiScale::SetupData( const std::string& FileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SetupData( FileName ) \n";
#endif

    super::SetupData( FileName );

    // Useful variables
    UInt id;

    modelsTypes model;
    couplingsTypes coupling;
    boost::array< Real, NDIM > geometryScale, geometryRotate, geometryTranslate;

    std::map< UInt, UInt > modelsIDMap;

    std::vector< UInt > modelsIDVector;
    std::vector< UInt > flagsIDVector;

    UInt modelsColumnsNumber    = 3.0;
    UInt couplingsColumnsNumber = 5.0;
    UInt geometryColumnsNumber  = 10.0;

    GetPot DataFile( FileName );

    UInt modelsLinesNumber    = DataFile.vector_variable_size( "Problem/models" ) / modelsColumnsNumber;
    UInt couplingsLinesNumber = DataFile.vector_variable_size( "Problem/couplings" ) / couplingsColumnsNumber;
    UInt geometryLinesNumber  = DataFile.vector_variable_size( "Problem/offset" ) / geometryColumnsNumber;

    M_modelsList.resize( modelsLinesNumber );
    M_couplingsList.resize( couplingsLinesNumber );

    // Load Models
    std::string path = DataFile( "Problem/modelsPath", "./" );
    for ( UInt i( 0 ); i < modelsLinesNumber; ++i )
    {
        id    = DataFile( "Problem/models", 0, i * modelsColumnsNumber );
        modelsIDMap[id] = i;
        model = MS_modelsMap[DataFile( "Problem/models", "undefined", i * modelsColumnsNumber + 1 )];

        // Set Geometry
        geometryScale[0]     = M_geometryScale[0];
        geometryScale[1]     = M_geometryScale[1];
        geometryScale[2]     = M_geometryScale[2];

        geometryRotate[0]    = M_geometryRotate[0];
        geometryRotate[1]    = M_geometryRotate[1];
        geometryRotate[2]    = M_geometryRotate[2];

        geometryTranslate[0] = M_geometryTranslate[0];
        geometryTranslate[1] = M_geometryTranslate[1];
        geometryTranslate[2] = M_geometryTranslate[2];

        for ( UInt j( 0 ); j < geometryLinesNumber; ++j )
            if ( id == DataFile( "Problem/offset", 1., j * geometryColumnsNumber ) )
            {
                geometryScale[0]     *= DataFile( "Problem/offset", 1., j * geometryColumnsNumber + 1 );
                geometryScale[1]     *= DataFile( "Problem/offset", 1., j * geometryColumnsNumber + 2 );
                geometryScale[2]     *= DataFile( "Problem/offset", 1., j * geometryColumnsNumber + 3 );

                geometryRotate[0]    += DataFile( "Problem/offset", 0., j * geometryColumnsNumber + 4 ) * Pi / 180;
                geometryRotate[1]    += DataFile( "Problem/offset", 0., j * geometryColumnsNumber + 5 ) * Pi / 180;
                geometryRotate[2]    += DataFile( "Problem/offset", 0., j * geometryColumnsNumber + 6 ) * Pi / 180;

                geometryTranslate[0] += DataFile( "Problem/offset", 0., j * geometryColumnsNumber + 7 );
                geometryTranslate[1] += DataFile( "Problem/offset", 0., j * geometryColumnsNumber + 8 );
                geometryTranslate[2] += DataFile( "Problem/offset", 0., j * geometryColumnsNumber + 9 );
            }

        M_modelsList[i] = MS_Model_PtrType( MS_Model_Factory::instance().createObject( model ) );
        M_modelsList[i]->SetCommunicator( M_comm );
        M_modelsList[i]->SetGeometry( geometryScale, geometryRotate, geometryTranslate );
        M_modelsList[i]->SetGlobalData( M_globalData );
        M_modelsList[i]->SetupData( path + Enum2String( model, MS_modelsMap ) + "/"
                                    + DataFile( "Problem/models", "undefined", i * modelsColumnsNumber + 2 ) + ".dat" );
    }

    // Load couplings
    path = DataFile( "Problem/couplingsPath", "./" );
    for ( UInt i( 0 ); i < couplingsLinesNumber; ++i )
    {
        coupling = MS_couplingsMap[DataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 1 )];

        M_couplingsList[i] = MS_Coupling_PtrType( MS_Coupling_Factory::instance().createObject( coupling ) );
        M_couplingsList[i]->SetCommunicator( M_comm );
        M_couplingsList[i]->SetGlobalData( M_globalData );
        M_couplingsList[i]->SetupData( path + Enum2String( coupling, MS_couplingsMap ) + "/"
                                       + DataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 2 ) + ".dat" );

        string2numbersVector< UInt > ( DataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 3 ), modelsIDVector );
        string2numbersVector< UInt > ( DataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 4 ), flagsIDVector );
        for ( UInt j( 0 ); j < static_cast< UInt > ( modelsIDVector.size() ); ++j )
        {
            M_couplingsList[i]->AddModel( M_modelsList[modelsIDMap[modelsIDVector[j]]] );
            M_couplingsList[i]->AddFlagID( flagsIDVector[j] );
            M_modelsList[modelsIDMap[modelsIDVector[j]]]->AddCoupling( M_couplingsList[i] );
        }
        modelsIDVector.clear();
        flagsIDVector.clear();
    }
}

void
MS_Model_MultiScale::SetupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SetupModel() \n";
#endif

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->SetupCoupling();

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->SetupModel();
}

void
MS_Model_MultiScale::BuildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::BuildSystem() \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->BuildSystem();

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->InitializeCouplingVariables();
}

void
MS_Model_MultiScale::UpdateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::UpdateSystem() \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->UpdateSystem();
}

void
MS_Model_MultiScale::SolveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SolveSystem() \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->SolveSystem();
}

void
MS_Model_MultiScale::SaveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SaveSolution() \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->SaveSolution();

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
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

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->ShowMe();

    if ( M_displayer->isLeader() )
        std::cout << "=================== Couplings Information ===================" << std::endl << std::endl;

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->ShowMe();
}

// ===================================================
// Methods
// ===================================================
void
MS_Model_MultiScale::CreateCouplingMap( EpetraMap& couplingMap )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::CreateCouplingMap( couplingMap ) \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->CreateCouplingMap( couplingMap );

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->CreateCouplingMap( couplingMap );
}

void
MS_Model_MultiScale::InitializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::InitializeCouplingVariables() \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->InitializeCouplingVariables();

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->InitializeCouplingVariables();
}

void
MS_Model_MultiScale::ExtrapolateCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ExtrapolateCouplingVariables() \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ExtrapolateCouplingVariables();

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->ExtrapolateCouplingVariables();
}

void
MS_Model_MultiScale::ImportCouplingVariables( const MS_Vector_Type& CouplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ImportCouplingVariables( CouplingVariables ) \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ImportCouplingVariables( CouplingVariables );

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->ImportCouplingVariables( CouplingVariables );
}

void
MS_Model_MultiScale::ExportCouplingVariables( MS_Vector_Type& CouplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ExportCouplingVariables( CouplingVariables ) \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ExportCouplingVariables( CouplingVariables);

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->ExportCouplingVariables( CouplingVariables );
}

void
MS_Model_MultiScale::ExportCouplingResiduals( MS_Vector_Type& CouplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ExportCouplingResiduals( CouplingResiduals ) \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ExportCouplingResiduals( CouplingResiduals );

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->ExportCouplingResiduals( CouplingResiduals );
}

void
MS_Model_MultiScale::ExportJacobian( MS_Matrix_Type& Jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::ExportJacobian() \n";
#endif

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->ExportJacobian( Jacobian );

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->ExportJacobian( Jacobian );
}

// ===================================================
// Get Methods
// ===================================================
UInt
MS_Model_MultiScale::GetCouplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::GetCouplingVariablesNumber() \n";
#endif

    UInt CouplingVariablesNumber = 0;

    for ( MS_ModelsVector_ConstIterator i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->GetType() == MultiScale )
            CouplingVariablesNumber += ( MS_DynamicCast< MS_Model_MultiScale > ( *i ) )->GetCouplingVariablesNumber();

    for ( MS_CouplingsVector_ConstIterator i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        CouplingVariablesNumber += ( *i )->GetCouplingVariablesNumber();

    return CouplingVariablesNumber;
}

} // Namespace LifeV
