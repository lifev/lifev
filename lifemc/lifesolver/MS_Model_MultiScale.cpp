/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-03-12

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Model_MultiScale.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-03-12
 */

#include <lifemc/lifesolver/MS_Model_MultiScale.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
MS_Model_MultiScale::MS_Model_MultiScale() :
    M_models            (),
    M_couplings         (),
    M_modelsNumber      (),
    M_couplingsNumber   ()
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::MS_Model_MultiScale() \n";
#endif

    M_type = MultiScale;
}

MS_Model_MultiScale::MS_Model_MultiScale( const MS_Model_MultiScale& multiscale ) :
    super               ( multiscale ),
    M_models            ( multiscale.M_models ),
    M_couplings         ( multiscale.M_couplings ),
    M_modelsNumber      ( multiscale.M_modelsNumber ),
    M_couplingsNumber   ( multiscale.M_couplingsNumber )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::MS_Model_MultiScale( multiscale ) \n";
#endif

}

// ===================================================
//! Methods
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
        M_models          = multiscale.M_models;
        M_couplings       = multiscale.M_couplings;
        M_modelsNumber    = multiscale.M_modelsNumber;
        M_couplingsNumber = multiscale.M_couplingsNumber;
    }
    return *this;
}

void
MS_Model_MultiScale::SetupData( void )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SetupData() \n";
#endif

    // Load
    loadModels();
    loadCouplings();
    loadGeometry();

    // SetData & SetCommunicator
    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->SetData( M_dataPhysics, M_dataTime );

    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->SetCommunicator( M_comm );

    for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
        M_couplings[i]->SetCommunicator( M_comm );

    // SetupData
    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->SetupData();

    for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
        M_couplings[i]->SetupData();
}

void
MS_Model_MultiScale::SetupModel( void )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SetupModel() \n";
#endif

    for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
        M_couplings[i]->SetupCoupling();

    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->SetupModel();
}

void
MS_Model_MultiScale::SetupImplicitCoupling( ContainerOfVectors< EpetraVector >& couplingVariables, ContainerOfVectors< EpetraVector >& couplingResiduals )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SetupImplicitCoupling( couplingVariables, couplingResiduals ) \n";
#endif

    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->SetupImplicitCoupling( couplingVariables, couplingResiduals );

    for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
        M_couplings[i]->SetupImplicitCoupling( couplingVariables, couplingResiduals );
}

void
MS_Model_MultiScale::BuildSystem( void )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::BuildSystem() \n";
#endif

    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->BuildSystem();

    for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
        M_couplings[i]->UpdateCouplingVariables();
}

void
MS_Model_MultiScale::UpdateSystem( void )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::UpdateSystem() \n";
#endif

    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->UpdateSystem();

    for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
        M_couplings[i]->UpdateCouplingVariables();
}

void
MS_Model_MultiScale::SolveSystem( void )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SolveSystem() \n";
#endif

    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->SolveSystem();

    for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
        M_couplings[i]->UpdateCouplingResiduals();
}

void
MS_Model_MultiScale::SaveSolution( void )
{

#ifdef DEBUG
    Debug( 8110 ) << "MS_Model_MultiScale::SaveSolution() \n";
#endif

    for ( UInt i( 0 ); i < M_modelsNumber; ++i )
        M_models[i]->SaveSolution();
}

void
MS_Model_MultiScale::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Models number      = " << M_modelsNumber << std::endl
                  << "Couplings number   = " << M_couplingsNumber << std::endl << std::endl;

        std::cout << "==================== Models Information =====================" << std::endl << std::endl;

        for ( UInt i( 0 ); i < M_modelsNumber; ++i )
            M_models[i]->ShowMe();

        std::cout << "=================== Couplings Information ===================" << std::endl << std::endl;

        for ( UInt i( 0 ); i < M_couplingsNumber; ++i )
            M_couplings[i]->ShowMe();
    }
}

// ===================================================
//! Private Methods
// ===================================================
inline void
MS_Model_MultiScale::loadModels()
{
    UInt id;
    modelsTypes model;

    UInt columnNumber = 3.0;
    UInt linesNumber = M_dataFile.vector_variable_size( "Problem/models" ) / columnNumber;

    std::string path = M_dataFile( "Problem/modelsPath", "./" );
    for ( UInt i( 0 ); i < linesNumber; ++i )
    {
        id    =           M_dataFile( "Problem/models", 0, i * columnNumber );
        model = modelsMap[M_dataFile( "Problem/models", "undefined", i * columnNumber + 1 )];

        M_models[id] = PhysicalModel_ptr( FactoryModels::instance().createObject( model ) );
        M_models[id]->SetDataFile( path + M_dataFile( "Problem/models", "undefined", i * columnNumber + 2 ) + ".dat" );
    }

    M_modelsNumber = static_cast< UInt > ( M_models.size() );
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
    for ( UInt i( 0 ); i < linesNumber; ++i )
    {
        id       =              M_dataFile( "Problem/couplings", 0, i * columnNumber );
        coupling = couplingsMap[M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 1 )];

        M_couplings[id] = PhysicalCoupling_ptr( FactoryCouplings::instance().createObject( coupling ) );
        M_couplings[id]->SetDataFile( path + M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 2 ) + ".dat" );

        modelsIDVector = string2numVect< UInt > ( M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 3 ) );
        flagsIDVector  = string2numVect< UInt > ( M_dataFile( "Problem/couplings", "undefined", i * columnNumber + 4 ) );
        for ( UInt j( 0 ); j < static_cast< UInt > ( modelsIDVector.size() ); ++j )
        {
            M_couplings[id]->AddModel( M_models[modelsIDVector[j]] );
            M_couplings[id]->AddFlagID( flagsIDVector[j] );
        }
    }

    M_couplingsNumber = static_cast< UInt > ( M_couplings.size() );
}

inline void
MS_Model_MultiScale::loadGeometry()
{
    UInt id;

    boost::array< Real, nDimensions > geometryScale, geometryRotate, geometryTranslate;

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

        M_models[id]->SetGeometry( geometryScale, geometryRotate, geometryTranslate );
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
