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
 *  @brief File containing the MultiScale Model MultiScale
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleModelMultiscale.hpp>

namespace LifeV
{
namespace multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelMultiscale::MultiscaleModelMultiscale() :
        multiscaleModel_Type       (),
        M_modelsList               (),
        M_couplingsList            ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::MultiscaleModelMultiscale() \n";
#endif

    M_type = MultiScale;
}

MultiscaleModelMultiscale::~MultiscaleModelMultiscale()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::~MultiscaleModelMultiscale( ) \n";
#endif

    // Disconnect models and couplings to allow their destruction
    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->clearCouplingsList();

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->clearModelsList();
}

// ===================================================
// MultiScale Physical Model
// ===================================================
void
MultiscaleModelMultiscale::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::SetupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData( fileName );

    // Useful variables
    UInt id;

    //!\todo pass a std::string to the factories
    //std::string model;
    //std::string coupling;

    models_Type model;
    couplings_Type coupling;

    boost::array< Real, NDIM > geometryScale, geometryRotate, geometryTranslate;

    std::map< UInt, UInt > modelsIDMap;

    std::vector< UInt > modelsIDVector;
    std::vector< UInt > flagsIDVector;

    UInt modelsColumnsNumber    = 3.0;
    UInt couplingsColumnsNumber = 5.0;
    UInt geometryColumnsNumber  = 10.0;

    GetPot dataFile( fileName );

    UInt modelsLinesNumber    = dataFile.vector_variable_size( "Problem/models" ) / modelsColumnsNumber;
    UInt couplingsLinesNumber = dataFile.vector_variable_size( "Problem/couplings" ) / couplingsColumnsNumber;
    UInt geometryLinesNumber  = dataFile.vector_variable_size( "Problem/offset" ) / geometryColumnsNumber;

    M_modelsList.resize( modelsLinesNumber );
    M_couplingsList.resize( couplingsLinesNumber );

    // Load Models
    std::string path = dataFile( "Problem/modelsPath", "./" );
    for ( UInt i( 0 ); i < modelsLinesNumber; ++i )
    {
        id    = dataFile( "Problem/models", 0, i * modelsColumnsNumber );
        modelsIDMap[id] = i;

        //!\todo pass a std::string to the factories
        //model = dataFile( "Problem/models", "undefined", i * modelsColumnsNumber + 1 );
        model = multiscaleModelsMap[dataFile( "Problem/models", "undefined", i * modelsColumnsNumber + 1 )];

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
            if ( id == dataFile( "Problem/offset", 1., j * geometryColumnsNumber ) )
            {
                geometryScale[0]     *= dataFile( "Problem/offset", 1., j * geometryColumnsNumber + 1 );
                geometryScale[1]     *= dataFile( "Problem/offset", 1., j * geometryColumnsNumber + 2 );
                geometryScale[2]     *= dataFile( "Problem/offset", 1., j * geometryColumnsNumber + 3 );

                geometryRotate[0]    += dataFile( "Problem/offset", 0., j * geometryColumnsNumber + 4 ) * Pi / 180;
                geometryRotate[1]    += dataFile( "Problem/offset", 0., j * geometryColumnsNumber + 5 ) * Pi / 180;
                geometryRotate[2]    += dataFile( "Problem/offset", 0., j * geometryColumnsNumber + 6 ) * Pi / 180;

                geometryTranslate[0] += dataFile( "Problem/offset", 0., j * geometryColumnsNumber + 7 );
                geometryTranslate[1] += dataFile( "Problem/offset", 0., j * geometryColumnsNumber + 8 );
                geometryTranslate[2] += dataFile( "Problem/offset", 0., j * geometryColumnsNumber + 9 );
            }

        M_modelsList[i] = multiscaleModelPtr_Type( multiscaleModelFactory_Type::instance().createObject( model ) );
        M_modelsList[i]->setCommunicator( M_comm );
        M_modelsList[i]->setGeometry( geometryScale, geometryRotate, geometryTranslate );
        M_modelsList[i]->setGlobalData( M_globalData );

//!\todo pass a std::string to the factories
//         M_modelsList[i]->setupData( path + model + "/"
//                                     + dataFile( "Problem/models", "undefined", i * modelsColumnsNumber + 2 ) + ".dat" );
        M_modelsList[i]->setupData( path + Enum2String( model, multiscaleModelsMap ) + "/"
                                    + dataFile( "Problem/models", "undefined", i * modelsColumnsNumber + 2 ) + ".dat" );

    }

    // Load couplings
    path = dataFile( "Problem/couplingsPath", "./" );
    for ( UInt i( 0 ); i < couplingsLinesNumber; ++i )
    {
        //!\todo pass a std::string to the factories
        //coupling = dataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 1 );
        coupling = multiscaleCouplingsMap[dataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 1 )];

        M_couplingsList[i] = multiscaleCouplingPtr_Type( multiscaleCouplingFactory_Type::instance().createObject( coupling ) );
        M_couplingsList[i]->setCommunicator( M_comm );
        M_couplingsList[i]->setGlobalData( M_globalData );
        //!\todo pass a std::string to the factories
//         M_couplingsList[i]->setupData( path + coupling + "/"
//                                        + dataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 2 ) + ".dat" );
        M_couplingsList[i]->setupData( path + Enum2String( coupling, multiscaleCouplingsMap ) + "/"
                                       + dataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 2 ) + ".dat" );

        string2numbersVector< UInt > ( dataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 3 ), modelsIDVector );
        string2numbersVector< UInt > ( dataFile( "Problem/couplings", "undefined", i * couplingsColumnsNumber + 4 ), flagsIDVector );
        for ( UInt j( 0 ); j < static_cast< UInt > ( modelsIDVector.size() ); ++j )
        {
            M_couplingsList[i]->addModel( M_modelsList[modelsIDMap[modelsIDVector[j]]] );
            M_couplingsList[i]->addFlagID( flagsIDVector[j] );
            M_modelsList[modelsIDMap[modelsIDVector[j]]]->addCoupling( M_couplingsList[i] );
        }
        modelsIDVector.clear();
        flagsIDVector.clear();
    }
}

void
MultiscaleModelMultiscale::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::SetupModel() \n";
#endif

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->setupCoupling();

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->setupModel();
}

void
MultiscaleModelMultiscale::buildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::BuildSystem() \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->buildSystem();

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->initializeCouplingVariables();
}

void
MultiscaleModelMultiscale::updateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::UpdateSystem() \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->updateSystem();
}

void
MultiscaleModelMultiscale::solveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::SolveSystem() \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->solveSystem();
}

void
MultiscaleModelMultiscale::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::SaveSolution() \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->saveSolution();

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->saveSolution();
}

void
MultiscaleModelMultiscale::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleModel_Type::showMe();
        std::cout << "Models number       = " << M_modelsList.size() << std::endl
                  << "Couplings number    = " << M_couplingsList.size() << std::endl << std::endl;

        std::cout << "==================== Models Information =====================" << std::endl << std::endl;
    }

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->showMe();

    if ( M_displayer->isLeader() )
        std::cout << "=================== Couplings Information ===================" << std::endl << std::endl;

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->showMe();
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleModelMultiscale::createCouplingMap( EpetraMap& couplingMap )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::CreateCouplingMap( couplingMap ) \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->createCouplingMap( couplingMap );

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->createCouplingMap( couplingMap );
}

void
MultiscaleModelMultiscale::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::InitializeCouplingVariables() \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->initializeCouplingVariables();

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->initializeCouplingVariables();
}

void
MultiscaleModelMultiscale::extrapolateCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::ExtrapolateCouplingVariables() \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->extrapolateCouplingVariables();

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->extrapolateCouplingVariables();
}

void
MultiscaleModelMultiscale::importCouplingVariables( const multiscaleVector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::ImportCouplingVariables( couplingVariables ) \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->importCouplingVariables( couplingVariables );

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->importCouplingVariables( couplingVariables );
}

void
MultiscaleModelMultiscale::exportCouplingVariables( multiscaleVector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::ExportCouplingVariables( couplingVariables ) \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportCouplingVariables( couplingVariables);

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->exportCouplingVariables( couplingVariables );
}

void
MultiscaleModelMultiscale::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::ExportCouplingResiduals( couplingResiduals ) \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportCouplingResiduals( couplingResiduals );

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->exportCouplingResiduals( couplingResiduals );
}

void
MultiscaleModelMultiscale::exportJacobian( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::ExportJacobian() \n";
#endif

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportJacobian( jacobian );

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->exportJacobian( jacobian );
}

// ===================================================
// Get Methods
// ===================================================
UInt
MultiscaleModelMultiscale::couplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::GetCouplingVariablesNumber() \n";
#endif

    UInt couplingVariablesNumber = 0;

    for ( multiscaleModelsVectorConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == MultiScale )
            couplingVariablesNumber += ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->couplingVariablesNumber();

    for ( multiscaleCouplingsVectorConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        couplingVariablesNumber += ( *i )->couplingVariablesNumber();

    return couplingVariablesNumber;
}

} // Namespace multiscale
} // Namespace LifeV
