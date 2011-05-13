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
 *  @brief File containing the Multiscale Model Multiscale
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleModelMultiscale.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelMultiscale::MultiscaleModelMultiscale() :
        multiscaleModel_Type       (),
        M_commManager              (),
        M_modelsList               (),
        M_couplingsList            ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::MultiscaleModelMultiscale() \n";
#endif

    M_type = Multiscale;
}

MultiscaleModelMultiscale::~MultiscaleModelMultiscale()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::~MultiscaleModelMultiscale( ) \n";
#endif

    // Disconnect models and couplings to allow their destruction
    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->clearCouplingsList();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->clearModelsList();
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelMultiscale::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData( fileName );

    // Definitions
    UInt fileID(0);
    UInt myIDCounter(0);
    std::map< UInt, UInt > modelsIDMap;

    Real load(0);

    models_Type model;
    couplings_Type coupling;

    boost::array< Real, NDIM > geometryScale, geometryRotate, geometryTranslate;

    std::vector< UInt > modelsIDVector;
    std::vector< UInt > flagsIDVector;

    // Open file
    GetPot dataFile( fileName );

    UInt modelsColumnsNumber    = 3.0;
    UInt couplingsColumnsNumber = 5.0;
    UInt mpiGroupsColumnsNumber = 3.0;
    UInt geometryColumnsNumber  = 10.0;

    UInt modelsLinesNumber    = dataFile.vector_variable_size( "Problem/models" ) / modelsColumnsNumber;
    UInt couplingsLinesNumber = dataFile.vector_variable_size( "Problem/couplings" ) / couplingsColumnsNumber;
    UInt mpiGroupsLinesNumber = dataFile.vector_variable_size( "Problem/mpiGroups" ) / mpiGroupsColumnsNumber;
    UInt geometryLinesNumber  = dataFile.vector_variable_size( "Problem/offset" ) / geometryColumnsNumber;

    // Load MPI groups
    M_commManager.setCommunicator( M_comm );

    for ( UInt fileMPIGroupsLine( 0 ); fileMPIGroupsLine < mpiGroupsLinesNumber; ++fileMPIGroupsLine )
    {
        load = dataFile( "Problem/mpiGroups", 0.0, fileMPIGroupsLine * mpiGroupsColumnsNumber + 1 );
        string2numbersVector< UInt > ( dataFile( "Problem/mpiGroups", "undefined", fileMPIGroupsLine * mpiGroupsColumnsNumber + 2 ), modelsIDVector );

        M_commManager.addGroup( load, modelsIDVector );
        modelsIDVector.clear();
    }

    // Split the communicator
    M_commManager.splitCommunicator();

    // Load Models
    std::string path = dataFile( "Problem/modelsPath", "./" );
    M_modelsList.resize( M_commManager.myModelsNumber() );
    for ( UInt fileModelsLine( 0 ); fileModelsLine < modelsLinesNumber; ++fileModelsLine )
    {
        fileID = dataFile( "Problem/models", 0, fileModelsLine * modelsColumnsNumber );
        if ( M_commManager.myModel( fileID ) )
        {
            modelsIDMap[fileID] = myIDCounter;
            model = multiscaleModelsMap[dataFile( "Problem/models", "undefined", fileModelsLine * modelsColumnsNumber + 1 )];

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

            for ( UInt fileGeometryLine( 0 ); fileGeometryLine < geometryLinesNumber; ++fileGeometryLine )
                if ( fileID == dataFile( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber ) )
                {
                    geometryScale[0]     *= dataFile( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber + 1 );
                    geometryScale[1]     *= dataFile( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber + 2 );
                    geometryScale[2]     *= dataFile( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber + 3 );

                    geometryRotate[0]    += dataFile( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 4 ) * M_PI / 180;
                    geometryRotate[1]    += dataFile( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 5 ) * M_PI / 180;
                    geometryRotate[2]    += dataFile( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 6 ) * M_PI / 180;

                    geometryTranslate[0] += dataFile( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 7 );
                    geometryTranslate[1] += dataFile( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 8 );
                    geometryTranslate[2] += dataFile( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 9 );
                }

            M_modelsList[myIDCounter] = multiscaleModelPtr_Type( multiscaleModelFactory_Type::instance().createObject( model, multiscaleModelsMap ) );
            M_modelsList[myIDCounter]->setID( fileModelsLine + 1 );
            M_modelsList[myIDCounter]->setCommunicator( M_commManager.modelCommunicator( fileID ) );
            M_modelsList[myIDCounter]->setGeometry( geometryScale, geometryRotate, geometryTranslate );
            M_modelsList[myIDCounter]->setGlobalData( M_globalData );
            M_modelsList[myIDCounter]->setupData( path + enum2String( model, multiscaleModelsMap ) + "/"
                                                  + dataFile( "Problem/models", "undefined", fileModelsLine * modelsColumnsNumber + 2 ) + ".dat" );

            // Increment my counter
            ++myIDCounter;
        }
    }

    // Reset my counter
    myIDCounter = 0;

    // Load Couplings
    path = dataFile( "Problem/couplingsPath", "./" );
    for ( UInt fileCouplingsLine( 0 ); fileCouplingsLine < couplingsLinesNumber; ++fileCouplingsLine )
    {
        string2numbersVector< UInt > ( dataFile( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 3 ), modelsIDVector );
        string2numbersVector< UInt > ( dataFile( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 4 ), flagsIDVector );
        for ( UInt j( 0 ); j < modelsIDVector.size(); ++j )
        {
            if ( M_commManager.myModel( modelsIDVector[j] ) )
            {
                if ( M_couplingsList.size() == myIDCounter )
                {
                    //id = dataFile( "Problem/couplings", 0, i * couplingsColumnsNumber );
                    coupling = multiscaleCouplingsMap[dataFile( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 1 )];
                    M_couplingsList.reserve( ++myIDCounter );
                    M_couplingsList.push_back( multiscaleCouplingPtr_Type( multiscaleCouplingFactory_Type::instance().createObject( coupling, multiscaleCouplingsMap ) ) );
                    M_couplingsList.back()->setID( fileCouplingsLine );
                    M_couplingsList.back()->setCommunicator( M_commManager.modelCommunicator( modelsIDVector[j] ) );
                    M_couplingsList.back()->setGlobalData( M_globalData );
                    M_couplingsList.back()->setupData( path + enum2String( coupling, multiscaleCouplingsMap ) + "/"
                                                       + dataFile( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 2 ) + ".dat" );
                }

                M_couplingsList.back()->addModel( M_modelsList[modelsIDMap[modelsIDVector[j]]] );
                M_couplingsList.back()->addFlagID( flagsIDVector[j] );
                M_modelsList[modelsIDMap[modelsIDVector[j]]]->addCoupling( M_couplingsList.back() );
            }
        }
        modelsIDVector.clear();
        flagsIDVector.clear();
    }
}

void
MultiscaleModelMultiscale::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::setupModel() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->setupCoupling();

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->setupModel();
}

void
MultiscaleModelMultiscale::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::buildModel() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->buildModel();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->initializeCouplingVariables();
}

void
MultiscaleModelMultiscale::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::updateModel() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->updateModel();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->updateCoupling();
}

void
MultiscaleModelMultiscale::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::solveModel() \n";
#endif

    displayModelStatus( "Solve" );
    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->solveModel();
}

void
MultiscaleModelMultiscale::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::saveSolution() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->saveSolution();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->saveSolution();
}

void
MultiscaleModelMultiscale::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();
        std::cout << "Models number       = " << M_modelsList.size() << std::endl
                  << "Couplings number    = " << M_couplingsList.size() << std::endl << std::endl;

        std::cout << "==================== Models Information =====================" << std::endl << std::endl;
    }

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        ( *i )->showMe();

    if ( M_comm->MyPID() == 0 )
        std::cout << "=================== Couplings Information ===================" << std::endl << std::endl;

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->showMe();

    if ( M_comm->MyPID() == 0 )
        std::cout << "================= Communicators Information =================" << std::endl << std::endl;

    M_commManager.showMe();
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleModelMultiscale::createCouplingMap( MapEpetra& couplingMap )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::createCouplingMap( couplingMap ) \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->createCouplingMap( couplingMap );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->createCouplingMap( couplingMap );
}

void
MultiscaleModelMultiscale::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::initializeCouplingVariables() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->initializeCouplingVariables();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->initializeCouplingVariables();
}

void
MultiscaleModelMultiscale::extrapolateCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::extrapolateCouplingVariables() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->extrapolateCouplingVariables();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->extrapolateCouplingVariables();
}

void
MultiscaleModelMultiscale::importCouplingVariables( const multiscaleVector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::importCouplingVariables( couplingVariables ) \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->importCouplingVariables( couplingVariables );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->importCouplingVariables( couplingVariables );
}

void
MultiscaleModelMultiscale::exportCouplingVariables( multiscaleVector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::exportCouplingVariables( couplingVariables ) \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportCouplingVariables( couplingVariables);

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->exportCouplingVariables( couplingVariables );
}

void
MultiscaleModelMultiscale::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::exportCouplingResiduals( couplingResiduals ) \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportCouplingResiduals( couplingResiduals );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->exportCouplingResiduals( couplingResiduals );
}

void
MultiscaleModelMultiscale::exportJacobian( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::exportJacobian() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportJacobian( jacobian );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        ( *i )->exportJacobian( jacobian );
}

bool
MultiscaleModelMultiscale::topologyChange()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::topologyChange() \n";
#endif

    bool topologyChange( false );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        topologyChange = topologyChange || ( *i )->topologyChange();

    return topologyChange;
}

// ===================================================
// Get Methods
// ===================================================
UInt
MultiscaleModelMultiscale::couplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8110 ) << "MultiscaleModelMultiscale::couplingVariablesNumber() \n";
#endif

    UInt couplingVariablesNumber = 0;

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
        if ( ( *i )->type() == Multiscale )
            couplingVariablesNumber += ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->couplingVariablesNumber();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
        couplingVariablesNumber += ( *i )->couplingVariablesNumber();

    return couplingVariablesNumber;
}

} // Namespace multiscale
} // Namespace LifeV
