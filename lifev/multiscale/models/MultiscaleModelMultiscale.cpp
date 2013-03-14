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

#include <lifev/multiscale/solver/MultiscaleModelMultiscale.hpp>

#include <lifev/multiscale/solver/MultiscaleAlgorithmAitken.hpp>
#include <lifev/multiscale/solver/MultiscaleAlgorithmBroyden.hpp>
#include <lifev/multiscale/solver/MultiscaleAlgorithmExplicit.hpp>
#include <lifev/multiscale/solver/MultiscaleAlgorithmNewton.hpp>

#include <lifev/multiscale/solver/MultiscaleCouplingBoundaryCondition.hpp>
#include <lifev/multiscale/solver/MultiscaleCouplingMeanNormalStress.hpp>
#include <lifev/multiscale/solver/MultiscaleCouplingMeanNormalStressValve.hpp>
#include <lifev/multiscale/solver/MultiscaleCouplingMeanTotalNormalStress.hpp>

#if defined(LIFEV_HAS_ONEDFSI) && defined(LIFEV_HAS_FSI)
#include <lifev/multiscale/solver/MultiscaleCouplingMeanNormalStressArea.hpp>
#include <lifev/multiscale/solver/MultiscaleCouplingMeanTotalNormalStressArea.hpp>
#endif

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
    M_couplingsList            (),
    M_algorithm                ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::MultiscaleModelMultiscale() \n";
#endif

    M_type = Multiscale;

    multiscaleAlgorithmFactory_Type::instance().registerProduct ( Aitken,                    &createMultiscaleAlgorithmAitken );
    multiscaleAlgorithmFactory_Type::instance().registerProduct ( Broyden,                   &createMultiscaleAlgorithmBroyden );
    multiscaleAlgorithmFactory_Type::instance().registerProduct ( Explicit,                  &createMultiscaleAlgorithmExplicit );
    multiscaleAlgorithmFactory_Type::instance().registerProduct ( Newton,                    &createMultiscaleAlgorithmNewton );

    multiscaleCouplingFactory_Type::instance().registerProduct (  BoundaryCondition,         &createMultiscaleCouplingBoundaryCondition );
    multiscaleCouplingFactory_Type::instance().registerProduct (  MeanNormalStress,          &createMultiscaleCouplingMeanNormalStress );
    multiscaleCouplingFactory_Type::instance().registerProduct (  MeanNormalStressValve,     &createMultiscaleCouplingMeanNormalStressValve );
    multiscaleCouplingFactory_Type::instance().registerProduct (  MeanTotalNormalStress,     &createMultiscaleCouplingMeanTotalNormalStress );

#if defined(LIFEV_HAS_ONEDFSI) && defined(LIFEV_HAS_FSI)
    multiscaleCouplingFactory_Type::instance().registerProduct (  MeanNormalStressArea,      &createMultiscaleCouplingMeanNormalStressArea );
    multiscaleCouplingFactory_Type::instance().registerProduct (  MeanTotalNormalStressArea, &createMultiscaleCouplingMeanTotalNormalStressArea );
#endif
}

MultiscaleModelMultiscale::~MultiscaleModelMultiscale()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::~MultiscaleModelMultiscale() \n";
#endif

    // Disconnect models and couplings to allow their destruction
    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    {
        ( *i )->clearCouplingsList();
    }

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->clearModelsList();
    }
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelMultiscale::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData ( fileName );

    // Definitions
    UInt fileID (0);
    UInt myIDCounter (0);
    std::map< UInt, UInt > modelsIDMap;

    Real load (0);

    models_Type model;
    couplings_Type coupling;
    algorithms_Type algorithm;

    boost::array< Real, NDIM > geometryScale, geometryRotate, geometryTranslate;

    multiscaleIDContainer_Type modelsIDVector;
    multiscaleIDContainer_Type boundaryIDVector;

    // Open file
    GetPot dataFile ( fileName );

    UInt modelsColumnsNumber    = 3.0;
    UInt couplingsColumnsNumber = 5.0;
    UInt mpiGroupsColumnsNumber = 3.0;
    UInt geometryColumnsNumber  = 10.0;

    UInt modelsLinesNumber    = dataFile.vector_variable_size ( "Problem/models" ) / modelsColumnsNumber;
    UInt couplingsLinesNumber = dataFile.vector_variable_size ( "Problem/couplings" ) / couplingsColumnsNumber;
    UInt mpiGroupsLinesNumber = dataFile.vector_variable_size ( "Problem/mpiGroups" ) / mpiGroupsColumnsNumber;
    UInt geometryLinesNumber  = dataFile.vector_variable_size ( "Problem/offset" ) / geometryColumnsNumber;

    // Load MPI groups
    M_commManager.setCommunicator ( M_comm );

    for ( UInt fileMPIGroupsLine ( 0 ); fileMPIGroupsLine < mpiGroupsLinesNumber; ++fileMPIGroupsLine )
    {
        load = dataFile ( "Problem/mpiGroups", 0.0, fileMPIGroupsLine * mpiGroupsColumnsNumber + 1 );
        string2numbersVector< UInt > ( dataFile ( "Problem/mpiGroups", "undefined", fileMPIGroupsLine * mpiGroupsColumnsNumber + 2 ), modelsIDVector );

        M_commManager.addGroup ( load, modelsIDVector );
        modelsIDVector.clear();
    }

    // Split the communicator
    M_commManager.splitCommunicator();

    // Load Models
    std::string path = dataFile ( "Problem/modelsPath", "./" );
    M_modelsList.resize ( M_commManager.myModelsNumber() );
    for ( UInt fileModelsLine ( 0 ); fileModelsLine < modelsLinesNumber; ++fileModelsLine )
    {
        fileID = dataFile ( "Problem/models", 0, fileModelsLine * modelsColumnsNumber );
        if ( M_commManager.myModel ( fileID ) )
        {
            modelsIDMap[fileID] = myIDCounter;
            model = multiscaleModelsMap[dataFile ( "Problem/models", "undefined", fileModelsLine * modelsColumnsNumber + 1 )];

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

            for ( UInt fileGeometryLine ( 0 ); fileGeometryLine < geometryLinesNumber; ++fileGeometryLine )
                if ( fileID == dataFile ( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber ) )
                {
                    geometryScale[0]     *= dataFile ( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber + 1 );
                    geometryScale[1]     *= dataFile ( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber + 2 );
                    geometryScale[2]     *= dataFile ( "Problem/offset", 1., fileGeometryLine * geometryColumnsNumber + 3 );

                    geometryRotate[0]    += dataFile ( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 4 ) * M_PI / 180;
                    geometryRotate[1]    += dataFile ( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 5 ) * M_PI / 180;
                    geometryRotate[2]    += dataFile ( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 6 ) * M_PI / 180;

                    geometryTranslate[0] += dataFile ( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 7 );
                    geometryTranslate[1] += dataFile ( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 8 );
                    geometryTranslate[2] += dataFile ( "Problem/offset", 0., fileGeometryLine * geometryColumnsNumber + 9 );
                }

            M_modelsList[myIDCounter] = multiscaleModelPtr_Type ( multiscaleModelFactory_Type::instance().createObject ( model, multiscaleModelsMap ) );
            M_modelsList[myIDCounter]->setID ( fileModelsLine + 1 );
            M_modelsList[myIDCounter]->setCommunicator ( M_commManager.modelCommunicator ( fileID ) );
            M_modelsList[myIDCounter]->setGeometry ( geometryScale, geometryRotate, geometryTranslate );
            M_modelsList[myIDCounter]->setGlobalData ( M_globalData );
            M_modelsList[myIDCounter]->setupData ( path + enum2String ( model, multiscaleModelsMap ) + "/"
                                                   + dataFile ( "Problem/models", "undefined", fileModelsLine * modelsColumnsNumber + 2 ) + ".dat" );

            // Increment my counter
            ++myIDCounter;
        }
    }

    // Load Couplings
    M_couplingsList.resize ( couplingsLinesNumber );
    path = dataFile ( "Problem/couplingsPath", "./" );
    for ( UInt fileCouplingsLine ( 0 ); fileCouplingsLine < couplingsLinesNumber; ++fileCouplingsLine )
    {
        //id = dataFile( "Problem/couplings", 0, i * couplingsColumnsNumber );
        coupling = multiscaleCouplingsMap[dataFile ( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 1 )];
        M_couplingsList.reserve ( ++myIDCounter );
        M_couplingsList[fileCouplingsLine] = multiscaleCouplingPtr_Type ( multiscaleCouplingFactory_Type::instance().createObject ( coupling, multiscaleCouplingsMap ) );
        M_couplingsList[fileCouplingsLine]->setID ( fileCouplingsLine );
        M_couplingsList[fileCouplingsLine]->setCommunicator ( M_comm );

        string2numbersVector< UInt > ( dataFile ( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 3 ), modelsIDVector );
        string2numbersVector< UInt > ( dataFile ( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 4 ), boundaryIDVector );

        M_couplingsList[fileCouplingsLine]->setModelsNumber ( modelsIDVector.size() );
        for ( UInt j ( 0 ); j < modelsIDVector.size(); ++j )
            if ( M_commManager.myModel ( modelsIDVector[j] ) )
            {
                M_couplingsList[fileCouplingsLine]->setModel ( j, M_modelsList[modelsIDMap[modelsIDVector[j]]] );
                M_couplingsList[fileCouplingsLine]->setBoundaryID ( j, boundaryIDVector[j] );
                M_modelsList[modelsIDMap[modelsIDVector[j]]]->addCoupling ( M_couplingsList[fileCouplingsLine] );
            }
        modelsIDVector.clear();
        boundaryIDVector.clear();

        M_couplingsList[fileCouplingsLine]->setGlobalData ( M_globalData );
        M_couplingsList[fileCouplingsLine]->setupData ( path + enum2String ( coupling, multiscaleCouplingsMap ) + "/"
                                                        + dataFile ( "Problem/couplings", "undefined", fileCouplingsLine * couplingsColumnsNumber + 2 ) + ".dat" );
    }

    // Load Algorithm
    path = dataFile ( "Problem/algorithmsPath", "./" );
    algorithm = multiscaleAlgorithmsMap[ dataFile ( "Problem/algorithm", "Explicit", 0 ) ];
    M_algorithm = multiscaleAlgorithmPtr_Type ( multiscaleAlgorithmFactory_Type::instance().createObject ( algorithm, multiscaleAlgorithmsMap ) );
    M_algorithm->setCommunicator ( M_comm );
    M_algorithm->setupData ( path + enum2String ( algorithm, multiscaleAlgorithmsMap ) + "/" + dataFile ( "Problem/algorithm", "undefined", 1 ) + ".xml" );
}

void
MultiscaleModelMultiscale::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::setupModel() \n";
#endif

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->setupCoupling();
    }

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    {
        ( *i )->setupModel();
    }

    M_algorithm->setMultiscaleModel ( this );
    M_algorithm->setupAlgorithm();
}

void
MultiscaleModelMultiscale::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::buildModel() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    {
        ( *i )->buildModel();
    }

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->initializeCouplingVariables();

        // Necessary to perform a correct interpolation at the first time step
        if ( M_algorithm->type() != Explicit )
        {
            ( *i )->extrapolateCouplingVariables();
        }
    }
}

void
MultiscaleModelMultiscale::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::updateModel() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    {
        ( *i )->updateModel();
    }

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->updateCoupling();

        if ( M_algorithm->type() != Explicit )
        {
            ( *i )->extrapolateCouplingVariables();
        }
        else
        {
            ( *i )->initializeCouplingVariables();
        }
    }
}

void
MultiscaleModelMultiscale::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::solveModel() \n";
#endif

    M_algorithm->subIterate();
}

void
MultiscaleModelMultiscale::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::updateSolution() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    {
        ( *i )->updateSolution();
    }
}

void
MultiscaleModelMultiscale::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::saveSolution() \n";
#endif

    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    {
        ( *i )->saveSolution();
    }

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->saveSolution();
    }

    // Save the framework numbering
    if ( M_globalData->dataTime()->isFirstTimeStep() )
    {
        // File name
        std::string filename = multiscaleProblemFolder + multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) + "_" + number2string ( multiscaleProblemStep ) + ".mfile";

        // Remove the old file if present
        if ( M_comm->MyPID() == 0 )
        {
            std::remove ( filename.c_str() );
        }

        M_comm->Barrier();

        // Open the file
        std::ofstream output;
        output.open ( filename.c_str(), std::ios::app );

        // Models list
        for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
            if ( ( *i )->communicator()->MyPID() == 0 )
            {
                output << "Model ID: " << ( *i )->ID() << ", Name: " << ( *i )->modelName() << std::endl;
            }

        M_comm->Barrier();

        // Couplings list
        if ( M_comm->MyPID() == 0 )
            for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
            {
                output << "Coupling ID: " << ( *i )->ID() << ", Name: " << ( *i )->couplingName() << std::endl;
            }

        // Close the file
        output.close();
    }
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
    {
        ( *i )->showMe();
    }

    if ( M_comm->MyPID() == 0 )
    {
        std::cout << "=================== Couplings Information ===================" << std::endl << std::endl;
    }

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->showMe();
    }

    if ( M_comm->MyPID() == 0 )
    {
        std::cout << "================= Communicators Information =================" << std::endl << std::endl;
    }

    M_commManager.showMe();

    if ( M_comm->MyPID() == 0 )
    {
        std::cout << "=================== Algorithm Information ===================" << std::endl << std::endl;
    }

    M_algorithm->showMe();
}

Real
MultiscaleModelMultiscale::checkSolution() const
{
    return M_algorithm->couplingVariables()->norm2();
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleModelMultiscale::createCouplingMap ( MapEpetra& couplingMap )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::createCouplingMap( couplingMap ) \n";
#endif

    //    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    //        if ( ( *i )->type() == Multiscale )
    //            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->createCouplingMap( couplingMap );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->createCouplingMap ( couplingMap );
    }
}

void
MultiscaleModelMultiscale::importCouplingVariables ( const multiscaleVector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::importCouplingVariables( couplingVariables ) \n";
#endif

    //    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    //        if ( ( *i )->type() == Multiscale )
    //           ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->importCouplingVariables( couplingVariables );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->importCouplingVariables ( couplingVariables );
    }
}

void
MultiscaleModelMultiscale::exportCouplingVariables ( multiscaleVector_Type& couplingVariables )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::exportCouplingVariables( couplingVariables ) \n";
#endif

    //    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    //        if ( ( *i )->type() == Multiscale )
    //            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportCouplingVariables( couplingVariables);

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->exportCouplingVariables ( couplingVariables );
    }
}

void
MultiscaleModelMultiscale::computeCouplingResiduals()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::computeCouplingResiduals() \n";
#endif

    displayModelStatus ( "Solve" );
    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    {
        ( *i )->solveModel();
    }

    //    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    //        if ( ( *i )->type() == Multiscale )
    //            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->computeCouplingResiduals();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->computeCouplingResiduals();
    }
}

void
MultiscaleModelMultiscale::exportCouplingResiduals ( multiscaleVector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::exportCouplingResiduals( couplingResiduals ) \n";
#endif

    //    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    //        if ( ( *i )->type() == Multiscale )
    //            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportCouplingResiduals( couplingResiduals );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->exportCouplingResiduals ( couplingResiduals );
    }
}

void
MultiscaleModelMultiscale::exportJacobian ( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::exportJacobian() \n";
#endif

    //    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    //        if ( ( *i )->type() == Multiscale )
    //            ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->exportJacobian( jacobian );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        ( *i )->exportJacobian ( jacobian );
    }
}

bool
MultiscaleModelMultiscale::topologyChange()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::topologyChange() \n";
#endif

    bool topologyChange ( false );

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        topologyChange = topologyChange || ( *i )->topologyChange();
    }

    return topologyChange;
}

// ===================================================
// Get Methods
// ===================================================
UInt
MultiscaleModelMultiscale::couplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8110 ) << "MultiscaleModelMultiscale::couplingVariablesNumber() \n";
#endif

    UInt couplingVariablesNumber = 0;

    //    for ( multiscaleModelsContainerConstIterator_Type i = M_modelsList.begin(); i != M_modelsList.end(); ++i )
    //        if ( ( *i )->type() == Multiscale )
    //            couplingVariablesNumber += ( multiscaleDynamicCast< MultiscaleModelMultiscale > ( *i ) )->couplingVariablesNumber();

    for ( multiscaleCouplingsContainerConstIterator_Type i = M_couplingsList.begin(); i != M_couplingsList.end(); ++i )
    {
        couplingVariablesNumber += ( *i )->couplingVariablesNumber();
    }

    return couplingVariablesNumber;
}

} // Namespace multiscale
} // Namespace LifeV
