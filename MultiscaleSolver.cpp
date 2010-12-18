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
//@HEADERR

/*!
 *  @file
 *  @brief File containing the MultiScale Solver
 *
 *  @date 28-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleSolver.hpp>

namespace LifeV
{
namespace Multiscale
{

std::string multiscaleProblemFolder = "./";
UInt        multiscaleProblemStep   = 0;
bool        multiscaleExitFlag      = EXIT_SUCCESS;

// ===================================================
// Constructors
// ===================================================
MultiscaleSolver::MultiscaleSolver() :
        M_model             (),
        M_algorithm         (),
        M_globalData        ( new multiscaleData_Type() ),
        M_comm              (),
        M_displayer         (),
        M_chrono            ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::MultiscaleSolver() \n";
#endif

    //Define the maps of MS objects
    multiscaleMapsDefinition();

    //Register the objects
    //!\todo pass a std::string to the factories
//     multiscaleModelFactory_Type::instance().registerProduct   (  "MultiScale",          &createMultiscaleModelMultiscale );
//     multiscaleModelFactory_Type::instance().registerProduct   (  "Fluid3D",             &createMultiscaleModelFluid3D );
//     multiscaleModelFactory_Type::instance().registerProduct   (  "OneDimensional",      &createMultiscaleModelOneDimensional );
//     multiscaleModelFactory_Type::instance().registerProduct   (  "FSI3D",               &createMultiscaleModelFSI3D );

//     multiscaleCouplingFactory_Type::instance().registerProduct(  "Stress",              &createMultiscaleCouplingStress );
//     multiscaleCouplingFactory_Type::instance().registerProduct(  "FlowRateStress",      &createMultiscaleCouplingFlowRateStress );
//     multiscaleCouplingFactory_Type::instance().registerProduct(  "BoundaryCondition",   &createMultiscaleCouplingBoundaryCondition );

//     multiscaleAlgorithmFactory_Type::instance().registerProduct( "Aitken",              &createMultiscaleAlgorithmAitken );
//     multiscaleAlgorithmFactory_Type::instance().registerProduct( "Explicit",            &createMultiscaleAlgorithmExplicit );
//     multiscaleAlgorithmFactory_Type::instance().registerProduct( "Newton",              &createMultiscaleAlgorithmNewton );

    multiscaleModelFactory_Type::instance().registerProduct   (  MultiScale,          &createMultiscaleModelMultiscale );
    multiscaleModelFactory_Type::instance().registerProduct   (  Fluid3D,             &createMultiscaleModelFluid3D );
    multiscaleModelFactory_Type::instance().registerProduct   (  OneDimensional,      &createMultiscaleModelOneDimensional );
    multiscaleModelFactory_Type::instance().registerProduct   (  FSI3D,               &createMultiscaleModelFSI3D );

    multiscaleCouplingFactory_Type::instance().registerProduct(  Stress,              &createMultiscaleCouplingStress );
    multiscaleCouplingFactory_Type::instance().registerProduct(  FlowRateStress,      &createMultiscaleCouplingFlowRateStress );
    multiscaleCouplingFactory_Type::instance().registerProduct(  BoundaryCondition,   &createMultiscaleCouplingBoundaryCondition );

    multiscaleAlgorithmFactory_Type::instance().registerProduct( Aitken,              &createMultiscaleAlgorithmAitken );
    multiscaleAlgorithmFactory_Type::instance().registerProduct( Explicit,            &createMultiscaleAlgorithmExplicit );
    multiscaleAlgorithmFactory_Type::instance().registerProduct( Newton,              &createMultiscaleAlgorithmNewton );
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleSolver::setCommunicator( const multiscaleCommPtr_Type& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

void
MultiscaleSolver::setupProblem( const std::string& fileName, const std::string& problemFolder )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::SetupData( fileName, problemFolder ) \n";
#endif

    // Load data file
    GetPot dataFile( fileName );

    // Define the folder containing the problem
    multiscaleProblemFolder = problemFolder;

    // Define the step of the problem
    if ( dataFile( "Solver/Restart/Restart", false ) )
        multiscaleProblemStep = dataFile( "Solver/Restart/RestartFromStepNumber", 0 ) + 1;

    // Create the main model and set the communicator
    //!\todo pass a std::string to the factories
//     M_model = multiscaleModelPtr_Type( multiscaleModelFactory_Type::instance().createObject(  dataFile( "Problem/ProblemType", "MultiScale" ) ) );

    M_model = multiscaleModelPtr_Type( multiscaleModelFactory_Type::instance().createObject( multiscaleModelsMap[ dataFile( "Problem/ProblemType", "MultiScale" ) ] ) );

    M_model->setCommunicator( M_comm );

    // Setup data
    M_globalData->readData( dataFile );
    M_model->setGlobalData( M_globalData );
    M_model->setupData( dataFile( "Problem/ProblemFile", "./MultiScaleData/Models/Model.dat" ) );

    // Setup Models
    M_model->setupModel();

    // Algorithm parameters
    if ( M_model->type() == MultiScale )
    {
        //!\todo pass a std::string to the factories
        //        M_algorithm = multiscaleAlgorithmPtr_Type( multiscaleAlgorithmFactory_Type::instance().createObject( dataFile( "Solver/Algorithm/AlgorithmType", "Newton" ) ) );
        M_algorithm = multiscaleAlgorithmPtr_Type( multiscaleAlgorithmFactory_Type::instance().createObject( multiscaleAlgorithmsMap[ dataFile( "Solver/Algorithm/AlgorithmType", "Newton" ) ] ) );

        M_algorithm->setCommunicator( M_comm );
        M_algorithm->setModel( M_model );
        M_algorithm->setupData( fileName );
        M_algorithm->initializeCouplingVariables();
    }
}

bool
MultiscaleSolver::solveProblem( const Real& externalResidual )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::SolveProblem() \n";
#endif

    // save initial solution if it is the very first time step
    if ( !multiscaleProblemStep )
        M_model->saveSolution();

    // Move to the "true" first time-step
    M_globalData->dataTime()->updateTime();
    M_globalData->dataTime()->setInitialTime( M_globalData->dataTime()->time() );

    for ( ; M_globalData->dataTime()->canAdvance(); M_globalData->dataTime()->updateTime() )
    {
        M_chrono.start();

        if ( M_displayer->isLeader() )
        {
            std::cout << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
            std::cout << "                    MULTISCALE SIMULATION" << std::endl;
            std::cout << "             time = " << M_globalData->dataTime()->time() << " s; "  <<
                      "time step number = " << M_globalData->dataTime()->timeStepNumber()  << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
        }

        // Build or Update System
        if ( M_globalData->dataTime()->isFirstTimeStep() )
            M_model->buildSystem();
        else
        {
            if ( M_model->type() == MultiScale )
                M_algorithm->updateCouplingVariables();
            M_model->updateSystem();
        }

        // solveSystem
        M_model->solveSystem();

        // If it is a MultiScale model, call algorithms for subiterations
        if ( M_model->type() == MultiScale )
            M_algorithm->subIterate();

        // saveSolution
        M_model->saveSolution();

        M_chrono.stop();

        if ( M_displayer->isLeader() )
            std::cout << " MS-  Total iteration time:                    " << M_chrono.diff() << " s" << std::endl;
    }

    // Redisual check
    Real algorithmResidual( M_algorithm->computeResidual() );
    if ( externalResidual >= 0. && std::abs( externalResidual - algorithmResidual ) > 1e-8 )
        multiscaleErrorCheck( Residual, "Algorithm Residual: " + number2string( algorithmResidual ) +
                       " (External Residual: " + number2string( externalResidual ) + ")\n" );

    return multiscaleExitFlag;
}

void
MultiscaleSolver::showMe()
{
    if ( M_displayer->isLeader() )
    {
        std::cout << std::endl << std::endl
                  << "=============== MultiScale Solver Information ===============" << std::endl << std::endl;

        std::cout << "Problem folder                = " << multiscaleProblemFolder << std::endl
                  << "Problem step                  = " << multiscaleProblemStep << std::endl << std::endl;

        M_globalData->showMe();

        std::cout << std::endl << std::endl;
    }

    M_model->showMe();
    if ( M_model->type() == MultiScale )
        M_algorithm->showMe();

    if ( M_displayer->isLeader() )
        std::cout << "=============================================================" << std::endl << std::endl;
}

} // Namespace multiscale
} // Namespace LifeV
