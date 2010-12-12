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

#include <lifemc/lifesolver/MS_Solver.hpp>

namespace LifeV
{

std::string MS_ProblemFolder = "./";
UInt        MS_ProblemStep   = 0;
bool        MS_ExitFlag      = EXIT_SUCCESS;

// ===================================================
// Constructors
// ===================================================
MS_Solver::MS_Solver() :
        M_model             (),
        M_algorithm         (),
        M_globalData        ( new MS_GlobalDataContainer_Type() ),
        M_comm              (),
        M_displayer         (),
        M_chrono            ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MS_Solver::MS_Solver() \n";
#endif

    //Define the maps of MS objects
    MS_MapsDefinition();

    //Register the objects
    MS_Model_Factory::instance().registerProduct   (  MultiScale,          &createMultiscaleModelMultiscale );
    MS_Model_Factory::instance().registerProduct   (  Fluid3D,             &createMultiscaleModelFluid3D );
    MS_Model_Factory::instance().registerProduct   (  OneDimensional,      &createMultiscaleModelOneDimensional );
    MS_Model_Factory::instance().registerProduct   (  FSI3D,               &createMultiscaleModelFSI3D );

    MS_Coupling_Factory::instance().registerProduct(  Stress,              &createMultiscaleCouplingStress );
    MS_Coupling_Factory::instance().registerProduct(  FlowRateStress,      &createMultiscaleCouplingFlowRateStress );
    MS_Coupling_Factory::instance().registerProduct(  BoundaryCondition,   &createMultiscaleCouplingBoundaryCondition );

    MS_Algorithm_Factory::instance().registerProduct( Aitken,              &createMultiscaleAlgorithmAitken );
    MS_Algorithm_Factory::instance().registerProduct( Explicit,            &createMultiscaleAlgorithmExplicit );
    MS_Algorithm_Factory::instance().registerProduct( Newton,              &createMultiscaleAlgorithmNewton );
}

// ===================================================
// Methods
// ===================================================
void
MS_Solver::setCommunicator( const MS_Comm_PtrType& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MS_Solver::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

void
MS_Solver::setupProblem( const std::string& fileName, const std::string& problemFolder )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MS_Solver::SetupData( fileName, problemFolder ) \n";
#endif

    // Load data file
    GetPot dataFile( fileName );

    // Define the folder containing the problem
    MS_ProblemFolder = problemFolder;

    // Define the step of the problem
    if ( dataFile( "Solver/Restart/Restart", false ) )
        MS_ProblemStep = dataFile( "Solver/Restart/RestartFromStepNumber", 0 ) + 1;

    // Create the main model and set the communicator
    M_model = MS_Model_PtrType( MS_Model_Factory::instance().createObject( MS_modelsMap[ dataFile( "Problem/ProblemType", "MultiScale" ) ] ) );
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
        M_algorithm = MS_Algorithm_PtrType( MS_Algorithm_Factory::instance().createObject( MS_algorithmsMap[ dataFile( "Solver/Algorithm/AlgorithmType", "Newton" ) ] ) );
        M_algorithm->setCommunicator( M_comm );
        M_algorithm->setModel( M_model );
        M_algorithm->setupData( fileName );
        M_algorithm->initializeCouplingVariables();
    }
}

bool
MS_Solver::solveProblem( const Real& externalResidual )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MS_Solver::SolveProblem() \n";
#endif

    // save initial solution if it is the very first time step
    if ( !MS_ProblemStep )
        M_model->saveSolution();

    // Move to the "true" first time-step
    M_globalData->dataTime()->updateTime();
    M_globalData->dataTime()->setInitialTime( M_globalData->dataTime()->getTime() );

    for ( ; M_globalData->dataTime()->canAdvance(); M_globalData->dataTime()->updateTime() )
    {
        M_chrono.start();

        if ( M_displayer->isLeader() )
        {
            std::cout << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
            std::cout << "                    MULTISCALE SIMULATION" << std::endl;
            std::cout << "             time = " << M_globalData->dataTime()->getTime() << " s; "  <<
                      "time step number = " << M_globalData->dataTime()->getTimeStepNumber()  << std::endl;
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
        MS_ErrorCheck( MS_Residual,
                       "Algorithm Residual: " + number2string( algorithmResidual ) +
                       " (External Residual: " + number2string( externalResidual ) + ")\n" );

    return MS_ExitFlag;
}

void
MS_Solver::showMe()
{
    if ( M_displayer->isLeader() )
    {
        std::cout << std::endl << std::endl
                  << "=============== MultiScale Solver Information ===============" << std::endl << std::endl;

        std::cout << "Problem folder                = " << MS_ProblemFolder << std::endl
                  << "Problem step                  = " << MS_ProblemStep << std::endl << std::endl;

        M_globalData->showMe();

        std::cout << std::endl << std::endl;
    }

    M_model->showMe();
    if ( M_model->type() == MultiScale )
        M_algorithm->showMe();

    if ( M_displayer->isLeader() )
        std::cout << "=============================================================" << std::endl << std::endl;
}

} // Namespace LifeV
