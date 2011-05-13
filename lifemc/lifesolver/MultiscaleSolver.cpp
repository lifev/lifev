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
 *  @brief File containing the Multiscale Solver
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

UInt        multiscaleCoresPerNode  = 1;
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
        M_chrono            ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::MultiscaleSolver() \n";
#endif

    //Define the maps of MS objects
    multiscaleMapsDefinition();

    //Register the objects
    multiscaleModelFactory_Type::instance().registerProduct   (  Fluid3D,             &createMultiscaleModelFluid3D );
    multiscaleModelFactory_Type::instance().registerProduct   (  FSI3D,               &createMultiscaleModelFSI3D );
    multiscaleModelFactory_Type::instance().registerProduct   (  Multiscale,          &createMultiscaleModelMultiscale );
    multiscaleModelFactory_Type::instance().registerProduct   (  OneDimensional,      &createMultiscaleModelOneDimensional );
    multiscaleModelFactory_Type::instance().registerProduct   (  Windkessel0D,        &createMultiscaleModelWindkessel0D );

    multiscaleCouplingFactory_Type::instance().registerProduct(  BoundaryCondition,   &createMultiscaleCouplingBoundaryCondition );
    multiscaleCouplingFactory_Type::instance().registerProduct(  FlowRate,            &createMultiscaleCouplingFlowRate );
    multiscaleCouplingFactory_Type::instance().registerProduct(  FlowRateValve,       &createMultiscaleCouplingFlowRateValve );
    multiscaleCouplingFactory_Type::instance().registerProduct(  FlowRateStress,      &createMultiscaleCouplingFlowRateStress );
    multiscaleCouplingFactory_Type::instance().registerProduct(  Stress,              &createMultiscaleCouplingStress );

    multiscaleAlgorithmFactory_Type::instance().registerProduct( Aitken,              &createMultiscaleAlgorithmAitken );
    multiscaleAlgorithmFactory_Type::instance().registerProduct( Broyden,             &createMultiscaleAlgorithmBroyden );
    multiscaleAlgorithmFactory_Type::instance().registerProduct( Explicit,            &createMultiscaleAlgorithmExplicit );
    multiscaleAlgorithmFactory_Type::instance().registerProduct( Newton,              &createMultiscaleAlgorithmNewton );
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleSolver::setupProblem( const std::string& fileName, const std::string& problemFolder, const UInt& coresPerNode )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::setupData( fileName, problemFolder ) \n";
#endif

    // Load data file
    GetPot dataFile( fileName );

    // Define the folder containing the problem
    multiscaleProblemFolder = problemFolder;

    // Define the number of cores on each node for the machine
    multiscaleCoresPerNode = coresPerNode;

    // Define the step of the problem
    if ( dataFile( "Solver/Restart/Restart", false ) )
        multiscaleProblemStep = dataFile( "Solver/Restart/RestartFromStepNumber", 0 ) + 1;

    // Create the main model and set the communicator
    M_model = multiscaleModelPtr_Type( multiscaleModelFactory_Type::instance().createObject( multiscaleModelsMap[ dataFile( "Problem/ProblemType", "Multiscale" ) ], multiscaleModelsMap ) );

    M_model->setID( 0 );
    M_model->setCommunicator( M_comm );

    // Setup data
    M_globalData->readData( dataFile );
    M_model->setGlobalData( M_globalData );
    M_model->setupData( dataFile( "Problem/ProblemFile", "./MultiscaleData/Models/NoModel.dat" ) + ".dat" );

    // Setup Models
    M_model->setupModel();

    // Algorithm parameters
    if ( M_model->type() == Multiscale )
    {
        M_algorithm = multiscaleAlgorithmPtr_Type( multiscaleAlgorithmFactory_Type::instance().createObject( multiscaleAlgorithmsMap[ dataFile( "Solver/Algorithm/type", "Newton" ) ], multiscaleAlgorithmsMap ) );

        M_algorithm->setCommunicator( M_comm );
        M_algorithm->setModel( M_model );
        M_algorithm->setSubiterationsMaximumNumber( dataFile( "Solver/Algorithm/subiterationsMaximumNumber", 100 ) );
        M_algorithm->setTolerance( dataFile( "Solver/Algorithm/tolerance", 1e-2 ) );
        std::string path = "./MultiscaleDatabase/Algorithms/"; // TODO Add this to files
        M_algorithm->setupData( path + enum2String( M_algorithm->type(), multiscaleAlgorithmsMap ) + "/" + dataFile( "Solver/Algorithm/file", "undefined" ) + ".dat" );
    }
}

bool
MultiscaleSolver::solveProblem( const Real& referenceSolution )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::solveProblem() \n";
#endif

    // save initial solution if it is the very first time step
    if ( !multiscaleProblemStep )
        M_model->saveSolution();

    // Move to the "true" first time-step
    M_globalData->dataTime()->updateTime();
    M_globalData->dataTime()->setInitialTime( M_globalData->dataTime()->time() );

    Real totalSimulationTime(0);
    for ( ; M_globalData->dataTime()->canAdvance(); M_globalData->dataTime()->updateTime() )
    {
        M_chrono.start();

        if ( M_comm->MyPID() == 0 )
        {
            std::cout << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
            std::cout << "                    MULTISCALE FRAMEWORK" << std::endl;
            std::cout << "             time = " << M_globalData->dataTime()->time() << " s; "
                      << "time step number = " << M_globalData->dataTime()->timeStepNumber()  << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
        }

        // Build or Update System
        if ( M_globalData->dataTime()->isFirstTimeStep() )
        {
            if ( M_model->type() == Multiscale )
                M_algorithm->initializeCouplingVariables();
            M_model->buildModel();
        }
        else
        {
            if ( M_model->type() == Multiscale )
                M_algorithm->updateCouplingVariables();
            M_model->updateModel();
        }

        // Solve the model
        M_model->solveModel();

        // If it is a Multiscale model, call algorithm for subiterations
        if ( M_model->type() == Multiscale )
            M_algorithm->subIterate();

        // Save the solution
        M_model->saveSolution();

        // Chrono stop
        M_chrono.stop();

        // Updating total simulation time
        totalSimulationTime += M_chrono.diff();

        if ( M_comm->MyPID() == 0 )
            std::cout << " MS-  Total iteration time:                    " << M_chrono.diff() << " s" << std::endl;
    }

    if ( M_comm->MyPID() == 0 )
        std::cout << " MS-  Total simulation time:                   " << totalSimulationTime << " s" << std::endl;

    // Check on the last iteration
    if ( M_model->type() == Multiscale )
    {
        Real computedSolution( M_algorithm->couplingVariables()->norm2() );
        if ( referenceSolution >= 0. && std::abs( referenceSolution - computedSolution ) > 1e-8 )
            multiscaleErrorCheck( Solution, "Algorithm Solution: "  + number2string( computedSolution ) +
                                            " (External Residual: " + number2string( referenceSolution ) + ")\n", M_comm->MyPID() );
    }

    return multiscaleExitFlag;
}

void
MultiscaleSolver::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        std::cout << std::endl << std::endl
                  << "=============== Multiscale Solver Information ===============" << std::endl << std::endl;

        std::cout << "Cores per node                = " << multiscaleCoresPerNode << std::endl
                  << "Problem folder                = " << multiscaleProblemFolder << std::endl
                  << "Problem step                  = " << multiscaleProblemStep << std::endl << std::endl;

        M_globalData->showMe();

        std::cout << std::endl << std::endl;
    }

    M_model->showMe();
    if ( M_model->type() == Multiscale )
        M_algorithm->showMe();

    if ( M_comm->MyPID() == 0 )
        std::cout << "=============================================================" << std::endl << std::endl;
}

} // Namespace multiscale
} // Namespace LifeV
