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
        M_globalData        ( new multiscaleData_Type() ),
        M_comm              ()
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
    M_model->setupData( dataFile( "Problem/ProblemFile", "./MultiscaleDatabase/Models/NoModel" ) + ".dat" );

    // Setup Models
    M_model->setupModel();
}

bool
MultiscaleSolver::solveProblem( const Real& referenceSolution )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8000 ) << "MultiscaleSolver::solveProblem() \n";
#endif

    // Save the initial solution if it is the very first time step
    if ( !multiscaleProblemStep )
        M_model->saveSolution();

    // Move to the "true" first time-step
    M_globalData->dataTime()->updateTime();
    M_globalData->dataTime()->setInitialTime( M_globalData->dataTime()->time() );

    // Chrono definitions
    LifeChrono buildUpdateChrono;
    LifeChrono solveChrono;
    LifeChrono saveChrono;
    LifeChrono globalChrono;
    Real       totalSimulationTime(0);
    Real       timeStepTime(0);

    for ( ; M_globalData->dataTime()->canAdvance(); M_globalData->dataTime()->updateTime() )
    {
        // Global chrono start
        globalChrono.start();

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
        buildUpdateChrono.start();
        if ( M_globalData->dataTime()->isFirstTimeStep() )
            M_model->buildModel();
        else
            M_model->updateModel();
        buildUpdateChrono.stop();

        // Solve the model
        solveChrono.start();
        M_model->solveModel();
        solveChrono.stop();

        // Save the solution
        saveChrono.start();
        M_model->saveSolution();
        saveChrono.stop();

        // Global chrono stop
        globalChrono.stop();

        // Compute time step time
        saveCPUTime( buildUpdateChrono.globalDiff( *M_comm ), solveChrono.globalDiff( *M_comm ), saveChrono.globalDiff( *M_comm ) );
        timeStepTime = globalChrono.globalDiff( *M_comm );

        if ( M_comm->MyPID() == 0 )
            std::cout << " MS-  Total iteration time:                    " << timeStepTime << " s" << std::endl;

        // Updating total simulation time
        totalSimulationTime += timeStepTime;
    }

    if ( M_comm->MyPID() == 0 )
        std::cout << " MS-  Total simulation time:                   " << totalSimulationTime << " s" << std::endl;

    // Check on the last iteration
    Real computedSolution( M_model->checkSolution() );
    if ( referenceSolution >= 0. && std::abs( referenceSolution - computedSolution ) > 1e-8 )
        multiscaleErrorCheck( Solution, "Problem solution: "  + number2string( computedSolution ) +
                                        " (External solution: " + number2string( referenceSolution ) + ")\n", M_comm->MyPID() == 0 );

    return multiscaleExitFlag;
}

void
MultiscaleSolver::showMe() const
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

    if ( M_comm->MyPID() == 0 )
        std::cout << "=============================================================" << std::endl << std::endl;
}

void
MultiscaleSolver::saveCPUTime( const Real& buildUpdateCPUTime, const Real& solveCPUTime, const Real& saveCPUTime ) const
{
    if ( M_comm->MyPID() == 0 )
    {
        std::ofstream output;
        output << std::scientific << std::setprecision( 15 );

        std::string filename = multiscaleProblemFolder + "Step_" + number2string( multiscaleProblemStep )
                                                       + "_CPUTime.mfile";

        if ( M_globalData->dataTime()->isFirstTimeStep() )
        {
            output.open( filename.c_str(), std::ios::trunc );
            output << "% TIME                     TOTAL                    BUILD/UPDATE             SOLVE                    SAVE" << std::endl;
        }
        else
        {
            output.open( filename.c_str(), std::ios::app );
        }
        output << "  " << M_globalData->dataTime()->time() << "    " << buildUpdateCPUTime+solveCPUTime+saveCPUTime
               << "    " << buildUpdateCPUTime << "    " << solveCPUTime  << "    " << saveCPUTime << std::endl;
        output.close();
    }
}

} // Namespace multiscale
} // Namespace LifeV
