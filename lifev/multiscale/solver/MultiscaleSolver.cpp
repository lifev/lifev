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

#include <lifev/multiscale/solver/MultiscaleSolver.hpp>

namespace LifeV
{
namespace Multiscale
{

UInt        multiscaleCoresPerNode       = 1;
std::string multiscaleProblemFolder      = "./";
std::string multiscaleProblemPrefix      = "Multiscale";
UInt        multiscaleProblemStep        = 0;
UInt        multiscaleSaveEachNTimeSteps = 1;
bool        multiscaleExitFlag           = EXIT_SUCCESS;

// ===================================================
// Constructors
// ===================================================
MultiscaleSolver::MultiscaleSolver() :
    M_model                 (),
    M_globalData            ( new multiscaleData_Type() ),
    M_comm                  ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8000 ) << "MultiscaleSolver::MultiscaleSolver() \n";
#endif

    //Define the maps of Multiscale objects
    multiscaleMapsDefinition();

    //Register the available models
    multiscaleModelFactory_Type::instance().registerProduct ( Fluid3D,         &createMultiscaleModelFluid3D );
    multiscaleModelFactory_Type::instance().registerProduct ( FSI3D,           &createMultiscaleModelFSI3D );
    multiscaleModelFactory_Type::instance().registerProduct ( FSI1D,           &createMultiscaleModelFSI1D );
    multiscaleModelFactory_Type::instance().registerProduct ( Multiscale,      &createMultiscaleModelMultiscale );
    multiscaleModelFactory_Type::instance().registerProduct ( Windkessel0D,    &createMultiscaleModelWindkessel0D );
    multiscaleModelFactory_Type::instance().registerProduct ( ZeroDimensional, &createMultiscaleModelZeroDimensional );
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleSolver::setupProblem ( const std::string& fileName, const std::string& problemFolder, const UInt& coresPerNode )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8000 ) << "MultiscaleSolver::setupData( fileName, problemFolder ) \n";
#endif

    // Load data file
    GetPot dataFile ( fileName );

    // Define the folder containing the problem
    multiscaleProblemFolder = problemFolder;

    // Define the number of cores on each node for the machine
    multiscaleCoresPerNode = coresPerNode;

    // Define the step of the problem
    if ( dataFile ( "Solver/Restart/Restart", false ) )
    {
        multiscaleProblemStep = dataFile ( "Solver/Restart/RestartFromStepNumber", 0 ) + 1;
    }

    // Define the filename prefix for the multiscale output
    multiscaleProblemPrefix = dataFile ( "Solver/Output/ProblemPrefix", "Multiscale" );

    // Define how many time step between two calls of the saveSolution() method
    multiscaleSaveEachNTimeSteps = dataFile ( "Solver/Output/SaveEach", 1 );

    // Create the main model and set the communicator
    M_model = multiscaleModelPtr_Type ( multiscaleModelFactory_Type::instance().createObject ( multiscaleModelsMap[ dataFile ( "Problem/ProblemType", "Multiscale" ) ], multiscaleModelsMap ) );

    M_model->setID ( 0 );
    M_model->setCommunicator ( M_comm );

    // Setup data
    M_globalData->readData ( dataFile );
    M_model->setGlobalData ( M_globalData );
    M_model->setupData ( dataFile ( "Problem/ProblemFile", "./MultiscaleDatabase/Models/NoModel" ) + ".dat" );

    // Setup Models
    if ( multiscaleProblemStep )
    {
        importIterationNumber();
    }
    M_model->setupModel();

    // Save the initial solution if it is the very first time step
    if ( !multiscaleProblemStep )
    {
        M_model->saveSolution();
    }

    // Move to the "true" first time-step (needed to perform initializations of the different models/couplings/algorithms)
    M_globalData->dataTime()->updateTime();
    M_globalData->dataTime()->setInitialTime ( M_globalData->dataTime()->time() );
}

bool
MultiscaleSolver::solveProblem ( const Real& referenceSolution, const Real& tolerance )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8000 ) << "MultiscaleSolver::solveProblem() \n";
#endif

    // Chrono definitions
    LifeChrono buildUpdateChrono;
    LifeChrono solveChrono;
    LifeChrono saveChrono;
    LifeChrono updateSolutionChrono;
    LifeChrono globalChrono;
    Real       totalSimulationTime (0);
    Real       timeStepTime (0);

    for ( ; M_globalData->dataTime()->canAdvance(); M_globalData->dataTime()->updateTime() )
    {
        // Global chrono start
        globalChrono.start();

        if ( M_comm->MyPID() == 0 )
        {
            std::cout << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
            std::cout << "                    MULTISCALE FRAMEWORK" << std::endl;
            std::cout << "             time = " << M_globalData->dataTime()->time() << " s;"
                      << " time step number = " << M_globalData->dataTime()->timeStepNumber() << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
        }

        // Build or Update System
        buildUpdateChrono.start();
        if ( M_globalData->dataTime()->isFirstTimeStep() )
        {
            M_model->buildModel();
        }
        else
        {
            M_model->updateModel();
        }
        buildUpdateChrono.stop();

        // Solve the model
        solveChrono.start();
        M_model->solveModel();
        solveChrono.stop();

        // Update solution
        updateSolutionChrono.start();
        M_model->updateSolution();
        updateSolutionChrono.stop();

        // Save the solution
        saveChrono.start();
        if ( M_globalData->dataTime()->timeStepNumber() % multiscaleSaveEachNTimeSteps == 0 || M_globalData->dataTime()->isLastTimeStep() )
        {
            M_model->saveSolution();
        }
        saveChrono.stop();

        // Global chrono stop
        globalChrono.stop();

        // Compute time step time
        timeStepTime = globalChrono.globalDiff ( *M_comm );

        // Updating total simulation time
        totalSimulationTime += timeStepTime;

        if ( M_comm->MyPID() == 0 )
        {
            std::cout << " MS-  Total iteration time:                    " << timeStepTime << " s" << std::endl;
            std::cout << " MS-  Total simulation time:                   " << totalSimulationTime << " s" << std::endl;
        }

        // Save CPU time
        saveCPUTime ( timeStepTime, buildUpdateChrono.globalDiff ( *M_comm ), solveChrono.globalDiff ( *M_comm ),
                      updateSolutionChrono.globalDiff ( *M_comm ), saveChrono.globalDiff ( *M_comm ) );
    }

    // Numerical check of the last iteration solution (used for the night testsuite check)
    Real computedSolution ( M_model->checkSolution() );
    Real relativeError ( std::abs ( ( referenceSolution - computedSolution ) / referenceSolution ) );
    if ( referenceSolution >= 0. && relativeError > tolerance )
        multiscaleErrorCheck ( Solution, "Problem solution: "    + number2string ( computedSolution ) +
                               " (Reference solution: " + number2string ( referenceSolution ) +
                               "; Relative error: " + number2string ( relativeError ) + ")\n", M_comm->MyPID() == 0 );

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
                  << "Problem prefix                = " << multiscaleProblemPrefix << std::endl
                  << "Problem step                  = " << multiscaleProblemStep << std::endl
                  << "Save each                     = " << multiscaleSaveEachNTimeSteps << " time steps" << std::endl << std::endl;

        M_globalData->showMe();

        std::cout << std::endl << std::endl;
    }

    M_model->showMe();

    if ( M_comm->MyPID() == 0 )
    {
        std::cout << "=============================================================" << std::endl << std::endl;
    }
}


// ===================================================
// Private Methods
// ===================================================
void
MultiscaleSolver::saveCPUTime ( const Real& totalCPUTime, const Real& buildUpdateCPUTime,    const Real& solveCPUTime,
                                const Real& updateSolutionCPUTime, const Real& saveCPUTime ) const
{
    if ( M_comm->MyPID() == 0 )
    {
        std::ofstream outputFile;
        outputFile << std::scientific << std::setprecision ( 15 );

        std::string filename = multiscaleProblemFolder + multiscaleProblemPrefix + "_CPUTime_" + number2string ( multiscaleProblemStep ) + ".mfile";

        if ( M_globalData->dataTime()->isFirstTimeStep() )
        {
            outputFile.open ( filename.c_str(), std::ios::trunc );
            outputFile << "% ITERATION                TIME                     TOTAL                    BUILD/UPDATE             "
                       "SOLVE                    UPDATE SOLUTION          SAVE" << std::endl;
        }
        else
        {
            outputFile.open ( filename.c_str(), std::ios::app );
        }
        outputFile << "  " << number2string ( M_globalData->dataTime()->timeStepNumber() )
                   << "                        " << M_globalData->dataTime()->time()
                   << "    " << totalCPUTime << "    " << buildUpdateCPUTime << "    " << solveCPUTime
                   << "    " << updateSolutionCPUTime  << "    " << saveCPUTime << std::endl;
        outputFile.close();
    }
}

void
MultiscaleSolver::importIterationNumber()
{
    // Initialize the iteration number
    Int iterationNumber ( 0 );
    Real initialTime ( 0. );

    if ( M_comm->MyPID() == 0 )
    {
        std::string fileName = multiscaleProblemFolder + multiscaleProblemPrefix + "_CPUTime_" + number2string ( multiscaleProblemStep - 1 ) + ".mfile";

        std::ifstream inputFile;
        inputFile.open ( fileName.c_str(), std::ios::in );

        if ( inputFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;
            std::vector< std::pair< Int, Real > > iterationAndTime;
            std::pair< Int, Real > selectedIterationAndTime;

            // Read the first line with comments
            std::getline ( inputFile, line, '\n' );

            // Read one-by-one all the other lines of the file
            while ( std::getline ( inputFile, line, '\n' ) )
            {
                boost::split ( stringsVector, line, boost::is_any_of ( " " ), boost::token_compress_on );
                if ( string2number ( stringsVector[7] ) > 0 ) // check if we have saved the data
                {
                    iterationAndTime.push_back ( std::make_pair ( string2number ( stringsVector[1] ), string2number ( stringsVector[2] ) ) );
                }
            }

            // Close file
            inputFile.close();

            // Find the closest time step
            selectedIterationAndTime = iterationAndTime.front();
            for ( std::vector< std::pair< Int, Real > >::const_iterator i = iterationAndTime.begin(); i < iterationAndTime.end() ; ++i )
                if ( std::fabs ( selectedIterationAndTime.second - M_globalData->dataTime()->time() ) >= std::fabs ( (*i).second - M_globalData->dataTime()->time() ) )
                {
                    selectedIterationAndTime = *i;
                }

            // Select the iteration number
            iterationNumber = selectedIterationAndTime.first;
            initialTime = selectedIterationAndTime.second;
        }
        else
        {
            std::cerr << " !!! Error: cannot open file: " << fileName.c_str() << " !!!" << std::endl;
        }
    }

    // Share the values with the other processes
    M_comm->Broadcast ( &iterationNumber, 1, 0 );
    M_comm->Broadcast ( &initialTime, 1, 0 );

    // Set the iteration number and the initial time on the basis of the available saved data
    M_globalData->dataTime()->setTimeStepNumber ( iterationNumber );
    M_globalData->dataTime()->setInitialTime ( initialTime );
    M_globalData->dataTime()->setTime ( initialTime );
}

} // Namespace multiscale
} // Namespace LifeV
