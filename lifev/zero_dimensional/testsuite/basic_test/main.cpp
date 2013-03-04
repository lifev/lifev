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
 *  @brief File containing the ZeroDimensional test
 *
 *  @date 10-02-2012
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <sys/stat.h>
#include <sys/types.h>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/bc_interface/fem/BCInterface0D.hpp>
#include <lifev/zero_dimensional/solver/ZeroDimensionalData.hpp>
#include <lifev/zero_dimensional/solver/ZeroDimensionalSolver.hpp>

using namespace LifeV;

bool checkValue (const Real val, const Real test, const Real tol = 1.e-5, const bool verbose = true)
{
    Real norm = std::abs (val - test);

    if ( verbose )
    {
        std::cout << " value = " << val << " computed value = " << test << " diff = " << norm << std::endl;
    }

    return (norm < tol);
}

Int
main ( Int argc, char** argv )
{
    //Setup main communicator
    boost::shared_ptr< Epetra_Comm > comm;

    //Setup MPI variables
    Int numberOfProcesses (1);
    Int rank (0);

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );

    MPI_Comm_size ( MPI_COMM_WORLD, &numberOfProcesses );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
#endif

    if ( rank == 0 )
    {
        std::cout << std::endl;
        std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
        std::cout << " THE ZERO DIMENSIONAL SOLVER IS AN ALPHA VERSION UNDER STRONG DEVELOPMENT" << std::endl;
        std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;

        std::cout << "MPI Processes: " << numberOfProcesses << std::endl;
    }

#ifdef HAVE_MPI
    if ( numberOfProcesses > 1 )
    {
        if ( rank == 0 )
        {
            std::cout << "test_ZeroDimensional not enabled in parallel, failing gracefully." << std::endl;
            std::cout << "MPI Finalization" << std::endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }
#endif

#ifdef EPETRA_MPI
    if ( rank == 0 )
    {
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    }
    comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    std::cout << "SERIAL Epetra Initialization ... " << std::endl;
    comm.reset ( new Epetra_SerialComm() );
#endif

    bool exitFlag = EXIT_SUCCESS;

#if ( defined(HAVE_NOX_THYRA) && defined(HAVE_TRILINOS_RYTHMOS) )
    // Command line parameters
    GetPot commandLine ( argc, argv );
    const bool check = commandLine.search ( 2, "-c", "--check" );
    string fileName  = commandLine.follow ( "data", 2, "-f", "--file" );

    // SetupData
    GetPot dataFile ( fileName  + ".dat" );

    std::string circuitDataFile = dataFile ( "0D_Model/CircuitDataFile", "./inputFile.dat" );
    BCInterface0D< ZeroDimensionalBCHandler, ZeroDimensionalData >  zeroDimensionalBC;
    zeroDimensionalBC.createHandler();
    zeroDimensionalBC.fillHandler ( circuitDataFile, "Files" );

    boost::shared_ptr< ZeroDimensionalData > zeroDimensionalData ( new ZeroDimensionalData );
    zeroDimensionalData->setup ( dataFile, zeroDimensionalBC.handler() );

    boost::shared_ptr< ZeroDimensionalSolver > zeroDimensionalSolver ( new ZeroDimensionalSolver ( zeroDimensionalData->unknownCounter(), comm, zeroDimensionalData->circuitData() ) );
    zeroDimensionalSolver->setup ( zeroDimensionalData->solverData() );

    zeroDimensionalData->showMe();

    // SetupModel
    zeroDimensionalData->dataTime()->setInitialTime (0);
    zeroDimensionalData->initializeSolution();

    // Create output folder
    if ( comm->MyPID() == 0 )
    {
        mkdir ( "output", 0777 );
    }

    // Save initial solution
    zeroDimensionalData->saveSolution();

    zeroDimensionalData->dataTime()->updateTime();
    zeroDimensionalData->dataTime()->setInitialTime (zeroDimensionalData->dataTime()->time() );

    // Definitions for the time loop
    LifeChrono chronoTotal;
    LifeChrono chronoSystem;
    LifeChrono chronoIteration;

    Int count = 0;
    chronoTotal.start();

    for ( ; zeroDimensionalData->dataTime()->canAdvance() ; zeroDimensionalData->dataTime()->updateTime(), ++count )
    {
        std::cout << std::endl << "--------- Iteration " << count << " time = " << zeroDimensionalData->dataTime()->time() << std::endl;

        chronoIteration.start();
        chronoSystem.start();

        zeroDimensionalSolver->takeStep ( zeroDimensionalData->dataTime()->previousTime(), zeroDimensionalData->dataTime()->time() );

        chronoSystem.stop();

        //Save solution
        zeroDimensionalData->saveSolution();

        chronoIteration.stop();

        std::cout << " System solved in " << chronoSystem.diff() << " s, (total time " << chronoIteration.diff() << " s)." << std::endl;
    }

    chronoTotal.stop();
    std::cout << std::endl << " Simulation ended successfully in " << chronoTotal.diff()  << " s" << std::endl;

    // Final check
    if ( check )
    {
        bool ok = true;

        ok = ok && checkValue ( 0.001329039627, zeroDimensionalData->circuitData()->Nodes()->nodeListAt (1)->voltage() );
        ok = ok && checkValue ( 0.000787475119, zeroDimensionalData->circuitData()->Elements()->elementListAt (1)->current() );
        if (ok)
        {
            std::cout << " Test succesful" << std::endl;
            exitFlag = EXIT_SUCCESS;
        }
        else
        {
            std::cout << " Test unsuccesful" << std::endl;
            exitFlag = EXIT_FAILURE;
        }
    }
#else
    std::cout << "ZeroDimensional requires configuring Trilinos with Rythmos/NOX/Thyra. Skipping test." << std::endl;
    exitFlag = EXIT_SUCCESS;
#endif /* HAVE_NOX_THYRA && HAVE_TRILINOS_RYTHMOS */

    if (rank == 0)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

#ifdef HAVE_MPI
    if ( rank == 0 )
    {
        std::cout << "MPI Finalization" << std::endl;
    }
    MPI_Finalize();
#endif

    return exitFlag;
}
