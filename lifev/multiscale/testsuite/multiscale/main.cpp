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
 *  @brief File containing the Multiscale Test
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  This is a very general main file to run a Multiscale simulation.
 *
 *  Models available:
 *  <ol>
 *      <li> Fluid3D
 *      <li> FSI3D
 *      <li> OneDimensional
 *      <li> Multiscale
 *      <li> Windkessel0D
 *  </ol>
 *
 *  Couplings available:
 *  <ol>
 *      <li> BoundaryCondition
 *      <li> MeanNormalStress
 *      <li> MeanNormalStressArea
 *      <li> MeanNormalStressValve
 *      <li> MeanTotalNormalStress
 *      <li> MeanTotalNormalStressArea
 *  </ol>
 *
 *  Algorithms available:
 *  <ol>
 *      <li> Aitken
 *      <li> Broyden
 *      <li> Explicit
 *      <li> Newton
 *  </ol>
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

#include <lifev/core/LifeV.hpp>
#include <lifev/multiscale/framework/MultiscaleSolver.hpp>

using namespace LifeV;
using namespace Multiscale;

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
        std::cout << "MPI Processes: " << numberOfProcesses << std::endl;
    }

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

    // Setup Multiscale problem
    bool exitFlag = EXIT_SUCCESS;
    MultiscaleSolver multiscale;

    // Set the communicator
    multiscale.setCommunicator ( comm );

    // Command line parameters
    GetPot commandLine ( argc, argv );
    std::string dataFile      = commandLine.follow ( "./Multiscale.dat", 2, "-f", "--file" );
    bool verbose              = commandLine.follow ( true, 2, "-s", "--showme" );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    Real referenceSolution    = commandLine.follow ( -1., 2, "-c", "--check" );
    UInt coresPerNode         = commandLine.follow (  1, 2, "-ns", "--nodesize" );
    Real tolerance            = commandLine.follow (  1e-8, 2, "-t", "--tolerance" );

    if ( coresPerNode > static_cast<UInt> ( numberOfProcesses ) )
    {
        coresPerNode = numberOfProcesses;
    }

    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    // Setup the problem
    multiscale.setupProblem ( dataFile, problemFolder, coresPerNode );

    // Display problem information
    if ( verbose )
    {
        multiscale.showMe();
    }

    // Solve the problem
    exitFlag = multiscale.solveProblem ( referenceSolution, tolerance );

#ifdef HAVE_MPI
    if ( rank == 0 )
    {
        std::cout << "MPI Finalization" << std::endl;
    }
    MPI_Finalize();
#endif

    return exitFlag;
}
