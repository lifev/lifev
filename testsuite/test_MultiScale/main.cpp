//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 * @file
 * @brief MultiScale Test
 *
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 * @date 12-03-2009
 *
 *  This is a very general main file to run a MultiScale simulation.
 *
 *  Models available:
 *  <ol>
 *      <li> Fluid3D (Oseen)
 *      <li> MultiScale
 *  </ol>
 *
 *  Couplings available:
 *  <ol>
 *      <li> Stress
 *      <li> FluxStress
 *  </ol>
 *
 *  Algorithms available:
 *  <ol>
 *      <li> Aitken
 *      <li> Explicit
 *      <li> Newton
 *  </ol>
 *
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/life.hpp>

#include <lifemc/lifesolver/MS_Solver.hpp>

#include <sys/stat.h>
#include <sys/types.h>

using namespace LifeV;

int
main( int argc, char** argv )
{
    //Setup main communicator
    boost::shared_ptr< Epetra_Comm >	comm;

    //Setup MPI variables
    int nprocs(1);
    int rank(0);

#ifdef HAVE_MPI
    MPI_Init( &argc, &argv );

    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif

    if ( rank == 0 )
        std::cout << "MPI Processes: " << nprocs << std::endl;

#ifdef EPETRA_MPI
    if ( rank == 0 )
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
#else
    std::cout << "SERIAL Epetra Initialization ... " << std::endl;
    comm.reset( new Epetra_SerialComm() );
#endif

    // Setup MultiScale problem
    bool exitFlag = EXIT_SUCCESS;
    MS_Solver multiscale;

    // Set the communicator
    multiscale.setCommunicator( comm );

    // Command line parameters
    GetPot commandLine( argc, argv );
    std::string dataFile      = commandLine.follow( "./MultiScale.dat", 2, "-f", "--file" );
    bool verbose              = commandLine.follow( false, 2, "-s", "--showme" );
    std::string problemFolder = commandLine.follow( "MultiScale", 2, "-n", "--name" );
    Real externalResidual     = commandLine.follow( -1., 2, "-c", "--check" );

    // Create the problem folder
    if ( problemFolder.compare("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
            mkdir( problemFolder.c_str(), 0777 );
    }

    // Setup the problem
    multiscale.setupProblem( dataFile, problemFolder );

    // Display problem information
    if ( verbose )
        multiscale.showMe();

    // Solve the problem
    exitFlag = multiscale.solveProblem( externalResidual );

#ifdef HAVE_MPI
    if ( rank == 0 )
        std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    return exitFlag;
}
