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
 *      <li> Aitken (and fixpoint)
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

#include <life/lifecore/application.hpp>
#include <life/lifecore/life.hpp>

#include <lifemc/lifesolver/MS_Solver.hpp>

#include <sys/stat.h>
#include <sys/types.h>

LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "MultiScale" ,
                            "LifeV MultiScale Tester" ,
                            "0.1",
                            "Under development",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2009 EPFL");

    about.addAuthor("Cristiano Malossi", "Developer", "cristiano.malossi@epfl.ch", "");
    return about;
}

using namespace LifeV;

int
main( int argc, char** argv )
{
    //Setup main communicator
    boost::shared_ptr<Epetra_Comm>	comm;

#ifdef HAVE_MPI
        std::cout << "MPI Initialization" << std::endl;
        MPI_Init( &argc, &argv );
#endif

	//MPI Preprocessing
#ifdef EPETRA_MPI
	    int nprocs;
        int rank;

        MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );

        if ( rank == 0 )
        {    std::cout << "MPI processes: " << nprocs << std::endl;

            std::cout << "MPI Epetra Initialization ... " << std::endl;
        }
        comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );

        comm->Barrier();
#else
        std::cout << "MPI SERIAL Epetra Initialization ... " << std::endl;
        comm.reset( new Epetra_SerialComm() );
#endif

    // Setup MultiScale problem
    bool exitFlag = EXIT_SUCCESS;
    MS_Solver MS;

    // Set the communicator
    MS.SetCommunicator( comm );

    // Command line parameters
    GetPot commandLine( argc, argv );
    std::string dataFile    = commandLine.follow( "./MultiScale.dat", 2, "-f", "--file" );
    bool verbose            = commandLine.follow( false, 2, "-s", "--showme" );
    std::string problemFolder = commandLine.follow( "MultiScale", 2, "-n", "--name" );

    // Create the problem folder
    if ( problemFolder.compare("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
            mkdir( problemFolder.c_str(), 0777 );
    }

    // Setup the problem
    MS.SetupProblem( dataFile, problemFolder );

    // Display problem information
    if ( verbose )
        MS.ShowMe();

    // Solve the problem
    exitFlag = MS.SolveProblem();

#ifdef HAVE_MPI
        std::cout << "MPI Finalization" << std::endl;
        MPI_Finalize();
#endif

    return exitFlag;
}
