/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-03-12

  Copyright (C) 2009 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file main.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-03-12
 */

#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
	#include <mpi.h>
    #include "Epetra_MpiComm.h"
#else
    #include "Epetra_SerialComm.h"
#endif

#include <life/lifecore/application.hpp>
#include <life/lifecore/life.hpp>

#include <lifemc/lifesolver/MS_Solver.hpp>



LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "MultiScale" ,
                            "LifeV MultiScale Assembler" ,
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

	int	nprocs;
	int rank;

	#ifdef HAVE_MPI
		std::cout << "MPI Initialization" << std::endl;
		MPI_Init( &argc, &argv );
	#endif

	//MPI Preprocessing
	#ifdef EPETRA_MPI

		MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
		MPI_Comm_rank( MPI_COMM_WORLD, &rank );

		if ( rank == 0 )
		{	std::cout << "MPI processes: " << nprocs << std::endl;

			std::cout << "MPI Epetra Initialization ... " << std::endl;
		}
		comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );

		comm->Barrier();

	#else

		std::cout << "MPI SERIAL Epetra Initialization ... " << std::endl;
		comm.reset( new Epetra_SerialComm() );

	#endif

	// Setup MultiScale problem
	MS_Solver MS;

	//Command line parameters
	GetPot commandLine( argc, argv );
	std::string dataFile = commandLine.follow( "./MultiScale.dat", 2, "-f", "--file" );
	bool verbose = commandLine.follow( false, 2, "-s", "--showme" );

	MS.SetCommunicator( comm );

	MS.SetupProblem( dataFile );

	if ( verbose )
		MS.ShowMe();

	MS.SolveProblem();

	#ifdef HAVE_MPI
		std::cout << "MPI Finalization" << std::endl;
		MPI_Finalize();
	#endif

	return( EXIT_SUCCESS );
}
