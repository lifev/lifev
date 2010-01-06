//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
    @file
    @brief A short description of the test content

    @author Name Surname <name.surname@email.org>
    @date 00-00-0000

    Here write a long and detailed description of the test.

    Please note that the test should be quick (maximum 5 minutes)
    and should not be an application.

    Remember to add  and configure a testsuite.at file to allow
    night execution of the test inside the testsuite.
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

LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "Name of the application" ,
                            "Name of the test" ,
                            "Test version 0.0",
                            "Short description",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2009 EPFL");

    about.addAuthor("Name Surname", "Developer", "name.surname@epfl.ch", "");
    about.addAuthor("Name Surname", "Developer", "name.surname@epfl.ch", "");

    return about;
}

using namespace LifeV;

int
main( int argc, char** argv )
{
    //MPI communicator initialization
    boost::shared_ptr<Epetra_Comm> comm;

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
    {
        std::cout << "MPI processes: " << nprocs << std::endl;
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    }
    comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );

    comm->Barrier();

#else

    std::cout << "MPI SERIAL Epetra Initialization ... " << std::endl;
    comm.reset( new Epetra_SerialComm() );

#endif

    // ----- Test calls -----
    //
    // The test should return a result (for tolerance comparison)
    // or a flag (0: EXIT_SUCCESS, 1:EXIT_FAILURE).
    //

    // The test must verify if tolerance is satisfied!
    if ( result > tolerance)
        return EXIT_FAILURE;

    // ----- End of test calls -----

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    // If everything runs correctly we return EXIT_SUCCESS flag
    return EXIT_SUCCESS;
}
