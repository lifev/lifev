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
    @file
    @brief test for converting from VTK to Ensight

    @author Tiziano Passerini <tiziano@mathcs.emory.edu>
    @contributor
    @maintainer

    @date 13-1-2011

 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#include <Epetra_Comm.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include "../importExport/TestImportExport.hpp"

using namespace LifeV;

int
main ( int argc, char** argv )
{
    //MPI communicator initialization
    boost::shared_ptr<Epetra_Comm> commPtr;

#ifdef HAVE_MPI
    std::cout << "MPI Initialization" << std::endl;
    MPI_Init ( &argc, &argv );
#endif

    //MPI Preprocessing
#ifdef EPETRA_MPI

    int nprocs;
    int rank;

    MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );

    if ( rank == 0 )
    {
        std::cout << "MPI processes: " << nprocs << std::endl;
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    }
    commPtr.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );

    commPtr->Barrier();

#else

    std::cout << "MPI SERIAL Epetra Initialization ... " << std::endl;
    commPtr.reset ( new Epetra_SerialComm() );

#endif

    GetPot command_line (argc, argv);
    TestImportExport testImportExport ( commPtr );

    bool passed (false);

    typedef ExporterVTK<mesh_Type> exporter1_Type;
    typedef ExporterEnsight<mesh_Type> exporter2_Type;
    passed = testImportExport.run<exporter1_Type, exporter2_Type > ( command_line, "import" );

    // ----- End of test calls -----

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    if (passed)
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}
