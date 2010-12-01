/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2008-08-11

  Copyright (C) 2008 EPFL

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
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2008-08-11
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/life.hpp>

#include "ensightToHdf5.hpp"
#include <mpi.h>


using namespace LifeV;


int
main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
        cout << "% using MPI" << endl;
#else
    Epetra_SerialComm Comm;
    cout << "% using serial Version" << endl;
#endif

//**************** cylinder
//    MPI_Init(&argc,&argv);

    EnsightToHdf5 es( argc, argv );
    es.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}


