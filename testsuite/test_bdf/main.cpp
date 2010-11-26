/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  Umberto Villa <uvilla@emory.edu>
       Date: 2010-04-14

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
/* ========================================================

  Simple Fourier test with Dirichlet Boundary condition

  Solve the problem

           \partial_t u - \nu(t) \Delta u + sigma(t) u = f(t)

                    u = g on the boundary
                        u(t=0) = u0 initial condition

            on a cube
\nu, \sigma and \source can be function of time
(which implies that the matrix needs to be reassembled each time)

 Purpose: Test BDF of different order
*/


/**
   \file main.cpp
   \author U. Villa <uvilla@emory.edu>
   \date 2010-04-14
 */


// ===================================================
//! Includes
// ===================================================
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
    #include <mpi.h>
	#include <Epetra_MpiComm.h>
#else
	#include <Epetra_SerialComm.h>
#endif

#include <boost/program_options.hpp>

#include <life/lifecore/life.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include "test_bdf.hpp"


// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;
//Register products in the preconditioner factory, so we can use Ifpack and ML as preconditioners.
namespace
{
	static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
	static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}





// ===================================================
//! Main
// ===================================================
int main(int argc, char** argv)
{

	#ifdef HAVE_MPI
		MPI_Init(&argc, &argv);
		std::cout << "MPI Initialization" << std::endl;
	#endif

    test_bdf bdf_t( argc, argv );
    bdf_t.run();


	#ifdef HAVE_MPI
		MPI_Finalize();
		std::cout << "MPI Finalization" << std::endl;
	#endif

    return( EXIT_SUCCESS );
}
