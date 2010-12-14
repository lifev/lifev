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
    @brief Simple Fourier test with Dirichlet Boundary condition

    @author Umberto Villa <uvilla@emory.edu>
    @date 14-04-2010

      Solve the problem

           \partial_t u - \nu(t) \Delta u + sigma(t) u = f(t)

                    u = g on the boundary
                        u(t=0) = u0 initial condition

            on a cube
       \nu, \sigma and \source can be function of time
      (which implies that the matrix needs to be reassembled each time)

     Purpose: Test BDF of different order
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
    bool check(true);
    check = bdf_t.check();


#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout << "MPI Finalization" << std::endl;
#endif

    if(check)
    	return( EXIT_SUCCESS );
    else
    	return( EXIT_FAILURE );
}
