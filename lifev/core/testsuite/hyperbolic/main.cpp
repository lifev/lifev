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
/* ========================================================

Simple hyperbolic test with Dirichlet, Neumann and Robin boundary conditions

*/


/**
   @file
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @author M. Kern <michel.kern@inria.fr>
   @date 2010-07-29

*/


// ===================================================
//! Includes
// ===================================================

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include "hyperbolic.hpp"



// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;

// ===================================================
//! Main
// ===================================================
int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::cout << "MPI Initialization" << std::endl;
#endif


    // Error of the problem
    LifeV::Real error (0);
    // Error known
    const LifeV::Real errorKnown ( 0.357740920454703 );
    // Tollerance between the error and the errorKnown
    const LifeV::Real tolerance ( 1e-8 );

    hyperbolic Hyperbolic ( argc, argv );

    // Error of the problem
    error = Hyperbolic.run();

    bool success = std::fabs ( error - errorKnown ) <= tolerance;

#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout << "MPI Finalization" << std::endl;
#endif

    // For tribits handling of success/failure
    //! @todo Add verbose to avoid all processes printing this stuff
    if (!success)
    {
        std::cout << "End Result: TEST NOT PASSED" << std::endl;
        return ( EXIT_FAILURE );
    }
    else
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

    return ( EXIT_SUCCESS );
}
