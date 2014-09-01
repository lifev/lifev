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

Simple multirate 2d test

*/


/**
   @file
   @author L. Oldani <luca.oldani@mail.polimi.it>
   @date 2013-12-27
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

#include "multirate.hpp"

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
    MPI_Init(&argc, &argv);
#endif

    // Error known
    const LifeV::Real errorKnown ( 0.0366812221943312 );

    // Tollerance between the error and the errorKnown
    const LifeV::Real tolerance ( 1e-10 );

    multirate_2d Multirate ( argc, argv );

    // Error of the problem
    const LifeV::Real error = Multirate.run();
    const bool unsuccess = std::fabs( error - errorKnown ) > tolerance;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    // For tribits handling of success/failure
    //! @todo Add verbose to avoid all processes printing this stuff
    if ( unsuccess )
    {
        return ( EXIT_FAILURE );
    }
    else
    {
        return ( EXIT_SUCCESS );
    }
}
