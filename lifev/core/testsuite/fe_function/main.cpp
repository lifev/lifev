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

Simple test for FEField and FEFunction

*/


/**
   @file main.hpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2011-05-02
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


#include <lifev/core/LifeV.hpp>

#include "fefct.hpp"

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

    MPI_Init ( &argc, &argv );

    std::cout << "MPI Initialization" << std::endl;

#endif

    // Tolerance between the error and the errorKnown
    const LifeV::Real tolerance ( 1e-8 );

    fefct Fefct ( argc, argv );

    // Error of the problem
    const LifeV::Real error = Fefct.run();

#ifdef HAVE_MPI

    MPI_Finalize();

    std::cout << "MPI Finalization" << std::endl;

#endif

    if ( error > tolerance )
    {
        return ( EXIT_FAILURE );
    }
    else
    {
        return ( EXIT_SUCCESS );
    }
}

