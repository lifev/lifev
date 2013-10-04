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

Simple Darcy test with Dirichlet, Neumann and Robin boundary conditions

*/


/**
   @file main.hpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2012-06-13
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

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include "darcy.hpp"


// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}

// ===================================================
//! Main
// ===================================================
int main (int argc, char* argv[])
{

#ifdef HAVE_MPI

    MPI_Init ( &argc, &argv );

#endif

    // Error known
    const LifeV::Real errorKnown ( 0.9554685918458008 );

    // Tolerance between the error and the error known
    const LifeV::Real tolerance ( 1e-10 );

    darcy_nonlinear Darcy ( argc, argv );

    // Error of the problem
    const LifeV::Real error = Darcy.run();
    const bool unsuccess = std::fabs ( error - errorKnown ) > tolerance;

#ifdef HAVE_MPI

    MPI_Finalize();

#endif

    if ( unsuccess )
    {
        return ( EXIT_FAILURE );
    }
    else
    {
        return ( EXIT_SUCCESS );
    }
}
