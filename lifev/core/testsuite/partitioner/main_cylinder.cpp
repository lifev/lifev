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
    @brief

    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
            Gilles Fourestey <gilles.fourestey@epfl.ch>
    @date 16-04-2005
 */

#ifdef TWODIM
#error test_cylinder cannot be compiled in 2D
#endif

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

#include <life/lifecore/LifeV.hpp>
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>

#ifdef HAVE_HDF5
#include "cylinder.hpp"
#endif

using namespace LifeV;

namespace LifeV
{
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}
}

int
main( int argc, char** argv )
{

#ifdef HAVE_HDF5

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

    int rank;
    int return_value = EXIT_SUCCESS;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//**************** cylinder
//    MPI_Init(&argc,&argv);

    Cylinder cyl( argc, argv  );
    bool success = cyl.run();

    if (rank == 0) {
        if (success) {
            std::cout << "\nTEST PARTITION WAS SUCCESSFUL.\n\n";
        } else {
            std::cout << "\nTEST PARTITION FAILED.\n\n";
            return_value = EXIT_FAILURE;
        }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return return_value;

#endif // HAVE_HDF5

    return 0;

}

