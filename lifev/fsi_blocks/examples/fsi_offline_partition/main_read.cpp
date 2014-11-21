//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010, 2011, 2012 EPFL, Politecnico di Milano, Emory University

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
    @brief Example for DOFInterfaceIO class - read

    @author Radu Popescu <radu.popescu@epfl.ch>
    @maintainer Radu Popescu <radu.popescu@epfl.ch>

    @date 2013-04-30

 */

#include <lifev/core/LifeV.hpp>

#include <map>

#ifdef LIFEV_HAS_HDF5

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "Epetra_config.h"

#ifdef HAVE_MPI

#include <mpi.h>

#include <Epetra_MpiComm.h>

#include <boost/shared_ptr.hpp>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/fsi_blocks/filter/DOFInterfaceIO.hpp>


using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

int
main ( int argc, char** argv )
{
#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );

    const bool verbose (comm->MyPID() == 0);

    // Read first the data needed

    if (verbose)
    {
        std::cout << " -- Reading the data ... " << std::flush;
    }
    GetPot dataFile ( "data" );

    const std::string fluidHdf5File (dataFile ("test/fluid_hdf5_file_name",
                                               "fluid.h5") );
    const std::string solidHdf5File (dataFile ("test/solid_hdf5_file_name",
                                               "fluid.h5") );
    const std::string interfaceHdf5File (dataFile ("test/interface_hdf5_file_name",
                                                   "interface.h5") );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    boost::shared_ptr<mesh_Type> fluidMesh;
    boost::shared_ptr<mesh_Type> solidMesh;
    boost::shared_ptr<std::map<UInt, UInt> > interfaceMap;

    {
        // Load fluid mesh part from HDF5
        if (verbose)
        {
            std::cout << " -- Reading the fluid mesh part ... " << std::endl;
            std::cout << fluidHdf5File << std::endl;
        }
        PartitionIO<RegionMesh<LinearTetra> > partitionIO (fluidHdf5File, comm);
        partitionIO.read (fluidMesh);

        fluidMesh->showMe();

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }
    }
    {
        // Load solid mesh part from HDF5
        if (verbose)
        {
            std::cout << " -- Reading the solid mesh part ... " << std::endl;
            std::cout << solidHdf5File << std::endl;
        }
        PartitionIO<RegionMesh<LinearTetra> > partitionIO (solidHdf5File, comm);
        partitionIO.read (solidMesh);

        solidMesh->showMe();

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }
    }
    {
        // Load interface from HDF5
        if (verbose)
        {
            std::cout << " -- Reading the interface ... " << std::endl;
            std::cout << fluidHdf5File << std::endl;
        }
        DOFInterfaceIO interfaceIO (interfaceHdf5File, comm);
        interfaceIO.read (interfaceMap);

        std::cout << "Interface " << comm->MyPID()
                  << " size: " << interfaceMap->size() << std::endl;

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }
    }

    MPI_Finalize();

#else
    std::cout << "This test needs MPI to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* HAVE_MPI */
#else
    std::cout << "This test needs HDF5 to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* LIFEV_HAS_HDF5 */

    return ( EXIT_SUCCESS );
}


