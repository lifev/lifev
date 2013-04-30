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
    @brief Example for the DOFInterfaceIO class - write

    @author Radu Popescu <radu.popescu@epfl.ch>
    @maintainer Radu Popescu <radu.popescu@epfl.ch>
    @date 2013-04-30

 */

#include <lifev/core/LifeV.hpp>

#ifdef LIFEV_HAS_HDF5

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "Epetra_config.h"

#ifdef HAVE_MPI

#include <Epetra_MpiComm.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/filter/ExporterHDF5Mesh3D.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/fsi/mesh/MeshPartitionerOfflineFSI.hpp>
#include <lifev/fsi/filter/DOFInterfaceIO.hpp>

using namespace LifeV;

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

int main (int argc, char** argv)
{
#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

    typedef RegionMesh<LinearTetra> mesh_Type;

    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );

    if (comm->NumProc() != 1)
    {
        std::cout << "This test needs to be run "
                  << "with a single process. Aborting."
                  << std::endl;
        return (EXIT_FAILURE);
    }

    GetPot commandLine (argc, argv);
    string dataFileName = commandLine.follow ("data", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

    const UInt numParts (dataFile ("test/num_parts", 4) );
    const std::string fluidPartsFileName (dataFile ("test/fluid_hdf5_file_name",
    												"fluid.h5") );
    const std::string solidPartsFileName (dataFile ("test/solid_hdf5_file_name",
    												"solid.h5") );
    const std::string interfacePartsFileName (dataFile ("test/interface_hdf5_file_name",
    													"solid.h5") );

    std::cout << "Number of parts: " << numParts << std::endl;
    std::cout << "Name of fluid HDF5 container: "
    		  << fluidPartsFileName << std::endl;
    std::cout << "Name of solid HDF5 container: "
    		  << solidPartsFileName << std::endl;
    std::cout << "Name of interface HDF5 container: "
    		  << interfacePartsFileName << std::endl;

    boost::shared_ptr<mesh_Type> fluidMeshPtr (new mesh_Type ( comm ) );
    boost::shared_ptr<mesh_Type> solidMeshPtr (new mesh_Type ( comm ) );

	//Fluid
	MeshData fluidMeshData(dataFile, "fluid_mesh");
	readMesh(*fluidMeshPtr, fluidMeshData);

	//Solid
	MeshData solidMeshData(dataFile, "solid_mesh");
	readMesh(*solidMeshPtr, solidMeshData);

	// Create the FSI partitioner
	MeshPartitionerOfflineFSI<mesh_Type> fsiPartitioner(
		fluidMeshPtr, solidMeshPtr, numParts, numParts, "P1", "P1",
		1, 1, 0, 0,	0, comm);

	// Release the original mesh from the MeshPartitioner object and
	// delete the RegionMesh object
	fluidMeshPtr.reset();
	solidMeshPtr.reset();

	// Write fluid, solid, and interface parts
	boost::shared_ptr<Epetra_MpiComm> mpiComm =
		boost::dynamic_pointer_cast<Epetra_MpiComm>(comm);

	PartitionIO<mesh_Type> fluidPartitionIO (fluidPartsFileName, mpiComm);
	fluidPartitionIO.write (fsiPartitioner.fluidPartitions() );
	PartitionIO<mesh_Type> solidPartitionIO (solidPartsFileName, mpiComm);
	solidPartitionIO.write (fsiPartitioner.solidPartitions() );
	DOFInterfaceIO interfaceIO (interfacePartsFileName, mpiComm);
	interfaceIO.write(fsiPartitioner.dofStructureToHarmonicExtension());

    MPI_Finalize();

#else
    std::cout << "This test needs MPI to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* HAVE_MPI */
#else
    std::cout << "This test needs HDF5 to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* LIFEV_HAS_HDF5 */

    return (EXIT_SUCCESS);
}
