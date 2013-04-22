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
    @brief Test for PartitionIO class - cut and write

    @author Radu Popescu <radu.popescu@epfl.ch>
    @maintainer Radu Popescu <radu.popescu@epfl.ch>
    @date 10-05-2012

    Partition a mesh using a single (MPI) process and save mesh parts
    to an HDF5 file.
 */

#include <lifev/core/LifeV.hpp>

#include <iostream>
#include <string>

#include "Epetra_config.h"

#ifdef HAVE_HDF5
#ifdef HAVE_MPI

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <mpi.h>

#include <Epetra_MpiComm.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

using namespace LifeV;

#endif /* HAVE_MPI */
#endif /* HAVE_HDF5 */

int main (int argc, char** argv)
{
#ifdef HAVE_HDF5
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

    const UInt numParts (dataFile ("offlinePartition/num_parts", 4) );
    std::string stringFileName (dataFile ("offlinePartition/nameFile", "NO_DEFAULT_VALE") );

    stringFileName += "Partitioned.h5";

    std::cout << "Number of parts: " << numParts << std::endl;
    std::cout << "Name of HDF5 container: " << stringFileName << std::endl;

    //Creation mesh data object to read the whole mesh
    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    //Creation of the object whole mesh
    boost::shared_ptr<mesh_Type> fullMeshPtr (new mesh_Type ( comm ) );
    readMesh (*fullMeshPtr, meshData);

    fullMeshPtr->showMe();

    //Creation object mesh partitioner
    MeshPartitioner<mesh_Type> meshPart;
    meshPart.setup (numParts, comm);

    meshPart.attachUnpartitionedMesh (fullMeshPtr);
    meshPart.doPartitionGraph();
    meshPart.doPartitionMesh();

    // Release the original mesh from the MeshPartitioner object and
    // delete the RegionMesh object
    meshPart.releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    // Write mesh parts to HDF5 container
    PartitionIO<mesh_Type> partitionIO (stringFileName, comm);
    partitionIO.write (meshPart.meshPartitions() );

    MPI_Finalize();

#else
    std::cout << "This test needs MPI to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* HAVE_MPI */
#else
    std::cout << "This test needs HDF5 to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* HAVE_HDF5 */

    return (EXIT_SUCCESS);
}
