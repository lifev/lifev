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

#ifdef LIFEV_HAS_HDF5


#include "Epetra_config.h"

#ifdef HAVE_MPI

#include <Epetra_MpiComm.h>


#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshPartitionTool.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>

using namespace LifeV;

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

typedef MeshPartitionTool < RegionMesh<LinearTetra> > meshCutter_Type;

int main (int argc, char** argv)
{
#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

    typedef RegionMesh<LinearTetra>           mesh_Type;
    typedef std::shared_ptr<mesh_Type>      meshPtr_Type;
    typedef std::vector<meshPtr_Type>         meshParts_Type;
    typedef std::shared_ptr<meshParts_Type> meshPartsPtr_Type;

    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );

    if (comm->NumProc() != 1)
    {
        std::cout << "This test needs to be run "
                  << "with a single process. Aborting."
                  << std::endl;
        return EXIT_FAILURE;
    }

    GetPot cl (argc, argv);
    const UInt numElements = cl.follow (9, "--num-elem");
    const Int numParts = cl.follow (3, "--num-parts");
    // partitionerType should be MeshPartitioner, MeshPartitionTool_ParMETIS or
    // MeshPartitionTool_Zoltan
    const std::string partitionerType = cl.follow ("MeshPartitioner",
                                                   "--partitioner-type");
    std::string partsFile;
    partsFile.reserve (50);
    partsFile += "cube_";
    partsFile += partitionerType;
    partsFile += ".h5";

    std::cout << "Number of elements in mesh: " << numElements << std::endl;
    std::cout << "Number of parts: " << numParts << std::endl;
    std::cout << "Mesh partitioner type: " << partitionerType << std::endl;
    std::cout << "Name of HDF5 container: " << partsFile << std::endl;

    meshPtr_Type fullMeshPtr (new mesh_Type ( comm ) );
    meshPartsPtr_Type meshPartPtr;
    regularMesh3D (*fullMeshPtr, 1, numElements, numElements, numElements,
                   false, 2.0, 2.0, 2.0, -1.0, -1.0, -1.0);

    if (partitionerType == "MeshPartitioner")
    {
        // Using old MeshPartitioner class
        MeshPartitioner<mesh_Type> meshCutter;
        meshCutter.setup (numParts, comm);
        meshCutter.attachUnpartitionedMesh (fullMeshPtr);
        meshCutter.doPartitionGraph();
        meshCutter.fillEntityPID();
        meshCutter.doPartitionMesh();
        meshCutter.releaseUnpartitionedMesh();
        meshPartPtr = meshCutter.meshPartitions();
    }
    else if (partitionerType == "MeshPartitionTool_ParMETIS")
    {
        // Using new MeshPartitionTool class with ParMETIS
        Teuchos::ParameterList meshParameters;
        meshParameters.set ("num-parts", numParts, "");
        meshParameters.set ("offline-mode", true, "");
        meshParameters.set ("graph-lib", "parmetis", "");
        meshCutter_Type meshCutter (fullMeshPtr, comm, meshParameters);
        if (! meshCutter.success() )
        {
            std::cout << "Mesh partition failed.";
            return EXIT_FAILURE;
        }
        meshPartPtr = meshCutter.allMeshParts();
    }
    else if (partitionerType == "MeshPartitionTool_Zoltan")
    {
        // Using new MeshPartitionTool class with Zoltan
        Teuchos::ParameterList meshParameters;
        meshParameters.set ("num-parts", numParts, "");
        meshParameters.set ("offline-mode", true, "");
        meshParameters.set ("graph-lib", "zoltan", "");
        meshCutter_Type meshCutter (fullMeshPtr, comm, meshParameters);
        if (! meshCutter.success() )
        {
            std::cout << "Mesh partition failed.";
            return EXIT_FAILURE;
        }
        meshPartPtr = meshCutter.allMeshParts();
    }

    // delete the RegionMesh object
    fullMeshPtr.reset();

    // Write mesh parts to HDF5 container
    {
        std::shared_ptr<Epetra_MpiComm> mpiComm =
            std::dynamic_pointer_cast<Epetra_MpiComm> (comm);
        PartitionIO<mesh_Type> partitionIO (partsFile, mpiComm);
        partitionIO.write (meshPartPtr);
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

    return EXIT_SUCCESS;
}
