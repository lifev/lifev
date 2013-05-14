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

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "Epetra_config.h"

#ifdef HAVE_MPI

#include <Epetra_MpiComm.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshPartitionTool.hpp>
#include <lifev/core/mesh/GraphCutterParMETIS.hpp>
#include <lifev/core/mesh/MeshPartBuilder.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>

using namespace LifeV;

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

typedef MeshPartitionTool < RegionMesh<LinearTetra>,
        GraphCutterParMETIS,
        MeshPartBuilder > meshCutterParMETIS_Type;

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
        return EXIT_FAILURE;
    }

    GetPot commandLine (argc, argv);
    string dataFileName = commandLine.follow ("data", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

    const UInt numElements (dataFile ("mesh/nelements", 10) );
    const Int numParts (dataFile ("test/num_parts", 3) );
    const bool hierarchical (dataFile ("test/hierarchical", false) );
    const std::string topology (dataFile ("test/topology", "1") );

    const std::string partsFileName (dataFile ("test/hdf5_file_name", "cube.h5") );

    std::cout << "Number of elements in mesh: " << numElements << std::endl;
    std::cout << "Number of parts: " << numParts << std::endl;
    std::cout << "Name of HDF5 container: " << partsFileName << std::endl;

    boost::shared_ptr<mesh_Type> fullMeshPtr (new mesh_Type ( comm ) );
    regularMesh3D (*fullMeshPtr, 1, numElements, numElements, numElements,
                   false, 2.0, 2.0, 2.0, -1.0, -1.0, -1.0);

    Teuchos::ParameterList meshParameters;
    meshParameters.set ("num_parts", numParts, "");
    meshParameters.set ("offline_mode", true, "");
    meshParameters.set ("hierarchical", hierarchical, "");
    meshParameters.set ("topology", topology, "");
    meshCutterParMETIS_Type meshCutter (fullMeshPtr, comm, meshParameters);
    if (! meshCutter.success() )
    {
        std::cout << "Mesh partition failed.";
        return EXIT_FAILURE;
    }

    // delete the RegionMesh object
    fullMeshPtr.reset();

    // Write mesh parts to HDF5 container
    {
        boost::shared_ptr<Epetra_MpiComm> mpiComm =
            boost::dynamic_pointer_cast<Epetra_MpiComm> (comm);
        PartitionIO<mesh_Type> partitionIO (partsFileName, mpiComm);
        partitionIO.write (meshCutter.allMeshParts() );
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
