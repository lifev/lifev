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
/**
   \file test_ghosthandler.cpp
   \author D. Forti <davide.forti@epfl.ch>
   \date 2012-28-12
 */


// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/GhostHandler.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/GetPot.hpp>

#include <lifev/core/filter/Exporter.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

using namespace LifeV;
using namespace std;

int main ( int argc, char* argv[] )
{
    boost::shared_ptr<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    Comm.reset (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    comm.reset ( new Epetra_SerialComm() );
#endif

    GetPot command_line (argc, argv);
    GetPot dataFile ( command_line.follow ("data_neighborsCircle", 2, "-f", "--file" ) );

    typedef LinearTriangle                     geoElement_Type;
    typedef RegionMesh < geoElement_Type >     mesh_Type;
    typedef boost::shared_ptr < mesh_Type >    meshPtr_Type;
    typedef VectorEpetra                       vector_Type;
    typedef boost::shared_ptr<vector_Type>     vectorPtr_Type;

    MeshData meshData;
    meshData.setup (dataFile, "space_discretization");

    meshPtr_Type fullMeshPtr ( new mesh_Type ( Comm ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner<mesh_Type> meshPart;
    meshPtr_Type localMeshPtr;

    // Partitioning the mesh with a number of overlapping regions equal to leveloverlap
    int levelOverlap = 0;
    meshPart.setPartitionOverlap (levelOverlap);
    meshPart.doPartition (fullMeshPtr, Comm);
    localMeshPtr = meshPart.meshPartition();

    // Creating and setting up a GhostHandler object
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > FESpaceP1 (new FESpace<mesh_Type, MapEpetra> (localMeshPtr, "P1", 1, Comm) );
    GhostHandler<mesh_Type> ghostObj ( fullMeshPtr, localMeshPtr, FESpaceP1->mapPtr(), Comm );
    ghostObj.setUpNeighbors();

    // Remark: by calling the setup method, such a class identifies the closest neighbors of grid each node, as it is shown below:
    //
    //                     n6
    // o ------ o ------ o ------ o ------ o
    // |        |        |        |        |       LEGEND: o   node
    // |        |        |        |        |               --- connectivity between nodes
    // |        | n13    | n2     | n7     |
    // o ------ o ------ o ------ o ------ o
    // |        |        |        |        |
    // |        |        |        |        |
    // | n12    | n5     | n1     | n3     | n8
    // o ------ o ------ o ------ o ------ o
    // |        |        |        |        |
    // |        |        |        |        |
    // |        | n11    | n4     | n9     |
    // o ------ o ------ o ------ o ------ o
    // |        |        |        |        |
    // |        |        |        |        |
    // |        |        | n10    |        |
    // o ------ o ------ o ------ o ------ o
    //
    // As an example, concerning node n1, the setUp() will identify as neighbors nodes n2, n3, n4 and n5.
    // The aim of this test is to extend the possibility of finding all the neighbors which are within
    // a user-defined number of circles. Regarding the example above, if one chooses as number of circles equal to 1,
    // we identify as additional neighbors the nodes n6, n7, n8, n9, n10, n11, n12 and n13.

    UInt ID_trial = 40;
    neighbors_Type Neighbors;
    UInt nc = 2;

    Neighbors = ghostObj.circleNeighbors ( fullMeshPtr->point (ID_trial).id(), nc );
    Neighbors.insert (fullMeshPtr->point (ID_trial).id() );

    // createCircleNodeNodeNeighborsMap takes both the number of circles (nc) where to find the neighbors, and the
    // globalID of the grid node localMeshPtr->point(i).id()() for which we are looking for its neighbors.
    // As output the method returns a set Neighbors[i] containing the globalID of each neighbor.

    // EXPORTING THE RESULT FOR ONE NODE (with globalID: ID_trial), IN ORDER TO VISUALIZE THE RESULT WITH PARAVIEW
    vectorPtr_Type TrialOutput (new vector_Type (FESpaceP1->map(), Unique) );

    for (neighbors_Type::iterator ii = Neighbors.begin(); ii != Neighbors.end(); ++ii)
        if (TrialOutput->blockMap().LID (static_cast<int> (*ii) ) != -1)
        {
            if (*ii == ID_trial)
            {
                (*TrialOutput) [*ii] = -1;
            }
            else
            {
                (*TrialOutput) [*ii] =  1;
            }
        }

    ExporterHDF5<mesh_Type> exporter (dataFile, localMeshPtr, "Output_test_neighborsCircle", Comm->MyPID() );
    exporter.setMeshProcId (localMeshPtr, Comm->MyPID() );
    exporter.exportPID (localMeshPtr, Comm, true );
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Neighbors, red color", FESpaceP1, TrialOutput, UInt (0) );
    exporter.postProcess (0);
    exporter.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;

}
