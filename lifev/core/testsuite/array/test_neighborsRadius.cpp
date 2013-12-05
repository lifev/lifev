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

int main ( int argc, char* argv[] )
{

    typedef LinearTetra                        geoElement_Type;
    typedef RegionMesh < geoElement_Type >     mesh_Type;
    typedef boost::shared_ptr < mesh_Type >    meshPtr_Type;
    typedef VectorEpetra                       vector_Type;
    typedef boost::shared_ptr<vector_Type>     vectorPtr_Type;

    boost::shared_ptr<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    Comm.reset (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    comm.reset ( new Epetra_SerialComm() );
#endif

    // Loading data
    GetPot command_line (argc, argv);
    GetPot dataFile ( command_line.follow ("data_neighborsRadius", 2, "-f", "--file" ) );

    // Loading mesh
    MeshData meshData;
    meshData.setup (dataFile, "space_discretization");

    meshPtr_Type fullMeshPtr ( new mesh_Type ( Comm ) );
    readMesh (*fullMeshPtr, meshData);

    // Mesh partitioning
    MeshPartitioner<mesh_Type> meshPart;
    meshPtr_Type localMeshPtr;

    // Partitioning the mesh with a number of overlapping regions equal to leveloverlap
    int levelOverlap = 5;
    meshPart.setPartitionOverlap (levelOverlap);
    meshPart.doPartition (fullMeshPtr, Comm);
    localMeshPtr = meshPart.meshPartition();

    // Creating and setting up a GhostHandler object
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > FESpaceP1 (new FESpace<mesh_Type, MapEpetra> (localMeshPtr, "P1", 1, Comm) );
    GhostHandler<mesh_Type> ghostObj ( fullMeshPtr, localMeshPtr, FESpaceP1->mapPtr(), Comm );

    // Creating node-node map over the interface
    std::vector<markerID_Type> interfaceMarkers (2);
    interfaceMarkers[0] = 20;
    interfaceMarkers[1] = 1;
    ghostObj.setUp ( interfaceMarkers );

    UInt ID_trial = 314;
    std::set<ID> Neighbors;
    double radius = 4 * (double) MeshUtility::MeshStatistics::computeSize (*fullMeshPtr).maxH;

    Neighbors = ghostObj.neighborsWithinRadius ( fullMeshPtr->point (ID_trial).id(), radius );
    Neighbors.insert (fullMeshPtr->point (ID_trial).id() );

    // Exporting the solution
    vectorPtr_Type TrialOutput (new vector_Type (FESpaceP1->map(), Unique) );

    for (std::set<ID>::iterator ii = Neighbors.begin(); ii != Neighbors.end(); ++ii)
        if (TrialOutput->blockMap().LID (*ii) != -1)
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

    ExporterHDF5<mesh_Type> exporter (dataFile, localMeshPtr, "Output_test_neighborsRadius", Comm->MyPID() );
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
