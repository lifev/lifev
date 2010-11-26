/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author: Radu Popescu <radu.popescu@epfl.ch>
  Copyright (C) 2010 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file main.cpp
   \author Radu Popescu <radu.popescu@epfl.ch>
   \date 2010-08-24
 */

/*
  Test of the offline partitioning for FSI simulations
*/

#include <iostream>
#include <string>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <boost/shared_ptr.hpp>
#include <mpi.h>

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/FSIOfflinePartitioner.hpp>
#include <life/lifefilters/HDF5Filter3DMesh.hpp>

using namespace LifeV;

typedef RegionMesh3D<LinearTetra> Mesh;

int main( int argc, char** argv )
{
    boost::shared_ptr<Epetra_Comm> comm;

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    comm.reset(new Epetra_MpiComm(MPI_COMM_WORLD) );
    std::cout << "% using MPI version" << std::endl;
#else
    comm.reset( new Epetra_SerialComm() );
    std::cout << "% using serial version" << std::end;
#endif

    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    int fluidInterfaceFlag = dataFile("interface/fluid_flag", 1);
    int solidInterfaceFlag = dataFile("interface/structure_flag", fluidInterfaceFlag);
    Real interfaceTolerance = dataFile("interface/tolerance", 0.0);

    int fluidInterfaceVertexFlag;
    int solidInterfaceVertexFlag;
    int vertexFlag;
    vertexFlag = dataFile("interface/edgeFlag",      -1);
    vertexFlag = dataFile("interface/fluid_vertex_flag", vertexFlag);

    fluidInterfaceVertexFlag = vertexFlag;

    solidInterfaceVertexFlag = dataFile("interface/structure_vertex_flag", -1);

    std::string fluidOrder(dataFile("fluid/space_discretization/vel_order", "P1"));
    std::string solidOrder(dataFile("solid/space_discretization/order", "P1"));

    boost::shared_ptr<DataMesh> fluidDataMesh(new DataMesh);
    fluidDataMesh->setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<Mesh> uncutFluidMesh(new Mesh);
    readMesh(*uncutFluidMesh, *fluidDataMesh);

    boost::shared_ptr<DataMesh> solidDataMesh(new DataMesh);
    solidDataMesh->setup(dataFile, "solid/space_discretization");

    boost::shared_ptr<Mesh> uncutSolidMesh(new Mesh);
    readMesh(*uncutSolidMesh, *solidDataMesh);

    boost::shared_ptr<FSIOfflinePartitioner<Mesh> >
        cutter(new FSIOfflinePartitioner<Mesh>);
    cutter->setup(uncutFluidMesh, uncutSolidMesh, 4, 4, fluidOrder, solidOrder,
                  fluidInterfaceFlag, solidInterfaceFlag, interfaceTolerance,
                  fluidInterfaceVertexFlag, solidInterfaceVertexFlag, comm);
    cutter->showMe();

    cutter->execute();

    HDF5Filter3DMesh<Mesh> fluidOutput(dataFile, uncutFluidMesh, "FSIFluidPartitions",
                                       comm->MyPID());
    fluidOutput.addPartitionGraph(cutter->fluidGraph(), comm);
    fluidOutput.addMeshPartitionAll(cutter->fluidPartitions(), comm);
    fluidOutput.addDOFInterface(cutter->dofStructureToHarmonicExtension(),
                                "dofStructureToHarmonicExtension",
                                cutter->fluidInterfaceFlag(),
                                cutter->solidInterfaceFlag(),
                                comm);
    fluidOutput.postProcess(0);
    fluidOutput.CloseFile();

    HDF5Filter3DMesh<Mesh> solidOutput(dataFile, uncutSolidMesh, "FSISolidPartitions",
                                       comm->MyPID());
    solidOutput.addPartitionGraph(cutter->solidGraph(), comm);
    solidOutput.addMeshPartitionAll(cutter->solidPartitions(), comm);
    solidOutput.addDOFInterface(cutter->dofStructureToHarmonicExtension(),
                                "dofStructureToHarmonicExtension",
                                cutter->fluidInterfaceFlag(),
                                cutter->solidInterfaceFlag(),
                                comm);
    solidOutput.postProcess(0);
    solidOutput.CloseFile();


#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return(EXIT_SUCCESS);
}

