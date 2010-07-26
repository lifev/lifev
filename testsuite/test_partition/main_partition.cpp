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
   \file main_partition.cpp
   \author Radu Popescu <radu.popescu@epfl.ch>
   \date 2010-07-02
 */

/*
    Test of the serial mesh partitioning - part 1
    Partition a mesh using a single (MPI) process and save mesh partitions
    to a HDF5 file.
*/

#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefilters/HDF5Filter3DMesh.hpp>

#include <iostream>
#include <string>
#include <mpi.h>

LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "life_partitioning" ,
                            "life_partitioning" ,
                            "0.1",
                            "serial mesh partitioning test",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2010 EPFL");

    about.addAuthor("Radu Popescu", "developer", "radu.popescu@epfl.ch", "");
    return about;

}

using namespace LifeV;

int main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    std::cout << "% using MPI version" << std::endl;
#else
    Epetra_SerialComm comm;
    std::cout << "% using serial version" << std::end;
#endif

    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    DataNavierStokes dataNavierStokes;
    dataNavierStokes.setup(dataFile);

    DataMesh dataMesh;
    dataMesh.setup(dataFile, "fluid/space_discretization");

    RegionMesh3D<LinearTetra> mesh;
    readMesh(mesh, dataMesh);

    partitionMesh<RegionMesh3D<LinearTetra> > meshPart;
    meshPart.setup(4, comm);

    meshPart.attachUnpartitionedMesh(mesh);
    meshPart.doPartitionGraph();
    meshPart.doPartitionMesh();

    HDF5Filter3DMesh<RegionMesh3D<LinearTetra> > HDF5Output(dataFile, meshPart.mesh(), "cylinderPart",
                                                            comm.MyPID());
    HDF5Output.addPartitionGraph(meshPart.graph(), &comm);
    HDF5Output.addMeshPartitionAll(meshPart.meshAllPartitions(), &comm);
    HDF5Output.postProcess(0);
    HDF5Output.CloseFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return(EXIT_SUCCESS);
}

