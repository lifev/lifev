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
    @brief Test of the serial mesh partitioning - part 1

    @author Radu Popescu <radu.popescu@epfl.ch>
    @date 02-07-2010

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

#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefilters/HDF5Filter3DMesh.hpp>

#include <iostream>
#include <string>
#include <mpi.h>


using namespace LifeV;

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

    DataNavierStokes dataNavierStokes;
    dataNavierStokes.setup(dataFile);

    DataMesh dataMesh;
    dataMesh.setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);
    readMesh(*fullMeshPtr, dataMesh);

    partitionMesh<RegionMesh3D<LinearTetra> > meshPart;
    meshPart.setup(4, comm);

    meshPart.attachUnpartitionedMesh(fullMeshPtr);
    meshPart.doPartitionGraph();
    meshPart.doPartitionMesh();

    // Release the original mesh from the partitionMesh object and delete the RegionMesh3D object
    meshPart.releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    HDF5Filter3DMesh<RegionMesh3D<LinearTetra> > HDF5Output(dataFile, meshPart.meshPartition(), "cylinderPart",
                                                            comm->MyPID());
    HDF5Output.addPartitionGraph(meshPart.elementDomains(), comm);
    HDF5Output.addMeshPartitionAll(meshPart.meshPartitions(), comm);
    HDF5Output.postProcess(0);
    HDF5Output.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return(EXIT_SUCCESS);
}

