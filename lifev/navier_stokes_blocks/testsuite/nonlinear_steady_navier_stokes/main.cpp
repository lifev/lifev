/* -*- mode: c++ -*-

  This file is part of the LifeV library.
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
   \author Davide Forti <davide.forti@epfl.ch>
   \date 2014-02-06
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes_blocks/solver/NavierStokesSolverBlocks.hpp>
#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

#include "boundaryConditions.hpp"

using namespace LifeV;

int
main ( int argc, char** argv )
{
    bool verbose (false);
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        verbose = true;
    }
#else
    std::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm () );
    verbose = true;
#endif

    Real normTwo_Velo;
    Real normTwo_Pres;

    {

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef VectorEpetra vector_Type;
    typedef std::shared_ptr<vector_Type> vectorPtr_Type;

    // Reading the dataFile
    const std::string defaultDataName = "data";
    GetPot command_line (argc, argv);
    std::string data_file_name = command_line.follow (defaultDataName.c_str(), 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    // Mesh
    std::shared_ptr<mesh_Type > localMeshPtr ( new mesh_Type ( Comm ) );

    if ( dataFile ( "offline_partitioner/useOfflinePartitionedMesh", false) )
    {
    	std::shared_ptr<Epetra_MpiComm> comm = std::dynamic_pointer_cast<Epetra_MpiComm>(Comm);
    	const std::string partsFileName (dataFile ("offline_partitioner/hdf5_file_name", "name.h5") );
    	PartitionIO<mesh_Type > partitionIO (partsFileName, comm);
    	partitionIO.read (localMeshPtr);
    }
    else
    {
    	// reading the mesh
    	std::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );
    	MeshData meshData;
    	meshData.setup (dataFile, "fluid/space_discretization");
    	readMesh (*fullMeshPtr, meshData);

    	// mesh partitioner
    	MeshPartitioner< mesh_Type >  meshPart (fullMeshPtr, Comm);
    	localMeshPtr = meshPart.meshPartition();
    	fullMeshPtr.reset();
    }

    // create the solver
    NavierStokesSolverBlocks ns( dataFile, Comm);
    ns.setup(localMeshPtr);
    ns.setParameters();
    ns.buildSystem();

    // Exporter
    std::string outputName = dataFile ( "exporter/filename", "result");
    std::shared_ptr< Exporter<mesh_Type > > exporter;
    std::string const exporterType =  dataFile ( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
    	exporter.reset ( new ExporterHDF5<mesh_Type > ( dataFile, outputName ) );
    	exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
    	exporter->setMeshProcId ( localMeshPtr, Comm->MyPID() );
    }
#endif
    else if(exporterType.compare ("vtk") == 0)
    {
    	exporter.reset ( new ExporterVTK<mesh_Type > ( dataFile, outputName ) );
    	exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
    	exporter->setMeshProcId ( localMeshPtr, Comm->MyPID() );
    }

    // Vectors post-processing
    vectorPtr_Type velocity( new vector_Type(ns.uFESpace()->map(), exporter->mapType() ) );
    vectorPtr_Type pressure( new vector_Type(ns.pFESpace()->map(), exporter->mapType() ) );
    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", ns.uFESpace(), velocity, UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", ns.pFESpace(), pressure, UInt (0) );
    exporter->postProcess ( 0.0 );

    // Boundary conditions
    std::shared_ptr<BCHandler> bc ( new BCHandler (*BCh_fluid ()) );
    
    // Set boundary conditions
    ns.setBoundaryConditions( bc );
    
    // Solve problem
    ns.iterate_steady( );

    // Get the velocity
    ns.updateVelocity(velocity);
    
    // Get the pressure
    ns.updatePressure(pressure);

    // Do post-processing
    exporter->postProcess ( 1.0 );

    // Close file post-processing
    exporter->closeFile();

    normTwo_Velo  = velocity->norm2();
    normTwo_Pres  = pressure->norm2();

	}

#ifdef HAVE_MPI
    if (verbose)
    {
        std::cout << "\nMPI Finalization" << std::endl;
    }
    MPI_Finalize();
#endif

    if ( std::abs(normTwo_Velo - 108.783606 ) < 1.0e-4 && std::abs(normTwo_Pres - 129.597065 ) < 1.0e-4 )
    {
    	return ( EXIT_SUCCESS );
    }
    else
    {
    	return ( EXIT_FAILURE );
    }
}


