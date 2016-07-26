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

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/fsi_blocks/solver/FSIHandler.hpp>

#include "boundaryConditions.hpp"

using namespace LifeV;

int
main ( int argc, char** argv )
{
    bool verbose (false);
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        verbose = true;
    }
#else
    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm () );
    verbose = true;
#endif
    
    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    
    {
    // ---------------------------//
    // Reading the input datafile //
    // ---------------------------//
    
    const std::string defaultDataName = "data";
    GetPot command_line (argc, argv);
    std::string data_file_name = command_line.follow (defaultDataName.c_str(), 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    
    // --------------------------------------------------//
    // Initializing the handler and setting the datafile //
    // --------------------------------------------------//
    
    FSIHandler fsi ( Comm );
    
    fsi.setDatafile ( dataFile );
    
    // check about the use of meshes: using partitioned meshes or do we have to partition them online?

    bool usePartitionedMeshes = dataFile("offlinePartioner/readPartitionedMeshes", false);

    if ( usePartitionedMeshes )
    {
    	if ( verbose )
    		std::cout << "\n[PreProcessing] - Using partitioned meshes\n";

    	int numParts = dataFile("offlinePartioner/parts", 2);

    	if( numParts != Comm->NumProc())
    	{
    		if ( verbose )
    			std::cout << "\n[PreProcessing] - Using the wrong number of processes\n";

    		return EXIT_FAILURE;
    	}

    	// -------------------------------------------------------//
    	// Loading the partitioned fluid and the structure meshes //
    	// -------------------------------------------------------//

    	fsi.readPartitionedMeshes ( );
    }
    else
    {
    	if ( verbose )
    		std::cout << "\n[PreProcessing] - The meshes will be partioned during this simulation\n\n";

    	// -------------------------------------------//
    	// Loading the fluid and the structure meshes //
    	// -------------------------------------------//

    	fsi.readMeshes ( );

    	// ------------------------------------------------//
    	// Partitioning the fluid and the structure meshes //
    	// ------------------------------------------------//

    	fsi.partitionMeshes ( );
    }


    // --------------------------------//
    // Getting the boundary conditions //
    // --------------------------------//

    boost::shared_ptr<BCHandler> fluidBC ( new BCHandler (*BCh_fluid () ) );
    boost::shared_ptr<BCHandler> fluidBC_residual ( new BCHandler (*BCh_fluid_residual () ) );
    boost::shared_ptr<BCHandler> structureBC ( new BCHandler (*BCh_structure () ) );
    boost::shared_ptr<BCHandler> structureBC_residual ( new BCHandler (*BCh_structure_residual () ) );
    boost::shared_ptr<BCHandler> aleBC ( new BCHandler (*BCh_ale () ) );
    boost::shared_ptr<BCHandler> aleBC_residual ( new BCHandler (*BCh_ale_residual () ) );
    boost::shared_ptr<BCHandler> pcdBC ( new BCHandler (*BCh_PCD () ) );
    boost::shared_ptr<BCHandler> interfaceFluidBC ( new BCHandler (*BCh_interfaceFluid ( ) ) );

    fsi.setBoundaryConditions(fluidBC, fluidBC_residual, structureBC, structureBC_residual, aleBC, aleBC_residual);

    fsi.setFluidInterfaceBoundaryConditions(interfaceFluidBC);

    std::string preconditioner = dataFile("fluid/preconditionerType","none");

    if ( preconditioner.compare("PCD") == 0 )
    	fsi.setBoundaryConditionsPCD(pcdBC);

    // ---------------------------------------------------------------------------------------------//
    // Reading the physical informations for the fluid and the structure and initialize the solvers //
    // ---------------------------------------------------------------------------------------------//
    
    fsi.setup();

    // ---------------------//
    // Create inteface maps //
    // ---------------------//

    fsi.buildInterfaceMaps ( );

    // -------------------------------//
    // Create blocks for the coupling //
    // -------------------------------//

    if ( !dataFile( "interface/nonconforming", false ) )
    	fsi.assembleCoupling ( );

    // ----------//
    // Time Loop //
    // ----------//
        
    fsi.solveFSIproblem ( );

    }

#ifdef HAVE_MPI
    if (verbose)
    {
        std::cout << "\nMPI Finalization\n" << std::endl;
    }
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}
