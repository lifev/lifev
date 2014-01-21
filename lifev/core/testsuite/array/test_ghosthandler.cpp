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
   \author A. Cervone <ant.cervone@gmail.com>
   \date 2011-11-03
 */

/*!
todo
*/

// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/array/GhostHandler.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/NeighborMarker.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/filter/Exporter.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

int main ( int argc, char* argv[] )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::cout << "MPI Initialization" << std::endl;
#endif

    // this brace is important to destroy the Epetra_Comm object before calling MPI_Finalize
    {

        typedef RegionMesh<LinearTriangle, neighborMarkerCommon_Type> mesh_Type;
        typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
        typedef boost::shared_ptr< feSpace_Type > feSpacePtr_Type;

        LifeChrono chronoTotal;
        LifeChrono chronoMesh;
        LifeChrono chronoGhost;

        // Start chronoTotal for measure the total time for the computation
        chronoTotal.start();

        GetPot command_line (argc, argv);
        const std::string data_file_name = command_line.follow ("data_ghosthandler", 2, "-f", "--file");
        GetPot dataFile ( data_file_name );

#ifdef EPETRA_MPI
        std::cout << "Epetra Initialization" << std::endl;
        boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
        boost::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm() );
#endif

        // Create the leader process, i.e. the process with MyPID equal to zero
        bool isLeader = ( comm->MyPID() == 0 );

#ifdef HAVE_LIFEV_DEBUG
        // Create a stream that is different for each process
        std::ofstream fileOut ( ( "gh." + boost::lexical_cast<std::string> ( comm->MyPID() ) + ".out" ).c_str() );
#else
        // discard file output in opt mode
        std::ofstream fileOut ( "/dev/null" );
#endif

        if ( isLeader )
        {
            std::cout << "GhostHendler test" << std::endl;
        }

        // Start chronoMesh for measure the total time for the creation of the local meshes
        chronoMesh.start();

        // Create the mesh file handler
        MeshData meshData;

        // Set up the mesh file
        meshData.setup ( dataFile, "space_discretization" );

        // Create the mesh
        boost::shared_ptr<mesh_Type> fullMeshPtr ( new mesh_Type ( comm ) );

        // Set up the mesh
        readMesh ( *fullMeshPtr, meshData );

        // Partition the mesh using ParMetis
        boost::shared_ptr<mesh_Type> meshPtr ( new mesh_Type ( comm ) );
        {
            MeshPartitioner< mesh_Type >  meshPart ( fullMeshPtr, comm );
            meshPtr = meshPart.meshPartition();
        }

        // Stop chronoReadAndPartitionMesh
        chronoMesh.stop();

        // The leader process print chronoReadAndPartitionMesh
        if ( isLeader )
        {
            std::cout << "  C- Time for mesh " << chronoMesh.diff() << std::endl;
        }

        // Start chronoGhost for measure the total time for GhostHandler routines
        chronoGhost.start();

        feSpacePtr_Type feSpaceP1 ( new feSpace_Type ( meshPtr,
                                                       feTriaP1,
                                                       quadRuleTria4pt,
                                                       quadRuleSeg2pt,
                                                       1,
                                                       comm ) );

        fileOut << "=================== P1" << std::endl;
        fileOut << "9---8---7---6---5" << std::endl;
        fileOut << "|  /|  /|  /|  /|" << std::endl;
        fileOut << "| / | / | / | / |" << std::endl;
        fileOut << "|/  |/  |/  |/  |" << std::endl;
        fileOut << "0---1---2---3---4" << std::endl;

        GhostHandler<mesh_Type> ghostP1 ( fullMeshPtr, meshPtr, feSpaceP1->mapPtr(), comm );

        ghostP1.setUpNeighbors();

        MapEpetra mapP1 ( feSpaceP1->map() );
        MapEpetra mapP1Overlap ( ghostP1.ghostMapOnPoints ( dataFile ( "ghost/overlap", 2 ) ) );

        fileOut << "=================== mapP1 Unique" << std::endl;
        fileOut << "9---8---7---+---+    +---+---+---6---5" << std::endl;
        fileOut << "|  /|  /|  /|  /|    |  /|  /|  /|  /|" << std::endl;
        fileOut << "| / | / | / | / |    | / | / | / | / |" << std::endl;
        fileOut << "|/  |/  |/  |/  |    |/  |/  |/  |/  |" << std::endl;
        fileOut << "0---1---2---+---+    +---+---+---3---4" << std::endl;
        fileOut << *mapP1.map ( Unique );

        fileOut << "=================== mapP1 Repeated" << std::endl;
        fileOut << "9---8---7---+---+    +---+---7---6---5" << std::endl;
        fileOut << "|  /|  /|  /|  /|    |  /|  /|  /|  /|" << std::endl;
        fileOut << "| / | / | / | / |    | / | / | / | / |" << std::endl;
        fileOut << "|/  |/  |/  |/  |    |/  |/  |/  |/  |" << std::endl;
        fileOut << "0---1---2---+---+    +---+---2---3---4" << std::endl;
        fileOut << *mapP1.map ( Repeated );

        fileOut << "=================== mapP1 Repeated overlap 2" << std::endl;
        fileOut << "9---8---7---6---5    9---8---7---6---5" << std::endl;
        fileOut << "|  /|  /|  /|  /|    |  /|  /|  /|  /|" << std::endl;
        fileOut << "| / | / | / | / |    | / | / | / | / |" << std::endl;
        fileOut << "|/  |/  |/  |/  |    |/  |/  |/  |/  |" << std::endl;
        fileOut << "0---1---2---3---4    0---1---2---3---4" << std::endl;
        fileOut << *mapP1Overlap.map ( Repeated );

#ifdef HAVE_HDF5

        ghostP1.exportToHDF5();

        ghostP1.importFromHDF5();

#endif // HAVE_HDF5

        ghostP1.clean();

        boost::shared_ptr<VectorEpetra> vP1 ( new VectorEpetra ( mapP1Overlap, Unique ) );

        // get all elements from the repeated map
        Int* pointer ( mapP1Overlap.map ( Repeated )->MyGlobalElements() );
        for ( Int ii = 0; ii < mapP1Overlap.map ( Repeated )->NumMyElements(); ++ii, ++pointer )
        {
            vP1->sumIntoGlobalValues ( *pointer, 1 );
        }

        vP1->globalAssemble();

        // check that the overlapping map has all the points in the mesh
        if ( mapP1Overlap.map ( Repeated )->NumMyElements() != 10 )
        {
            return 1;
        }

        feSpacePtr_Type feSpaceP0 ( new feSpace_Type ( meshPtr,
                                                       feTriaP0,
                                                       quadRuleTria1pt,
                                                       quadRuleSeg1pt,
                                                       1,
                                                       comm ) );

        fileOut << "=================== P0" << std::endl;
        fileOut << "+---+---+---+---+" << std::endl;
        fileOut << "|1 /|3 /|5 /|7 /|" << std::endl;
        fileOut << "| / | / | / | / |" << std::endl;
        fileOut << "|/ 0|/ 2|/ 4|/ 6|" << std::endl;
        fileOut << "+---+---+---+---+" << std::endl;

        GhostHandler<mesh_Type> ghostP0 ( fullMeshPtr, meshPtr, feSpaceP0->mapPtr(), comm );

        ghostP0.setUpNeighbors();

        MapEpetra mapP0 ( feSpaceP0->map() );
        MapEpetra mapP0P0 ( ghostP0.ghostMapOnElementsFV() );
        MapEpetra mapP0P1 ( ghostP0.ghostMapOnElementsFE ( dataFile ( "ghost/overlap", 2 ) ) );

        fileOut << "=================== mapP0 Unique" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << "|1 /|3 /|  /|  /|   |  /|  /|5 /|7 /|" << std::endl;
        fileOut << "| / | / | / | / |   | / | / | / | / |" << std::endl;
        fileOut << "|/ 0|/ 2|/  |/  |   |/  |/  |/ 4|/ 6|" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << *mapP0.map ( Unique );

        fileOut << "=================== mapP0 Repeated" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << "|1 /|3 /|  /|  /|   |  /|  /|5 /|7 /|" << std::endl;
        fileOut << "| / | / | / | / |   | / | / | / | / |" << std::endl;
        fileOut << "|/ 0|/ 2|/  |/  |   |/  |/  |/ 4|/ 6|" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << *mapP0.map ( Repeated );

        fileOut << "=================== mapP0 Repeated face neighbors" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << "|1 /|3 /|5 /|  /|   |  /|  /|5 /|7 /|" << std::endl;
        fileOut << "| / | / | / | / |   | / | / | / | / |" << std::endl;
        fileOut << "|/ 0|/ 2|/  |/  |   |/  |/ 2|/ 4|/ 6|" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << *mapP0P0.map ( Repeated );

        fileOut << "=================== mapP0 Repeated point neighbors overlap 2" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << "|1 /|3 /|5 /|7 /|   |1 /|3 /|5 /|7 /|" << std::endl;
        fileOut << "| / | / | / | / |   | / | / | / | / |" << std::endl;
        fileOut << "|/ 0|/ 2|/ 4|/ 6|   |/ 0|/ 2|/ 4|/ 6|" << std::endl;
        fileOut << "+---+---+---+---+   +---+---+---+---+" << std::endl;
        fileOut << *mapP0P1.map ( Repeated );

        ghostP0.showMe ( true, fileOut );

        ghostP0.clean();

        boost::shared_ptr<VectorEpetra> vP0 ( new VectorEpetra ( mapP0P0, Unique ) );

        // get all elements from the repeated map
        pointer = mapP0P1.map ( Repeated )->MyGlobalElements();
        if ( isLeader )
            for ( Int ii = 0; ii < mapP0P1.map ( Repeated )->NumMyElements(); ++ii, ++pointer )
            {
                vP0->sumIntoGlobalValues ( *pointer, 1 );
            }

        vP0->globalAssemble();

        fileOut << "=================== vector" << std::endl;
        fileOut << vP0->epetraVector();

        // check that the overlapping map has all the elements in the mesh
        if ( mapP0P1.map ( Repeated )->NumMyElements() != 8 )
        {
            return 1;
        }

        // Stop chronoGhost
        chronoGhost.stop();

        // The leader process print chronoGhost
        if ( isLeader )
        {
            std::cout << "  C- Time for ghost " << chronoGhost.diff() << std::endl;
        }

        boost::shared_ptr< Exporter< mesh_Type > > exporter;

        // Type of the exporter
        std::string const exporterType =  dataFile ( "exporter/type", "ensight" );

        // Choose the exporter
#ifdef HAVE_HDF5
        if ( exporterType.compare ( "hdf5" ) == 0 )
        {
            exporter.reset ( new ExporterHDF5< mesh_Type > ( dataFile, dataFile ( "exporter/file_name", "GH" ) ) );
        }
        else
#endif
        {
            if ( exporterType.compare ("none") == 0 )
            {
                exporter.reset ( new ExporterEmpty< mesh_Type > ( dataFile, dataFile ( "exporter/file_name", "GH" ) ) );
            }
            else
            {
                exporter.reset ( new ExporterEnsight< mesh_Type > ( dataFile, dataFile ( "exporter/file_name", "GH" ) ) );
            }
        }

        // Set directory where to save the solution
        exporter->setPostDir ( dataFile ( "exporter/folder", "./" ) );
        exporter->setMeshProcId ( meshPtr, comm->MyPID() );

        // Export the partitioning
        exporter->exportPID ( meshPtr, comm );

        // Add the solution to the exporter
        exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,
                                "MapP1", feSpaceP1,
                                vP1,
                                static_cast<UInt> ( 0 ),
                                ExporterData<mesh_Type>::SteadyRegime,
                                ExporterData<mesh_Type>::Node );

        // Add the solution to the exporter
        exporter->addVariable ( ExporterData<mesh_Type>::ScalarField,
                                "MapP0", feSpaceP0,
                                vP0,
                                static_cast<UInt> ( 0 ),
                                ExporterData<mesh_Type>::SteadyRegime,
                                ExporterData<mesh_Type>::Cell );

        // Save the initial solution into the exporter
        exporter->postProcess ( 0 );

        // Stop chronoTotal
        chronoTotal.stop();

        // The leader process print chronoTotal
        if ( isLeader )
        {
            std::cout << "  C- Time total " << chronoTotal.diff() << std::endl;
        }

    }

#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout << "MPI Finalization" << std::endl;
#endif

    return 0;

}
