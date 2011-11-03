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

#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifearray/GhostHandler.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifefilters/Exporter.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/ExporterHDF5.hpp>
#endif
#include <life/lifefilters/ExporterEmpty.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

int main( int argc, char* argv[] )
{
    typedef RegionMesh3D<LinearTetra,neighborMarkerCommon_Type> RegionMesh;
    typedef FESpace< RegionMesh, MapEpetra >            feSpace_Type;
    typedef boost::shared_ptr< feSpace_Type >           feSpacePtr_Type;

    LifeChrono chronoTotal;
    LifeChrono chronoMesh;
    LifeChrono chronoGhost;

    // Start chronoTotal for measure the total time for the computation
    chronoTotal.start();

    GetPot command_line(argc, argv);
    const std::string data_file_name = command_line.follow("data_ghosthandler", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    std::cout << "MPI Initialization" << std::endl;
#endif

#ifdef EPETRA_MPI
    std::cout << "Epetra Initialization" << std::endl;
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm() );
#endif

    // Create the leader process, i.e. the process with MyPID equal to zero
    bool isLeader = ( comm->MyPID() == 0 );

    // Create a stream that is different for each process
    std::ofstream fileOut ( ( "gh." + boost::lexical_cast<std::string>( comm->MyPID() ) + ".out" ).c_str() );

    if ( isLeader ) std::cout << "GhostHendler test" << std::endl;

    // Start chronoMesh for measure the total time for the creation of the local meshes
    chronoMesh.start();

    // Create the mesh file handler
    MeshData meshData;

    // Set up the mesh file
    meshData.setup( dataFile, "space_discretization" );

    // Create the mesh
    boost::shared_ptr<RegionMesh> fullMeshPtr( new RegionMesh );

    // Set up the mesh
    readMesh( *fullMeshPtr, meshData );

    // create node neighbors
//    createNodeNeighbors ( *fullMeshPtr );

    // Partition the mesh using ParMetis
    MeshPartitioner< RegionMesh >  meshPart( fullMeshPtr, comm );

    // Stop chronoReadAndPartitionMesh
    chronoMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    if ( isLeader ) std::cout << "C - Time for mesh " << chronoMesh.diff() << std::endl;

    // Start chronoGhost for measure the total time for GhostHandler routines
    chronoGhost.start();

    feSpacePtr_Type feSpaceP1( new feSpace_Type( meshPart,
                                                 feTetraP1,
                                                 quadRuleTetra15pt,
                                                 quadRuleTria1pt,
                                                 1,
                                                 comm ) );

    GhostHandler<RegionMesh> ghostP1 ( fullMeshPtr, meshPart.meshPartition(), feSpaceP1->map(), comm );

//    ghostP1.setUp();

    MapEpetra mapP1 ( feSpaceP1->map() );
    MapEpetra mapP1Overlap1 ( ghostP1.ghostMapOnNodes() );
    MapEpetra mapP1Overlap1bis ( ghostP1.ghostMapOnNodes( dataFile( "ghost/overlap", 1 ) ) );

    fileOut << "=================== mapP1" << std::endl;
    fileOut << *mapP1.map( Unique );
    fileOut << "=================== mapP1" << std::endl;
    fileOut << *mapP1.map( Repeated );
    fileOut << "=================== mapP1Overlap1" << std::endl;
    fileOut << *mapP1Overlap1.map( Repeated );
    fileOut << "=================== mapP1Overlap1bis" << std::endl;
    fileOut << *mapP1Overlap1bis.map( Repeated );

    ghostP1.clean();

    feSpacePtr_Type feSpaceP0( new feSpace_Type( meshPart,
                                                 feTetraP0,
                                                 quadRuleTetra15pt,
                                                 quadRuleTria1pt,
                                                 1,
                                                 comm ) );

    GhostHandler<RegionMesh> ghostP0 ( fullMeshPtr, meshPart.meshPartition(), feSpaceP0->map(), comm );

    ghostP0.setUp();

    MapEpetra mapP0 ( feSpaceP0->map() );
    MapEpetra mapP0P0 ( ghostP0.ghostMapOnElementsP0() );
    MapEpetra mapP0P1 ( ghostP0.ghostMapOnElementsP1() );

    fileOut << "=================== mapP0" << std::endl;
    fileOut << *mapP0.map( Unique );
    fileOut << "=================== mapP0" << std::endl;
    fileOut << *mapP0.map( Repeated );
    fileOut << "=================== mapP0P0" << std::endl;
    fileOut << *mapP0P0.map( Repeated );
    fileOut << "=================== mapP0P1" << std::endl;
    fileOut << *mapP0P1.map( Repeated );

    ghostP0.showMe( true, fileOut );

    ghostP0.clean();

    // Stop chronoGhost
    chronoGhost.stop();

    // The leader process print chronoGhost
    if ( isLeader ) std::cout << "C - Time for ghost " << chronoGhost.diff() << std::endl;

    boost::shared_ptr< Exporter< RegionMesh > > exporter;

    // Type of the exporter
    std::string const exporterType =  dataFile( "exporter/type", "ensight" );

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare( "hdf5" ) == 0 )
        exporter.reset( new ExporterHDF5< RegionMesh > ( dataFile, dataFile( "exporter/file_name", "Concentration" ) ) );
    else
#endif
    {
        if ( exporterType.compare("none") == 0 )
        {
            exporter.reset( new ExporterEmpty< RegionMesh > ( dataFile, dataFile( "exporter/file_name", "Concentration" ) ) );
        }
        else
        {
            exporter.reset( new ExporterEnsight< RegionMesh > ( dataFile, dataFile( "exporter/file_name", "Concentration" ) ) );
        }
    }

    // Set directory where to save the solution
    exporter->setPostDir( dataFile( "exporter/folder", "./" ) );
    exporter->setMeshProcId( meshPart.meshPartition(), comm->MyPID() );

    // Export the partitioning
    exporter->exportPID( meshPart );

    // Set the exporter solution
//    exporterSolution.reset( new vector_type ( *hyperbolicSolver.solution(),exporter->mapType() ) );

    // Add the solution to the exporter
//    exporter->addVariable( ExporterData<RegionMesh>::ScalarField,
//                           "Concentration", feSpacePtr,
//                           exporterSolution,
//                           static_cast<UInt>( 0 ),
//                           ExporterData<RegionMesh>::UnsteadyRegime,
//                           ExporterData<RegionMesh>::Cell );

    // Copy the initial solution to the exporter
//    *exporterSolution = *hyperbolicSolver.solution();

    // Save the initial solution into the exporter
//    exporter->postProcess( dataHyperbolic.dataTime()->initialTime() );

    // Stop chronoTotal
    chronoTotal.stop();

    // The leader process print chronoTotal
    if ( isLeader ) std::cout << "Time total " << chronoTotal.diff() << std::endl;


#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout << "MPI Finalization" << std::endl;
#endif

    return 0;

}
