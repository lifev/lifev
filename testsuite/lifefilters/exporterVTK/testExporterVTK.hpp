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

/*
 * testExporterVTK.hpp
 *
 *  Created on: Jan 13, 2011
 *      Author: tiziano
 */

#ifndef TESTEXPORTERVTK_HPP_
#define TESTEXPORTERVTK_HPP_ 1

// LifeV definition files
#include <life/lifecore/Displayer.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifecore/LifeV.hpp>
#include <life/lifefem/TimeData.hpp>
#include <life/lifefilters/ExporterVTK.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifefunctions/Womersley.hpp>
#include <life/lifefunctions/RossEthierSteinmanInc.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// Trilinos-MPI communication definitions
#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// Object type definitions
typedef LifeV::RegionMesh<LifeV::LinearTetra>       mesh_Type;
typedef LifeV::RossEthierSteinmanUnsteadyInc          problem_Type;
typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;
typedef boost::shared_ptr<Epetra_Comm>                commPtr_Type;


class TestExporterVTK
{
public:
    TestExporterVTK( const commPtr_Type& commPtr );

    bool testExport( GetPot& commandLine );
    bool testImport( GetPot& commandLine );

private:
    void loadData( GetPot& commandLine );
    void buildMesh();
    void buildFESpaces();
    void buildExporter( boost::shared_ptr< LifeV::ExporterVTK< mesh_Type > >& exporterPtr );

    commPtr_Type                                             M_commPtr;
    LifeV::Displayer                                         M_displayer;
    GetPot                                                   M_dataFile;
    boost::shared_ptr< LifeV::MeshPartitioner< mesh_Type > > M_meshPartPtr;

    feSpacePtr_Type                                          M_velFESpacePtr;
    feSpacePtr_Type                                          M_pressFESpacePtr;
};


TestExporterVTK::TestExporterVTK( const commPtr_Type& commPtr ) :
    M_commPtr  ( commPtr ),
    M_displayer( commPtr )
{}


void
TestExporterVTK::loadData( GetPot& commandLine )
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Loading the data]\n" );
    chrono.start();

    const std::string problemName = commandLine.follow("", 2, "-p", "--problem");
    const std::string defaultDataFileName("data"+problemName);
    const std::string dataFileName = commandLine.follow(defaultDataFileName.c_str(), 2, "-f","--file");
    M_dataFile = GetPot(dataFileName);

    problem_Type::setParamsFromGetPot(M_dataFile);

    chrono.stop();
    M_displayer.leaderPrint( "[...done in ", chrono.diff(), "s]\n" );
}


void
TestExporterVTK::buildMesh()
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |              Building the mesh                |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Building the mesh]\n" );
    chrono.start();

    boost::shared_ptr< mesh_Type > fullMeshPtr(new mesh_Type);
    UInt nEl(M_dataFile("space_discretization/dimension", 1));
    regularMesh3D(*fullMeshPtr, 0, nEl, nEl, nEl);
    //  The following is when reading from file
    /*
    MeshData meshData;
    meshData.setup(M_dataFile, "fluid/space_discretization");
    if (verbose) std::cout << "Mesh file: " << dataMesh.meshDir() << dataMesh.meshFile() << std::endl;
    readMesh(*fullMeshPtr, dataMesh);
    */
    // Split the mesh between processors
    M_meshPartPtr.reset( new MeshPartitioner< mesh_Type >( fullMeshPtr, M_commPtr ) );
    // Release the original mesh from the MeshPartitioner object and delete the RegionMesh object
    M_meshPartPtr->releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    chrono.stop();
    M_displayer.leaderPrint( "[...done in ", chrono.diff(), "s]\n" );
}


void
TestExporterVTK::buildFESpaces()
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Creating the FE spaces...]\n" );
    chrono.start();

    const std::string velFE =  M_dataFile( "space_discretization/velFE", "P2");
    M_displayer.leaderPrint( "\t-o FE for the velocity: ", velFE, "\n" );

    M_displayer.leaderPrint( "\t-o Building the velocity FE space...\n" );

    M_velFESpacePtr.reset( new feSpace_Type(*M_meshPartPtr, velFE, nDimensions, M_commPtr) );

    M_displayer.leaderPrint( "\t\t...ok.\n" );

    const std::string pressFE =  M_dataFile( "space_discretization/pressFE", "P1");
    M_displayer.leaderPrint( "\t-o FE for the pressure: ", pressFE, "\n" );

    M_displayer.leaderPrint( "\t-o Building the pressure FE space...\n" );

    M_pressFESpacePtr.reset( new feSpace_Type(*M_meshPartPtr, pressFE, 1, M_commPtr) );

    M_displayer.leaderPrint( "\t\t...ok.\n" );

    // Total degrees of freedom (elements of matrix)
    UInt velTotalDof   = M_velFESpacePtr->map().map(Unique)->NumGlobalElements();
    UInt pressTotalDof = M_pressFESpacePtr->map().map(Unique)->NumGlobalElements();

    M_displayer.leaderPrint( "\t-o Total Dof for the velocity: ", velTotalDof, "\n" );
    M_displayer.leaderPrint( "\t-o Total Dof for the pressure: ", pressTotalDof, "\n" );

    chrono.stop();
    M_displayer.leaderPrint( "[...done in ", chrono.diff(), "s]\n" );

}


void
TestExporterVTK::buildExporter( boost::shared_ptr< LifeV::ExporterVTK< mesh_Type > >& exporterPtr )
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |            Creating the exporter              |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Creating the exporter...]\n" );
    chrono.start();

    exporterPtr.reset( new ExporterVTK<mesh_Type > ( M_dataFile, "testExporterVTK" ) );
    exporterPtr->setPostDir( "./" );
    exporterPtr->setMeshProcId( M_meshPartPtr->meshPartition(), M_commPtr->MyPID() );

    chrono.stop();
    M_displayer.leaderPrint( "[...done in ", chrono.diff(), "s]\n" );

}


bool
TestExporterVTK::testExport( GetPot& commandLine )
{
    using namespace LifeV;

    loadData( commandLine );
    buildMesh();
    buildFESpaces();

    boost::shared_ptr< LifeV::ExporterVTK< mesh_Type > > exporterPtr;

    buildExporter( exporterPtr );

    // Chronometer
    LifeChrono globalChrono;

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Entering the time loop]\n" );
    globalChrono.start();

    Exporter<mesh_Type >::vectorPtr_Type velInterpolantPtr(
        new Exporter<mesh_Type >::vector_Type   ( M_velFESpacePtr->map(), Repeated ) );
    Exporter<mesh_Type >::vectorPtr_Type pressInterpolantPtr(
        new Exporter<mesh_Type >::vector_Type   ( M_pressFESpacePtr->map(), Repeated ) );

    exporterPtr->addVariable( ExporterData<mesh_Type>::VectorField, "velocity",
                           M_velFESpacePtr, velInterpolantPtr, UInt(0) );

    exporterPtr->addVariable( ExporterData<mesh_Type>::ScalarField, "pressure",
                           M_pressFESpacePtr, pressInterpolantPtr, UInt(0) );
    exporterPtr->postProcess( 0 );

    TimeData timeData;
    timeData.setup( M_dataFile, "time_discretization" );

    for ( ; timeData.canAdvance(); timeData.updateTime() )
    {

        M_displayer.leaderPrint( "[t = ", timeData.time(), " s.]\n" );

        // Computation of the interpolation
        M_velFESpacePtr->interpolate( problem_Type::uexact, *velInterpolantPtr, timeData.time() );
        M_pressFESpacePtr->interpolate( problem_Type::pexact, *pressInterpolantPtr, timeData.time() );

        // Exporting the solution
        exporterPtr->postProcess( timeData.time() );

        MPI_Barrier(MPI_COMM_WORLD);

    }

    globalChrono.stop();
    M_displayer.leaderPrint( "[Time loop, elapsed time:  ", globalChrono.diff(), " s.]\n" );

    exporterPtr->closeFile();

    return 0;
}


bool
TestExporterVTK::testImport( GetPot& commandLine )
{
    using namespace LifeV;

    loadData( commandLine );
    buildMesh();
    buildFESpaces();

    boost::shared_ptr< LifeV::ExporterVTK< mesh_Type > > importerPtr, exporterPtr;
    buildExporter( importerPtr );
    buildExporter( exporterPtr );

    // Chronometer
    LifeChrono globalChrono;

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Entering the time loop]\n" );
    globalChrono.start();

    Exporter<mesh_Type >::vectorPtr_Type velImportedPtr(
        new Exporter<mesh_Type >::vector_Type   ( M_velFESpacePtr->map(), Repeated ) );
    Exporter<mesh_Type >::vectorPtr_Type pressImportedPtr(
        new Exporter<mesh_Type >::vector_Type   ( M_pressFESpacePtr->map(), Repeated ) );

    importerPtr->addVariable( ExporterData<mesh_Type>::VectorField, "velocity",
                           M_velFESpacePtr, velImportedPtr, UInt(0) );

    importerPtr->addVariable( ExporterData<mesh_Type>::ScalarField, "pressure",
                           M_pressFESpacePtr, pressImportedPtr, UInt(0) );

    exporterPtr->addVariable( ExporterData<mesh_Type>::VectorField, "importedVelocity",
                           M_velFESpacePtr, velImportedPtr, UInt(0) );

    exporterPtr->addVariable( ExporterData<mesh_Type>::ScalarField, "importedPressure",
                           M_pressFESpacePtr, pressImportedPtr, UInt(0) );

    TimeData timeData;
    timeData.setup( M_dataFile, "time_discretization" );

    M_displayer.leaderPrint( "\t[t = ", timeData.time(), " s.]\n" );

    // Import test
    importerPtr->setTimeIndex( 1 );
    importerPtr->import( timeData.time() );

    exporterPtr->postProcess( 0 );

    globalChrono.stop();
    M_displayer.leaderPrint( "[Time loop, elapsed time:  ", globalChrono.diff(), " s.\n" );

    return 0;
}

#endif /* TESTEXPORTERVTK_HPP_ */
