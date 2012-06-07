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
 * testImportExport.hpp
 *
 *  Created on: Jan 13, 2011
 *      Author: tiziano
 */

#ifndef TESTIMPORTEXPORT_HPP_
#define TESTIMPORTEXPORT_HPP_ 1

// LifeV definition files
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/LifeChrono.hpp>
//#include <life/lifecore/LifeV.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes/function/Womersley.hpp>
#include <lifev/navier_stokes/function/RossEthierSteinmanInc.hpp>

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
typedef LifeV::RegionMesh<LifeV::LinearTetra>         mesh_Type;
typedef LifeV::RossEthierSteinmanUnsteadyInc          problem_Type;
typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;
typedef boost::shared_ptr<Epetra_Comm>                commPtr_Type;
typedef LifeV::Exporter<mesh_Type >::vectorPtr_Type   vectorPtr_Type;

class TestImportExport
{
public:
    TestImportExport( const commPtr_Type& commPtr );

    template<typename ImporterType, typename ExporterType>
    bool run( GetPot& commandLine, const std::string& testString = "import" );

    template<typename ImporterType, typename ExporterType>
    bool importLoop( const boost::shared_ptr< ImporterType > & importerPtr,
					 const boost::shared_ptr< ImporterType > & importerMeshPtr,
                     const boost::shared_ptr< ExporterType > & exporterPtr );

    template<typename ImporterType, typename ExporterType>
    bool exportLoop( const boost::shared_ptr< ImporterType > & importerPtr,
                     const boost::shared_ptr< ExporterType > & exporterPtr );

private:
    void loadData( GetPot& commandLine );
    void buildMesh();
    void buildFESpaces();
    template<typename ExporterType>
    void buildExporter( boost::shared_ptr< ExporterType >& exporterPtr,
                        const std::string& prefix );
    //template<typename ImporterType>
    //void buildImporter( boost::shared_ptr< ImporterType >& exporterPtr  );

    commPtr_Type                                             M_commPtr;
    LifeV::Displayer                                         M_displayer;
    GetPot                                                   M_dataFile;
    boost::shared_ptr< LifeV::MeshPartitioner< mesh_Type > > M_meshPartPtr;
    LifeV::TimeData                                          M_timeData;

    feSpacePtr_Type                                          M_velFESpacePtr;
    feSpacePtr_Type                                          M_pressFESpacePtr;
    vectorPtr_Type                                           M_velInterpolantPtr;
    vectorPtr_Type                                           M_pressInterpolantPtr;
    vectorPtr_Type                                           M_velImportedPtr;
    vectorPtr_Type                                           M_pressImportedPtr;
    vectorPtr_Type                                           M_meshDispImportedPtr;
    vectorPtr_Type                                           M_meshVelImportedPtr;
    bool													 M_moving;
    bool													 M_loadFluid;
    std::string 											M_velocityName;
    std::string 											M_pressureName;
    std::string 											M_meshVelocityName;
    std::string 											M_meshDisplacementName;
};


TestImportExport::TestImportExport( const commPtr_Type& commPtr ) :
                    M_commPtr  ( commPtr ),
                    M_displayer( commPtr )
{}


void
TestImportExport::loadData( GetPot& commandLine )
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
    M_moving = M_dataFile("importer/moving", 0);
    M_loadFluid = M_dataFile("importer/loadFluid", 1);

    M_velocityName = M_dataFile("importer/velocityName", "");
    M_pressureName = M_dataFile("importer/pressureName", "");
    M_meshVelocityName = M_dataFile("importer/meshVelocityName", "meshVelocity");
    M_meshDisplacementName = M_dataFile("importer/meshDisplacementName", "meshDisplacement");

    problem_Type::setParamsFromGetPot(M_dataFile);

    chrono.stop();
    M_displayer.leaderPrint( "[...done in ", chrono.diff(), "s]\n" );
}


void
TestImportExport::buildMesh()
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

    if( M_dataFile("space_discretization/mesh_from_file", false) )
    {
        //  The following is when reading from file
        MeshData meshData;
        meshData.setup(M_dataFile, "space_discretization");
        //if (verbose) std::cout << "Mesh file: " << meshData.meshDir() << meshData.meshFile() << std::endl;
        readMesh(*fullMeshPtr, meshData);
    }
    else
    {
        UInt nEl(M_dataFile("space_discretization/dimension", 1));
        regularMesh3D(*fullMeshPtr, 0, nEl, nEl, nEl);
    }
    // Split the mesh between processors
    M_meshPartPtr.reset( new MeshPartitioner< mesh_Type >( fullMeshPtr, M_commPtr ) );
    // Release the original mesh from the MeshPartitioner object and delete the RegionMesh3D object
    M_meshPartPtr->releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    chrono.stop();
    M_displayer.leaderPrint( "[...done in ", chrono.diff(), "s]\n" );
}


void
TestImportExport::buildFESpaces()
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Creating the FE spaces...]\n" );
    chrono.start();

    const std::string velFE =  M_dataFile( "space_discretization/velocity_fespace", "P2");
    M_displayer.leaderPrint( "\t-o FE for the velocity: ", velFE, "\n" );

    M_displayer.leaderPrint( "\t-o Building the velocity FE space...\n" );

    M_velFESpacePtr.reset( new feSpace_Type(*M_meshPartPtr, velFE, nDimensions, M_commPtr) );

    M_displayer.leaderPrint( "\t\t...ok.\n" );

    const std::string pressFE =  M_dataFile( "space_discretization/pressure_fespace", "P1");
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


template<typename ExporterType>
void
TestImportExport::buildExporter( boost::shared_ptr< ExporterType >& exporterPtr,
                                 const std::string& prefix )
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |            Creating the exporter              |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Creating the exporter...]\n" );
    chrono.start();

    exporterPtr.reset( new ExporterType ( M_dataFile, prefix ) );

    //exporterPtr->setPostDir( "./" );
    exporterPtr->setMeshProcId( M_meshPartPtr->meshPartition(), M_commPtr->MyPID() );
    chrono.stop();
    M_displayer.leaderPrint( "[...done in ", chrono.diff(), "s]\n" );

}


template<typename ImporterType, typename ExporterType>
bool
TestImportExport::run( GetPot& commandLine, const std::string& testString )
{
    using namespace LifeV;

    bool passed( true );

    loadData( commandLine );
    buildMesh();
    buildFESpaces();

    boost::shared_ptr< ExporterType > exporterPtr;
    boost::shared_ptr< ImporterType > importerPtr;
    boost::shared_ptr< ImporterType > importerMeshPtr;

    buildExporter ( importerPtr, M_dataFile("importer/testName", "testImporter" ) );
    buildExporter ( importerMeshPtr, M_dataFile("importer/testMeshName", "testImporterMesh" ) );
    buildExporter ( exporterPtr, M_dataFile("exporter/testName", "testExporter" ) );
    importerPtr->setDataFromGetPot(M_dataFile, "importer");
    importerMeshPtr->setDataFromGetPot(M_dataFile, "importer");
    exporterPtr->setDataFromGetPot(M_dataFile, "exporter");

    // Chronometer
    LifeChrono globalChrono;

    // Set up the IMPORTER
    M_velImportedPtr.reset(
                    new Exporter<mesh_Type >::vector_Type   ( M_velFESpacePtr->map(), Repeated ) );
    M_pressImportedPtr.reset(
                    new Exporter<mesh_Type >::vector_Type   ( M_pressFESpacePtr->map(), Repeated ) );

    ExporterData<mesh_Type>::WhereEnum whereVelocity =
                    ( M_velFESpacePtr->fe().refFE().type() == FE_P0_3D )
                    ?
                     ExporterData<mesh_Type>::Cell : ExporterData<mesh_Type>::Node;

    ExporterData<mesh_Type>::WhereEnum wherePressure =
                    ( M_pressFESpacePtr->fe().refFE().type() == FE_P0_3D )
                    ?
                     ExporterData<mesh_Type>::Cell : ExporterData<mesh_Type>::Node;

    if(M_loadFluid)
    {
        importerPtr->addVariable( ExporterData<mesh_Type>::VectorField, M_velocityName,
                                  M_velFESpacePtr, M_velImportedPtr, UInt(0),
                                  ExporterData<mesh_Type>::UnsteadyRegime,
                                  whereVelocity );
        importerPtr->addVariable( ExporterData<mesh_Type>::ScalarField, M_pressureName,
                                  M_pressFESpacePtr, M_pressImportedPtr, UInt(0),
                                  ExporterData<mesh_Type>::UnsteadyRegime,
                                  wherePressure );
    }

    if(M_moving)
    {
        M_meshDispImportedPtr.reset(
                        new Exporter<mesh_Type >::vector_Type   ( M_velFESpacePtr->map(), Repeated ) );
        M_meshVelImportedPtr.reset(
                        new Exporter<mesh_Type >::vector_Type   ( M_velFESpacePtr->map(), Repeated ) );

        importerMeshPtr->addVariable( ExporterData<mesh_Type>::VectorField, M_meshDisplacementName,
                                      M_velFESpacePtr, M_meshDispImportedPtr, UInt(0) );

        importerMeshPtr->addVariable( ExporterData<mesh_Type>::VectorField, M_meshVelocityName,
                                        M_velFESpacePtr, M_meshVelImportedPtr, UInt(0) );
    }

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    M_displayer.leaderPrint( "[Entering the time loop]\n" );
    globalChrono.start();

    M_timeData.setup( M_dataFile, "time_discretization" );

    M_timeData.updateTime();
    M_displayer.leaderPrint( "\t[t = ", M_timeData.time(), " s.]\n" );

    if( testString.compare( "import" ) == 0 )
        passed = passed && importLoop(importerPtr, importerMeshPtr, exporterPtr);
    else
        passed = passed && exportLoop(importerPtr, exporterPtr);

    globalChrono.stop();
    M_displayer.leaderPrint( "[Time loop, elapsed time:  ", globalChrono.diff(), " s.]\n" );

    exporterPtr->closeFile();

    return passed;
}


template<typename ImporterType, typename ExporterType>
bool
TestImportExport::exportLoop( const boost::shared_ptr< ImporterType > & importerPtr,
                              const boost::shared_ptr< ExporterType > & exporterPtr )
{
    using namespace LifeV;

    bool passed(true);

   // Set up the EXPORTER
    M_velInterpolantPtr.reset(
        new Exporter<mesh_Type >::vector_Type   ( M_velFESpacePtr->map(), Repeated ) );
    M_pressInterpolantPtr.reset(
        new Exporter<mesh_Type >::vector_Type   ( M_pressFESpacePtr->map(), Repeated ) );

    exporterPtr->addVariable( ExporterData<mesh_Type>::VectorField, M_velocityName,
                           M_velFESpacePtr, M_velInterpolantPtr, UInt(0) );

    exporterPtr->addVariable( ExporterData<mesh_Type>::ScalarField, M_pressureName,
                           M_pressFESpacePtr, M_pressInterpolantPtr, UInt(0) );

    //exporterPtr->postProcess( 0 );

    Exporter<mesh_Type >::vector_Type velDiff( M_velFESpacePtr->map(), Repeated );
    Exporter<mesh_Type >::vector_Type pressDiff( M_pressFESpacePtr->map(), Repeated );

    for ( ; M_timeData.canAdvance(); M_timeData.updateTime() )
    {

        M_displayer.leaderPrint( "[t = ", M_timeData.time(), " s.]\n" );

        // Computation of the interpolation
        M_velFESpacePtr->interpolate( problem_Type::uexact, *M_velInterpolantPtr, M_timeData.time() );
        M_pressFESpacePtr->interpolate( problem_Type::pexact, *M_pressInterpolantPtr, M_timeData.time() );

        // Exporting the solution
        exporterPtr->postProcess( M_timeData.time() );

    }
    MPI_Barrier(MPI_COMM_WORLD);

    M_timeData.setup( M_dataFile, "time_discretization" );
    M_timeData.updateTime();

    for ( ; M_timeData.canAdvance(); M_timeData.updateTime() )
    {
        // Computation of the interpolation
        M_velFESpacePtr->interpolate( problem_Type::uexact, *M_velInterpolantPtr, M_timeData.time() );
        M_pressFESpacePtr->interpolate( problem_Type::pexact, *M_pressInterpolantPtr, M_timeData.time() );

        // Importing the solution
        importerPtr->import( M_timeData.time() );

        velDiff = *M_velInterpolantPtr;
        velDiff += (-*M_velImportedPtr);
        Real maxDiff( velDiff.normInf() );

        M_velInterpolantPtr->spy("interpolant");
        M_velImportedPtr->spy("imported");

        M_displayer.leaderPrint( "[velDiff.normInf() = ", velDiff.normInf(), "]\n" );

        passed = passed && (maxDiff < 1.e-6);

    }

    return passed;
}


template<typename ImporterType, typename ExporterType>
bool
TestImportExport::importLoop( const boost::shared_ptr< ImporterType > & importerPtr, const boost::shared_ptr< ImporterType > & importerMeshPtr,
                              const boost::shared_ptr< ExporterType > & exporterPtr )
{
    using namespace LifeV;

    bool passed(true);

    if(M_loadFluid)
    {
    	exporterPtr->addVariable( ExporterData<mesh_Type>::VectorField, M_velocityName,
                              M_velFESpacePtr, M_velImportedPtr, UInt(0) );

    	exporterPtr->addVariable( ExporterData<mesh_Type>::ScalarField, M_pressureName,
                              M_pressFESpacePtr, M_pressImportedPtr, UInt(0) );
    }
    if(M_moving)
    {
        exporterPtr->addVariable( ExporterData<mesh_Type>::VectorField, M_meshDisplacementName,
									  M_velFESpacePtr, M_meshDispImportedPtr, UInt(0) );

        exporterPtr->addVariable( ExporterData<mesh_Type>::VectorField, M_meshVelocityName,
									M_velFESpacePtr, M_meshVelImportedPtr, UInt(0) );
    }

    for( ; M_timeData.canAdvance(); M_timeData.updateTime())
    {
        M_displayer.leaderPrint( "\t[t = ", M_timeData.time(), " s.]\n" );

        importerPtr->import( M_timeData.time() );
        importerMeshPtr->import( M_timeData.time() );

        exporterPtr->postProcess( M_timeData.time() );
    }

    return passed;
}

#endif /* TESTIMPORTEXPORT_HPP_ */
