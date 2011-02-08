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
typedef LifeV::RegionMesh3D<LifeV::LinearTetra>       mesh_Type;
typedef LifeV::RossEthierSteinmanUnsteadyInc          problem_Type;
typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;

template<typename ProblemType>
bool testExporterVTK(/*const*/ boost::shared_ptr<Epetra_Comm> & comm, GetPot & commandLine)
{
    using namespace LifeV;

    // Chronometer
    LifeChrono globalChrono;
    LifeChrono initChrono;
    LifeChrono iterChrono;
    globalChrono.start();
    initChrono.start();

    Displayer displayer(comm);

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    displayer.leaderPrint( "[Loading the data]\n" );
    const std::string problemName = commandLine.follow("", 2, "-p", "--problem");
    const std::string defaultDataFileName("data"+problemName);
    const std::string dataFileName = commandLine.follow(defaultDataFileName.c_str(), 2, "-f","--file");
    GetPot dataFile(dataFileName);

    problem_Type::setParamsFromGetPot(dataFile);

    // +-----------------------------------------------+
    // |              Building the mesh                |
    // +-----------------------------------------------+
    displayer.leaderPrint( "[Building the mesh]\n" );

    boost::shared_ptr< mesh_Type > fullMeshPtr(new mesh_Type);
    UInt nEl(dataFile("space_discretization/dimension", 1));
    regularMesh3D(*fullMeshPtr, 0, nEl, nEl, nEl);
    //  The following is when reading from file
    /*
    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");
    if (verbose) std::cout << "Mesh file: " << dataMesh.meshDir() << dataMesh.meshFile() << std::endl;
    readMesh(*fullMeshPtr, dataMesh);
    */
    // Split the mesh between processors
    MeshPartitioner< mesh_Type > meshPart( fullMeshPtr, comm );
    // Release the original mesh from the MeshPartitioner object and delete the RegionMesh3D object
    meshPart.releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    displayer.leaderPrint( "[Creating the FE spaces]\n" );

    const std::string velFE =  dataFile( "space_discretization/velFE", "P2");
    displayer.leaderPrint( "\tFE for the velocity: ", velFE );

    displayer.leaderPrint( "\tBuilding the velocity FE space... " );

    feSpacePtr_Type velFESpacePtr( new feSpace_Type(meshPart, velFE, nDimensions, comm) );

    displayer.leaderPrint( "\tok." );

    const std::string pressFE =  dataFile( "space_discretization/pressFE", "P1");
    displayer.leaderPrint( "\tFE for the pressure: ", velFE );

    displayer.leaderPrint( "\tBuilding the pressure FE space... " );

    feSpacePtr_Type pressFESpacePtr( new feSpace_Type(meshPart, pressFE, nDimensions, comm) );

    displayer.leaderPrint( "\tok." );

    // Total degrees of freedom (elements of matrix)
    UInt velTotalDof   = velFESpacePtr->map().map(Unique)->NumGlobalElements();
    UInt pressTotalDof = pressFESpacePtr->map().map(Unique)->NumGlobalElements();

    displayer.leaderPrint( "\tTotal Dof for the velocity: ", velTotalDof );
    displayer.leaderPrint( "\tTotal Dof for the pressure: ", pressTotalDof );

    // +-----------------------------------------------+
    // |            Creating the exporter              |
    // +-----------------------------------------------+
    Exporter<mesh_Type >::vectorPtr_Type velInterpolantPtr(
        new Exporter<mesh_Type >::vector_Type   ( velFESpacePtr->map(), Repeated ) );
    Exporter<mesh_Type >::vectorPtr_Type pressInterpolantPtr(
        new Exporter<mesh_Type >::vector_Type   ( pressFESpacePtr->map(), Repeated ) );

    boost::shared_ptr< ExporterVTK<mesh_Type > > exporter;
    exporter.reset( new ExporterVTK<mesh_Type > ( dataFile, "testExporterVTK" ) );
    exporter->setPostDir( "./" );
    exporter->setMeshProcId( meshPart.meshPartition(), comm->MyPID() );

    exporter->addVariable( ExporterData<mesh_Type>::VectorField, "velocity",
                           velFESpacePtr, velInterpolantPtr, UInt(0) );

    exporter->addVariable( ExporterData<mesh_Type>::ScalarField, "pressure",
                           pressFESpacePtr, pressInterpolantPtr, UInt(0) );
    exporter->postProcess( 0 );

    initChrono.stop();
    displayer.leaderPrint( "[Initialization time]", initChrono.diff() );

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    displayer.leaderPrint( "[Interpolating the analytic solution]\n" );

    TimeData timeData;
    timeData.setup( dataFile, "time_discretization" );

    for ( ; timeData.canAdvance(); timeData.updateTime() )
    {

        displayer.leaderPrint( "\t[t = ", timeData.time(), " s.]\n" );

        // Computation of the interpolation
        velFESpacePtr->interpolate( problem_Type::uexact, *velInterpolantPtr, timeData.time() );
        pressFESpacePtr->interpolate( problem_Type::pexact, *pressInterpolantPtr, timeData.time() );

        // Exporting the solution
        exporter->postProcess( timeData.time() );

        MPI_Barrier(MPI_COMM_WORLD);

    }

    globalChrono.stop();
    displayer.leaderPrint( "Total simulation time:  ", globalChrono.diff(), " s.\n" );

    exporter->closeFile();

    return 0;
}



#endif /* TESTEXPORTERVTK_HPP_ */
