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
    @brief TestImportExport

    @author Tiziano Passerini <tiziano@mathcs.emory.edu>
    @contributor
    @maintainer

    @date 05-2011

 */

#ifndef TESTIMPORTEXPORT_HPP_
#define TESTIMPORTEXPORT_HPP_ 1

// LifeV definition files
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes/function/Womersley.hpp>
#include <lifev/navier_stokes/function/RossEthierSteinmanDec.hpp>


// Trilinos-MPI communication definitions
#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


class TestImportExport
{
public:

  // Object type definitions
  typedef LifeV::RegionMesh<LifeV::LinearTetra>         mesh_Type;
  typedef LifeV::RossEthierSteinmanUnsteadyDec          problem_Type;
  typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
  typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;
  typedef boost::shared_ptr<Epetra_Comm>                commPtr_Type;
  typedef LifeV::Exporter<mesh_Type >::vectorPtr_Type   vectorPtr_Type;
  typedef boost::function < LifeV::Real ( LifeV::Real const&,
					  LifeV::Real const&,
					  LifeV::Real const&,
					  LifeV::Real const&,
					  LifeV::UInt const& ) > function_Type;

    TestImportExport ( const commPtr_Type& commPtr );

    template<typename ImporterType, typename ExporterType>
    bool run ( GetPot& commandLine, const std::string& testString = "import" );

    template<typename ImporterType, typename ExporterType>
    bool importLoop ( const boost::shared_ptr< ImporterType >& importerPtr,
                      const boost::shared_ptr< ExporterType >& exporterPtr );

    template<typename ImporterType, typename ExporterType>
    bool exportLoop ( const boost::shared_ptr< ImporterType >& importerPtr,
                      const boost::shared_ptr< ExporterType >& exporterPtr );

private:
    void loadData ( GetPot& commandLine );
    void buildMesh();
    void buildFESpaces();
    template<typename ExporterType>
    void buildExporter ( boost::shared_ptr< ExporterType >& exporterPtr,
                         const std::string& prefix );

    commPtr_Type                                             M_commPtr;
    LifeV::Displayer                                         M_displayer;
    GetPot                                                   M_dataFile;
    boost::shared_ptr< mesh_Type >                           M_meshPtr;
    LifeV::TimeData                                          M_timeData;

    feSpacePtr_Type                                          M_vectorFESpacePtr;
    feSpacePtr_Type                                          M_scalarFESpacePtr;

    std::vector< vectorPtr_Type >                            M_vectorInterpolantPtr;
    std::vector< vectorPtr_Type >                            M_scalarInterpolantPtr;

    std::vector< vectorPtr_Type >                            M_vectorImportedPtr;
    std::vector< vectorPtr_Type >                            M_scalarImportedPtr;

    std::vector< std::string >                               M_vectorName;
    std::vector< std::string >                               M_scalarName;
};


TestImportExport::TestImportExport ( const commPtr_Type& commPtr ) :
    M_commPtr  ( commPtr ),
    M_displayer ( commPtr )
{}


void
TestImportExport::loadData ( GetPot& commandLine )
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    M_displayer.leaderPrint ( "[Loading the data]\n" );
    chrono.start();

    const std::string problemName = commandLine.follow ("", 2, "-p", "--problem");
    const std::string defaultDataFileName ("data" + problemName);
    const std::string dataFileName = commandLine.follow (defaultDataFileName.c_str(), 2, "-f", "--file");
    M_dataFile = GetPot (dataFileName);

    UInt numVectors = M_dataFile ("importer/numVectors", 0);
    UInt numScalars = M_dataFile ("importer/numScalars", 0);

    M_vectorInterpolantPtr.resize (numVectors);
    M_scalarInterpolantPtr.resize (numScalars);

    M_vectorImportedPtr.resize (numVectors);
    M_scalarImportedPtr.resize (numScalars);

    M_vectorName.resize (numVectors);
    M_scalarName.resize (numScalars);;

    for ( UInt iVec = 0; iVec < numVectors; ++iVec )
    {
        std::stringstream ss;
        ss.str ("");
        ss << "importer/vector" << iVec << "Name";
        M_vectorName[iVec] = M_dataFile (ss.str().c_str(), "");
    }
    for ( UInt iScal = 0; iScal < numScalars; ++iScal )
    {
        std::stringstream ss;
        ss.str ("");
        ss << "importer/scalar" << iScal << "Name";
        M_scalarName[iScal] = M_dataFile (ss.str().c_str(), "");
    }

    problem_Type::setParamsFromGetPot (M_dataFile);

    chrono.stop();
    M_displayer.leaderPrint ( "[...done in ", chrono.diff(), "s]\n" );
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
    M_displayer.leaderPrint ( "[Building the mesh]\n" );
    chrono.start();
    boost::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( M_commPtr ) );

    if ( M_dataFile ("space_discretization/mesh_from_file", false) )
    {
        //  The following is when reading from file
        MeshData meshData;
        meshData.setup (M_dataFile, "space_discretization");
        //if (verbose) std::cout << "Mesh file: " << meshData.meshDir() << meshData.meshFile() << std::endl;
        readMesh (*fullMeshPtr, meshData);
    }
    else
    {
        UInt nEl (M_dataFile ("space_discretization/dimension", 1) );
        regularMesh3D (*fullMeshPtr, 0, nEl, nEl, nEl);
    }
    // Split the mesh between processors
    MeshPartitioner<mesh_Type> meshPart ( fullMeshPtr, M_commPtr );
    // Get the mesh for the current partition
    M_meshPtr = meshPart.meshPartition();
    // Release the original mesh from the MeshPartitioner object and delete the RegionMesh3D object
    meshPart.releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    chrono.stop();
    M_displayer.leaderPrint ( "[...done in ", chrono.diff(), "s]\n" );
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
    M_displayer.leaderPrint ( "[Creating the FE spaces...]\n" );
    chrono.start();

    const std::string vectorFE =  M_dataFile ( "space_discretization/vector_fespace", "P2");
    M_displayer.leaderPrint ( "\t-o FE for the vector: ", vectorFE, "\n" );

    M_displayer.leaderPrint ( "\t-o Building the vector FE space...\n" );

    M_vectorFESpacePtr.reset ( new feSpace_Type ( M_meshPtr, vectorFE, nDimensions, M_commPtr ) );

    M_displayer.leaderPrint ( "\t\t...ok.\n" );

    const std::string scalarFE =  M_dataFile ( "space_discretization/scalar_fespace", "P1");
    M_displayer.leaderPrint ( "\t-o FE for the scalar: ", scalarFE, "\n" );

    M_displayer.leaderPrint ( "\t-o Building the scalar FE space...\n" );

    M_scalarFESpacePtr.reset ( new feSpace_Type ( M_meshPtr, scalarFE, 1, M_commPtr ) );

    M_displayer.leaderPrint ( "\t\t...ok.\n" );

    // Total degrees of freedom (elements of matrix)
    UInt vectorTotalDof   = M_vectorFESpacePtr->map().map (Unique)->NumGlobalElements();
    UInt scalarTotalDof = M_scalarFESpacePtr->map().map (Unique)->NumGlobalElements();

    M_displayer.leaderPrint ( "\t-o Total Dof for the vector: ", vectorTotalDof, "\n" );
    M_displayer.leaderPrint ( "\t-o Total Dof for the scalar: ", scalarTotalDof, "\n" );

    chrono.stop();
    M_displayer.leaderPrint ( "[...done in ", chrono.diff(), "s]\n" );

}


template<typename ExporterType>
void
TestImportExport::buildExporter ( boost::shared_ptr< ExporterType >& exporterPtr,
                                  const std::string& prefix )
{
    using namespace LifeV;
    // Chronometer
    LifeChrono chrono;

    // +-----------------------------------------------+
    // |            Creating the exporter              |
    // +-----------------------------------------------+
    M_displayer.leaderPrint ( "[Creating the exporter...]\n" );
    chrono.start();

    exporterPtr.reset ( new ExporterType ( M_dataFile, prefix ) );

    //exporterPtr->setPostDir( "./" );
    exporterPtr->setMeshProcId ( M_meshPtr, M_commPtr->MyPID() );
    chrono.stop();
    M_displayer.leaderPrint ( "[...done in ", chrono.diff(), "s]\n" );

}


template<typename ImporterType, typename ExporterType>
bool
TestImportExport::run ( GetPot& commandLine, const std::string& testString )
{
    using namespace LifeV;

    bool passed ( true );

    loadData ( commandLine );
    buildMesh();
    buildFESpaces();

    boost::shared_ptr< ExporterType > exporterPtr;
    boost::shared_ptr< ImporterType > importerPtr;

    buildExporter ( importerPtr, M_dataFile ("importer/prefix", "testImporter" ) );
    buildExporter ( exporterPtr, M_dataFile ("exporter/prefix", "testExporter" ) );

    importerPtr->setDataFromGetPot (M_dataFile, "importer");
    exporterPtr->setDataFromGetPot (M_dataFile, "exporter");

    // Chronometer
    LifeChrono globalChrono;

    // Set up the IMPORTER
    ExporterData<mesh_Type>::WhereEnum whereVector =
        ( M_vectorFESpacePtr->fe().refFE().type() == FE_P0_3D )
        ?
        ExporterData<mesh_Type>::Cell : ExporterData<mesh_Type>::Node;

    ExporterData<mesh_Type>::WhereEnum whereScalar =
        ( M_scalarFESpacePtr->fe().refFE().type() == FE_P0_3D )
        ?
        ExporterData<mesh_Type>::Cell : ExporterData<mesh_Type>::Node;

    for ( UInt iVec (0); iVec < M_vectorImportedPtr.size(); ++iVec )
    {
        M_vectorImportedPtr[iVec].reset (
            new Exporter<mesh_Type >::vector_Type   ( M_vectorFESpacePtr->map(), ExporterType::MapType ) );

        importerPtr->addVariable ( ExporterData<mesh_Type>::VectorField, M_vectorName[iVec],
                                   M_vectorFESpacePtr, M_vectorImportedPtr[iVec], UInt (0),
                                   ExporterData<mesh_Type>::UnsteadyRegime,
                                   whereVector );
    }
    for ( UInt iScal (0); iScal < M_scalarImportedPtr.size(); ++iScal )
    {
        M_scalarImportedPtr[iScal].reset (
            new Exporter<mesh_Type >::vector_Type   ( M_scalarFESpacePtr->map(), ExporterType::MapType ) );
        importerPtr->addVariable ( ExporterData<mesh_Type>::ScalarField, M_scalarName[iScal],
                                   M_scalarFESpacePtr, M_scalarImportedPtr[iScal], UInt (0),
                                   ExporterData<mesh_Type>::UnsteadyRegime,
                                   whereScalar );
    }

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    M_displayer.leaderPrint ( "[Entering the time loop]\n" );
    globalChrono.start();

    M_timeData.setup ( M_dataFile, "time_discretization" );

    M_timeData.updateTime();
    M_displayer.leaderPrint ( "\t[t = ", M_timeData.time(), " s.]\n" );


    // Exporting
    if ( testString.compare ( "export" ) == 0 )
    {
    passed = passed && exportLoop (importerPtr, exporterPtr);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    M_displayer.leaderPrint ( "Exporting finished. Now importing\n" );

    // Exporting and checking
    passed = passed && importLoop (importerPtr, exporterPtr);

    globalChrono.stop();
    M_displayer.leaderPrint ( "[Time loop, elapsed time:  ", globalChrono.diff(), " s.]\n" );

    return passed;
}


template<typename ImporterType, typename ExporterType>
bool
TestImportExport::exportLoop ( const boost::shared_ptr< ImporterType >& importerPtr,
                               const boost::shared_ptr< ExporterType >& exporterPtr )
{
    using namespace LifeV;

    bool passed (true);

    ASSERT ( M_vectorImportedPtr.size() + M_scalarImportedPtr.size(), "There's no data on which to work!" )

    const UInt vectorImportedPtrSize = M_vectorImportedPtr.size();
    const UInt scalarImportedPtrSize = M_scalarImportedPtr.size();

    // Set up the EXPORTER
    for ( UInt iVec (0); iVec < vectorImportedPtrSize; ++iVec )
    {
        M_vectorInterpolantPtr[iVec].reset (
            new Exporter<mesh_Type >::vector_Type   ( M_vectorFESpacePtr->map(), ExporterType::MapType ) );
        exporterPtr->addVariable ( ExporterData<mesh_Type>::VectorField, M_vectorName[iVec],
                                   M_vectorFESpacePtr, M_vectorInterpolantPtr[iVec], UInt (0) );
    }
    for ( UInt iScal (0); iScal < scalarImportedPtrSize; ++iScal )
    {
        M_scalarInterpolantPtr[iScal].reset (
            new Exporter<mesh_Type >::vector_Type   ( M_scalarFESpacePtr->map(), ExporterType::MapType ) );
        exporterPtr->addVariable ( ExporterData<mesh_Type>::ScalarField, M_scalarName[iScal],
                                   M_scalarFESpacePtr, M_scalarInterpolantPtr[iScal], UInt (0) );
    }

    //exporterPtr->postProcess( 0 );

    for ( ; M_timeData.canAdvance(); M_timeData.updateTime() )
    {

        M_displayer.leaderPrint ( "[t = ", M_timeData.time(), " s.]\n" );

        // Computation of the interpolation
        if ( vectorImportedPtrSize )
        {
            M_vectorFESpacePtr->interpolate ( static_cast<function_Type> ( problem_Type::uexact ), *M_vectorInterpolantPtr[0], M_timeData.time() );
        }
        if ( scalarImportedPtrSize )
        {
            M_scalarFESpacePtr->interpolate ( static_cast<function_Type> ( problem_Type::pexact ), *M_scalarInterpolantPtr[0], M_timeData.time() );
        }

        // Exporting the solution
        exporterPtr->postProcess ( M_timeData.time() );

    }

    exporterPtr->closeFile();

    return passed;

}


template<typename ImporterType, typename ExporterType>
bool
TestImportExport::importLoop ( const boost::shared_ptr< ImporterType >& importerPtr,
                               const boost::shared_ptr< ExporterType >& exporterPtr )
{
    using namespace LifeV;

    bool passed (true);

    ASSERT ( M_vectorImportedPtr.size() + M_scalarImportedPtr.size(), "There's no data on which to work!" )

    const UInt vectorImportedPtrSize = M_vectorImportedPtr.size();
    const UInt scalarImportedPtrSize = M_scalarImportedPtr.size();

    Exporter<mesh_Type >::vector_Type vectorDiff ( M_vectorFESpacePtr->map(), ExporterType::MapType );
    Exporter<mesh_Type >::vector_Type scalarDiff ( M_scalarFESpacePtr->map(), ExporterType::MapType );

    // Set up the EXPORTER
    for ( UInt iVec (0); iVec < vectorImportedPtrSize; ++iVec )
    {
        M_vectorInterpolantPtr[iVec].reset (
            new Exporter<mesh_Type >::vector_Type   ( M_vectorFESpacePtr->map(), ExporterType::MapType ) );
    }
    for ( UInt iScal (0); iScal < scalarImportedPtrSize; ++iScal )
    {
        M_scalarInterpolantPtr[iScal].reset (
            new Exporter<mesh_Type >::vector_Type   ( M_scalarFESpacePtr->map(), ExporterType::MapType ) );
    }

    M_timeData.setup ( M_dataFile, "time_discretization" );
    M_timeData.updateTime();

    for ( ; M_timeData.canAdvance(); M_timeData.updateTime() )
    {
        // Computation of the interpolation
        if ( vectorImportedPtrSize )
        {
            M_vectorFESpacePtr->interpolate ( static_cast<function_Type> ( problem_Type::uexact ), *M_vectorInterpolantPtr[0], M_timeData.time() );
        }
        if ( scalarImportedPtrSize )
        {
            M_scalarFESpacePtr->interpolate ( static_cast<function_Type> ( problem_Type::pexact ), *M_scalarInterpolantPtr[0], M_timeData.time() );
        }

        // Importing the solution
        importerPtr->import ( M_timeData.time() );

        Real maxDiff (1.e6);
        if ( vectorImportedPtrSize )
        {
            vectorDiff = *M_vectorInterpolantPtr[0];
            vectorDiff += (-*M_vectorImportedPtr[0]);
            maxDiff = vectorDiff.normInf();

            M_vectorInterpolantPtr[0]->spy ("interpolant");
            M_vectorImportedPtr[0]->spy ("imported");

            M_displayer.leaderPrint ( "[vectorDiff.normInf() = ", vectorDiff.normInf(), "]\n" );
        }
        if ( scalarImportedPtrSize )
        {
            scalarDiff = *M_scalarInterpolantPtr[0];
            scalarDiff += (-*M_scalarImportedPtr[0]);
            maxDiff = std::max (scalarDiff.normInf(), maxDiff);

            M_scalarInterpolantPtr[0]->spy ("interpolant");
            M_scalarImportedPtr[0]->spy ("imported");

            M_displayer.leaderPrint ( "[scalarDiff.normInf() = ", scalarDiff.normInf(), "]\n" );
        }
        passed = passed && (maxDiff < 1.e-4);

    }

    return passed;
}

#endif /* TESTIMPORTEXPORT_HPP_ */
