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
    @brief

    @author Antonio Cervone <ant.cervone@gmail.com>
    @date 05-2012
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

#include <lifev/core/util/LifeChronoManager.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

using namespace LifeV;


Real epsilon (1);

Real exactSolution ( const Real& /* t */, const Real& x, const Real& y, const Real& /* z */, const ID& /* i */ )
{
    return  sin ( x ) + y * y / 2.;
}

Real fRhs ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
    return  sin ( x ) - 1.;
}


typedef RegionMesh<LinearTriangle> mesh_Type;
typedef boost::shared_ptr<Epetra_Comm> commPtr_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
typedef VectorEpetra vector_Type;
typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif

    int returnValue = 0;

    // introducing a local scope in order to properly destroy all objects
    // before calling MPI_Finalize()
    {

#ifdef HAVE_MPI
        commPtr_Type comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
        commPtr_Type comm ( new Epetra_SerialComm );
#endif
        LifeChronoManager<> chronoMgr( comm );

        LifeChrono initTime;
        chronoMgr.add( "Initialization Time", &initTime );
        initTime.start();

        GetPot dataFile ( "data_2d" );
        const bool isLeader ( comm->MyPID() == 0 );
        const bool verbose ( dataFile ( "miscellaneous/verbose", 0 ) && isLeader );

#ifdef HAVE_LIFEV_DEBUG
        std::ofstream debugOut (
            ( "rm." +
              ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( comm->MyPID() ) : "s" ) +
              ".out" ).c_str() );
#else
        std::ofstream debugOut ( "/dev/null" );
#endif

        initTime.stop();

        // Build and partition the mesh

        LifeChrono meshTime;
        chronoMgr.add( "Mesh reading/creation Time", &initTime );
        meshTime.start();

        if ( verbose )
        {
            std::cout << " -- Reading the mesh ... " << std::flush;
        }
        meshPtr_Type fullMeshPtr (new mesh_Type() );
        if ( dataFile ( "mesh/mesh_type", "structured" ) == "structured" )
        {
            regularMesh2D ( *fullMeshPtr, 0,
                            dataFile ( "mesh/nx", 20 ), dataFile ( "mesh/ny", 20 ),
                            dataFile ( "mesh/verbose", false ),
                            dataFile ( "mesh/lx", 1. ), dataFile ( "mesh/ly", 1. ) );
        }
        else
        {
            MeshData meshData (dataFile, "mesh");
            readMesh (*fullMeshPtr, meshData);
        }
        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }
        if ( isLeader )
        {
            std::cout << "mesh elements = " << fullMeshPtr->numElements() << "\n"
                      << "mesh points   = " << fullMeshPtr->numPoints() << std::endl;
        }

        meshTime.stop();

        if ( verbose )
        {
            std::cout << " -- Partitioning the mesh ... " << std::flush;
        }

        LifeChrono partTime;
        chronoMgr.add( "Partition Time", &partTime );
        partTime.start();
        meshPtr_Type localMesh;
        {
            MeshPartitioner< mesh_Type >   meshPart;
            meshPart.doPartition ( fullMeshPtr, comm );
            localMesh = meshPart.meshPartition();
        }
        partTime.stop();
        //if ( isLeader )
        //{
        //    std::cout << "partitioning time  = " << partTime.diff() << std::endl;
        //}

        Int localMeshNum[ 2 ];
        localMeshNum[ 0 ] = localMesh->numElements();
        localMeshNum[ 1 ] = localMesh->numPoints();
        Int maxMeshNum[ 2 ] = { 0, 0 };
        comm->MaxAll( localMeshNum, maxMeshNum, 2 );

        if ( isLeader )
        {
            std::cout << "part mesh elements = " << maxMeshNum[ 0 ] << "\n"
                      << "part mesh points   = " << maxMeshNum[ 1 ] << std::endl;
        }

        LifeChrono partTimeR;
        chronoMgr.add( "Partition Time (R)", &partTimeR );
        partTimeR.start();
        meshPtr_Type localMeshR;
        {
            MeshPartitioner< mesh_Type >   meshPartR;
            meshPartR.setPartitionOverlap ( 1 );
            meshPartR.doPartition ( fullMeshPtr, comm );
            localMeshR = meshPartR.meshPartition();
        }
        partTimeR.stop();
        //if ( isLeader )
        //{
        //    std::cout << "partitioningR time = " << partTimeR.diff() << std::endl;
        //}

        localMeshNum[ 0 ] = localMeshR->numElements();
        localMeshNum[ 1 ] = localMeshR->numPoints();
        maxMeshNum[ 0 ] = 0;
        maxMeshNum[ 1 ] = 0;
        comm->MaxAll( localMeshNum, maxMeshNum, 2 );

        if ( isLeader )
        {
            std::cout << "part mesh elements (R) = " << maxMeshNum[ 0 ] << "\n"
                      << "part mesh points   (R) = " << maxMeshNum[ 1 ] << std::endl;
        }

        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }

        if ( verbose )
        {
            std::cout << " -- Freeing the global mesh ... " << std::flush;
        }
        fullMeshPtr.reset();
#ifdef HAVE_LIFEV_DEBUG
        ASSERT( fullMeshPtr.use_count() == 0, "full mesh not properly freed." );
#endif
        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }

        // Build the FESpaces

        if ( verbose )
        {
            std::cout << " -- Building FESpaces ... " << std::flush;
        }
        std::string uOrder ( dataFile ( "fe/type", "P1" ) );
        LifeChrono feSpaceTime;
        chronoMgr.add( "FESpace creation Time", &feSpaceTime );
        feSpaceTime.start();
        feSpacePtr_Type uFESpace ( new feSpace_Type ( localMesh, uOrder, 1, comm ) );
        feSpaceTime.stop();

        LifeChrono feSpaceTimeR;
        chronoMgr.add( "FESpace creation Time (R)", &feSpaceTimeR );
        feSpaceTimeR.start();
        feSpacePtr_Type uFESpaceR ( new feSpace_Type ( localMeshR, uOrder, 1, comm ) );
        feSpaceTimeR.stop();

        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }
        if ( verbose )
        {
            std::cout << " ---> Dofs: " << uFESpaceR->dof().numTotalDof() << std::endl;
        }

        // Build the assembler and the matrices

        if ( verbose )
        {
            std::cout << " -- Building assembler ... " << std::flush;
        }

        LifeChrono matTime;
        chronoMgr.add( "Matrix initialization Time", &matTime );
        matTime.start();
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        adrAssembler.setFespace ( uFESpace );
        matrixPtr_Type systemMatrix ( new matrix_Type ( uFESpace->map() ) );
        matTime.stop();

        LifeChrono matTimeR;
        chronoMgr.add( "Matrix initialization Time (R)", &matTimeR );
        matTimeR.start();
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssemblerR;
        adrAssemblerR.setFespace ( uFESpaceR );
        matrixPtr_Type systemMatrixR ( new matrix_Type ( uFESpaceR->map(), 50, true ) );
        matTimeR.stop();

        if ( verbose )
        {
            std::cout << " done! " << std::endl;
        }

        // Perform the assembly of the matrix

        if ( verbose )
        {
            std::cout << " -- Adding the diffusion ... " << std::flush;
        }

        if ( verbose )
        {
            std::cout << " done! " << std::endl;
        }

        UInt timeSteps = dataFile ( "miscellaneous/timesteps", 1 );

        LifeChrono assemblyTime;
        chronoMgr.add( "Assembly Time", &assemblyTime );
        assemblyTime.start();
        for( UInt n = 0; n < timeSteps; n++ )
        {
            systemMatrix->zero();
            adrAssembler.addDiffusion ( systemMatrix, epsilon );
            adrAssembler.addMass ( systemMatrix, epsilon );
            systemMatrix->globalAssemble();
        }
        assemblyTime.stop();
        if ( isLeader )
        {
            std::cout << "assembly time  = " << assemblyTime.diff() << std::endl;
        }

        LifeChrono assemblyTimeR;
        chronoMgr.add( "Assembly Time (R)", &assemblyTimeR );
        assemblyTimeR.start();
        for( UInt n = 0; n < timeSteps; n++ )
        {
            systemMatrixR->zero();
            adrAssemblerR.addDiffusion ( systemMatrixR, epsilon );
            adrAssemblerR.addMass ( systemMatrixR, epsilon );
            systemMatrixR->fillComplete();
        }
        assemblyTimeR.stop();
        if ( isLeader )
        {
            std::cout << "assemblyR time = " << assemblyTimeR.diff() << std::endl;
        }

        //if ( verbose )
        //{
        //    std::cout << " Time needed : " << adrAssembler.diffusionAssemblyChrono().diffCumul() << std::endl;
        //}
        //if ( verbose )
        //{
        //    std::cout << " Time needed : " << adrAssemblerR.diffusionAssemblyChrono().diffCumul() << std::endl;
        //}

        // SPY
        //systemMatrix->spy("matrixNoBC");
        //systemMatrixR->spy("matrixNoBCR");

        // check that the assembled matrices are the same

        LifeChrono checkMatTime;
        chronoMgr.add( "Check (Matrix) Time", &checkMatTime );
        checkMatTime.start();
        matrix_Type matrixDiff ( *systemMatrix );
        matrixDiff -= *systemMatrixR;

        Real diff = matrixDiff.normInf();
        checkMatTime.stop();

        if ( isLeader )
        {
            std::cout << "Norm of the difference between the 2 matrices = " << diff << std::endl;
        }

        if ( verbose )
        {
            std::cout << " -- Building the RHS ... " << std::flush;
        }

        LifeChrono rhsTime;
        chronoMgr.add( "Rhs build Time", &rhsTime );
        rhsTime.start();
        vector_Type rhs ( uFESpace->map(), Unique );
        adrAssembler.addMassRhs ( rhs, fRhs, 0. );
        rhs.globalAssemble();
        rhsTime.stop();

        LifeChrono rhsTimeR;
        chronoMgr.add( "Rhs build Time (R)", &rhsTimeR );
        rhsTimeR.start();
        vector_Type rhsR ( uFESpaceR->map(), Unique, Zero );
        adrAssemblerR.addMassRhs ( rhsR, fRhs, 0. );
        rhsR.globalAssemble ();
        rhsTimeR.stop();

        LifeChrono checkVecTime;
        chronoMgr.add( "Check (Vector) Time", &checkVecTime );
        checkVecTime.start();
        vector_Type vectorDiff ( rhs );
        vectorDiff -= rhsR;

        Real diffNormV = vectorDiff.normInf();
        checkVecTime.stop();

        if ( isLeader )
        {
            std::cout << "Norm of the difference between the 2 vectors = " << diffNormV << std::endl;
        }

        diff += diffNormV;

        vector_Type rhs2 ( uFESpace->map(), Unique );
        vector_Type f ( uFESpace->map(), Repeated );
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), f, 0.0 );
        adrAssembler.addMassRhs ( rhs2, f );
        rhs2.globalAssemble();

        vector_Type rhs2R ( uFESpaceR->map(), Unique );
        vector_Type fR ( uFESpaceR->map(), Repeated );
        uFESpaceR->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), fR, 0.0 );
        adrAssemblerR.addMassRhs ( rhs2R, fR );
        rhs2R.globalAssemble ( Zero );

        vector_Type vectorDiff2 ( rhs2 );
        vectorDiff2 -= rhs2R;

        Real diffNormV2 = vectorDiff2.normInf();

        if ( isLeader )
        {
            std::cout << "Norm of the difference between the 2 vectors = " << diffNormV2 << std::endl;
        }

        diff += diffNormV2;

        // test exporting of a repeated mesh
#ifdef HAVE_HDF5
        ExporterHDF5<mesh_Type> exporter ( dataFile, localMeshR, "pid_2d", comm->MyPID() );
        exporter.exportPID ( localMeshR, comm, true );
        exporter.postProcess ( 0. );
#endif

        vector_Type rhsCopy ( rhs, Repeated );
        vector_Type rhsCopyR ( rhsR, Repeated, Add );
        Real l2Error  = uFESpace->l2Error ( fRhs, rhsCopy, 0.0 );
        Real l2ErrorR = uFESpaceR->l2Error ( fRhs, rhsCopyR, 0.0 );
        Real diffL2Error = std::fabs ( l2Error - l2ErrorR );

        if ( isLeader )
        {
            std::cout << "difference in l2error  = " << diffL2Error << std::endl;
        }

        diff += diffL2Error;

        if ( diff < 1.e-14 )
        {
            returnValue = EXIT_SUCCESS;
            if ( isLeader )
            {
                std::cout << "End Result: TEST PASSED" << std::endl;
            }
        }
        else
        {
            returnValue = EXIT_FAILURE;
        }

        // print out times
        chronoMgr.print();
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return returnValue;
}
