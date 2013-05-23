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
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif

    // introducing a local scope in order to properly destroy all objects
    // before calling MPI_Finalize()
    {

#ifdef HAVE_MPI
        boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
        boost::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm );
#endif

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

        // Build and partition the mesh

        if ( verbose )
        {
            std::cout << " -- Reading the mesh ... " << std::flush;
        }
        boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type() );
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

        if ( verbose )
        {
            std::cout << " -- Partitioning the mesh ... " << std::flush;
        }

        LifeChrono partTime;
        partTime.start();
        boost::shared_ptr< mesh_Type > localMesh;
        {
            MeshPartitioner< mesh_Type >   meshPart;
            meshPart.doPartition ( fullMeshPtr, comm );
            localMesh = meshPart.meshPartition();
        }
        partTime.stop();
        if ( isLeader )
        {
            std::cout << "partitioning time  = " << partTime.diff() << std::endl;
        }
        if ( isLeader )
        {
            std::cout << "part mesh elements = " << localMesh->numElements() << "\n"
                      << "part mesh points   = " << localMesh->numPoints() << std::endl;
        }

        // localMesh->mesh_Type::showMe( true, debugOut );

        LifeChrono partTimeR;
        partTimeR.start();
        boost::shared_ptr< mesh_Type > localMeshR;
        {
            MeshPartitioner< mesh_Type >   meshPartR;
            meshPartR.setPartitionOverlap ( 1 );
            meshPartR.doPartition ( fullMeshPtr, comm );
            localMeshR = meshPartR.meshPartition();
        }
        partTimeR.stop();
        if ( isLeader )
        {
            std::cout << "partitioningR time = " << partTimeR.diff() << std::endl;
        }

        // debugOut << "============================" << std::endl;
        // localMeshR->mesh_Type::showMe( true, debugOut );
        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }

        if ( verbose )
        {
            std::cout << " -- Freeing the global mesh ... " << std::flush;
        }
        fullMeshPtr.reset();
        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }

        // Build the FESpaces

        if ( verbose )
        {
            std::cout << " -- Building FESpaces ... " << std::flush;
        }
        std::string uOrder ( "P1" );
        feSpacePtr_Type uFESpace ( new feSpace_Type ( localMesh, uOrder, 1, comm ) );
        feSpacePtr_Type uFESpaceR ( new feSpace_Type ( localMeshR, uOrder, 1, comm ) );
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
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssemblerR;
        if ( verbose )
        {
            std::cout << " done! " << std::endl;
        }

        if ( verbose )
        {
            std::cout << " -- Setting up assembler ... " << std::flush;
        }
        adrAssembler.setFespace ( uFESpace );
        adrAssemblerR.setFespace ( uFESpaceR );
        if ( verbose )
        {
            std::cout << " done! " << std::endl;
        }

        if ( verbose )
        {
            std::cout << " -- Defining the matrix ... " << std::flush;
        }
        boost::shared_ptr<matrix_Type> systemMatrix ( new matrix_Type ( uFESpace->map() ) );
        boost::shared_ptr<matrix_Type> systemMatrixR ( new matrix_Type ( uFESpaceR->map(), 50, true ) );
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

        LifeChrono assemblyTime;
        assemblyTime.start();
        adrAssembler.addDiffusion ( systemMatrix, epsilon );

        if ( verbose )
        {
            std::cout << " -- Closing the matrix ... " << std::flush;
        }
        systemMatrix->globalAssemble();
        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }
        assemblyTime.stop();
        if ( isLeader )
        {
            std::cout << "assembly time  = " << assemblyTime.diff() << std::endl;
        }
        LifeChrono assemblyTimeR;
        assemblyTimeR.start();
        adrAssemblerR.addDiffusion ( systemMatrixR, epsilon );
        systemMatrixR->fillComplete();
        assemblyTimeR.stop();
        if ( isLeader )
        {
            std::cout << "assemblyR time = " << assemblyTimeR.diff() << std::endl;
        }

        if ( verbose )
        {
            std::cout << " Time needed : " << adrAssembler.diffusionAssemblyChrono().diffCumul() << std::endl;
        }
        if ( verbose )
        {
            std::cout << " Time needed : " << adrAssemblerR.diffusionAssemblyChrono().diffCumul() << std::endl;
        }

        // SPY
        //systemMatrix->spy("matrixNoBC");
        //systemMatrixR->spy("matrixNoBCR");

        // check that the assembled matrices are the same

        matrix_Type matrixDiff ( *systemMatrix );
        matrixDiff -= *systemMatrixR;

        Real diff = matrixDiff.normInf();

        if ( isLeader )
        {
            std::cout << "Norm of the difference between the 2 matrices = " << diff << std::endl;
        }

        if ( verbose )
        {
            std::cout << " -- Building the RHS ... " << std::flush;
        }

        vector_Type rhs ( uFESpace->map(), Unique );
        vector_Type rhsR ( uFESpaceR->map(), Unique, Zero );

        adrAssembler.addMassRhs ( rhs, fRhs, 0. );
        rhs.globalAssemble();

        adrAssemblerR.addMassRhs ( rhsR, fRhs, 0. );
        rhsR.globalAssemble ();

        vector_Type vectorDiff ( rhs );
        vectorDiff -= rhsR;

        Real diffNormV = vectorDiff.normInf();

        if ( isLeader )
        {
            std::cout << "Norm of the difference between the 2 vectors = " << diffNormV << std::endl;
        }

        diff += diffNormV;

        vector_Type rhs2 ( uFESpace->map(), Unique );
        vector_Type rhs2R ( uFESpaceR->map(), Unique );
        vector_Type f ( uFESpace->map(), Repeated );
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), f, 0.0 );
        adrAssembler.addMassRhs ( rhs2, f );
        rhs2.globalAssemble();
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
        ExporterHDF5<mesh_Type> exporter ( dataFile, localMeshR, "pid", comm->MyPID() );
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
            if ( isLeader )
            {
                std::cout << "End Result: TEST PASSED" << std::endl;
            }
        }
        else
        {
            return EXIT_FAILURE;
        }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}
