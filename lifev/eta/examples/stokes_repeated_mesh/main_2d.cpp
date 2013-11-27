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
    @date 2013-07-15
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/LifeChronoManager.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

using namespace LifeV;

Real exactSolution ( const Real& /* t */, const Real& x, const Real& y, const Real& /* z */, const ID& /* i */ )
{
    return  sin ( x ) + y * y / 2.;
}

Real fRhs ( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */ , const ID& i )
{
    switch ( i )
    {
        case 0:
            return 0.;
            break;
        case 1:
            return 1.;
            break;
        default:
            ERROR_MSG ( "component not available!" );
    }

    return 0.;
}


typedef RegionMesh<LinearTriangle> mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

typedef boost::shared_ptr<Epetra_Comm> commPtr_Type;

typedef MatrixEpetraStructured<Real> matrix_Type;
typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
typedef VectorEpetra vectorStd_Type;
typedef VectorEpetraStructured vector_Type;

typedef FESpace<mesh_Type, MapEpetra> uSpaceStd_Type;
typedef boost::shared_ptr<uSpaceStd_Type> uSpaceStdPtr_Type;
typedef ETFESpace< mesh_Type, MapEpetra, 2, 2 > uSpace_Type;
typedef boost::shared_ptr<uSpace_Type> uSpacePtr_Type;
typedef ETFESpace< mesh_Type, MapEpetra, 2, 1 > pSpace_Type;
typedef boost::shared_ptr<pSpace_Type> pSpacePtr_Type;

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
        LifeChronoManager<> chronoMgr ( comm );

        LifeChrono initTime;
        chronoMgr.add ( "Initialization Time", &initTime );
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
        chronoMgr.add ( "Mesh reading/creation Time", &initTime );
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
        chronoMgr.add ( "Partition Time", &partTime );
        partTime.start();
        meshPtr_Type localMesh;
        {
            MeshPartitioner< mesh_Type >   meshPart;
            meshPart.doPartition ( fullMeshPtr, comm );
            localMesh = meshPart.meshPartition();
        }
        partTime.stop();

        Int localMeshNum[ 2 ];
        localMeshNum[ 0 ] = localMesh->numElements();
        localMeshNum[ 1 ] = localMesh->numPoints();
        Int maxMeshNum[ 2 ] = { 0, 0 };
        comm->MaxAll ( localMeshNum, maxMeshNum, 2 );

        if ( isLeader )
        {
            std::cout << "part mesh elements = " << maxMeshNum[ 0 ] << "\n"
                      << "part mesh points   = " << maxMeshNum[ 1 ] << std::endl;
        }

        LifeChrono partTimeR;
        chronoMgr.add ( "Partition Time (R)", &partTimeR );
        partTimeR.start();
        meshPtr_Type localMeshR;
        {
            MeshPartitioner< mesh_Type >   meshPartR;
            meshPartR.setPartitionOverlap ( 1 );
            meshPartR.doPartition ( fullMeshPtr, comm );
            localMeshR = meshPartR.meshPartition();
        }
        partTimeR.stop();

        localMeshNum[ 0 ] = localMeshR->numElements();
        localMeshNum[ 1 ] = localMeshR->numPoints();
        maxMeshNum[ 0 ] = 0;
        maxMeshNum[ 1 ] = 0;
        comm->MaxAll ( localMeshNum, maxMeshNum, 2 );

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
        ASSERT ( fullMeshPtr.use_count() == 0, "full mesh not properly freed." );
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
        LifeChrono feSpaceTime;
        chronoMgr.add ( "FESpace creation Time", &feSpaceTime );
        feSpaceTime.start();
        uSpaceStdPtr_Type uSpaceStd ( new uSpaceStd_Type ( localMesh, "P2", 2, comm ) );
        uSpacePtr_Type uSpace ( new uSpace_Type ( localMesh, &feTriaP2, & (uSpaceStd->fe().geoMap() ), comm) );
        pSpacePtr_Type pSpace ( new pSpace_Type ( localMesh, &feTriaP1, & (uSpaceStd->fe().geoMap() ), comm) );
        feSpaceTime.stop();

        LifeChrono feSpaceTimeR;
        chronoMgr.add ( "FESpace creation Time (R)", &feSpaceTimeR );
        feSpaceTimeR.start();
        uSpaceStdPtr_Type uSpaceStdR ( new uSpaceStd_Type ( localMeshR, "P2", 2, comm ) );
        uSpacePtr_Type uSpaceR ( new uSpace_Type ( localMeshR, &feTriaP2, & (uSpaceStdR->fe().geoMap() ), comm) );
        pSpacePtr_Type pSpaceR ( new pSpace_Type ( localMeshR, &feTriaP1, & (uSpaceStdR->fe().geoMap() ), comm) );
        feSpaceTimeR.stop();

        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }
        if ( verbose )
        {
            std::cout << " ---> Dofs: " << uSpace->dof().numTotalDof() << std::endl;
        }

        // Build the assembler and the matrices

        if ( verbose )
        {
            std::cout << " -- Defining and filling the matrix ... " << std::flush;
        }

        LifeChrono matTime;
        chronoMgr.add ( "Matrix creation Time", &matTime );
        matTime.start();
        matrixPtr_Type systemMatrix ( new matrix_Type ( uSpace->map() | pSpace->map() ) );
        matTime.stop();

        LifeChrono assemblyTime;
        chronoMgr.add ( "Assembly Time", &assemblyTime );
        assemblyTime.start();
        {
            using namespace ExpressionAssembly;

            integrate ( elements ( localMesh ),
                        quadRuleTria3pt,
                        uSpace,
                        uSpace,

                        dot ( grad (phi_i) , grad (phi_j) )

                      )
                    >> systemMatrix->block (0, 0);

            integrate ( elements ( localMesh ),
                        quadRuleTria3pt,
                        uSpace,
                        pSpace,

                        phi_j * div (phi_i)

                      )
                    >> systemMatrix->block (0, 1);

            integrate ( elements ( localMesh ),
                        quadRuleTria3pt,
                        pSpace,
                        uSpace,

                        phi_i * div (phi_j)

                      )
                    >> systemMatrix->block (1, 0);
        }
        systemMatrix->globalAssemble();
        assemblyTime.stop();

        LifeChrono matTimeR;
        chronoMgr.add ( "Matrix creation Time (R)", &matTimeR );
        matTimeR.start();
        matrixPtr_Type systemMatrixR ( new matrix_Type ( uSpaceR->map() | pSpaceR->map(), 50, true ) );
        matTimeR.stop();

        LifeChrono assemblyTimeR;
        chronoMgr.add ( "Assembly Time (R)", &assemblyTimeR );
        assemblyTimeR.start();
        {
            using namespace ExpressionAssembly;

            integrate ( elements ( localMeshR ),
                        quadRuleTria3pt,
                        uSpaceR,
                        uSpaceR,

                        dot ( grad (phi_i) , grad (phi_j) )

                      )
                    >> systemMatrixR->block (0, 0);

            integrate ( elements ( localMeshR ),
                        quadRuleTria3pt,
                        uSpaceR,
                        pSpaceR,

                        phi_j * div (phi_i)

                      )
                    >> systemMatrixR->block (0, 1);

            integrate ( elements ( localMeshR ),
                        quadRuleTria3pt,
                        pSpaceR,
                        uSpaceR,

                        phi_i * div (phi_j)

                      )
                    >> systemMatrixR->block (1, 0);
        }
        systemMatrixR->fillComplete();
        assemblyTime.stop();

        if ( verbose )
        {
            std::cout << " done! " << std::endl;
        }

        // check that the assembled matrices are the same

        LifeChrono checkMatTime;
        chronoMgr.add ( "Check (Matrix) Time", &checkMatTime );
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
        chronoMgr.add ( "Rhs build Time", &rhsTime );
        rhsTime.start();
        vector_Type rhs ( uSpaceStd->map() | pSpace->map(), Unique );

        vectorStd_Type fInterpolated ( uSpace->map(), Repeated );
        uSpaceStd->interpolate ( static_cast<uSpaceStd_Type::function_Type> (fRhs), fInterpolated, 0.0 );

        {
            using namespace ExpressionAssembly;

            integrate ( elements ( localMesh ),
                        quadRuleTria3pt,
                        uSpace,
                        dot (value (uSpace, fInterpolated), phi_i)
                      )
                    >> rhs.block (0);
        }
        rhs.globalAssemble();
        rhsTime.stop();

        LifeChrono rhsTimeR;
        chronoMgr.add ( "Rhs build Time (R)", &rhsTimeR );
        rhsTimeR.start();
        vector_Type rhsR ( uSpaceStdR->map() | pSpaceR->map(), Unique, Zero );

        vectorStd_Type fInterpolatedR ( uSpaceR->map(), Repeated );
        uSpaceStdR->interpolate ( static_cast<uSpaceStd_Type::function_Type> (fRhs), fInterpolatedR, 0.0 );

        {
            using namespace ExpressionAssembly;

            integrate ( elements ( localMeshR ),
                        quadRuleTria3pt,
                        uSpaceR,
                        dot (value (uSpaceR, fInterpolatedR), phi_i)
                      )
                    >> rhsR.block (0);
        }
        rhsR.globalAssemble();
        rhsTimeR.stop();

        LifeChrono checkVecTime;
        chronoMgr.add ( "Check (Vector) Time", &checkVecTime );
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
