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


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/LifeChrono.hpp>

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
    switch( i )
    {
    case 0:
        return 0.;
        break;
    case 1:
        return 1.;
        break;
    default:
        ERROR_MSG( "component not available!" );
    }

    return 0.;
}


typedef RegionMesh<LinearTriangle> mesh_Type;

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
        uSpaceStdPtr_Type uSpaceStd ( new uSpaceStd_Type ( localMesh, "P2", 2, comm ) );
        uSpaceStdPtr_Type uSpaceStdR ( new uSpaceStd_Type ( localMeshR, "P2", 2, comm ) );

        uSpacePtr_Type uSpace ( new uSpace_Type ( localMesh, &feTriaP2, & (uSpaceStd->fe().geoMap() ), comm) );
        uSpacePtr_Type uSpaceR ( new uSpace_Type ( localMeshR, &feTriaP2, & (uSpaceStdR->fe().geoMap() ), comm) );

        pSpacePtr_Type pSpace ( new pSpace_Type ( localMesh, &feTriaP1, & (uSpaceStd->fe().geoMap() ), comm) );
        pSpacePtr_Type pSpaceR ( new pSpace_Type ( localMeshR, &feTriaP1, & (uSpaceStdR->fe().geoMap() ), comm) );

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
            std::cout << " -- Defining the matrix ... " << std::flush;
        }
        matrixPtr_Type systemMatrix ( new matrix_Type ( uSpace->map() | pSpace->map() ) );
        matrixPtr_Type systemMatrixR ( new matrix_Type ( uSpaceR->map() | pSpaceR->map(), 50, true ) );
        if ( verbose )
        {
            std::cout << " done! " << std::endl;
        }

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

        vector_Type rhs ( uSpaceStd->map() | pSpace->map(), Unique );

        vectorStd_Type fInterpolated ( uSpace->map(), Repeated );
        uSpaceStd->interpolate ( fRhs, fInterpolated, 0.0 );

        {
            using namespace ExpressionAssembly;

            integrate ( elements ( localMesh ),
                        quadRuleTria3pt,
                        uSpace,
                        dot (value (uSpace, fInterpolated), phi_i)
                      )
                    >> rhs.block (0);
        }

        debugOut << "noGA\n" << rhs.epetraVector();

        rhs.globalAssemble();

        debugOut << "\n" << rhs.epetraVector();

        vector_Type rhsR ( uSpaceStdR->map() | pSpaceR->map(), Unique, Zero );

        vectorStd_Type fInterpolatedR ( uSpaceR->map(), Repeated );
        uSpaceStdR->interpolate ( fRhs, fInterpolatedR, 0.0 );

        {
            using namespace ExpressionAssembly;

            integrate ( elements ( localMeshR ),
                        quadRuleTria3pt,
                        uSpaceR,
                        dot (value (uSpaceR, fInterpolatedR), phi_i)
                      )
                    >> rhsR.block (0);
        }

        debugOut << "RnoGA\n" << rhsR.epetraVector();

        rhsR.globalAssemble();

        debugOut << "R\n" << rhsR.epetraVector();

        vector_Type vectorDiff ( rhs );
        vectorDiff -= rhsR;

        Real diffNormV = vectorDiff.normInf();

        if ( isLeader )
        {
            std::cout << "Norm of the difference between the 2 vectors = " << diffNormV << std::endl;
        }

        diff += diffNormV;

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
