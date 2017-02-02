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
    @file main.cpp
    @brief Test of the MatrixblockView framework

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 04-01-2012

    In this test we create a reference block matrix. It is used
    to create a second matrix using block manipulation.
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

using namespace LifeV;

int
main ( int argc, char** argv )
{
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
    std::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    int nprocs;
    MPI_Comm_size ( MPI_COMM_WORLD, &nprocs );
#else
    std::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm );
#endif

    Displayer displayer ( comm );

    displayer.leaderPrint ( " +-----------------------------------------------+\n" );
    displayer.leaderPrint ( " |    MatrixEpetraStructured Framework test      |\n" );
    displayer.leaderPrint ( " +-----------------------------------------------+\n\n" );
    displayer.leaderPrint ( " +-----------------------------------------------+\n" );
    displayer.leaderPrint ( " |           Author: Gwenol Grandperrin          |\n" );
    displayer.leaderPrint ( " |             Date: 2012-01-04                  |\n" );
    displayer.leaderPrint ( " +-----------------------------------------------+\n" );

    displayer.leaderPrint ( "\n[Initilization of MPI]\n" );
#ifdef HAVE_MPI
    displayer.leaderPrint ( "Using MPI (", nprocs, " proc.)\n" );
#else
    displayer.leaderPrint ( "Using serial version\n" );
#endif

    // +-----------------------------------------------+
    // |        Initialization of the test             |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "\n[Initilization of the test]\n" );

    int numFailed ( 0 );
    int problemSize ( 10 );

    // Creating the map
    displayer.leaderPrint ( "Creating the map... " );
    int numGlobalElement ( problemSize );
    int indexBase ( 0 );
    MapEpetra map ( numGlobalElement, indexBase, comm );
    displayer.leaderPrint ( "done\n" );

    /*
     * Creating the reference matrix, say A
     *   /1 1 2 2\
     *   |1 1 2 2|
     * A=|3 3 4 4|
     *   \3 3 4 4/
     */
    displayer.leaderPrint ( "Creating the source matrix... " );
    int numEntries ( 10 );
    std::vector<UInt> blockNumRows   ( 2, numEntries / 2 );
    std::vector<UInt> blockNumColumns ( 2, numEntries / 2 );
    std::shared_ptr< MatrixEpetraStructured<double> > A;
    A.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    A->setBlockStructure ( blockNumRows, blockNumColumns );
    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            for ( int j ( 0 ); j < problemSize / 2; ++j )
            {
                A->addToCoefficient ( i, j, 1 );
                A->addToCoefficient ( i, j + problemSize / 2, 2 );
            }
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            for ( int j ( 0 ); j < problemSize / 2; ++j )
            {
                A->addToCoefficient ( i + problemSize / 2, j, 3 );
                A->addToCoefficient ( i + problemSize / 2, j + problemSize / 2, 4 );
            }
        }
    }
    A->globalAssemble();
    displayer.leaderPrint ( "done\n" );

    // Getting the block structure of A
    MatrixEpetraStructuredView<double> A11, A12, A21, A22;
    A->blockView ( 0, 0, A11 );
    A->blockView ( 0, 1, A12 );
    A->blockView ( 1, 0, A21 );
    A->blockView ( 1, 1, A22 );

    // Creating a pointer for B
    std::shared_ptr< MatrixEpetraStructured<double> > B;

    // Defining the view
    MatrixEpetraStructuredView<double> B11, B12, B21, B22;

    // Creating a pointer for the test matrix
    std::shared_ptr< MatrixEpetraStructured<double> > testMatrix;

    // +-----------------------------------------------+
    // |        Core of the test                       |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "\n[Test of the framework]\n" );

    // +-----------------------------------------------+
    // |           TEST 1: Copy of blocks              |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "1) Testing blocks copies... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::copyBlock ( A11, B12 );
    MatrixEpetraStructuredUtility::copyBlock ( A12, B22 );
    MatrixEpetraStructuredUtility::copyBlock ( A22, B21 );
    MatrixEpetraStructuredUtility::copyBlock ( A21, B11 );

    B->globalAssemble();

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );
    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            for ( int j ( 0 ); j < problemSize / 2; ++j )
            {
                testMatrix->addToCoefficient ( i, j, -3 );
                testMatrix->addToCoefficient ( i, j + problemSize / 2, -1 );
            }
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            for ( int j ( 0 ); j < problemSize / 2; ++j )
            {
                testMatrix->addToCoefficient ( i + problemSize / 2, j, -4 );
                testMatrix->addToCoefficient ( i + problemSize / 2, j + problemSize / 2, -2 );
            }
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |           TEST 2: Identity blocks             |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "2) Creating identity blocks... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::createIdentityBlock ( B11 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B12 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B21 );
    MatrixEpetraStructuredUtility::createIdentityBlock ( B22 );

    B->globalAssemble();

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );

    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            testMatrix->addToCoefficient ( i, i, -1 );
            testMatrix->addToCoefficient ( i, i + problemSize / 2, -1 );
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            testMatrix->addToCoefficient ( i + problemSize / 2, i, -1 );
            testMatrix->addToCoefficient ( i + problemSize / 2, i + problemSize / 2, -1 );
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |           TEST 3: Diagonal blocks             |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "3) Creating diagonal blocks... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::createDiagBlock ( A11, B11 );
    MatrixEpetraStructuredUtility::createDiagBlock ( A12, B12 );
    MatrixEpetraStructuredUtility::createDiagBlock ( A21, B21 );
    MatrixEpetraStructuredUtility::createDiagBlock ( A22, B22 );

    B->globalAssemble();

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );

    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            testMatrix->addToCoefficient ( i, i, -1 );
            testMatrix->addToCoefficient ( i, i + problemSize / 2, -2 );
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            testMatrix->addToCoefficient ( i + problemSize / 2, i, -3 );
            testMatrix->addToCoefficient ( i + problemSize / 2, i + problemSize / 2, -4 );
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |           TEST 4: Inverse diagonal blocks     |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "4) Creating inverse diagonal blocks... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::createInvDiagBlock ( A11, B11 );
    MatrixEpetraStructuredUtility::createInvDiagBlock ( A12, B12 );
    MatrixEpetraStructuredUtility::createInvDiagBlock ( A21, B21 );
    MatrixEpetraStructuredUtility::createInvDiagBlock ( A22, B22 );

    B->globalAssemble();

    Real minusOneThird ( -1. / 3. );

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );
    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            testMatrix->addToCoefficient ( i, i, -1 );
            testMatrix->addToCoefficient ( i, i + problemSize / 2, -0.5 );
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            testMatrix->addToCoefficient ( i + problemSize / 2, i, minusOneThird );
            testMatrix->addToCoefficient ( i + problemSize / 2, i + problemSize / 2, -0.25 );
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |           TEST 5: Upper triangular blocks     |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "5) Creating upper triangular blocks... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::createUpperTriangularBlock ( A11, B11 );
    MatrixEpetraStructuredUtility::createUpperTriangularBlock ( A12, B12 );
    MatrixEpetraStructuredUtility::createUpperTriangularBlock ( A21, B21 );
    MatrixEpetraStructuredUtility::createUpperTriangularBlock ( A22, B22 );

    B->globalAssemble();

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );
    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            for ( int j ( i ); j < problemSize / 2; ++j )
            {
                testMatrix->addToCoefficient ( i, j, -1 );
                testMatrix->addToCoefficient ( i, j + problemSize / 2, -2 );
            }
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            for ( int j ( i ); j < problemSize / 2; ++j )
            {
                testMatrix->addToCoefficient ( i + problemSize / 2, j, -3 );
                testMatrix->addToCoefficient ( i + problemSize / 2, j + problemSize / 2, -4 );
            }
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |           TEST 6: Lower triangular blocks     |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "6) Lower triangular blocks... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::createLowerTriangularBlock ( A11, B11 );
    MatrixEpetraStructuredUtility::createLowerTriangularBlock ( A12, B12 );
    MatrixEpetraStructuredUtility::createLowerTriangularBlock ( A21, B21 );
    MatrixEpetraStructuredUtility::createLowerTriangularBlock ( A22, B22 );

    B->globalAssemble();

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );
    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            for ( int j ( 0 ); j <= i; ++j )
            {
                testMatrix->addToCoefficient ( i, j, -1 );
                testMatrix->addToCoefficient ( i, j + problemSize / 2, -2 );
            }
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            for ( int j ( 0 ); j <= i; ++j )
            {
                testMatrix->addToCoefficient ( i + problemSize / 2, j, -3 );
                testMatrix->addToCoefficient ( i + problemSize / 2, j + problemSize / 2, -4 );
            }
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |           TEST 7: Lumped blocks               |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "7) Creating lumped blocks... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::createLumpedBlock ( A11, B11 );
    MatrixEpetraStructuredUtility::createLumpedBlock ( A12, B12 );
    MatrixEpetraStructuredUtility::createLumpedBlock ( A21, B21 );
    MatrixEpetraStructuredUtility::createLumpedBlock ( A22, B22 );

    B->globalAssemble();

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );
    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            testMatrix->addToCoefficient ( i, i, -problemSize / 2 );
            testMatrix->addToCoefficient ( i, i + problemSize / 2, -problemSize );
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            testMatrix->addToCoefficient ( i + problemSize / 2, i, -problemSize * 1.5 );
            testMatrix->addToCoefficient ( i + problemSize / 2, i + problemSize / 2, -problemSize * 2 );
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |           TEST 8: Inverse lumped blocks       |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "8) Creating inverse lumped blocks... " );
    B.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    B->setBlockStructure ( blockNumRows, blockNumColumns );
    B->blockView ( 0, 0, B11 );
    B->blockView ( 0, 1, B12 );
    B->blockView ( 1, 0, B21 );
    B->blockView ( 1, 1, B22 );

    MatrixEpetraStructuredUtility::createInvLumpedBlock ( A11, B11 );
    MatrixEpetraStructuredUtility::createInvLumpedBlock ( A12, B12 );
    MatrixEpetraStructuredUtility::createInvLumpedBlock ( A21, B21 );
    MatrixEpetraStructuredUtility::createInvLumpedBlock ( A22, B22 );

    B->globalAssemble();

    testMatrix.reset ( new MatrixEpetraStructured<double> ( map, numEntries )  );
    testMatrix->setBlockStructure ( blockNumRows, blockNumColumns );
    for ( int i ( 0 ); i < problemSize / 2; ++i )
    {
        if ( A->matrixPtr()->MyGRID ( i ) )
        {
            testMatrix->addToCoefficient ( i, i, -2. / problemSize );
            testMatrix->addToCoefficient ( i, i + problemSize / 2, -1. / problemSize );
        }

        if ( A->matrixPtr()->MyGRID ( i + problemSize / 2 ) )
        {
            testMatrix->addToCoefficient ( i + problemSize / 2, i, -1. / (problemSize * 1.5) );
            testMatrix->addToCoefficient ( i + problemSize / 2, i + problemSize / 2, -1. / (problemSize * 2) );
        }
    }
    testMatrix->globalAssemble();
    *B += *testMatrix;
    if ( B->normInf() == 0)
    {
        displayer.leaderPrint ( "PASSED\n" );
    }
    else
    {
        numFailed += 1;
        displayer.leaderPrint ( "FAILED\n" );
    }


    // +-----------------------------------------------+
    // |        Results of the test                    |
    // +-----------------------------------------------+
    displayer.leaderPrint ( "\n[Test of the blockview framework]\n" );
    // The test must verify if tolerance is satisfied!
    if ( numFailed > 0 )
    {
        displayer.leaderPrint ( numFailed, "/8 tests failed\n" );
        return EXIT_FAILURE;
    }
    else
    {
        displayer.leaderPrint ( "The framework is working as expected\n" );
    }
    // ----- End of test calls -----

#ifdef HAVE_MPI
    displayer.leaderPrint ( "\n[MPI Finalization]\n" );
    MPI_Finalize();
#endif

    // If everything runs correctly we return EXIT_SUCCESS flag
    return EXIT_SUCCESS;
}
