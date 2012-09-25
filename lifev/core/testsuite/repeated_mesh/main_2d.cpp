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

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-10-2010
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

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}

//#define TEST_MASS
//#define TEST_ADVECTION
#define TEST_RHS


#ifdef TEST_MASS
Real epsilon(1);

Real exactSolution( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
    Real seps(sqrt(epsilon));
    return  exp(seps*x)/(exp(seps)-exp(-seps));
}
#endif

#ifdef TEST_ADVECTION
Real epsilon(1);

Real exactSolution( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */, const ID& /* i */ )
{
    return  (exp(x/epsilon) - 1 )/( exp(1/epsilon) - 1);
}

Real betaFct( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& i )
{
    return (i == 0);
}
#endif

#ifdef TEST_RHS
Real epsilon(1);

Real exactSolution( const Real& /* t */, const Real& x, const Real& y, const Real& /* z */, const ID& /* i */ )
{
    return  sin( x ) + y * y / 2.;
}


Real fRhs( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
    return  sin( x ) - 1.;
}
#endif


typedef RegionMesh<LinearTriangle> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

int
main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
#endif

    // introducing a local scope in order to properly destroy all objects
    // before calling MPI_Finalize()
    {

#ifdef HAVE_MPI
        boost::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
        boost::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm);
#endif

        GetPot dataFile( "data_2d" );
        const bool isLeader( comm->MyPID() == 0 );
        const bool verbose( dataFile( "miscellaneous/verbose", 0 ) && isLeader );

#ifdef HAVE_LIFEV_DEBUG
        std::ofstream debugOut(
                    ( "rm." +
                      ( comm->NumProc() > 1 ? boost::lexical_cast<std::string>( comm->MyPID() ) : "s" ) +
                      ".out" ).c_str() );
#else
        std::ofstream debugOut( "/dev/null" );
#endif

        // Build and partition the mesh

        if ( verbose ) std::cout << " -- Reading the mesh ... " << std::flush;
        boost::shared_ptr< mesh_Type > fullMeshPtr(new mesh_Type());
        if( dataFile( "mesh/mesh_type", "structured" ) == "structured" )
        {
            regularMesh2D( *fullMeshPtr, 0,
                           dataFile( "mesh/nx", 20 ), dataFile( "mesh/ny", 20 ),
                           dataFile( "mesh/verbose", false ),
                           dataFile( "mesh/lx", 1. ), dataFile( "mesh/ly", 1. ) );
        }
        else
        {
            MeshData meshData(dataFile, "mesh");
            readMesh(*fullMeshPtr,meshData);
        }
        if ( verbose ) std::cout << " done ! " << std::endl;
        if( isLeader )
        {
            std::cout << "mesh elements = " << fullMeshPtr->numElements() << "\n"
                      << "mesh points   = " << fullMeshPtr->numPoints() << std::endl;
        }

        if ( verbose ) std::cout << " -- Partitioning the mesh ... " << std::flush;

        LifeChrono partTime;
        partTime.start();
        boost::shared_ptr< mesh_Type > localMesh;
        {
            MeshPartitioner< mesh_Type >   meshPart;
            meshPart.doPartition( fullMeshPtr, comm );
            localMesh = meshPart.meshPartition();
        }
        partTime.stop();
        if( isLeader ) std::cout << "partitioning time  = " << partTime.diff() << std::endl;
        if( isLeader )
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
            meshPartR.setPartitionOverlap( 1 );
            meshPartR.doPartition( fullMeshPtr, comm );
            localMeshR = meshPartR.meshPartition();
        }
        partTimeR.stop();
        if( isLeader ) std::cout << "partitioningR time = " << partTimeR.diff() << std::endl;

        // debugOut << "============================" << std::endl;
        // localMeshR->mesh_Type::showMe( true, debugOut );
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Freeing the global mesh ... " << std::flush;
        fullMeshPtr.reset();
        if ( verbose ) std::cout << " done ! " << std::endl;

        // Build the FESpaces

        if ( verbose ) std::cout << " -- Building FESpaces ... " << std::flush;
        std::string uOrder("P1");
        std::string bOrder("P1");
        feSpacePtr_Type uFESpace( new feSpace_Type( localMesh, uOrder, 1, comm ) );
        feSpacePtr_Type uFESpaceR( new feSpace_Type( localMeshR, uOrder, 1, comm ) );
        feSpacePtr_Type betaFESpace( new feSpace_Type( localMesh, bOrder, 2, comm ) );
        feSpacePtr_Type betaFESpaceR( new feSpace_Type( localMeshR, bOrder, 2, comm ) );
        if ( verbose ) std::cout << " done ! " << std::endl;
        if ( verbose ) std::cout << " ---> Dofs: " << uFESpaceR->dof().numTotalDof() << std::endl;

        // Build the assembler and the matrices

        if ( verbose ) std::cout << " -- Building assembler ... " << std::flush;
        ADRAssembler<mesh_Type,matrix_Type,vector_Type> adrAssembler;
        ADRAssembler<mesh_Type,matrix_Type,vector_Type> adrAssemblerR;
        if ( verbose ) std::cout << " done! " << std::endl;

        if ( verbose ) std::cout << " -- Setting up assembler ... " << std::flush;
        adrAssembler.setup( uFESpace, betaFESpace);
        adrAssemblerR.setup( uFESpaceR, betaFESpaceR);
        if ( verbose ) std::cout << " done! " << std::endl;

        if ( verbose ) std::cout << " -- Defining the matrix ... " << std::flush;
        boost::shared_ptr<matrix_Type> systemMatrix(new matrix_Type( uFESpace->map() ));
        boost::shared_ptr<matrix_Type> systemMatrixR(new matrix_Type( uFESpaceR->map(), 50, true ));
        if ( verbose ) std::cout << " done! " << std::endl;

        // Perform the assembly of the matrix

        if ( verbose ) std::cout << " -- Adding the diffusion ... " << std::flush;
        if ( verbose ) std::cout << " done! " << std::endl;

        LifeChrono assemblyTime;
        assemblyTime.start();
        adrAssembler.addDiffusion(systemMatrix,epsilon);

//#ifdef TEST_ADVECTION
//        if ( verbose ) std::cout << " -- Adding the advection ... " << std::flush;
//        vector_type beta(betaFESpace->map(),Repeated);
//        betaFESpace->interpolate(betaFct,beta,0.0);
//        adrAssembler.addAdvection(systemMatrix,beta);
//        if ( verbose ) std::cout << " done! " << std::endl;
//#endif
//#ifdef TEST_MASS
//        if ( verbose ) std::cout << " -- Adding the mass ... " << std::flush;
//        adrAssembler.addMass(systemMatrix,1.0);
//        if ( verbose ) std::cout << " done! " << std::endl;
//#endif

        if ( verbose ) std::cout << " -- Closing the matrix ... " << std::flush;
        systemMatrix->globalAssemble();
        if ( verbose ) std::cout << " done ! " << std::endl;
        assemblyTime.stop();
        if( isLeader ) std::cout << "assembly time  = " << assemblyTime.diff() << std::endl;
        LifeChrono assemblyTimeR;
        assemblyTimeR.start();
        adrAssemblerR.addDiffusion(systemMatrixR,epsilon);
        systemMatrixR->matrixPtr()->FillComplete();
        assemblyTimeR.stop();
        if( isLeader ) std::cout << "assemblyR time = " << assemblyTimeR.diff() << std::endl;

        if ( verbose ) std::cout << " Time needed : " << adrAssembler.diffusionAssemblyChrono().diffCumul() << std::endl;
        if ( verbose ) std::cout << " Time needed : " << adrAssemblerR.diffusionAssemblyChrono().diffCumul() << std::endl;

        // SPY
//        systemMatrix->spy("matrixNoBC");
//        systemMatrixR->spy("matrixNoBCR");

#ifdef TEST_RHS
        Real matrixNorm(systemMatrix->norm1());
        if ( verbose ) std::cout << " ---> Norm 1 : " << matrixNorm << std::endl;
        if ( std::fabs(matrixNorm - 8 ) > 1e-3)
        {
            std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }
#endif

        // Definition and assembly of the RHS

        if ( verbose ) std::cout << " -- Building the RHS ... " << std::flush;
        //vector_type rhs(uFESpace->map(),Unique);
        vector_Type rhs(uFESpace->map(),Repeated);
        vector_Type rhsR(uFESpaceR->map(),Repeated);

#ifdef TEST_RHS
        adrAssembler.addMassRhs( rhs, fRhs, 0. );
        rhs.globalAssemble();
        adrAssemblerR.addMassRhs( rhsR, fRhs, 0. );
        rhsR.globalAssemble( Zero );
#endif

        if ( verbose ) std::cout << " done ! " << std::endl;

        // Definition and application of the BCs

        if ( verbose ) std::cout << " -- Building the BCHandler ... " << std::flush;
        BCHandler bchandler;
        BCHandler bchandlerR;
        BCFunctionBase BCu( exactSolution );
        for ( UInt side = 1; side <= 4; side++ )
        {
            bchandler.addBC ( "Dirichlet", side, Essential, Full, BCu, 1 );
            bchandlerR.addBC( "Dirichlet", side, Essential, Full, BCu, 1 );
        }
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Updating the BCs ... " << std::flush;
        bchandler.bcUpdate(  *uFESpace->mesh(),  uFESpace->feBd(),  uFESpace->dof() );
        bchandlerR.bcUpdate( *uFESpaceR->mesh(), uFESpaceR->feBd(), uFESpaceR->dof() );
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Applying the BCs ... " << std::flush;
        vector_Type rhsBC(rhs,Unique);
        vector_Type rhsBCR(rhsR,Unique);
        bcManage(*systemMatrix,rhsBC,*uFESpace->mesh(),uFESpace->dof(),bchandler,uFESpace->feBd(),1.0,0.0);
        bcManage(*systemMatrixR,rhsBCR,*uFESpaceR->mesh(),uFESpaceR->dof(),bchandlerR,uFESpaceR->feBd(),1.0,0.0);
        rhs = rhsBC;
        rhsR = rhsBCR;
        if ( verbose ) std::cout << " done ! " << std::endl;

        // SPY
//        systemMatrix->spy("matrix");
//        systemMatrixR->spy("matrixR");
//        rhs.spy("rhs");
//        rhsR.spy("rhsR");

        // Definition of the solver

        if ( verbose ) std::cout << " -- Building the solver ... " << std::flush;
        SolverAztecOO linearSolver;
        SolverAztecOO linearSolverR;
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Setting up the solver ... " << std::flush;
        linearSolver.setDataFromGetPot(    dataFile, "solver" );
        linearSolver.setupPreconditioner(  dataFile, "prec" );
        linearSolverR.setDataFromGetPot(   dataFile, "solver" );
        linearSolverR.setupPreconditioner( dataFile, "prec" );
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Setting matrix in the solver ... " << std::flush;
        linearSolver.setMatrix(  *systemMatrix  );
        linearSolverR.setMatrix( *systemMatrixR );
        if ( verbose ) std::cout << " done ! " << std::endl;

        linearSolver.setCommunicator( comm );
        linearSolverR.setCommunicator( comm );

        // Definition of the solution

        if ( verbose ) std::cout << " -- Defining the solution ... " << std::flush;
        vector_Type solution( uFESpace->map(), Unique );
        vector_Type solutionR( uFESpaceR->map(),Unique );
        if ( verbose ) std::cout << " done ! " << std::endl;

        // Solve the solution

        if ( verbose ) std::cout << " -- Solving the system ... " << std::flush;
        linearSolver.solveSystem( rhsBC, solution, systemMatrix );
        linearSolverR.solveSystem( rhsBCR, solutionR, systemMatrixR );
        if ( verbose ) std::cout << " done ! " << std::endl;

        // SPY
//        solution.spy("solution");
//        solutionR.spy("solutionR");

        // Error computation

        if ( verbose ) std::cout << " -- Computing the error ... " << std::flush;
        vector_Type error( solution );
        error -= solutionR;
        if ( verbose ) std::cout << " done ! " << std::endl;

        // Exporter definition and use

        if ( verbose ) std::cout << " -- Defining the exporter ... " << std::flush;
#ifdef HAVE_HDF5
        ExporterHDF5<mesh_Type> exporter ( dataFile, localMeshR, "solution", comm->MyPID()) ;
#else
        // TODO: ExporterEnsight with repeated mesh not tested!
        ExporterEnsight<mesh_type> exporter ( dataFile, localMesh, "solution", Comm->MyPID()) ;
#endif
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Defining the exported quantities ... " << std::flush;
        boost::shared_ptr<vector_Type> solutionPtr( new vector_Type( solution, Repeated ) );
        boost::shared_ptr<vector_Type> errorPtr( new vector_Type( error, Repeated ) );
        boost::shared_ptr<vector_Type> solutionPtrR( new vector_Type( solutionR, Repeated ) );
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Updating the exporter ... " << std::flush;
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "solution", uFESpace, solutionPtr, UInt(0) );
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "solutionR", uFESpaceR, solutionPtrR, UInt(0) );
        exporter.addVariable( ExporterData<mesh_Type>::ScalarField, "error", uFESpace, errorPtr, UInt(0) );
        if ( verbose ) std::cout << " done ! " << std::endl;

        if ( verbose ) std::cout << " -- Exporting ... " << std::flush;
        exporter.exportPID( localMesh, comm, false );
        exporter.exportPID( localMeshR, comm, true );
        exporter.postProcess( 0 );
        if ( verbose ) std::cout << " done ! " << std::endl;


        if( isLeader ) std::cout << "End Result: TEST PASSED" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return( EXIT_SUCCESS );
}
