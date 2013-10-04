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

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 30-03-2011
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/algorithm/PreconditionerLinearSolver.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/function/Laplacian.hpp>

#define TEST_TOLERANCE 1e-13

using namespace LifeV;

namespace
{
typedef RegionMesh<LinearTetra>           mesh_Type;
typedef boost::shared_ptr<mesh_Type>      meshPtr_Type;
typedef MatrixEpetra<Real>                matrix_Type;
typedef VectorEpetra                      vector_Type;
typedef boost::shared_ptr<VectorEpetra>   vectorPtr_Type;
typedef FESpace< mesh_Type, MapEpetra >   fespace_Type;
typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;

typedef LifeV::Preconditioner             basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
typedef LifeV::PreconditionerLinearSolver prec_Type;
typedef boost::shared_ptr<prec_Type>      precPtr_Type;
typedef boost::function < Real ( Real const&,
                                 Real const&,
                                 Real const&,
                                 Real const&,
                                 UInt const& ) > function_Type;
}

void printErrors ( const vector_Type& solution, fespacePtr_Type uFESpace, bool verbose )
{
    vector_Type velocity ( solution, Repeated );
    Real uRelativeError, uL2Error;
    uL2Error = uFESpace->l2Error (Laplacian::uexact, velocity, 0, &uRelativeError );
    if ( verbose )
    {
        std::cout << "Velocity" << std::endl;
    }
    if ( verbose )
    {
        std::cout << "  L2 error      : " << uL2Error << std::endl;
    }
    if ( verbose )
    {
        std::cout << "  Relative error: " << uRelativeError << std::endl;
    }
}


int
main ( int argc, char** argv )
{
    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif
    int exit_status;
    bool verbose ( true );

    {

#ifdef HAVE_MPI
        boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
        boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif

        verbose = (Comm->MyPID() == 0);

        if ( verbose )
        {
            std::cout
                    << " +-----------------------------------------------+" << std::endl
                    << " |       LinearSolverPreconditioner Test         |" << std::endl
                    << " +-----------------------------------------------+" << std::endl
                    << std::endl
                    << " +-----------------------------------------------+" << std::endl
                    << " |           Author: Gwenol Grandperrin          |" << std::endl
                    << " |             Date: 2013-01-16                  |" << std::endl
                    << " +-----------------------------------------------+" << std::endl
                    << std::endl;

            std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
            std::cout << "Using MPI (" << Comm->NumProc() << " proc.)" << std::endl;
#else
            std::cout << "Using serial version" << std::endl;
#endif
        }

        // +-----------------------------------------------+
        // |               Loading the data                |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Loading the data]" << std::endl;
        }
        LifeChrono globalChrono;

        globalChrono.start();

        // ********** GetPot **********
        GetPot command_line ( argc, argv );
        const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
        GetPot dataFile ( dataFileName );
        // ****************************

        // Space discretization
        const UInt numMeshElem    = dataFile ( "mesh/num_elements", 10);
        const std::string uOrder  = dataFile ( "finite_element/velocity", "P1" );

        // Solution initialization
        Laplacian::setModes ( 1, 1, 1 );

        // +-----------------------------------------------+
        // |               Loading the mesh                |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Loading the mesh]" << std::endl;
        }

        meshPtr_Type fullMeshPtr ( new mesh_Type );
        const UInt& geoDim = static_cast<UInt> ( mesh_Type::S_geoDimensions );


        // Building the mesh from the source

        regularMesh3D ( *fullMeshPtr,
                        1,
                        numMeshElem, numMeshElem, numMeshElem,
                        false,
                        1.0, 1.0, 1.0,
                        0.0, 0.0, 0.0 );

        if ( verbose ) std::cout << "Mesh source: regular mesh("
                                     << numMeshElem << "x" << numMeshElem << "x" << numMeshElem << ")" << std::endl;

        if ( verbose )
        {
            std::cout << "Mesh size  : " << MeshUtility::MeshStatistics::computeSize ( *fullMeshPtr ).maxH << std::endl;
        }
        if ( verbose )
        {
            std::cout << "Partitioning the mesh ... " << std::endl;
        }
        MeshPartitioner< mesh_Type >   meshPart ( fullMeshPtr, Comm );
        fullMeshPtr.reset(); //Freeing the global mesh to save memory

        // +-----------------------------------------------+
        // |            Creating the FE spaces             |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
        }
        if ( verbose )
        {
            std::cout << "FE for the velocity: " << uOrder << std::endl;
        }

        if ( verbose )
        {
            std::cout << "Building the velocity FE space ... " << std::flush;
        }
        fespacePtr_Type uFESpace ( new fespace_Type ( meshPart, uOrder, geoDim, Comm ) );
        if ( verbose )
        {
            std::cout << "ok." << std::endl;
        }

        // Pressure offset in the vector
        UInt numDofs = geoDim * uFESpace->dof().numTotalDof();

        if ( verbose )
        {
            std::cout << "Total Velocity Dof: " << numDofs << std::endl;
        }

        // +-----------------------------------------------+
        // |             Boundary conditions               |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Boundary conditions]" << std::endl;
        }
        BCHandler bcHandler;
        BCFunctionBase uExact ( Laplacian::uexact );
        BCFunctionBase fRHS ( Laplacian::f );

        if ( verbose )
        {
            std::cout << "Setting Dirichlet BC... " << std::flush;
        }
        for ( UInt iDirichlet ( 1 ); iDirichlet <= 26; ++iDirichlet )
        {
            bcHandler.addBC ( "Wall", iDirichlet, Essential, Full, uExact, geoDim );
        }
        if ( verbose )
        {
            std::cout << "ok." << std::endl;
        }

        // Update the BCHandler (internal data related to FE)
        bcHandler.bcUpdate ( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

        // +-----------------------------------------------+
        // |              Matrices Assembly                |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Matrices Assembly]" << std::endl;
        }

        if ( verbose )
        {
            std::cout << "Setting up assembler... " << std::flush;
        }
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        adrAssembler.setup ( uFESpace, uFESpace );
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }

        if ( verbose )
        {
            std::cout << "Defining the matrices... " << std::flush;
        }
        boost::shared_ptr<matrix_Type> systemMatrix ( new matrix_Type ( uFESpace->map() ) );
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }

        if ( verbose )
        {
            std::cout << "Adding the viscous stress... " << std::flush;
        }
        adrAssembler.addDiffusion ( systemMatrix, 1.0 );
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }

        // +-----------------------------------------------+
        // |            Solver initialization              |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Solvers initialization]" << std::endl;
        }
        prec_Type* precRawPtr;
        basePrecPtr_Type precPtr;
        precRawPtr = new prec_Type;
        precRawPtr->setDataFromGetPot ( dataFile, "prec/LinearSolver" );
        precPtr.reset ( precRawPtr );

        if ( verbose )
        {
            std::cout << "Setting up LinearSolver (Belos)... " << std::flush;
        }
        Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
        belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

        LinearSolver linearSolver;
        linearSolver.setCommunicator ( Comm );
        linearSolver.setParameters ( *belosList );
        linearSolver.setPreconditioner ( precPtr );
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }
        linearSolver.showMe();

        // +-----------------------------------------------+
        // |                   Simulation                  |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Initialization of the simulation]" << std::endl;
        }
        if ( verbose )
        {
            std::cout << "Creation of vectors... " << std::flush;
        }

        boost::shared_ptr<vector_Type> rhs;
        rhs.reset ( new vector_Type ( uFESpace->map(), Repeated ) );

        vector_Type fInterpolated ( uFESpace->map(), Repeated );
        adrAssembler.addMassRhs ( *rhs, fRHS, 0.0 );
        rhs->globalAssemble();
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }

        fInterpolated = 0.0;
        uFESpace->interpolate ( uExact, fInterpolated, 0.0 );


        if ( verbose )
        {
            std::cout << "Applying BC... " << std::flush;
        }
        systemMatrix->globalAssemble();
        boost::shared_ptr<vector_Type> rhsBC;
        rhsBC.reset ( new vector_Type ( *rhs, Unique ) );
        bcManage ( *systemMatrix, *rhsBC, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }

        if ( verbose )
        {
            std::cout << std::endl << "Solving the system with LinearSolver (Belos)... " << std::endl;
        }
        boost::shared_ptr<vector_Type> solution;
        solution.reset ( new vector_Type ( uFESpace->map(), Unique ) );
        linearSolver.setOperator ( systemMatrix );
        linearSolver.setRightHandSide ( rhsBC );
        linearSolver.solve ( solution );

        // +-----------------------------------------------+
        // |             Computing the error               |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Errors computation]" << std::endl;
        }

        vector_Type solution2Err ( *solution );
        solution2Err *= 0.0;
        uFESpace->interpolate ( static_cast<function_Type> ( Laplacian::uexact ), solution2Err, 0.0 );
        solution2Err -= *solution;
        solution2Err.abs();

        if ( verbose )
        {
            std::cout << "Linear solver Belos" << std::endl;
        }
        printErrors ( *solution, uFESpace, verbose );

        // +-----------------------------------------------+
        // |            Ending the simulation              |
        // +-----------------------------------------------+
        globalChrono.stop();
        if ( verbose )
        {
            std::cout << std::endl << "Total simulation time: " << globalChrono.diff() << " s." << std::endl;
        }
        if ( verbose )
        {
            std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;
        }

        exit_status = (linearSolver.numIterations() != 3 || precRawPtr->solverPtr()->numIterations() != 5);

    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( exit_status )
    {
        if ( verbose )
        {
            std::cout << "The number of iterations has changed." << std::endl;
        }
        if ( verbose )
        {
            std::cout << "Test status: FAILED" << std::endl;
        }
        return ( EXIT_FAILURE );
    }


    if ( verbose )
    {
        std::cout << "Test status: SUCCESS" << std::endl;
    }

    return ( EXIT_SUCCESS );
}
