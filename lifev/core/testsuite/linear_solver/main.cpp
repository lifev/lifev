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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include<math.h>
#include <string>

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
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/function/Laplacian.hpp>

using namespace LifeV;

namespace
{
typedef RegionMesh<LinearTetra>           mesh_type;
typedef MatrixEpetra<Real>                matrix_type;
typedef VectorEpetra                      vector_type;
typedef boost::shared_ptr<VectorEpetra>   vectorPtr_type;
typedef FESpace< mesh_type, MapEpetra >   fespace_type;
typedef boost::shared_ptr< fespace_type > fespacePtr_type;

typedef LifeV::Preconditioner             basePrec_type;
typedef boost::shared_ptr<basePrec_type>  basePrecPtr_type;
typedef LifeV::PreconditionerIfpack       prec_type;
typedef boost::shared_ptr<prec_type>      precPtr_type;
}

void printErrors( const vector_type& solution, fespacePtr_type uFESpace, bool verbose )
{
    vector_type velocity( solution, Repeated );
    Real uRelativeError, uL2Error;
    uL2Error = uFESpace->l2Error (Laplacian::uexact, velocity, 0, &uRelativeError );
    if( verbose ) std::cout << "Velocity" << std::endl;
    if( verbose ) std::cout << "  L2 error      : " << uL2Error << std::endl;
    if( verbose ) std::cout << "  Relative error: " << uRelativeError << std::endl;
}


int
main( int argc, char** argv )
{
    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int nproc;
    MPI_Comm_size( MPI_COMM_WORLD, &nproc );
#else
    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );
#endif

    const bool verbose( Comm->MyPID() == 0 );
    if( verbose ){
        std::cout
            << " +-----------------------------------------------+" << std::endl
            << " |             LinearSolver Test                 |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl
            << " +-----------------------------------------------+" << std::endl
            << " |           Author: Gwenol Grandperrin          |" << std::endl
            << " |             Date: 2010-08-03                  |" << std::endl
            << " +-----------------------------------------------+" << std::endl
            << std::endl;

        std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
        std::cout << "Using MPI (" << nproc << " proc.)" << std::endl;
#else
        std::cout << "Using serial version" << std::endl;
#endif
    }

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    if( verbose ) std::cout << std::endl << "[Loading the data]" << std::endl;
    LifeChrono globalChrono;

    globalChrono.start();

    // ********** GetPot **********
    GetPot command_line( argc, argv );
    const std::string dataFileName = command_line.follow( "data", 2, "-f", "--file" );
    GetPot dataFile( dataFileName );
    // ****************************

    // Space discretization
    const UInt numMeshElem    = dataFile( "mesh/num_elements", 10);
    const std::string uOrder  = dataFile( "finite_element/velocity", "P1" );

    // Solution initialization
    Laplacian::setModes( 1, 1, 1 );

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if( verbose ) std::cout << std::endl << "[Loading the mesh]" << std::endl;

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr( new RegionMesh<LinearTetra> );

    // Building the mesh from the source

    regularMesh3D( *fullMeshPtr,
                   1,
                   numMeshElem, numMeshElem, numMeshElem,
                   false,
                   1.0, 1.0, 1.0,
                   0.0, 0.0, 0.0 );

    if( verbose ) std::cout << "Mesh source: regular mesh("
                           << numMeshElem << "x" << numMeshElem << "x" << numMeshElem << ")" << std::endl;

    if( verbose ) std::cout << "Mesh size  : " << MeshUtility::MeshStatistics::computeSize( *fullMeshPtr ).maxH << std::endl;
    if( verbose ) std::cout << "Partitioning the mesh ... " << std::endl;
    MeshPartitioner< RegionMesh<LinearTetra> >   meshPart( fullMeshPtr, Comm );
    fullMeshPtr.reset(); //Freeing the global mesh to save memory

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if( verbose ) std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    if( verbose ) std::cout << "FE for the velocity: " << uOrder << std::endl;

    if( verbose ) std::cout << "Building the velocity FE space ... " << std::flush;
    fespacePtr_type uFESpace( new FESpace< mesh_type, MapEpetra >( meshPart, uOrder, nDimensions, Comm ) );
    if( verbose ) std::cout << "ok." << std::endl;

    // Pressure offset in the vector
    UInt numDofs = nDimensions * uFESpace->dof().numTotalDof();

    if( verbose ) std::cout << "Total Velocity Dof: " << numDofs << std::endl;

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if( verbose ) std::cout << std::endl << "[Boundary conditions]" << std::endl;
    BCHandler bcHandler;
    BCFunctionBase uExact( Laplacian::uexact );
    BCFunctionBase fRHS( Laplacian::f );

    if( verbose ) std::cout << "Setting Dirichlet BC... " << std::flush;
    for( UInt iDirichlet( 1 ); iDirichlet<=26; ++iDirichlet )
    {
        bcHandler.addBC( "Wall", iDirichlet, Essential, Full, uExact, 3 );
    }
    if( verbose ) std::cout << "ok." << std::endl;

    // Update the BCHandler (internal data related to FE)
    bcHandler.bcUpdate( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

    // +-----------------------------------------------+
    // |              Matrices Assembly                |
    // +-----------------------------------------------+
    if( verbose ) std::cout << std::endl << "[Matrices Assembly]" << std::endl;

    if( verbose ) std::cout << "Setting up assembler... " << std::flush;
    ADRAssembler<mesh_type,matrix_type,vector_type> adrAssembler;
    adrAssembler.setup( uFESpace, uFESpace );
    if( verbose ) std::cout << "done" << std::endl;

    if( verbose ) std::cout << "Defining the matrices... " << std::flush;
    boost::shared_ptr<matrix_type> systemMatrix( new matrix_type( uFESpace->map() ) );
    *systemMatrix *= 0.0;
    if( verbose ) std::cout << "done" << std::endl;

    if( verbose ) std::cout << "Adding the viscous stress... " << std::flush;
    adrAssembler.addDiffusion( systemMatrix, 1.0 );
    if( verbose ) std::cout << "done" << std::endl;

    // +-----------------------------------------------+
    // |            Solver initialization              |
    // +-----------------------------------------------+
    if( verbose ) std::cout << std::endl << "[Solvers initialization]" << std::endl;
    prec_type* precRawPtr;
    basePrecPtr_type precPtr;
    precRawPtr = new PreconditionerIfpack;
    precRawPtr->setDataFromGetPot( dataFile, "prec" );
    precPtr.reset( precRawPtr );

    if( verbose ) std::cout << "Setting up SolverAztecOO... " << std::flush;
    SolverAztecOO linearSolver1;
    linearSolver1.setCommunicator( Comm );
    linearSolver1.setDataFromGetPot( dataFile, "solver" );
    linearSolver1.setTolerance( 1e-10 );
    linearSolver1.setPreconditioner( precPtr );
    if( verbose ) std::cout << "done" << std::endl;

    if( verbose ) std::cout << "Setting up LinearSolver (Belos)... " << std::flush;
    Teuchos::RCP< Teuchos::ParameterList > belosList2 = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList2 = Teuchos::getParametersFromXmlFile( "SolverParamList2.xml" );

    LinearSolver linearSolver2;
    linearSolver2.setCommunicator( Comm );
    linearSolver2.setParameters( *belosList2 );
    linearSolver2.setPreconditioner( precPtr );
    if( verbose ) std::cout << "done" << std::endl;
    linearSolver2.showMe();

    if( verbose ) std::cout << "Setting up LinearSolver (AztecOO)... " << std::flush;
    Teuchos::RCP< Teuchos::ParameterList > belosList3 = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList3 = Teuchos::getParametersFromXmlFile( "SolverParamList3.xml" );

    LinearSolver linearSolver3;
    linearSolver3.setCommunicator( Comm );
    linearSolver3.setParameters( *belosList3 );
    linearSolver3.setPreconditioner( precPtr );
    if( verbose ) std::cout << "done" << std::endl;
    linearSolver3.showMe();

    // +-----------------------------------------------+
    // |                   Simulation                  |
    // +-----------------------------------------------+
    if( verbose ) std::cout<< std::endl << "[Initialization of the simulation]" << std::endl;
    if( verbose ) std::cout << "Creation of vectors... " << std::flush;

    boost::shared_ptr<vector_type> rhs;
    rhs.reset( new vector_type( uFESpace->map(), Repeated ) );
    *rhs *= 0.0;

    vector_type fInterpolated( uFESpace->map(), Repeated );
    fInterpolated *= 0.0;
    adrAssembler.addMassRhs( *rhs, fRHS, 0.0 );
    rhs->globalAssemble();
    if( verbose ) std::cout << "done" << std::endl;

    fInterpolated *= 0.0;
    uFESpace->interpolate( uExact, fInterpolated, 0.0 );


    if( verbose ) std::cout << "Applying BC... " << std::flush;
    systemMatrix->globalAssemble();
    boost::shared_ptr<vector_type> rhsBC;
    rhsBC.reset( new vector_type( *rhs, Unique ) );
    bcManage( *systemMatrix, *rhsBC, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    if( verbose ) std::cout << "done" << std::endl;

    if( verbose ) std::cout << std::endl << "Solving the system with SolverAztec00... " << std::endl;
    boost::shared_ptr<vector_type> solution;
    solution.reset( new vector_type( uFESpace->map(), Unique ) );
    *solution *= 0.0;
    linearSolver1.setMatrix( *systemMatrix );
    linearSolver1.solveSystem( *rhsBC, *solution, systemMatrix );

    if( verbose ) std::cout << std::endl << "Solving the system with LinearSolver (Belos)... " << std::endl;
    boost::shared_ptr<vector_type> solution2;
    solution2.reset( new vector_type( uFESpace->map(), Unique ) );
    *solution2 *= 0.0;
    linearSolver2.setOperator( systemMatrix );
    linearSolver2.setRightHandSide( rhsBC );
    linearSolver2.solve( solution2 );

    if( verbose ) std::cout << std::endl << "Solving the system with LinearSolver (AztecOO)... " << std::endl;
    boost::shared_ptr<vector_type> solution3;
    solution3.reset( new vector_type( uFESpace->map(), Unique ) );
    *solution3 *= 0.0;
    linearSolver3.setOperator( systemMatrix );
    linearSolver3.setRightHandSide( rhsBC );
    linearSolver3.solve( solution3 );

    // +-----------------------------------------------+
    // |             Computing the error               |
    // +-----------------------------------------------+
    if( verbose ) std::cout << std::endl << "[Errors computation]" << std::endl;
    vector_type solutionErr( *solution );
    solutionErr *= 0.0;
    uFESpace->interpolate( Laplacian::uexact, solutionErr, 0.0 );
    solutionErr -= *solution;
    solutionErr.abs();

    vector_type solution2Err( *solution2 );
    solution2Err *= 0.0;
    uFESpace->interpolate( Laplacian::uexact,solution2Err, 0.0 );
    solution2Err -= *solution2;
    solution2Err.abs();

    vector_type solution3Err( *solution3 );
    solution3Err *= 0.0;
    uFESpace->interpolate( Laplacian::uexact, solution3Err, 0.0 );
    solution3Err -= *solution3;
    solution3Err.abs();

    vector_type solutionsDiff( *solution2 );
    solutionsDiff -= *solution;
    Real solutionsDiffNorm = solutionsDiff.norm2();

    vector_type solutionsDiff2( *solution2 );
    solutionsDiff2 -= *solution3;
    Real solutionsDiffNorm2 = solutionsDiff2.norm2();

    if( verbose ) std::cout << "AztecOO solver" << std::endl;
    printErrors( *solution, uFESpace,verbose );

    if( verbose ) std::cout << "Linear solver Belos" << std::endl;
    printErrors( *solution2, uFESpace,verbose );

    if( verbose ) std::cout << "Linear solver AztecOO" << std::endl;
    printErrors( *solution3, uFESpace,verbose );

    if( verbose ) std::cout << "Difference between the Azteco and the Belos solutions: " << solutionsDiffNorm << std::endl;
    if( verbose ) std::cout << "Difference between the two AztecOO solvers solutions: " << solutionsDiffNorm2 << std::endl;

    // +-----------------------------------------------+
    // |            Ending the simulation              |
    // +-----------------------------------------------+
    globalChrono.stop();
    if( verbose ) std::cout << std::endl << "Total simulation time: " << globalChrono.diff() << " s." << std::endl;
    if( verbose ) std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}


