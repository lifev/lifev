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


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>


#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}

//#define TEST_MASS
//#define TEST_ADVECTION
#define TEST_RHS


#ifdef TEST_MASS
Real epsilon (1);

Real exactSolution ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
    Real seps (s td::sqrt (epsilon) );
    return  std::exp (seps * x) / (std::exp (seps) - std::exp (-seps) );
}
#endif

#ifdef TEST_ADVECTION
Real epsilon (1);

Real exactSolution ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */, const ID& /* i */ )
{
    return  ( std::exp (x / epsilon) - 1 ) / ( std::exp (1 / epsilon) - 1 );
}

Real betaFct ( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& i )
{
    return (i == 0);
}
#endif

#ifdef TEST_RHS
Real epsilon (1);

Real exactSolution ( const Real& /* t */, const Real& x, const Real& y, const Real& /* z */, const ID& /* i */ )
{
    return  std::sin (x) + y * y / 2;
}


Real fRhs ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
    return  std::sin (x) - 1;
}
#endif


typedef RegionMesh<LinearTriangle> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;

typedef LifeV::Preconditioner basePrec_Type;
typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
typedef LifeV::PreconditionerIfpack prec_Type;
typedef boost::shared_ptr<prec_Type> precPtr_Type;

typedef std::shared_ptr<feSpace_Type> feSpacePtr_Type;


int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);

    // Read first the data needed

    if (verbose)
    {
        std::cout << " -- Reading the data ... " << std::flush;
    }
    GetPot dataFile ( "data" );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // Build and partition the mesh

    if (verbose)
    {
        std::cout << " -- Reading the mesh ... " << std::flush;
    }
    MeshData meshData (dataFile, "mesh");
    std::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );

    // Select if the mesh is structured or not
    if ( meshData.meshType() != "structured" )
    {
        // Set up the mesh
        readMesh ( *fullMeshPtr, meshData );
    }
    else
    {
        // Set up the structured mesh
        regularMesh2D ( *fullMeshPtr, 0,
                        dataFile ( "mesh/nx", 20 ),
                        dataFile ( "mesh/ny", 20 ),
                        dataFile ( "mesh/verbose", false ),
                        dataFile ( "mesh/lx", 1. ),
                        dataFile ( "mesh/ly", 1. ) );
    }

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Partitioning the mesh ... " << std::flush;
    }
    std::shared_ptr< mesh_Type > meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
        meshPtr = meshPart.meshPartition();
    }
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Freeing the global mesh ... " << std::flush;
    }
    fullMeshPtr.reset();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    // Build the FESpaces

    if (verbose)
    {
        std::cout << " -- Building FESpaces ... " << std::flush;
    }
    std::string uOrder ("P1");
    std::string bOrder ("P1");
    feSpacePtr_Type uFESpace ( new feSpace_Type ( meshPtr, uOrder, 1, Comm ) );
    feSpacePtr_Type betaFESpace ( new feSpace_Type ( meshPtr, bOrder, 3, Comm ) );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << uFESpace->dof().numTotalDof() << std::endl;
    }

    // Build the assembler and the matrices

    if (verbose)
    {
        std::cout << " -- Building assembler ... " << std::flush;
    }
    ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Setting up assembler ... " << std::flush;
    }
    adrAssembler.setup (uFESpace, betaFESpace);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Defining the matrix ... " << std::flush;
    }
    std::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uFESpace->map() ) );
    *systemMatrix *= 0.0;
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    // Perform the assembly of the matrix

    if (verbose)
    {
        std::cout << " -- Adding the diffusion ... " << std::flush;
    }
    adrAssembler.addDiffusion (systemMatrix, epsilon);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time needed : " << adrAssembler.diffusionAssemblyChrono().diffCumul() << std::endl;
    }

#ifdef TEST_ADVECTION
    if (verbose)
    {
        std::cout << " -- Adding the advection ... " << std::flush;
    }
    vector_Type beta (betaFESpace->map(), Repeated);
    betaFESpace->interpolate (betaFct, beta, 0.0);
    adrAssembler.addAdvection (systemMatrix, beta);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
#endif
#ifdef TEST_MASS
    if (verbose)
    {
        std::cout << " -- Adding the mass ... " << std::flush;
    }
    adrAssembler.addMass (systemMatrix, 1.0);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
#endif

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }
    systemMatrix->globalAssemble();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

#ifdef TEST_RHS
    Real matrixNorm (systemMatrix->norm1() );
    if (verbose)
    {
        std::cout << " ---> Norm 1 : " << matrixNorm << std::endl;
    }
    if ( std::fabs (matrixNorm - 8 ) > 1e-3)
    {
        std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }
#endif

    // Definition and assembly of the RHS

    if (verbose)
    {
        std::cout << " -- Building the RHS ... " << std::flush;
    }
    //vector_Type rhs(uFESpace->map(),Unique);
    vector_Type rhs (uFESpace->map(), Repeated);
    rhs *= 0.0;

#ifdef TEST_RHS
    vector_Type fInterpolated (uFESpace->map(), Repeated);
    fInterpolated *= 0.0;
    uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), fInterpolated, 0.0 );
    adrAssembler.addMassRhs (rhs, fInterpolated);
    rhs.globalAssemble();
#endif

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    // Definition and application of the BCs

    if (verbose)
    {
        std::cout << " -- Building the BCHandler ... " << std::flush;
    }
    BCHandler bchandler;
    BCFunctionBase BCu ( static_cast<feSpace_Type::function_Type> ( exactSolution ) );
    bchandler.addBC ("Dirichlet", 1, Essential, Full, BCu, 1);
    for (UInt i (2); i <= 4; ++i)
    {
        bchandler.addBC ("Dirichlet", i, Essential, Full, BCu, 1);
    }
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Updating the BCs ... " << std::flush;
    }
    bchandler.bcUpdate (*uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Applying the BCs ... " << std::flush;
    }
    vectorPtr_Type rhsBC ( new vector_Type (rhs, Unique) );
    bcManage (*systemMatrix, *rhsBC, *uFESpace->mesh(), uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
    rhs = *rhsBC;
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    //************* SPY ***********
    //systemMatrix->spy("matrix");
    //rhs.spy("vector");
    //*****************************

    // Definition of the solver

    if (verbose)
    {
        std::cout << " -- Building the solver ... " << std::flush;
    }
    LinearSolver linearSolver;
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Setting up the solver ... " << std::flush;
    }

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

    linearSolver.setCommunicator ( Comm );
    linearSolver.setParameters ( *belosList );

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );
    linearSolver.setPreconditioner ( precPtr );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Setting matrix in the solver ... " << std::flush;
    }
    linearSolver.setOperator ( systemMatrix );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    // Definition of the solution

    if (verbose)
    {
        std::cout << " -- Defining the solution ... " << std::flush;
    }
    vectorPtr_Type solution ( new vector_Type (uFESpace->map(), Unique) );
    solution->zero();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    // Solve the solution

    if (verbose)
    {
        std::cout << " -- Solving the system ... " << std::flush;
    }
    linearSolver.setRightHandSide ( rhsBC );
    linearSolver.solve ( solution );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    //************* SPY ***********
    //solution.spy("solution");
    //*****************************

    // Error computation

    if (verbose)
    {
        std::cout << " -- Computing the error ... " << std::flush;
    }
    vector_Type solutionErr (*solution);
    solutionErr *= 0.0;
    uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactSolution ), solutionErr, 0.0 );
    solutionErr -= *solution;
    solutionErr.abs();
    Real l2error (uFESpace->l2Error (exactSolution, vector_Type (*solution, Repeated), 0.0) );
    if (verbose)
    {
        std::cout << " -- done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Norm L2  : " << l2error << std::endl;
    }
    Real linferror ( solutionErr.normInf() );
    if (verbose)
    {
        std::cout << " ---> Norm Inf : " << linferror << std::endl;
    }


    if (l2error > 0.0026)
    {
        std::cout << " <!> Solution has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }
    if (linferror > 0.000016)
    {
        std::cout << " <!> Solution has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }

    // Exporter definition and use

    if (verbose)
    {
        std::cout << " -- Defining the exporter ... " << std::flush;
    }
    ExporterEnsight<mesh_Type> exporter ( dataFile, meshPtr, "solution", Comm->MyPID() ) ;
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Defining the exported quantities ... " << std::flush;
    }
<<<<<<< HEAD
    boost::shared_ptr<vector_Type> solutionPtr (new vector_Type (*solution, Repeated) );
    boost::shared_ptr<vector_Type> solutionErrPtr (new vector_Type (solutionErr, Repeated) );
=======
    std::shared_ptr<vector_Type> solutionPtr (new vector_Type (solution, Repeated) );
    std::shared_ptr<vector_Type> solutionErrPtr (new vector_Type (solutionErr, Repeated) );
>>>>>>> 24ac07b... Versione c++11
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Updating the exporter ... " << std::flush;
    }
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "solution", uFESpace, solutionPtr, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "error", uFESpace, solutionErrPtr, UInt (0) );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Exporting ... " << std::flush;
    }
    exporter.postProcess (0);
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}


