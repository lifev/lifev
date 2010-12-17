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
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/life.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include <mpi.h>

#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifearray/EpetraMatrix.hpp>

#include <life/lifefilters/ensight.hpp>
//#include <life/lifefilters/hdf5exporter.hpp>

#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bcManage.hpp>

#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifemesh/structuredMesh3D.hpp>
#include <life/lifemesh/regionMesh3D.hpp>

#include <life/lifesolver/ADRAssembler.hpp>
#include <life/lifemesh/dataMesh.hpp>

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
    if (i == 1)
        return 1;
    return 0;
}
#endif

#ifdef TEST_RHS
Real epsilon(1);

Real exactSolution( const Real& /* t */, const Real& x, const Real& y, const Real& z , const ID& /* i */ )
{
    return  sin(x+y)+z*z/2;
}


Real fRhs( const Real& /* t */, const Real& x, const Real& y, const Real& /* z */ , const ID& /* i */ )
{
    return  2*sin(x+y)-1;
}
#endif


typedef RegionMesh3D<LinearTetra> mesh_type;
typedef EpetraMatrix<Real> matrix_type;
typedef EpetraVector vector_type;

int
main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_SerialComm);
#endif

    const bool verbose(Comm->MyPID()==0);

// Read first the data needed

    if (verbose) std::cout << " -- Reading the data ... " << std::flush;
    GetPot dataFile( "data" );
    if (verbose) std::cout << " done ! " << std::endl;

    const UInt Nelements(dataFile("mesh/nelements",10));
    if (verbose) std::cout << " ---> Number of elements : " << Nelements << std::endl;

// Build and partition the mesh

    if (verbose) std::cout << " -- Building the mesh ... " << std::flush;
    boost::shared_ptr< mesh_type > fullMeshPtr(new RegionMesh3D<LinearTetra>);
    regularMesh3D( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                   2.0,   2.0,   2.0,
                   -1.0,  -1.0,  -1.0);
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Partitioning the mesh ... " << std::flush;
    partitionMesh< mesh_type >   meshPart(fullMeshPtr, Comm);
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Freeing the global mesh ... " << std::flush;
    fullMeshPtr.reset();
    if (verbose) std::cout << " done ! " << std::endl;

// Build the FESpaces

    if (verbose) std::cout << " -- Building FESpaces ... " << std::flush;
    std::string uOrder("P1");
    std::string bOrder("P1");
    boost::shared_ptr<FESpace< mesh_type, EpetraMap > > uFESpace( new FESpace< mesh_type, EpetraMap >(meshPart,uOrder, 1, Comm));
    boost::shared_ptr<FESpace< mesh_type, EpetraMap > > betaFESpace( new FESpace< mesh_type, EpetraMap >(meshPart,bOrder, 3, Comm));
    if (verbose) std::cout << " done ! " << std::endl;
    if (verbose) std::cout << " ---> Dofs: " << uFESpace->dof().numTotalDof() << std::endl;

// Build the assembler and the matrices

    if (verbose) std::cout << " -- Building assembler ... " << std::flush;
    ADRAssembler<mesh_type,matrix_type,vector_type> adrAssembler;
    if (verbose) std::cout << " done! " << std::endl;

    if (verbose) std::cout << " -- Setting up assembler ... " << std::flush;
    adrAssembler.setup(uFESpace,betaFESpace);
    if (verbose) std::cout << " done! " << std::endl;

    if (verbose) std::cout << " -- Defining the matrix ... " << std::flush;
    boost::shared_ptr<matrix_type> systemMatrix(new matrix_type( uFESpace->map() ));
    *systemMatrix *=0.0;
    if (verbose) std::cout << " done! " << std::endl;

// Perform the assembly of the matrix

    if (verbose) std::cout << " -- Adding the diffusion ... " << std::flush;
    adrAssembler.addDiffusion(systemMatrix,epsilon);
    if (verbose) std::cout << " done! " << std::endl;
    if (verbose) std::cout << " Time needed : " << adrAssembler.diffusionAssemblyChrono().diffCumul() << std::endl;

#ifdef TEST_ADVECTION
    if (verbose) std::cout << " -- Adding the advection ... " << std::flush;
    vector_type beta(betaFESpace->map(),Repeated);
    betaFESpace->interpolate(betaFct,beta,0.0);
    adrAssembler.addAdvection(systemMatrix,beta);
    if (verbose) std::cout << " done! " << std::endl;
#endif
#ifdef TEST_MASS
    if (verbose) std::cout << " -- Adding the mass ... " << std::flush;
    adrAssembler.addMass(systemMatrix,1.0);
    if (verbose) std::cout << " done! " << std::endl;
#endif

    if (verbose) std::cout << " -- Closing the matrix ... " << std::flush;
    systemMatrix->GlobalAssemble();
    if (verbose) std::cout << " done ! " << std::endl;

#ifdef TEST_RHS
    Real matrixNorm(systemMatrix->NormOne());
    if (verbose) std::cout << " ---> Norm 1 : " << matrixNorm << std::endl;
    if ( abs(matrixNorm - 1.68421 ) > 1e-3)
    {
        std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }
#endif

// Definition and assembly of the RHS

    if (verbose) std::cout << " -- Building the RHS ... " << std::flush;
    //vector_type rhs(uFESpace->map(),Unique);
    vector_type rhs(uFESpace->map(),Repeated);
    rhs*=0.0;

#ifdef TEST_RHS
    vector_type fInterpolated(uFESpace->map(),Repeated);
    fInterpolated*=0.0;
    uFESpace->interpolate(fRhs,fInterpolated,0.0);
    adrAssembler.addMassRhs(rhs,fInterpolated);
    rhs.GlobalAssemble();
#endif

    if (verbose) std::cout << " done ! " << std::endl;

// Definition and application of the BCs

    if (verbose) std::cout << " -- Building the BCHandler ... " << std::flush;
    BCHandler bchandler(26);
    BCFunctionBase BCu( exactSolution );
    bchandler.addBC("Dirichlet",1,Essential,Full,BCu,1);
    for (UInt i(2); i<=26; ++i)
    {
        bchandler.addBC("Dirichlet",i,Essential,Full,BCu,1);
    }
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Updating the BCs ... " << std::flush;
    bchandler.bdUpdate(*uFESpace->mesh(),uFESpace->feBd(),uFESpace->dof());
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Applying the BCs ... " << std::flush;
    vector_type rhsBC(rhs,Unique);
    bcManage(*systemMatrix,rhsBC,*uFESpace->mesh(),uFESpace->dof(),bchandler,uFESpace->feBd(),1.0,0.0);
    rhs = rhsBC;
    if (verbose) std::cout << " done ! " << std::endl;

    //************* SPY ***********
    //systemMatrix->spy("matrix");
    //rhs.spy("vector");
    //*****************************

// Definition of the solver

    if (verbose) std::cout << " -- Building the solver ... " << std::flush;
    SolverTrilinos linearSolver;
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Setting up the solver ... " << std::flush;
    linearSolver.setDataFromGetPot(dataFile,"solver");
    linearSolver.setUpPrec(dataFile,"prec");
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Setting matrix in the solver ... " << std::flush;
    linearSolver.setMatrix(*systemMatrix);
    if (verbose) std::cout << " done ! " << std::endl;

    linearSolver.setCommunicator(Comm);

// Definition of the solution

    if (verbose) std::cout << " -- Defining the solution ... " << std::flush;
    vector_type solution(uFESpace->map(),Unique);
    solution*=0.0;
    if (verbose) std::cout << " done ! " << std::endl;

// Solve the solution

    if (verbose) std::cout << " -- Solving the system ... " << std::flush;
    linearSolver.solveSystem(rhsBC,solution,systemMatrix);
    if (verbose) std::cout << " done ! " << std::endl;

    //************* SPY ***********
    //solution.spy("solution");
    //*****************************

// Error computation

    if (verbose) std::cout << " -- Computing the error ... " << std::flush;
    vector_type solutionErr(solution);
    solutionErr*=0.0;
    uFESpace->interpolate(exactSolution,solutionErr,0.0);
    solutionErr-=solution;
    solutionErr.Abs();
    Real l2error(uFESpace->L2Error(exactSolution,vector_type(solution,Repeated),0.0));
    if (verbose) std::cout << " -- done ! " << std::endl;
    if (verbose) std::cout << " ---> Norm L2  : " << l2error << std::endl;
    Real linferror( solutionErr.NormInf());
    if (verbose) std::cout << " ---> Norm Inf : " << linferror << std::endl;


    if (l2error > 0.0055)
    {
        std::cout << " <!> Solution has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }
    if (linferror > 0.0046)
    {
        std::cout << " <!> Solution has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }

// Exporter definition and use

    if (verbose) std::cout << " -- Defining the exporter ... " << std::flush;
    Ensight<mesh_type> exporter ( dataFile, meshPart.meshPartition(), "solution", Comm->MyPID()) ;
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Defining the exported quantities ... " << std::flush;
    boost::shared_ptr<vector_type> solutionPtr (new vector_type(solution,Repeated));
    boost::shared_ptr<vector_type> solutionErrPtr (new vector_type(solutionErr,Repeated));
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Updating the exporter ... " << std::flush;
    exporter.addVariable( ExporterData::Scalar, "solution", solutionPtr, UInt(0), uFESpace->dof().numTotalDof() );
    exporter.addVariable( ExporterData::Scalar, "error", solutionErrPtr, UInt(0), uFESpace->dof().numTotalDof() );
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Exporting ... " << std::flush;
    exporter.postProcess(0);
    if (verbose) std::cout << " done ! " << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}


