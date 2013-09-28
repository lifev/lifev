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

    @author
    @date
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

#include <lifev/core/array/MatrixEpetra.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#else
#include <lifev/core/filter/ExporterVTK.hpp>
#endif

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}

Real exactSolution ( const Real& /* t */, const Real& x, const Real& /*y*/, const Real& /* z */, const ID& /* i */ )
{
    return std::sin ( M_PI * 0.5 * x );
}


Real fRhs ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
    return 1.25 * M_PI * M_PI * std::sin ( M_PI * 0.5 * x );
}


int
main ( int argc, char* argv[] )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif

    {
        // needed to properly destroy all objects inside before mpi finalize

#ifdef HAVE_MPI
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
        ASSERT ( Comm->NumProc() < 2, "The test does not run in parallel." );
#else
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

        typedef RegionMesh<LinearLine> mesh_Type;
        typedef MatrixEpetra<Real> matrix_Type;
        typedef VectorEpetra vector_Type;
        typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
        typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

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

        // Build the mesh

        if (verbose)
        {
            std::cout << " -- Reading the mesh ... " << std::flush;
        }
        MeshData meshData (dataFile, "mesh");
        boost::shared_ptr< mesh_Type > meshPtr ( new mesh_Type ( Comm ) );

        // Set up the structured mesh
        regularMesh1D ( *meshPtr, 0,
                        dataFile ( "mesh/n", 20 ),
                        dataFile ( "mesh/verbose", false ),
                        dataFile ( "mesh/length", 1. ),
                        dataFile ( "mesh/origin", 0. ) );

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Build the FESpaces
        if (verbose)
        {
            std::cout << " -- Building FESpaces ... " << std::flush;
        }
        feSpacePtr_Type uFESpace ( new feSpace_Type ( meshPtr, feSegP1, quadRuleSeg1pt, quadRuleNode1pt, 1, Comm ) );
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
        adrAssembler.setFespace (uFESpace);
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Defining the matrix ... " << std::flush;
        }
        boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uFESpace->map() ) );
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }

        // Perform the assembly of the matrix

        if (verbose)
        {
            std::cout << " -- Adding the diffusion ... " << std::flush;
        }
        adrAssembler.addDiffusion (systemMatrix, 1.);
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }
        if (verbose)
        {
            std::cout << " Time needed : " << adrAssembler.diffusionAssemblyChrono().diffCumul() << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Adding the mass ... " << std::flush;
        }
        adrAssembler.addMass (systemMatrix, M_PI * M_PI );
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Closing the matrix ... " << std::flush;
        }
        systemMatrix->globalAssemble();
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Definition and assembly of the RHS
        if (verbose)
        {
            std::cout << " -- Building the RHS ... " << std::flush;
        }
        vector_Type rhs (uFESpace->map(), Repeated);

        vector_Type fInterpolated (uFESpace->map(), Repeated);
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), fInterpolated, 0.0 );
        adrAssembler.addMassRhs (rhs, fInterpolated);
        rhs.globalAssemble();

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

        bchandler.addBC ("Dirichlet", Structured1DLabel::LEFT, Essential, Full, BCu, 1);
        bchandler.addBC ("Dirichlet", Structured1DLabel::RIGHT, Essential, Full, BCu, 1);

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
        vector_Type rhsBC (rhs, Unique);
        bcManage (*systemMatrix, rhsBC, *uFESpace->mesh(), uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
        rhs = rhsBC;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Definition of the solver
        if (verbose)
        {
            std::cout << " -- Building the solver ... " << std::flush;
        }
        SolverAztecOO linearSolver;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Setting up the solver ... " << std::flush;
        }
        linearSolver.setDataFromGetPot (dataFile, "solver");
        linearSolver.setupPreconditioner (dataFile, "prec");
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Setting matrix in the solver ... " << std::flush;
        }
        linearSolver.setMatrix (*systemMatrix);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        linearSolver.setCommunicator (Comm);

        // Definition of the solution
        if (verbose)
        {
            std::cout << " -- Defining the solution ... " << std::flush;
        }
        vector_Type solution (uFESpace->map(), Unique);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Solve the solution
        if (verbose)
        {
            std::cout << " -- Solving the system ... " << std::flush;
        }
        linearSolver.solveSystem (rhsBC, solution, systemMatrix);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Error computation
        if (verbose)
        {
            std::cout << " -- Computing the error ... " << std::flush;
        }
        vector_Type solutionErr (solution);
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactSolution ), solutionErr, 0.0 );
        solutionErr -= solution;
        solutionErr.abs();
        Real l2error (uFESpace->l2Error (exactSolution, vector_Type (solution, Repeated), 0.0) );
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

        // Exporter definition and use
        if (verbose)
        {
            std::cout << " -- Defining the exporter ... " << std::flush;
        }

#ifdef HAVE_HDF5
        ExporterHDF5<mesh_Type> exporter ( dataFile, "solution" );
#else
        ExporterVTK<mesh_Type> exporter ( dataFile, "solution" );
#endif
        exporter.setMeshProcId ( meshPtr, Comm->MyPID() ) ;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Defining the exported quantities ... " << std::flush;
        }
        boost::shared_ptr<vector_Type> solutionPtr (new vector_Type (solution, Repeated) );
        boost::shared_ptr<vector_Type> solutionErrPtr (new vector_Type (solutionErr, Repeated) );
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

        if ( std::fabs ( l2error - 0.0006169843149652788l ) > 1e-10 )
        {
            std::cout << " <!> Solution has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }
        if ( std::fabs ( linferror - 0.0001092814405985187l ) > 1e-10 )
        {
            std::cout << " <!> Solution has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }

        if (verbose)
        {
            std::cout << "End Result: TEST PASSED" << std::endl;
        }

    } // needed to properly destroy all objects inside before mpi finalize

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}
