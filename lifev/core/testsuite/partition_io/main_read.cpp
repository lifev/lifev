//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010, 2011, 2012 EPFL, Politecnico di Milano, Emory University

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
    @brief Test for PartitionIO class - read and solve

    @author Radu Popescu <radu.popescu@epfl.ch>
    @maintainer Radu Popescu <radu.popescu@epfl.ch>
    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 10-05-2012

    Loads mesh parts from HDF5 and solves a Laplacian problem.
    Based on the ADRAssembler class unit test.
 */

#include <lifev/core/LifeV.hpp>

#ifdef LIFEV_HAS_HDF5

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "Epetra_config.h"

#ifdef HAVE_MPI

#include <mpi.h>

#include <Epetra_MpiComm.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/ExporterHDF5Mesh3D.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>

using namespace LifeV;

Real exactSolution (const Real& /* t */,
                    const Real& x,
                    const Real& y,
                    const Real& z,
                    const ID& /* i */)
{
    return sin (x + y) + z * z / 2;
}
Real fRhs (const Real& /* t */,
           const Real& x,
           const Real& y,
           const Real& /* z */,
           const ID& /* i */ )
{
    return 2 * sin (x + y) - 1;
}

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef boost::function < Real ( Real const&,
                                 Real const&,
                                 Real const&,
                                 Real const&,
                                 UInt const& ) > function_Type;

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

int
main ( int argc, char** argv )
{
#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );

    const bool verbose (comm->MyPID() == 0);

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

    const UInt Nelements (dataFile ("mesh/nelements", 10) );
    if (verbose) std::cout << " ---> Number of elements : "
                               << Nelements << std::endl;

    // Load mesh part from HDF5
    const std::string partsFileName (dataFile ("test/hdf5_file_name", "cube.h5") );
    const std::string ioClass (dataFile ("test/io_class", "new") );

    boost::shared_ptr<mesh_Type> mesh;
    if (! ioClass.compare ("old") )
    {
        ExporterHDF5Mesh3D<mesh_Type> HDF5Input (dataFile, partsFileName);
        HDF5Input.setComm (comm);
        mesh = HDF5Input.getMeshPartition();
        HDF5Input.closeFile();
    }
    else
    {
        PartitionIO<RegionMesh<LinearTetra> > partitionIO (partsFileName, comm);
        partitionIO.read (mesh);
    }

    // Build the FESpaces

    if (verbose)
    {
        std::cout << " -- Building FESpaces ... " << std::flush;
    }
    std::string uOrder ("P1");
    std::string bOrder ("P1");
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> >
    uFESpace (new FESpace<mesh_Type, MapEpetra> (mesh, uOrder, 1, comm) );
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> >
    betaFESpace (new FESpace<mesh_Type, MapEpetra> (mesh, bOrder, 3, comm) );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose) std::cout << " ---> Dofs: "
                               << uFESpace->dof().numTotalDof() << std::endl;

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
    boost::shared_ptr<matrix_Type>
    systemMatrix (new matrix_Type (uFESpace->map() ) );
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
    adrAssembler.addDiffusion (systemMatrix, 1);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose) std::cout << " Time needed : "
                               << adrAssembler.diffusionAssemblyChrono().diffCumul()
                               << std::endl;

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }
    systemMatrix->globalAssemble();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    Real matrixNorm (systemMatrix->norm1() );
    if (verbose)
    {
        std::cout << " ---> Norm 1 : " << matrixNorm << std::endl;
    }
    if (std::fabs (matrixNorm - 1.68421) > 1e-3)
    {
        std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }

    // Definition and assembly of the RHS

    if (verbose)
    {
        std::cout << " -- Building the RHS ... " << std::flush;
    }
    vector_Type rhs (uFESpace->map(), Repeated);
    rhs *= 0.0;

    vector_Type fInterpolated (uFESpace->map(), Repeated);
    fInterpolated *= 0.0;
    uFESpace->interpolate ( static_cast<function_Type> ( fRhs ), fInterpolated, 0.0);
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
    BCFunctionBase BCu (exactSolution);
    bchandler.addBC ("Dirichlet", 1, Essential, Full, BCu, 1);
    for (UInt i (2); i <= 6; ++i)
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
    vector_Type rhsBC (rhs, Unique);
    bcManage (*systemMatrix, rhsBC, *uFESpace->mesh(),
              uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
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

    if (verbose) std::cout << " -- Setting matrix in the solver ... "
                               << std::flush;
    linearSolver.setMatrix (*systemMatrix);
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    linearSolver.setCommunicator (comm);

    // Definition of the solution

    if (verbose)
    {
        std::cout << " -- Defining the solution ... " << std::flush;
    }
    vector_Type solution (uFESpace->map(), Unique);
    solution *= 0.0;
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
    solutionErr *= 0.0;
    uFESpace->interpolate ( static_cast<function_Type> ( exactSolution ), solutionErr, 0.0);
    solutionErr -= solution;
    solutionErr.abs();
    Real l2error (uFESpace->l2Error (exactSolution,
                                     vector_Type (solution, Repeated), 0.0) );
    if (verbose)
    {
        std::cout << " -- done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Norm L2  : " << l2error << std::endl;
    }
    Real linferror (solutionErr.normInf() );
    if (verbose)
    {
        std::cout << " ---> Norm Inf : " << linferror << std::endl;
    }


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

    if (verbose)
    {
        std::cout << " -- Defining the exporter ... " << std::flush;
    }
    ExporterHDF5<mesh_Type> exporter (dataFile, mesh,
                                      "solution", comm->MyPID() ) ;
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose) std::cout << " -- Defining the exported quantities ... "
                               << std::flush;
    boost::shared_ptr<vector_Type>
    solutionPtr (new vector_Type (solution, Repeated) );
    boost::shared_ptr<vector_Type>
    solutionErrPtr (new vector_Type (solutionErr, Repeated) );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose) std::cout << " -- Updating the exporter ... "
                               << std::flush;
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,
                           "solution", uFESpace, solutionPtr, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,
                           "error", uFESpace, solutionErrPtr, UInt (0) );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Exporting ... " << std::flush;
    }
    exporter.postProcess (0);
    exporter.closeFile();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

    MPI_Finalize();

#else
    std::cout << "This test needs MPI to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* HAVE_MPI */
#else
    std::cout << "This test needs HDF5 to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* LIFEV_HAS_HDF5 */

    return ( EXIT_SUCCESS );
}


