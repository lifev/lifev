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
    @brief Test for the MeshPartitionTool class

    @author Radu Popescu <radu.popescu@epfl.ch>
    @date 14-11-2012
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

#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>

#include <lifev/core/mesh/MeshPartitionTool.hpp>
#include <lifev/core/mesh/GraphCutterZoltan.hpp>
#include <lifev/core/mesh/GraphCutterParMETIS.hpp>
#include <lifev/core/mesh/MeshPartBuilder.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

using namespace LifeV;

Real epsilon (1);

Real exactSolution (const Real& /* t */,
                    const Real& x,
                    const Real& y,
                    const Real& z,
                    const ID& /* i */)
{
    return  sin (x + y) + z * z / 2;
}


Real fRhs (const Real& /* t */,
           const Real& x,
           const Real& y,
           const Real& /* z */,
           const ID& /* i */)
{
    return  2 * sin (x + y) - 1;
}

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
typedef MeshPartitionTool < mesh_Type,
        GraphCutterZoltan,
        MeshPartBuilder > meshCutterZoltan_Type;
typedef MeshPartitionTool < mesh_Type,
        GraphCutterParMETIS,
        MeshPartBuilder > meshCutterParMETIS_Type;

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
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

    const UInt Nelements (dataFile ("mesh/nelements", 10) );
    if (verbose) std::cout << " ---> Number of elements : "
                               << Nelements << std::endl;

    const std::string zoltan_lb_method = dataFile ("mesh/zoltan_lb_method",
                                                   "GRAPH");
    const int zoltan_debug_level = dataFile ("mesh/zoltan_debug_level", 0);
    const int zoltan_hier_debug_level
        = dataFile ("mesh/zoltan_hier_debug_level", 0);

    // Build and partition the mesh

    if (verbose)
    {
        std::cout << " -- Building the mesh ... " << std::flush;
    }
    boost::shared_ptr< mesh_Type > fullMeshPtr (new RegionMesh<LinearTetra>);
    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Partitioning the mesh ... " << std::flush;
    }
    Teuchos::ParameterList meshParameters;
    meshParameters.set ("num_parts", Comm->NumProc(), "");
    meshParameters.set ("lb_method", zoltan_lb_method, "");
    meshParameters.set ("debug_level", zoltan_debug_level, "");
    meshParameters.set ("hier_debug_level", zoltan_hier_debug_level, "");
    if (verbose)
    {
        std::cout << " -- Using Zoltan ... " << std::flush;
    }
    meshCutterZoltan_Type meshPartZoltan (fullMeshPtr, Comm, meshParameters);
    if (! meshPartZoltan.success() )
    {
        if (verbose)
        {
            std::cout << "Zoltan partitioning failed." << std::endl;
        }
        return EXIT_FAILURE;
    }
    if (verbose)
    {
        std::cout << " -- Using ParMETIS ... " << std::flush;
    }
    meshCutterParMETIS_Type meshPartParMETIS (fullMeshPtr, Comm, meshParameters);
    if (! meshPartParMETIS.success() )
    {
        if (verbose)
        {
            std::cout << "ParMETIS partitioning failed." << std::endl;
        }
        return EXIT_FAILURE;
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

    {
        if (verbose)
        {
            std::cout << "\n=================================================\n"
                      << " -- Solving with the Zoltan mesh partition ...\n"
                      << "================================================="
                      << std::endl;
        }
        if (verbose)
        {
            std::cout << " -- Building FESpaces ... " << std::flush;
        }
        std::string uOrder ("P1");
        std::string bOrder ("P1");
        boost::shared_ptr < FESpace < mesh_Type,
              MapEpetra > >
              uFESpace (new FESpace < mesh_Type,
                        MapEpetra > (meshPartZoltan.meshPart(),
                                     uOrder,
                                     1,
                                     Comm) );

        boost::shared_ptr < FESpace < mesh_Type,
              MapEpetra > >
              betaFESpace (new FESpace < mesh_Type,
                           MapEpetra > (meshPartZoltan.meshPart(),
                                        bOrder,
                                        3,
                                        Comm) );

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
        adrAssembler.addDiffusion (systemMatrix, epsilon);
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }
        if (verbose)
        {
            std::cout << " Time needed : "
                      << adrAssembler.diffusionAssemblyChrono().diffCumul()
                      << std::endl;
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

        Real matrixNorm (systemMatrix->norm1() );
        if (verbose)
        {
            std::cout << " ---> Norm 1 : " << matrixNorm << std::endl;
        }
        if ( std::fabs (matrixNorm - 1.68421 ) > 1e-3)
        {
            std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }

        // Definition and assembly of the RHS

        if (verbose)
        {
            std::cout << " -- Building the RHS ... " << std::flush;
        }
        //vector_Type rhs(uFESpace->map(),Unique);
        vector_Type rhs (uFESpace->map(), Repeated);
        rhs *= 0.0;

        vector_Type fInterpolated (uFESpace->map(), Repeated);
        fInterpolated *= 0.0;
        uFESpace->interpolate (static_cast<feSpace_Type::function_Type> (fRhs),
                               fInterpolated, 0.0);
        adrAssembler.addMassRhs (rhs, fInterpolated);
        rhs.globalAssemble();

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Definition and application of the BCs

        if (verbose) std::cout << " -- Building the BCHandler ... "
                                   << std::flush;
        BCHandler bchandler;
        BCFunctionBase BCu ( exactSolution );
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
        boost::shared_ptr<vector_Type>
        rhsBC (new vector_Type (rhs, Unique) );
        bcManage (*systemMatrix, *rhsBC, *uFESpace->mesh(),
                  uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
        rhs = *rhsBC;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

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

        if (verbose) std::cout << " -- Setting up the solver ... "
                                   << std::flush;
        Teuchos::RCP<Teuchos::ParameterList>
        solverParam = Teuchos::rcp (new Teuchos::ParameterList);
        solverParam = Teuchos::getParametersFromXmlFile ("SolverParamList.xml");
        linearSolver.setCommunicator ( Comm );
        linearSolver.setParameters ( *solverParam );
        linearSolver.setPreconditionerFromGetPot (dataFile, "prec");
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose) std::cout << " -- Setting matrix in the solver ... "
                                   << std::flush;
        linearSolver.setOperator (systemMatrix);
        linearSolver.setRightHandSide (rhsBC);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        linearSolver.setCommunicator (Comm);

        // Definition of the solution

        if (verbose) std::cout << " -- Defining the solution ... "
                                   << std::flush;
        boost::shared_ptr<vector_Type>
        solution (new vector_Type (uFESpace->map(), Unique) );
        * (solution) *= 0.0;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Solve the solution

        if (verbose)
        {
            std::cout << " -- Solving the system ... " << std::flush;
        }
        linearSolver.solve (solution);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Error computation

        if (verbose)
        {
            std::cout << " -- Computing the error ... " << std::flush;
        }
        vector_Type solutionErr (*solution);
        solutionErr *= 0.0;
        uFESpace->interpolate (
            static_cast<feSpace_Type::function_Type> ( exactSolution ),
            solutionErr, 0.0);
        solutionErr -= *solution;
        solutionErr.abs();
        Real l2error (uFESpace->l2Error (exactSolution,
                                         vector_Type (*solution, Repeated), 0.0) );
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

        if (verbose) std::cout << " -- Defining the exporter ... "
                                   << std::flush;
        ExporterEnsight<mesh_Type>
        exporter (dataFile, meshPartZoltan.meshPart(),
                  "solution_zoltan", Comm->MyPID() ) ;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose) std::cout << " -- Defining the exported quantities ... "
                                   << std::flush;
        boost::shared_ptr<vector_Type>
        solutionPtr (new vector_Type (*solution, Repeated) );
        boost::shared_ptr<vector_Type>
        solutionErrPtr (new vector_Type (solutionErr, Repeated) );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose) std::cout << " -- Updating the exporter ... "
                                   << std::flush;
        exporter.addVariable (
            ExporterData<mesh_Type>::ScalarField, "solution",
            uFESpace, solutionPtr, UInt (0) );
        exporter.addVariable (
            ExporterData<mesh_Type>::ScalarField,
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
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Zoltan test successful! " << std::endl;
        }
    }

    {
        if (verbose)
        {
            std::cout << "\n=================================================\n"
                      << " -- Solving with the ParMETIS mesh partition ...\n"
                      << "=================================================\n"
                      << std::endl;
        }
        if (verbose)
        {
            std::cout << " -- Building FESpaces ... " << std::flush;
        }
        std::string uOrder ("P1");
        std::string bOrder ("P1");
        boost::shared_ptr < FESpace < mesh_Type,
              MapEpetra > >
              uFESpace (new FESpace < mesh_Type,
                        MapEpetra > (meshPartParMETIS.meshPart(),
                                     uOrder,
                                     1,
                                     Comm) );

        boost::shared_ptr < FESpace < mesh_Type,
              MapEpetra > >
              betaFESpace (new FESpace < mesh_Type,
                           MapEpetra > (meshPartParMETIS.meshPart(),
                                        bOrder,
                                        3,
                                        Comm) );

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
        adrAssembler.addDiffusion (systemMatrix, epsilon);
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }
        if (verbose)
        {
            std::cout << " Time needed : "
                      << adrAssembler.diffusionAssemblyChrono().diffCumul()
                      << std::endl;
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

        Real matrixNorm (systemMatrix->norm1() );
        if (verbose)
        {
            std::cout << " ---> Norm 1 : " << matrixNorm << std::endl;
        }
        if ( std::fabs (matrixNorm - 1.68421 ) > 1e-3)
        {
            std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }

        // Definition and assembly of the RHS

        if (verbose)
        {
            std::cout << " -- Building the RHS ... " << std::flush;
        }
        //vector_Type rhs(uFESpace->map(),Unique);
        vector_Type rhs (uFESpace->map(), Repeated);
        rhs *= 0.0;

        vector_Type fInterpolated (uFESpace->map(), Repeated);
        fInterpolated *= 0.0;
        uFESpace->interpolate (static_cast<feSpace_Type::function_Type> (fRhs),
                               fInterpolated, 0.0);
        adrAssembler.addMassRhs (rhs, fInterpolated);
        rhs.globalAssemble();

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Definition and application of the BCs

        if (verbose) std::cout << " -- Building the BCHandler ... "
                                   << std::flush;
        BCHandler bchandler;
        BCFunctionBase BCu ( exactSolution );
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
        boost::shared_ptr<vector_Type>
        rhsBC (new vector_Type (rhs, Unique) );
        bcManage (*systemMatrix, *rhsBC, *uFESpace->mesh(),
                  uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
        rhs = *rhsBC;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

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

        if (verbose) std::cout << " -- Setting up the solver ... "
                                   << std::flush;
        Teuchos::RCP<Teuchos::ParameterList>
        solverParam = Teuchos::rcp (new Teuchos::ParameterList);
        solverParam = Teuchos::getParametersFromXmlFile ("SolverParamList.xml");
        linearSolver.setCommunicator ( Comm );
        linearSolver.setParameters ( *solverParam );
        linearSolver.setPreconditionerFromGetPot (dataFile, "prec");
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose) std::cout << " -- Setting matrix in the solver ... "
                                   << std::flush;
        linearSolver.setOperator (systemMatrix);
        linearSolver.setRightHandSide (rhsBC);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        linearSolver.setCommunicator (Comm);

        // Definition of the solution

        if (verbose) std::cout << " -- Defining the solution ... "
                                   << std::flush;
        boost::shared_ptr<vector_Type>
        solution (new vector_Type (uFESpace->map(), Unique) );
        * (solution) *= 0.0;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Solve the solution

        if (verbose)
        {
            std::cout << " -- Solving the system ... " << std::flush;
        }
        linearSolver.solve (solution);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Error computation

        if (verbose)
        {
            std::cout << " -- Computing the error ... " << std::flush;
        }
        vector_Type solutionErr (*solution);
        solutionErr *= 0.0;
        uFESpace->interpolate (
            static_cast<feSpace_Type::function_Type> ( exactSolution ),
            solutionErr, 0.0);
        solutionErr -= *solution;
        solutionErr.abs();
        Real l2error (uFESpace->l2Error (exactSolution,
                                         vector_Type (*solution, Repeated), 0.0) );
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

        if (verbose) std::cout << " -- Defining the exporter ... "
                                   << std::flush;
        ExporterEnsight<mesh_Type>
        exporter (dataFile, meshPartParMETIS.meshPart(),
                  "solution_parmetis", Comm->MyPID() );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose) std::cout << " -- Defining the exported quantities ... "
                                   << std::flush;
        boost::shared_ptr<vector_Type>
        solutionPtr (new vector_Type (*solution, Repeated) );
        boost::shared_ptr<vector_Type>
        solutionErrPtr (new vector_Type (solutionErr, Repeated) );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose) std::cout << " -- Updating the exporter ... "
                                   << std::flush;
        exporter.addVariable (ExporterData<mesh_Type>::ScalarField,
                              "solution", uFESpace, solutionPtr, UInt (0) );
        exporter.addVariable (ExporterData<mesh_Type>::ScalarField,
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
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- ParMETIS test successful! " << std::endl;
        }
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


