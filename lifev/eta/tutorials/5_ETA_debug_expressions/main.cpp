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
    @brief Tutorial explaining how to debug a failing expression.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 29-06-2012

    The ETA framework relies heavily on template mechanisms to
    work. This provides a clean interface and good performances.
    A drawback is that the compilation errors that one can get
    when making a mistake in the expression are extremly
    complicated and nearly impossible to understand for non-expert
    developers.

    We explain here how to find out where is the problem and what
    are the exact criteria that makes an expression well-formed
    or ill-formed.

    Tutorials that should be read before: 1,2

 */

// ---------------------------------------------------------------
// We include the usual headers (see tutorial 1 for explanations)
// ---------------------------------------------------------------

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <boost/shared_ptr.hpp>


// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. For clarity, we also
// make two typedefs for the mesh type and matrix type. We define
// then the MPI communicator.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;


int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);


    // ---------------------------------------------------------------
    // The next step is to build the mesh. We use here a structured
    // cartesian mesh over the square domain (-1,1)x(-1,1)x(-1,1).
    // The mesh is the partitioned for the parallel computations and
    // the original mesh is deleted.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }

    const UInt Nelements (10);

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type ( Comm ) );

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);

    boost::shared_ptr< mesh_Type > meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
        meshPtr = meshPart.meshPartition();
    }

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We start the discussion with the scalar case, so we define a
    // scalar finite element space and the corresponding matrix.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the scalar ETFESpace ... " << std::flush;
    }

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > scalarSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, &feTetraP1, Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << scalarSpace->dof().numTotalDof() << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Defining the matrix ... " << std::flush;
    }

    if (verbose)
    {
        std::cout << " -- Defining the matrix ... " << std::flush;
    }

    boost::shared_ptr<matrix_Type> scalarMatrix (new matrix_Type ( scalarSpace->map() ) );

    *scalarMatrix *= 0.0;

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    // ---------------------------------------------------------------
    // We can now start the assembly. To understand whether an
    // expression is valid, the critical observation is that every
    // piece of expression is associated to a "fictitious" type (in
    // the sense that it is not the type of the expression in the
    // C++ sense), which matches the mathematical type.
    //
    // For example, in the case of a scalar finite element space, the
    // basis functions are scalar quantities, while their gradients
    // are vectorial quantities.
    //
    // With this in mind, we can now formulate the rules for an
    // expression to be valid:
    //
    // Rule A: The combinaisons between two expressions (through an
    // operator or a function) must be valid. For example, it is not
    // possible to sum a vectorial quantity and a scalar quantity.
    //
    // Rule B: The overall expression must be a scalar quantity. For
    // example, it not possible to integrate simply grad(phi_i)).
    //
    //
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // We can now start the assembly. To understand whether an
    // expression is valid, the critical observation is that every
    // piece of expression is associated to a "fictitious" type (in
    // the sense that it is not the type of the expression in the
    // C++ sense), which matches the mathematical type.
    //
    // For example, in the case of a scalar finite element space, the
    // basis functions are scalar quantities, while their gradients
    // are vectorial quantities.
    //
    // With this in mind, we can now formulate the rules for an
    // expression to be valid:
    //
    // Rule A: The combinaisons between two expressions (through an
    // operator or a function) must be valid. For example, it is not
    // possible to sum a vectorial quantity and a scalar quantity.
    //
    // Rule B: The overall expression must be a scalar quantity. For
    // example, it not possible to integrate simply grad(phi_i)).
    //
    //
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembling the scalar matrix ... " << std::flush;
    }

    {
        using namespace ExpressionAssembly;

        // This would NOT compile (Rule A)
        // Indeed,
        //  grad(phi_i) is a vectorial quantity
        //  grad(phi_j) is a vectorial quantity
        // but the product between two vectors is not
        // defined. If you mean the scalar product, you
        // should use the dot function.

        /*integrate(  elements(scalarSpace->mesh()),
                    quadRuleTetra4pt,
                    scalarSpace,
                    scalarSpace,

                    grad(phi_i) * grad(phi_j)
            )
            >> scalarMatrix;
        */

        // This would NOT compile (Rule B)
        // Indeed,
        //  grad(phi_i) is a vectorial quantity
        //  phi_j is a scalar quantity
        // the product between a scalar quantity and
        // a vectorial quantity is well defined and
        // yield a vectorial quantity.
        // However, the whole expression is a vectorial
        // quantity, therefore, it is not possible to
        // integrate it.

        /*integrate(  elements(scalarSpace->mesh()),
                    quadRuleTetra4pt,
                    scalarSpace,
                    scalarSpace,

                    grad(phi_i) * phi_j
            )
            >> scalarMatrix;
        */


        // This compile because Rule A and Rule B
        // are respected, and even if the expression
        // does not correspond to any real problem.

        VectorSmall<3> V1 (1.0, 0.0, 0.0);

        integrate (  elements (scalarSpace->mesh() ),
                     quadRuleTetra4pt,
                     scalarSpace,
                     scalarSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
                     + 0.0 *
                     (2.0 * phi_i / phi_j
                      - dot ( grad (phi_i), phi_i * V1) *phi_j
                      + 3.1415)
                  )
                >> scalarMatrix;
    }

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // The rule A and B are very simple, usually much more than the
    // compilation errors that can be issued. In case of problem, it
    // is then much easier to look at the expression with the rules A
    // and B to find where is the problem.
    //
    // To finish this tutorial, we compare the norm of the matrix with
    // the norm it should have.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }

    scalarMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    Real matrixNorm ( scalarMatrix->normInf() );

    if (verbose)
    {
        std::cout << " Matrix norm : " << matrixNorm << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    Real matrixNormDiff (std::abs (matrixNorm - 3.2) );

    if (verbose)
    {
        std::cout << " Error : " << matrixNormDiff << std::endl;
    }

    Real testTolerance (1e-10);

    if ( matrixNormDiff < testTolerance )
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );

}


