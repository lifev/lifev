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
    @brief Tutorial for the use of block structures with the ETA framework

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 02-07-2012

    In this tutorial, we explain how to assemble the linear systems associated
    with problems with several physical quantities. This is of particular
    interest for multiphysics problems.

    In this tutorial, we assemble a Stokes problem and a convection problem
    using the blocks.

    Tutorials that should be read before: 1,4
 */

// ---------------------------------------------------------------
// We still use the same files, see the previous tutorials for
// explanations of the different utilities.
//
// The only new files are MatrixEpetraStructured.hpp and
// VectorEpetraStructured.hpp, that provide the block versions
// of MatrixEpetra and VectorEpetra.
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
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <boost/shared_ptr.hpp>

// ---------------------------------------------------------------
// We work in the LifeV namespace and define the mesh, matrix and
// vector types that we will need several times. We add new
// typedefs for the block matrices and block vectors.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef MatrixEpetraStructured<Real> blockMatrix_Type;
typedef VectorEpetraStructured blockVector_Type;
typedef FESpace<mesh_Type, MapEpetra>::function_Type function_Type;

// ---------------------------------------------------------------
// We also define a function that is supposed to represent a
// force acting on the fluid.
// ---------------------------------------------------------------

Real forceFctRaw ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z , const ID& /*i*/)
{
    return z * z;
}
function_Type forceFct (forceFctRaw);

// ---------------------------------------------------------------
// As usual, we start by the MPI communicator, the definition of
// the mesh and its partitioning.
// ---------------------------------------------------------------

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);


    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }

    const UInt Nelements (10);

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type ( Comm ) );

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);

    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We define now the ETFESpaces. We need one space for the
    // velocity (vectorial, P2) and one space for the pressure
    // (scalar, P1).
    //
    // We also define velocity FESpace to use the interpolation.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the spaces ... " << std::flush;
    }

    std::string uOrder ("P2");
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPart, uOrder, 1, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPart, &feTetraP2, & (uSpace->fe().geoMap() ), Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETpSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, & (uSpace->fe().geoMap() ), Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Velocity dofs: " << ETuSpace->dof().numTotalDof() << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Pressure dofs: " << ETpSpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // We want to assemble a Stokes matrix. This matrix is divided
    // into four blocks (with the (1,1)-block empty). We define it
    // using a special constructor with "|" operators, which
    // separates the different blocks.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembly of the Stokes matrix ... " << std::flush;
    }

    boost::shared_ptr<blockMatrix_Type> ETsystemMatrix (new blockMatrix_Type ( ETuSpace->map() | ETpSpace->map() ) );
    *ETsystemMatrix *= 0.0;


    // ---------------------------------------------------------------
    // To fill them, we simply have to call one integrate for each
    // non-empty block and indicate in which block the contribution
    // has to go.
    //
    // Remark that, depending on the block assembled, the spaces for
    // the test and trial functions are different.
    // ---------------------------------------------------------------

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuSpace,
                    ETuSpace,

                    dot ( grad (phi_i) , grad (phi_j) )

                  )
                >> ETsystemMatrix->block (0, 0);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuSpace,
                    ETpSpace,

                    phi_j * div (phi_i)

                  )
                >> ETsystemMatrix->block (0, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuSpace,

                    phi_i * div (phi_j)

                  )
                >> ETsystemMatrix->block (1, 0);
    }

    ETsystemMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // For testing purposes, we also compute the norm of the matrix.
    // ---------------------------------------------------------------

    Real matrixNorm (ETsystemMatrix->normInf() );

    if (verbose)
    {
        std::cout << " Matrix norm " << matrixNorm << std::endl;
    }


    // ---------------------------------------------------------------
    // In a very similar way, we can assemble the block right hand
    // side. Here, we suppose that a force acts on  the fluid, so
    // only the block associated to the velocity is not empty.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the rhs ... " << std::flush;
    }

    vector_Type fInterpolated (ETuSpace->map(), Repeated);

    fInterpolated *= 0.0;

    uSpace->interpolate (forceFct, fInterpolated, 0.0);

    blockVector_Type ETrhs (ETuSpace->map() | ETpSpace->map() , Repeated);
    ETrhs *= 0.0;

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuSpace,
                    dot (value (ETuSpace, fInterpolated), phi_i)
                  )
                >> ETrhs.block (0);
    }

    ETrhs.globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // For testing purposes, we also compute the norm of the rhs. To
    // take the norm, we need a unique vector (globalAssemble is not
    // sufficient).
    // ---------------------------------------------------------------

    blockVector_Type ETrhsUnique (ETrhs, Unique);

    Real rhsNorm (ETrhsUnique.normInf() );

    if (verbose)
    {
        std::cout << " Rhs norm " << rhsNorm << std::endl;
    }


    // ---------------------------------------------------------------
    // We finalize the MPI if needed.
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    // ---------------------------------------------------------------
    // We finally compare the difference with the tolerance of the
    // test.
    // ---------------------------------------------------------------

    if (
        ( std::abs (matrixNorm - 3.856 ) < 1e-2) &&
        ( std::abs (rhsNorm - 0.00187384) < 1e-5)
    )
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );
}


