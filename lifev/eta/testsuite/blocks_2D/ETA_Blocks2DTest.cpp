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
/**
   @file ETA_Blocks2DTest.cpp
   @author L. Pasquale <luca.pasquale@mail.polimi.it>
   @date 2012-11-20
 */

// ===================================================
//! Includes
// ===================================================

#include "ETA_Blocks2DTest.hpp"

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>


// ===================================================
//! Namespaces & define
// ===================================================

// ---------------------------------------------------------------
// We work in the LifeV namespace and define the mesh, matrix and
// vector types that we will need several times.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTriangle> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef MatrixEpetraStructured<Real> blockMatrix_Type;
typedef VectorEpetraStructured blockVector_Type;

// ---------------------------------------------------------------
// We define then a function, which represents a velocity field,
// of magnitude 1 in the y direction (supposing an (x,y)
// reference frame).
// ---------------------------------------------------------------


// ===================================================
//!                   Functions
// ===================================================

// ---------------------------------------------------------------
// We define a function that is supposed to represent a
// force acting on the fluid.
// ---------------------------------------------------------------

Real forceFct ( const Real& /*t*/, const Real& /*x*/, const Real& y, const Real& /*z*/ , const ID& /*i*/)
{
    return y * y;
}

// ===================================================
//!                  Constructors
// ===================================================

ETA_Blocks2DTest::ETA_Blocks2DTest ()
{

#ifdef EPETRA_MPI
    M_comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    M_comm.reset ( new Epetra_SerialComm() );
#endif

}

// ===================================================
//!                      Methods
// ===================================================

std::vector<Real>
ETA_Blocks2DTest::run()
{
    bool verbose (M_comm->MyPID() == 0);
    // ---------------------------------------------------------------
    // We define the mesh and parition it. We use the domain
    // (-1,1)x(-1,1) and a structured mesh.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }

    const UInt Nelements (10);

    std::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh2D ( *fullMeshPtr, 0, Nelements, Nelements, false,
                    2.0,   2.0,
                    -0.0,  -0.0);

    std::shared_ptr< mesh_Type > meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_comm);
        meshPtr = meshPart.meshPartition();
    }

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
    std::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 1, M_comm) );

    std::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 2, 2 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 2, 2 > (meshPtr, &feTriaP2, & (uSpace->fe().geoMap() ), M_comm) );

    std::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 2, 1 > > ETpSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 2, 1 > (meshPtr, &feTriaP1, & (uSpace->fe().geoMap() ), M_comm) );

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

    std::shared_ptr<blockMatrix_Type> ETsystemMatrix (new blockMatrix_Type ( ETuSpace->map() | ETpSpace->map() ) );
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

    std::cout.precision (15);
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

    uSpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (forceFct), fInterpolated, 0.0);

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
    // Finally we return a vector containing both norms
    // ---------------------------------------------------------------

    std::vector<Real> errorNorms (2);
    errorNorms[0] = matrixNorm;
    errorNorms[1] = rhsNorm;

    return errorNorms;

} // run
