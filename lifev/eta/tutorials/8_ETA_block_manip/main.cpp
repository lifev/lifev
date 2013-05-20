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
    @brief Tutorial for a "better" usage of the blocks.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 02-07-2012

    In the previous tutorial, we saw how to use the ETA framework with
    block structured matrices. This tutorial aims at investigating
    that more in detail.

    In particular, 4 different strategies are proposed for the assembly
    of the Stokes matrix. The first one is the most simple strategy, but
    not the most efficient one. Two intermediate strategies are developed
    to show the influence of the different points. The last strategy
    represents the most efficient way of assembling the matrix.

    In practice, only the first (simplest) and last (fastest) strategy
    should be considered since intermediate strategies have no advantages.

    Tutorials that should be read before: 1,4,7

 */

// ---------------------------------------------------------------
// We reuse the same files as in the previous tutorial. We add the
// chrono to measure the timings and the file
// MatrixEpetraStructuredUtility to use the block copy features.
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
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/util/LifeChrono.hpp>

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
    // We define now the ETFESpaces. We need one space for the
    // velocity (vectorial, P2) and one space for the pressure
    // (scalar, P1).
    //
    // For reasons explained hereafter, we also need a scalar
    // ETFESpace for the velocity.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the spaces ... " << std::flush;
    }

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPtr, &feTetraP2, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETpSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, &feTetraP1, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETuCompSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, &feTetraP2, Comm) );

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
    // We recall here the implementation that we used for the
    // Stokes system in the previous tutorial. It consists in defining
    // a 2x2 structure on the matrix and assemble the different blocks
    // separately.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembly of the Stokes matrix (I) ... " << std::flush;
    }
    LifeChrono chronoI;
    chronoI.start();

    boost::shared_ptr<blockMatrix_Type> ETsystemMatrixI (new blockMatrix_Type ( ETuSpace->map() | ETpSpace->map() ) );
    *ETsystemMatrixI *= 0.0;

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuSpace,
                    ETuSpace,

                    dot ( grad (phi_i) , grad (phi_j) )

                  )
                >> ETsystemMatrixI->block (0, 0);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuSpace,
                    ETpSpace,

                    phi_j * div (phi_i)

                  )
                >> ETsystemMatrixI->block (0, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuSpace,

                    phi_i * div (phi_j)

                  )
                >> ETsystemMatrixI->block (1, 0);
    }

    chronoI.stop();

    ETsystemMatrixI->globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time: " << chronoI.diff() << std::endl;
    }


    // ---------------------------------------------------------------
    // This implementation is however "suboptimal". Indeed, the
    // (0,0)-block is itself block diagonal, but this is not taken
    // into account. So, many zeros are added to the matrix, which
    // might result in both a slower assembly and a longer time
    // required to compute the preconditioner.
    //
    // To avoid this bottlneck, a first option is to define the
    // structure to be 4x4 (3 velocity components and pressure).
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembly of the Stokes matrix (II) ... " << std::flush;
    }
    LifeChrono chronoII;
    chronoII.start();

    boost::shared_ptr<blockMatrix_Type> ETsystemMatrixII
    (new blockMatrix_Type ( ETuCompSpace->map() | ETuCompSpace->map() | ETuCompSpace->map() | ETpSpace->map() ) );
    *ETsystemMatrixII *= 0.0;

    {
        using namespace ExpressionAssembly;

        VectorSmall<3> e1 (1, 0, 0);
        VectorSmall<3> e2 (0, 1, 0);
        VectorSmall<3> e3 (0, 0, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETuCompSpace,

                    dot ( grad (phi_i) , grad (phi_j) )

                  )
                >> ETsystemMatrixII->block (0, 0);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETuCompSpace,

                    dot ( grad (phi_i) , grad (phi_j) )

                  )
                >> ETsystemMatrixII->block (1, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETuCompSpace,

                    dot ( grad (phi_i) , grad (phi_j) )

                  )
                >> ETsystemMatrixII->block (2, 2);



        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETpSpace,

                    phi_j * dot (grad (phi_i), e1)

                  )
                >> ETsystemMatrixII->block (0, 3);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETpSpace,

                    phi_j * dot (grad (phi_i), e2)

                  )
                >> ETsystemMatrixII->block (1, 3);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETpSpace,

                    phi_j * dot (grad (phi_i), e3)

                  )
                >> ETsystemMatrixII->block (2, 3);


        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuCompSpace,

                    phi_i * dot (grad (phi_j), e1)

                  )
                >> ETsystemMatrixII->block (3, 0);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuCompSpace,

                    phi_i * dot (grad (phi_j), e2)

                  )
                >> ETsystemMatrixII->block (3, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuCompSpace,

                    phi_i * dot (grad (phi_j), e3)

                  )
                >> ETsystemMatrixII->block (3, 2);
    }

    chronoII.stop();

    ETsystemMatrixII->globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time: " << chronoII.diff() << std::endl;
    }


    // ---------------------------------------------------------------
    // The difference in timings should be quite large between the
    // two first implementations. However, we can still remark that
    // the blocks (0,0), (1,1) and (2,2) (in the latest
    // implementation) are exact similar. So, instead of recomputing
    // this block 3 times, we can assemble it simply once and then
    // copy it where needed.
    //
    // The following code realizes this idea.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembly of the Stokes matrix (III) ... " << std::flush;
    }
    LifeChrono chronoIII;
    chronoIII.start();


    // ---------------------------------------------------------------
    // First, we build the small matrix, assemble it and close it.
    // ---------------------------------------------------------------

    boost::shared_ptr<blockMatrix_Type> ETcomponentMatrixIII
    (new blockMatrix_Type ( ETuCompSpace->map() ) );
    *ETcomponentMatrixIII *= 0.0;

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETuCompSpace,

                    dot ( grad (phi_i) , grad (phi_j) )

                  )
                >> ETcomponentMatrixIII->block (0, 0);
    }

    ETcomponentMatrixIII->globalAssemble();


    // ---------------------------------------------------------------
    // We define then the large matrix.
    // ---------------------------------------------------------------

    boost::shared_ptr<blockMatrix_Type> ETsystemMatrixIII
    (new blockMatrix_Type ( ETuCompSpace->map() | ETuCompSpace->map() | ETuCompSpace->map() | ETpSpace->map() ) );
    *ETsystemMatrixIII *= 0.0;


    // ---------------------------------------------------------------
    // We copy the three blocks
    // ---------------------------------------------------------------

    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMatrixIII->block (0, 0), *ETsystemMatrixIII->block (0, 0) );
    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMatrixIII->block (0, 0), *ETsystemMatrixIII->block (1, 1) );
    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMatrixIII->block (0, 0), *ETsystemMatrixIII->block (2, 2) );

    // ---------------------------------------------------------------
    // Assemble the remaing blocks
    // ---------------------------------------------------------------

    {
        using namespace ExpressionAssembly;

        VectorSmall<3> e1 (1, 0, 0);
        VectorSmall<3> e2 (0, 1, 0);
        VectorSmall<3> e3 (0, 0, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETpSpace,

                    phi_j * dot (grad (phi_i), e1)

                  )
                >> ETsystemMatrixIII->block (0, 3);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETpSpace,

                    phi_j * dot (grad (phi_i), e2)

                  )
                >> ETsystemMatrixIII->block (1, 3);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETpSpace,

                    phi_j * dot (grad (phi_i), e3)

                  )
                >> ETsystemMatrixIII->block (2, 3);


        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuCompSpace,

                    phi_i * dot (grad (phi_j), e1)

                  )
                >> ETsystemMatrixIII->block (3, 0);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuCompSpace,

                    phi_i * dot (grad (phi_j), e2)

                  )
                >> ETsystemMatrixIII->block (3, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuCompSpace,

                    phi_i * dot (grad (phi_j), e3)

                  )
                >> ETsystemMatrixIII->block (3, 2);
    }

    chronoIII.stop();

    ETsystemMatrixIII->globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time: " << chronoIII.diff() << std::endl;
    }


    // ---------------------------------------------------------------
    // The variant three is usually a bit faster than the variant II,
    // which shows that some more time can be gained.
    //
    // For the final variant, we exploit further the block structure
    // by changing the structure of the matrix on the fly. The
    // strategy consists in assembly the velocity-velocity block
    // as in the variant III (by copy of a small block), while
    // for the remain parts of the matrix, the structure is changed
    // and only two blocks are assembled.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembly of the Stokes matrix (IV) ... " << std::flush;
    }
    LifeChrono chronoIV;
    chronoIV.start();


    // ---------------------------------------------------------------
    // First, we build the small matrix, assemble it and close it.
    // ---------------------------------------------------------------

    boost::shared_ptr<blockMatrix_Type> ETcomponentMatrixIV
    (new blockMatrix_Type ( ETuCompSpace->map() ) );
    *ETcomponentMatrixIV *= 0.0;

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuCompSpace,
                    ETuCompSpace,

                    dot ( grad (phi_i) , grad (phi_j) )

                  )
                >> ETcomponentMatrixIV->block (0, 0);
    }

    ETcomponentMatrixIV->globalAssemble();


    // ---------------------------------------------------------------
    // We define then the large matrix.
    // ---------------------------------------------------------------

    boost::shared_ptr<blockMatrix_Type> ETsystemMatrixIV
    (new blockMatrix_Type ( ETuCompSpace->map() | ETuCompSpace->map() | ETuCompSpace->map() | ETpSpace->map() ) );
    *ETsystemMatrixIV *= 0.0;


    // ---------------------------------------------------------------
    // We copy the three blocks
    // ---------------------------------------------------------------

    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMatrixIV->block (0, 0), *ETsystemMatrixIV->block (0, 0) );
    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMatrixIV->block (0, 0), *ETsystemMatrixIV->block (1, 1) );
    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMatrixIV->block (0, 0), *ETsystemMatrixIV->block (2, 2) );

    // ---------------------------------------------------------------
    // Now change the structure of the matrix back to the 2x2
    // structure and assemble to two remaining parts.
    // ---------------------------------------------------------------

    ETsystemMatrixIV->setBlockStructure (ETuSpace->map() | ETpSpace->map() );

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETuSpace,
                    ETpSpace,

                    phi_j * div (phi_i)

                  )
                >> ETsystemMatrixIV->block (0, 1);

        integrate ( elements (ETuSpace->mesh() ),
                    quadRuleTetra4pt,
                    ETpSpace,
                    ETuSpace,

                    phi_i * div (phi_j)

                  )
                >> ETsystemMatrixIV->block (1, 0);
    }

    chronoIV.stop();

    ETsystemMatrixIV->globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time: " << chronoIV.diff() << std::endl;
    }


    // ---------------------------------------------------------------
    // We finally compare the different matrices obtained.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the error ... " << std::flush;
    }

    boost::shared_ptr<matrix_Type> checkMatrixIvsII (new matrix_Type ( ETuSpace->map() + ETpSpace->map() ) );
    *checkMatrixIvsII *= 0.0;

    *checkMatrixIvsII += *ETsystemMatrixI;
    *checkMatrixIvsII += (*ETsystemMatrixII) * (-1);

    checkMatrixIvsII->globalAssemble();

    Real errorNormIvsII ( checkMatrixIvsII->normInf() );


    boost::shared_ptr<matrix_Type> checkMatrixIvsIII (new matrix_Type ( ETuSpace->map() + ETpSpace->map() ) );
    *checkMatrixIvsIII *= 0.0;

    *checkMatrixIvsIII += *ETsystemMatrixI;
    *checkMatrixIvsIII += (*ETsystemMatrixIII) * (-1);

    checkMatrixIvsIII->globalAssemble();

    Real errorNormIvsIII ( checkMatrixIvsIII->normInf() );


    boost::shared_ptr<matrix_Type> checkMatrixIvsIV (new matrix_Type ( ETuSpace->map() + ETpSpace->map() ) );
    *checkMatrixIvsIV *= 0.0;

    *checkMatrixIvsIV += *ETsystemMatrixI;
    *checkMatrixIvsIV += (*ETsystemMatrixIV) * (-1);

    checkMatrixIvsIV->globalAssemble();

    Real errorNormIvsIV ( checkMatrixIvsIV->normInf() );



    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " I vs II : " << errorNormIvsII << std::endl;
    }
    if (verbose)
    {
        std::cout << " I vs III : " << errorNormIvsIII << std::endl;
    }
    if (verbose)
    {
        std::cout << " I vs IV : " << errorNormIvsIV << std::endl;
    }


#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    Real tolerance (1e-10);

    if ( (errorNormIvsII < tolerance)
            && (errorNormIvsIII < tolerance)
            && (errorNormIvsIV < tolerance)
       )
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );
}


