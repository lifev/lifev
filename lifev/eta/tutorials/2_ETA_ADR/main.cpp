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
    @brief Tutorial for the assembly of the ADR problem and comparison with classical assembly

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-10-2010

    This tutorial aims at presenting further expressions and compare the performances with
    the ones obtained with the classical assembly.

    Tutorials that should be read before: 1

 */

// ---------------------------------------------------------------
// We start by including the same files as in the tutorial 1, with
// the two "eta" files being specific to the ET assembly. We also
// need to use vectors so the VectorEpetra file is included as
// well.
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

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <boost/shared_ptr.hpp>


// ---------------------------------------------------------------
// We include two additional files which are useful to make the
// assembly using the "classical way". We also include a file to
// monitor the different timings.
// ---------------------------------------------------------------

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>

#include <lifev/core/util/LifeChrono.hpp>


// ---------------------------------------------------------------
// We work in the LifeV namespace and define the mesh, matrix and
// vector types that we will need several times.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef FESpace<mesh_Type, MapEpetra>::function_Type function_Type;

// ---------------------------------------------------------------
// We define then a function, which represents a velocity field,
// of magnitude 1 in the y direction (supposing an (x,y,z)
// reference frame).
// ---------------------------------------------------------------


Real betaFctRaw ( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& i )
{
    if (i == 1)
    {
        return 1;
    }
    return 0;
}
function_Type betaFct (betaFctRaw);



// ---------------------------------------------------------------
// As in the tutorial 1, we start with the definition of the
// MPI communicator and boolean for the outputs.
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


    // ---------------------------------------------------------------
    // We define the mesh and parition it. We use again the domain
    // (-1,1)x(-1,1)x(-1,1) and a structured mesh.
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
    // We start by defining the finite element spaces. We use the type
    // FESpace for the classical way and copy it in an ETFESpace for
    // the ET assembly.
    //
    // We also build similar spaces for the advection field. This
    // space does not need to be the same as the solution space, even
    // if this is the case here.
    //
    // We remark here two details:
    // 1. The spaces for the advection (betaSpace and ETbetaSpace) are
    //    vectorial.
    // 2. The constructor for the ETFESpace structures use an
    //    additional arguement, the geometric mapping. In the
    //    tutorial 1, this argument was omitted, so the geometric
    //    mapping was guessed from the mesh type.
    //
    // In the end, both solution spaces display their respective
    // number of degrees of freedom, which must be the same.
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // We start by defining the finite element spaces. We use the type
    // FESpace for the classical way and copy it in an ETFESpace for
    // the ET assembly.
    //
    // We also build similar spaces for the advection field. This
    // space does not need to be the same as the solution space, even
    // if this is the case here.
    //
    // We remark here two details:
    // 1. The spaces for the advection (betaSpace and ETbetaSpace) are
    //    vectorial.
    // 2. The constructor for the ETFESpace structures use an
    //    additional arguement, the geometric mapping. In the
    //    tutorial 1, this argument was omitted, so the geometric
    //    mapping was guessed from the mesh type.
    //
    // In the end, both solution spaces display their respective
    // number of degrees of freedom, which must be the same.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building FESpaces ... " << std::flush;
    }

    std::string uOrder ("P1");
    std::string bOrder ("P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 1, Comm) );

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > betaSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, bOrder, 3, Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << uSpace->dof().numTotalDof() << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Building ETFESpaces ... " << std::flush;
    }

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETbetaSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPtr, & (betaSpace->refFE() ), & (betaSpace->fe().geoMap() ), Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << ETuSpace->dof().numTotalDof() << std::endl;
    }

    // ---------------------------------------------------------------
    // We interpolate then the advection function of the mesh at hand.
    // This is performed with the classical FESpace only.
    //
    // Indeed, the interpolation has not yet been implemented for the
    // ETFESpace and vector of values for the FESpace and ETFESpace
    // are fully compatible (they use the same degrees of freedom
    // numbering). Therefore, they can be exchanged at will. This is
    // very important since some features have been implemented only
    // for regular FESpace but not yet for ETFESpace (e.g. output
    // functionalities).
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // We interpolate then the advection function of the mesh at hand.
    // This is performed with the classical FESpace only.
    //
    // Indeed, the interpolation has not yet been implemented for the
    // ETFESpace and vector of values for the FESpace and ETFESpace
    // are fully compatible (they use the same degrees of freedom
    // numbering). Therefore, they can be exchanged at will. This is
    // very important since some features have been implemented only
    // for regular FESpace but not yet for ETFESpace (e.g. output
    // functionalities).
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Interpolating the advection field ... " << std::flush;
    }

    vector_Type beta (betaSpace->map(), Repeated);
    betaSpace->interpolate (betaFct, beta, 0.0);

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We build now two matrices, one for each assembly procedure in
    // order to compare them in the end. There is no difference in the
    // construction of the matrices.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Defining the matrices ... " << std::flush;
    }

    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uSpace->map() ) );
    *systemMatrix *= 0.0;

    boost::shared_ptr<matrix_Type> ETsystemMatrix (new matrix_Type ( ETuSpace->map() ) );
    *ETsystemMatrix *= 0.0;

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We come now to the assembly itself. We start with the classical
    // way, which consists in instentiating an ADRAssembler and call
    // the different methods from that class.
    //
    // Using a LifeChrono, we monitor the time required for the
    // assembly.
    //
    // The terms assembled correspond, in the order of appearance, to
    // - a laplacian term
    // - an advection term
    // - a reaction term (with coefficient 2.0), aka mass term.
    // ---------------------------------------------------------------

    LifeChrono StdChrono;
    StdChrono.start();

    if (verbose)
    {
        std::cout << " -- Classical assembly ... " << std::flush;
    }

    ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;

    adrAssembler.setup (uSpace, betaSpace);

    adrAssembler.addDiffusion (systemMatrix, 1.0);

    adrAssembler.addAdvection (systemMatrix, beta);

    adrAssembler.addMass (systemMatrix, 2.0);

    StdChrono.stop();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time : " << StdChrono.diff() << std::endl;
    }


    // ---------------------------------------------------------------
    // We perform now the same assembly with the ET assembly and still
    // monitor the timings required.
    //
    // As in tutorial 1, we need to use the special namespace. The
    // arguments of the integrate function still have the same
    // meaning, but the expression is not a bit different.
    // Remark that to ensure a fair comparison, we use the quadrature
    // rule used for the classical way (stored in the uSpace).
    //
    // One of the differences between the two assembly is that with
    // the ET way, the weak formulation is immediately visible, while
    // it is not with the classical way.
    //
    // For the advective term, we remark that a new expression is used
    // for the interpolation of the velocity field. The value
    // function is used with, as first argument the ETFESpace in
    // which the velocity is given and in second argument the vector
    // of the values.
    //
    // Remark also that the Real constants and VectorSmall constants
    // can be used directly in the expressions to integrate.
    // ---------------------------------------------------------------

    LifeChrono ETChrono;
    ETChrono.start();

    if (verbose)
    {
        std::cout << " -- ET assembly ... " << std::flush;
    }


    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr(),
                    ETuSpace,
                    ETuSpace,

                    dot ( grad (phi_i) , grad (phi_j) )
                    + dot ( grad (phi_j) , value (ETbetaSpace, beta) ) *phi_i
                    + 2.0 * phi_i * phi_j

                  )
                >> ETsystemMatrix;
    }

    ETChrono.stop();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time : " << ETChrono.diff() << std::endl;
    }


    // ---------------------------------------------------------------
    // The timings displayed should indicate that the ET way is faster
    // than the classical way. Indeed, the classical way loops over
    // the elements for each term added (3 times here), assembles
    // different local contributions for each term and adds them for
    // each term. With the ET way, only one loop over the elements
    // is required, computations are reused and only one local
    // contribution is computed and added to the global matrix.
    //
    // In general, the longer is the expression, the better is the
    // performance of the ET with respect to the classical way.
    //
    // We finally need to check that both yield the same matrix. In
    // that aim, we need to finalize both matrices.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Closing the matrices ... " << std::flush;
    }

    systemMatrix->globalAssemble();
    ETsystemMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We compute now the matrix of the difference and finally the
    // norm of the difference. This should be very low if the two
    // matrices are identical.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the error ... " << std::flush;
    }

    boost::shared_ptr<matrix_Type> checkMatrix (new matrix_Type ( ETuSpace->map() ) );
    *checkMatrix *= 0.0;

    *checkMatrix += *systemMatrix;
    *checkMatrix += (*ETsystemMatrix) * (-1);

    checkMatrix->globalAssemble();

    Real errorNorm ( checkMatrix->normInf() );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We finalize the MPI if needed.
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    // ---------------------------------------------------------------
    // We finally display the error norm and compare with the
    // tolerance of the test.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " Matrix Error : " << errorNorm << std::endl;
    }

    Real testTolerance (1e-10);

    if (errorNorm < testTolerance)
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );
}


