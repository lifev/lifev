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
    @brief Tutorial for the use of functor in the ETA framework

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-10-2010

    In this tutorial, we show how the functors can be used
    to add functionalities to the ETA framework.

    We show, in particular, how to assemble a SUPG stabilized
    advection diffusion matrix. We also show how to integrate
    functions without interpolating them.

    Tutorials that should be read before: 1,2,3
 */

// ---------------------------------------------------------------
// We include the same files as in the previous tutorial. No new
// file is required to use the functor functionalities, they are
// already contained in the Integrate.hpp file.
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

#include <lifev/core/fem/FESpace.hpp>

#include <boost/shared_ptr.hpp>

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
// We define here a function for the sake of comparison for the
// right hand side that we will consider.
// ---------------------------------------------------------------

Real fFctRaw ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */, const ID& /*i*/ )
{
    return std::sin (x);
}
function_Type fFct (fFctRaw);

// ---------------------------------------------------------------
// To use the functors, we need to define new classes, which are
// the functors. First of all, we define such a class that is
// equivalent to the fFct function that we have just defined.
// Several elements must be present in this class:
//
// - A type, called return_Type, which represents the type of
//   the values returned by the class. If the evaluation of the
//   function is supposed to be scalar, it is Real, if it is
//   vectorial valued, the type is VectorSmall,... Here, the
//   represented function is scalar, so we define it as Real.
//
// - An operator() function, which takes in argument the wanted
//   type (with the same reasoning as the return_Type) and return
//   a return_Type. Here, we want the functor to return the sinus
//   of the first component of the position (which is a vectorial
//   quantity).
//
// - The copy constructor and the destructor are mandatory as
//   well (we explicit them here).
//
// The example that make in this tutorial are very simple, but a
// functor class can contain data, have much more methods,...
// There is no limitation, excepted the three requirements listed
// above.
// ---------------------------------------------------------------

class fFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<3>& position)
    {
        return std::sin ( position[0] );
    }

    fFunctor() {}
    fFunctor (const fFunctor&) {}
    ~fFunctor() {}
};


// ---------------------------------------------------------------
// We want also to implement a SUPG stabilized advection diffusion
// problem. To weight the stabilization term, we would like to
// have an expression that computes the norm of the advection
// field. As there is no expression already set up for it, we
// devise a functor.
//
// For the sake of exposition, we add a member to the functor
// representing the exponent of the norm.
//
// Following the same rules as before, we define the following
// class.
// ---------------------------------------------------------------

class normFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<3>& advection_field)
    {
        return std::pow ( std::pow ( std::abs (advection_field[0]), M_exponent)
                          + std::pow ( std::abs (advection_field[1]), M_exponent)
                          + std::pow ( std::abs (advection_field[2]), M_exponent)
                          , 1.0 / M_exponent);
    }

    normFunctor (const UInt& expo = 1)
        : M_exponent (expo) {}

    normFunctor (const normFunctor& nf)
        : M_exponent (nf.M_exponent) {}

    ~normFunctor() {}

    void setExponent (const UInt& expo)
    {
        M_exponent = expo;
    }

private:

    UInt M_exponent;
};


// ---------------------------------------------------------------
// As usual, we start by the MPI communicator, the definition of
// the mesh and its partitioning
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
    // We define then the ETFESpace, that we suppose scalar in this
    // case. We also need a standard FESpace to perform the
    // interpolation.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the ETFESpace ... " << std::flush;
    }

    std::string uOrder ("P1");
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace ( new FESpace< mesh_Type, MapEpetra > (meshPart, uOrder, 1, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETuSpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << ETuSpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // We start with the computation of the right hand side. First of
    // all, we reuse the procedure that we used in the tutorial 3,
    // i.e., we interpolate the function fFct before integration.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the rhs via interpolation ... " << std::flush;
    }

    vector_Type fInterpolated (ETuSpace->map(), Repeated);

    fInterpolated *= 0.0;

    uSpace->interpolate (fFct, fInterpolated, 0.0);

    vector_Type rhsET (ETuSpace->map(), Repeated);
    rhsET *= 0.0;

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr() ,
                    ETuSpace,
                    value (ETuSpace, fInterpolated) *phi_i
                  )
                >> rhsET;
    }

    rhsET.globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    // ---------------------------------------------------------------
    // We assemble now a similar right hand side but without using the
    // interpolation step, i.e. with a L2 projection.
    //
    // To this aim, we use the exprssion "X" which represents the
    // position (i.e. the position of the quadrature node).
    //
    // We also use the functor, with a shared pointer (to avoid copy
    // of functors that would contain large data). The eval keyword
    // is then used to evaluate the functor.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the rhs via functor ... " << std::flush;
    }

    vector_Type rhsETfunctor (ETuSpace->map(), Repeated);
    rhsETfunctor *= 0.0;

    {
        using namespace ExpressionAssembly;

        boost::shared_ptr<fFunctor> ffunctor (new fFunctor);

        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr() ,
                    ETuSpace,
                    eval (ffunctor, X) *phi_i
                  )
                >> rhsETfunctor;
    }

    rhsETfunctor.globalAssemble();

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We take now a look at the SUPG stabilization. The term that we
    // want to add is gamma (beta dot grad v)(beta dot grad u) where
    // gamma is h/2|beta|.
    //
    // To get the diameter of the element considered (h), the keyword
    // h_K is used.
    //
    // For the norm of the advection field (|beta|), we simply
    // evaluate the norm functor with beta as argument. Remark that
    // here, the norm of the advection field is known, so the functor
    // is not strictly necessary. However, in other cases, it can be
    // unknown, so the functor strategy must be used.
    //
    // Here is the resulting code:
    // ---------------------------------------------------------------


    if (verbose)
    {
        std::cout << " -- Assembly of the AD-SUPG ... " << std::flush;
    }

    boost::shared_ptr<matrix_Type> ETsystemMatrix (new matrix_Type ( ETuSpace->map() ) );
    *ETsystemMatrix *= 0.0;

    {
        using namespace ExpressionAssembly;

        VectorSmall<3> beta (0, 1.0, 2.0);

        boost::shared_ptr<normFunctor> norm ( new normFunctor);

        norm->setExponent (2);

        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr(),
                    ETuSpace,
                    ETuSpace,

                    dot ( grad (phi_i) , grad (phi_j) )
                    + dot ( grad (phi_j) , beta ) *phi_i
                    + h_K * dot (grad (phi_i), beta) *dot (grad (phi_j), beta) / (2.0 * eval (norm, beta) )

                  )
                >> ETsystemMatrix;
    }

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // Finally, we compute the difference between the two right hand
    // sides computed previously. The difference should not be too
    // large but is not supposed to be zero (because the interpolation
    // and the L2 projection are not the same operations).
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the error ... " << std::flush;
    }

    vector_Type checkRhs (ETuSpace->map(), Repeated);
    checkRhs *= 0.0;

    checkRhs += rhsET;
    checkRhs -= rhsETfunctor;

    checkRhs.globalAssemble();

    vector_Type checkRhsUnique (checkRhs, Unique);

    Real diffRhsNorm ( checkRhsUnique.normInf() );

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Rhs diff : " << diffRhsNorm << std::endl;
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

    Real testTolerance (1e-3);

    if (diffRhsNorm < testTolerance)
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );
}


