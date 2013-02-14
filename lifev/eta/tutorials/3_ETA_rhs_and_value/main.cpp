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
    @brief Tutorial for the assembly of a right hand side and the computation of
    benchmark values.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 29-06-2012

    In this tutorial, we show how to assemble a right hand side and compute a
    volume integral, using the ETA framework.

    Tutorials that should be read before: 1,2

 */

// ---------------------------------------------------------------
// We use the same files as for the tutorial 2, nothing else is
// required to compute a right hand side or a value.
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

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>


// ---------------------------------------------------------------
// All the work is done in the namespace LifeV. We also define
// the mesh type and the vector type. The matrix type is also
// required for the classical way of assembling.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef VectorEpetra vector_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef FESpace<mesh_Type, MapEpetra>::function_Type function_Type;


// ---------------------------------------------------------------
// The assembly of the right hand side by the ADRAssembler
// requires the definition of a free function. We take it here as
// a constant in space.
//
// We will also integrate a function over the whole domain. We
// chose the function g(x,y,z)=x so that an exact value for the
// integral is known.
// ---------------------------------------------------------------

Real fRhsRaw ( const Real& /* t */, const Real& /*x*/, const Real& /*y*/, const Real& /* z */ , const ID& /* i */ )
{
    return  2;
}
function_Type fRhs (fRhsRaw);

Real gRaw ( const Real& /* t */, const Real& x, const Real& /*y*/, const Real& /* z */ , const ID& /* i */ )
{
    return  x;
}
function_Type g (gRaw);


// ---------------------------------------------------------------
// As in the previous tutorials, we start the execution by
// constructing the MPI communicator.
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
    // We define the mesh: the domain is (-1,1)x(-1,1)x(-1,1) and the
    // mesh structured. We also partition the mesh.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }

    const UInt Nelements (10);

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

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
    // As in tutorial 2, we define now both FESpace and ETFESpace, for
    // the sake of comparison.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the FESpace ... " << std::flush;
    }

    std::string uOrder ("P1");
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace ( new FESpace< mesh_Type, MapEpetra > (meshPart, uOrder, 1, Comm) );

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
        std::cout << " -- Building the ETFESpace ... " << std::flush;
    }

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
    // We proceed now with the assembly of the right hand sides. For
    // the classical way, we need to interpolate the function and then
    // call the corresponding method in the ADRAssembler class.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Classical way ... " << std::flush;
    }

    vector_Type fInterpolated (uSpace->map(), Repeated);

    fInterpolated = 0.0;

    uSpace->interpolate (fRhs, fInterpolated, 0.0);

    vector_Type rhs (uSpace->map(), Repeated);

    rhs = 0.0;

    ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;

    adrAssembler.setup (uSpace, uSpace);

    adrAssembler.addMassRhs (rhs, fInterpolated);

    rhs.globalAssemble();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // Now we perform the same assembly with the ETA framework.
    //
    // The assembly looks very similar to those seen in the previous
    // tutorials. The main differences are:
    //
    // - There is one less argument in the integrate function.
    //   Indeed, the trial space (solution space) does not make sense
    //   when assembling the right hand side.
    //
    // - The expression to integrate is different. Remark that using
    //   the keyword phi_j in the expression to integrate would result
    //   in a compilation error. Indeed, phi_j refers to the trial
    //   space, which is not defined here.
    //
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- ETA way ... " << std::flush;
    }

    vector_Type rhsET (ETuSpace->map(), Repeated);
    rhsET = 0.0;

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
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We remark that in this particular case, since the function
    // is constant, we could directly use the value of the constant
    // in the expression, thus avoiding the interpolation step.
    //
    // The following code would then yield the same results.
    // ---------------------------------------------------------------
    /*
        if (verbose) std::cout << " -- ETA way ... " << std::flush;

        vector_Type rhsET(ETuSpace->map(),Repeated);
        rhsET=0.0;

        {
            using namespace ExpressionAssembly;

            integrate( elements(ETuSpace->mesh()),
                       uSpace->qr() ,
                       ETuSpace,
                       2.0*phi_i
                     )
                >> rhsET;
        }

        rhsET.globalAssemble();

        if (verbose) std::cout << " done ! " << std::endl;
    */

    // ---------------------------------------------------------------
    // Before checking that the two right hand sides that we assembled
    // are legal, we also integrate values on the domain.
    //
    // Integrating values is again slightly different than the other
    // integrations.
    //
    // We start by integrating the function g on the
    // domain. In this tutorial, we interpolate the function g to
    // integrate it. Another way would be to pass through functors
    // (which is explained in another tutorial).
    //
    // The result of the integration is put in a simple scalar value.
    //
    // The integrate function is now another argument less: no space
    // for the test function (nor for the trial functions) is given
    // in argument. The keyword phi_i (and the keyword phi_j) are
    // illegal in the expression and using it would result in a
    // compilation error.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Integrating g ... " << std::flush;
    }

    vector_Type gInterpolated (uSpace->map(), Repeated);

    gInterpolated = 0.0;

    uSpace->interpolate (g, gInterpolated, 0.0);

    Real localIntegralg (0.0);

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr() ,
                    value (ETuSpace, gInterpolated)
                  )
                >> localIntegralg;
    }


    // ---------------------------------------------------------------
    // When assembling a matrix or a vector, the computed local
    // contributions are communicated, if necessary, to the other
    // processes when the globalAssemble method is called. There is
    // no such method for scalar constants.

    // The communication must therefore be performed "by hand",
    // through the MPI SumAll instruction, which sums the
    // contributions of all the processes.
    //
    // Before performing the sum, when use a Barrier, which ensures
    // that all the processes have computed their contributions.
    // ---------------------------------------------------------------

    Real globalIntegralg (0.0);

    Comm->Barrier();
    Comm->SumAll (&localIntegralg, &globalIntegralg, 1);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We present finally a small example which corresponds to an
    // exception with respect to what has been said before: we stated
    // that the scalar constants can be used without adaptation in the
    // expressions. This is true as long as the expression does not
    // reduce to a scalar constant.
    //
    // For example, if we want to compute the volume of the domain,
    // which is the integral of 1 over it, we need to use the function
    // value (otherwise a compilation error is issued). Remark that
    // the value function can also be used for the constants in
    // all the expressions, but it is then only facultative.
    //
    // Again, to get the global volume, we need to sum the different
    // contributions.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the volume ... " << std::flush;
    }

    Real localVolume (0.0);

    {
        using namespace ExpressionAssembly;

        // This would create a compilation error
        /*
        integrate( elements(ETuSpace->mesh()),
                   uSpace->qr() ,
                   1.0
                 )
            >> localVolume;
        */

        // This is the right synthax
        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr() ,
                    value (1.0)
                  )
                >> localVolume;
    }

    Real globalVolume (0.0);

    Comm->Barrier();
    Comm->SumAll (&localVolume, &globalVolume, 1);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // Finally, we check that the two right hand sides are similar
    // with the two assembly procedures. We also check that the two
    // integrals correspond to their exact values.
    // ---------------------------------------------------------------

    vector_Type checkRhs (ETuSpace->map(), Repeated);
    checkRhs = 0.0;

    checkRhs += rhs;
    checkRhs -= rhsET;

    checkRhs.globalAssemble();

    vector_Type checkRhsUnique (checkRhs, Unique);

    Real errorRhsNorm ( checkRhsUnique.normInf() );

    Real errorIntegralg (std::abs ( globalIntegralg - 0.0) );

    Real errorVolume (std::abs ( globalVolume - 8.0) );

    if (verbose)
    {
        std::cout << " Rhs error : " << errorRhsNorm << std::endl;
    }
    if (verbose)
    {
        std::cout << " Integral error : " << errorIntegralg << std::endl;
    }
    if (verbose)
    {
        std::cout << " Volume error : " << errorVolume << std::endl;
    }


    // ---------------------------------------------------------------
    // As usual, finalize the MPI if needed.
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    // ---------------------------------------------------------------
    // Finally, we test the error with respect to the test tolerance.
    // ---------------------------------------------------------------

    Real testTolerance (1e-10);

    if ( (errorRhsNorm < testTolerance)
            && (errorIntegralg < testTolerance)
            && (errorVolume < testTolerance) )
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );
}


