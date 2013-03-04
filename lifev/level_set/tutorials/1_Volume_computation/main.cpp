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
    @brief Tutorial introducing the expression assembly specialized for the level set problem.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-08-2012

    In this simple tutorial, we show how to use a special integration scheme
    to compute the volume comprise in the domain {phi>0} where phi is a given
    level set function. To this aim, we define a mesh, interpolate the phi function
    on it and then use the special integration.

    Before reading this tutorial, the tutorials concerning the ETA module
    should be read.
 */

// ---------------------------------------------------------------
// We include here the MPI headers for the parallel computations.
// The specific "pragma" instructions are used to avoid warning
// coming from the MPI library, that are not useful to us.
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


// ---------------------------------------------------------------
// We include then the required headers from LifeV. First of all,
// the definition file and mesh related files.
// ---------------------------------------------------------------

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

// ---------------------------------------------------------------
// We will need the EpetraVector and the FESpace
// (for the interpolation).
// ---------------------------------------------------------------

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>


// ---------------------------------------------------------------
// We will also need the ETFESpace from the ETA framework
// ---------------------------------------------------------------

#include <lifev/eta/fem/ETFESpace.hpp>


// ---------------------------------------------------------------
// To have the ETA framework working here, we need to include
// this file, which defines the mechanism for the integration.
// ---------------------------------------------------------------

#include <lifev/eta/expression/Integrate.hpp>


// ---------------------------------------------------------------
// Finally, we include shared pointer from boost since we use
// them explicitly in this tutorial.
// ---------------------------------------------------------------

#include <boost/shared_ptr.hpp>


// ---------------------------------------------------------------
// We also use the integration scheme for the level set. To have
// it working, we need the following file:
// ---------------------------------------------------------------

#include <lifev/level_set/fem/LevelSetQRAdapter.hpp>


// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. For clarity, we also
// make two typedefs for the mesh type and vector type.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef VectorEpetra vector_Type;


// ---------------------------------------------------------------
// We define also the level set function. Here, we set it to
// phi(x,y,z) = (x-3y)/sqrt(10)
// ---------------------------------------------------------------

Real phiFct (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& /*i*/)
{
    return (x - 3 * y) / std::sqrt (10);
}


// ---------------------------------------------------------------
// Finally, we also need a functor, to represent the Heaviside
// function, that we will need to integrate.
// ---------------------------------------------------------------

class HeavisideFunctor
{
public:

    typedef Real return_Type;

    return_Type operator() (const Real& value)
    {
        if (value > 0)
        {
            return 1;
        }
        return 0;
    }

    HeavisideFunctor() {}

    HeavisideFunctor ( const HeavisideFunctor& /*f*/) {}

    ~HeavisideFunctor() {}
};


// ---------------------------------------------------------------
// We start the programm by the definition of the communicator
// (as usual) depending on whether MPI is available or not. We
// also define a boolean to allow only one process to display
// messages.
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

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);

    MeshPartitioner< mesh_Type >  meshPart (fullMeshPtr, Comm);

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We define now the finite element space for the level set
    // function. We use here piecewise linear elements.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building ETFESpace ... " << std::flush;
    }

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > phiSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << phiSpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // For the interpolation, we need the FESpace.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Interpolation ... " << std::flush;
    }

    std::string phiOrder ("P1");
    FESpace<mesh_Type, MapEpetra> interpolationSpace (meshPart, phiOrder, 1, Comm);

    vector_Type phiInterpolated (phiSpace->map(), Repeated);
    interpolationSpace.interpolate (phiFct, phiInterpolated, 0.0);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We compute now the volume taken by the phase {phi>0}, that is
    // the volume integral of the heaviside function of phi.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the volume " << std::flush;
    }

    Real localVolume (0.0);

    // ---------------------------------------------------------------
    // We can now perform the integration. The evaluation (3rd
    // argument) represents the heaviside function of phi.
    //
    // The interesting difference is the 2nd argument, which is
    // usually the quadrature rule. Here, it is replaced by this
    // call to the special function "adapt", which transforms the
    // quadrature rule depending on the element encountered.
    //
    // The arguments to the adapt are the following:
    // - phiSpace: the finite element space in which the level set
    //             values are defined
    // - phiInterpolated: the values of the level set in the FE space
    // - quadRuleTetra1pt: the quadrature to adapt.
    // ---------------------------------------------------------------

    {
        using namespace ExpressionAssembly;

        boost::shared_ptr<HeavisideFunctor> heaviside (new HeavisideFunctor);

        integrate (  elements (phiSpace->mesh() ),
                     adapt (phiSpace, phiInterpolated, quadRuleTetra1pt),
                     eval (heaviside, value (phiSpace, phiInterpolated) )

                  )
                >> localVolume;
    }

    Real globalVolume (0.0);

    Comm->Barrier();
    Comm->SumAll (&localVolume, &globalVolume, 1);

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Volume: " << globalVolume << std::endl;
    }


    // ---------------------------------------------------------------
    // We finalize the MPI session if MPI was used
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    // ---------------------------------------------------------------
    // Finally, we check the norm with respect to a previously
    // computed one to ensure that it has not changed.
    // ---------------------------------------------------------------

    Real tol (1e-10);
    Real absError (std::abs (globalVolume - 4) );

    if (verbose)
    {
        std::cout << " Error : " << absError << std::endl;
    }

    if ( absError < tol )
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );

}
