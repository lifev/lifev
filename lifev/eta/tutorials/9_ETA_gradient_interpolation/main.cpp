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
    @brief Tutorial for gradient interpolation.

    @author Luca Pasquale <luca.pasquale@mail.polimi.it>
    @date 2013-05

    In this tutorial we use two expressions available when
    integrating with Expression Templates:
    - interpoaltion of the gradient of a given function
    - multiplication by a MatrixSmall
    We do this by multiplying the gradient by the identity matrix
    (i.e. we're integrating the divergence)

    Tutorials that should be read before: 1,3

 */

// ---------------------------------------------------------------
// We reuse the same files as in the previous tutorial. We add the
// chrono to measure the timings.
// ---------------------------------------------------------------


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <boost/shared_ptr.hpp>


// ---------------------------------------------------------------
// We work in the LifeV namespace and define the mesh, matrix and
// vector types that we will need several times.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;


// ===================================================
//!                   Functions
// ===================================================

// ---------------------------------------------------------------
// We define a function whose gradient is
// (1 0 0)
// (0 5 0)
// (0 0 2)
// ---------------------------------------------------------------

Real uFct ( const Real& /* t */, const Real&  x , const Real&  y , const Real& z , const ID& i )
{
    if (i == 0)
    {
        return x;
    }
    else if (i == 1)
    {
        return 5 * y;
    }
    return 2 * z;
}

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

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);

    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
    boost::shared_ptr< mesh_Type > meshPtr (meshPart.meshPartition() );

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
        std::cout << " -- Building FESpaces ... " << std::flush;
    }

    std::string uOrder ("P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 3, Comm) );

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

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPtr, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << ETuSpace->dof().numTotalDof() << std::endl;
    }

    // ---------------------------------------------------------------
    // We interpolate then the function.
    // This can only be performed with the classical FESpace.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Interpolating the function ... " << std::flush;
    }

    vector_Type uInterpolated (uSpace->map(), Unique);
    uSpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (uFct), uInterpolated, 0.0);
    vector_Type uInterpolatedRepeated (uInterpolated, Repeated);

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We build the MatrixSmall used to extract the trace
    // and a variable used to store the result of the integration.
    // We only need to specify the diagonal elements of the MatrixSmall
    // because its default constructor sets all the components to 0
    // ---------------------------------------------------------------

    MatrixSmall<3, 3> identityMatrix;
    identityMatrix (0, 0) = 1;
    identityMatrix (1, 1) = 1;
    identityMatrix (2, 2) = 1;

    Real ETintegral (0);

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    // ---------------------------------------------------------------
    // We integrate on the domain the trace of the gradient
    // (1+5+2 in this case)
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

                    dot ( grad (ETuSpace, uInterpolatedRepeated) , value (identityMatrix) )

                  )
                >> ETintegral;

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
    // Integrals computed on each processor must be summed together
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Broadcasting and summing integrals across processes ... " << std::flush;
    }

    Real globalIntegral (0.0);

    Comm->Barrier();
    Comm->SumAll (&ETintegral, &globalIntegral, 1);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    // ---------------------------------------------------------------
    // We now compute the error as the difference between the integral
    // computed and the exact value expected ( (1+5+2)*8 = 64 )
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the error ... " << std::flush;
    }

    Real error = std::abs (globalIntegral - 64.0);

    if (verbose)
    {
        std::cout << "Error: " << error << std::endl;
    }

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    Real tolerance (1e-10);

    if (error < tolerance)
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );
}


