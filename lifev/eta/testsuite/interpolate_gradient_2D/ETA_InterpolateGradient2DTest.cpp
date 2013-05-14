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
   @file ETA_InterpolateGradient2DTest.cpp
   @author L. Pasquale <luca.pasquale@mail.polimi.it>
   @date 2012-11-20
 */

// ===================================================
//! Includes
// ===================================================

#include "ETA_InterpolateGradient2DTest.hpp"

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

// ===================================================
//!                   Functions
// ===================================================

// ---------------------------------------------------------------
// We define a function whose gradient is
// (1 0)
// (0 2)
// ---------------------------------------------------------------

Real uFct ( const Real& /* t */, const Real&  x , const Real&  y , const Real& /* z */, const ID& i )
{
    if (i == 0)
    {
        return x;
    }
    return 2 * y;
}

// ===================================================
//!                  Constructors
// ===================================================

ETA_InterpolateGradient2DTest::ETA_InterpolateGradient2DTest ()
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

Real
ETA_InterpolateGradient2DTest::run()
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

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh2D ( *fullMeshPtr, 0, Nelements, Nelements, false,
                    2.0,   2.0,
                    -1.0,  -1.0);

    boost::shared_ptr< mesh_Type > meshPtr;
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
    // We start by defining the finite element spaces. We still need
    // a FESpace because ETFESpace is still lacking some methods
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building FESpaces ... " << std::flush;
    }

    std::string uOrder ("P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 2, M_comm) );

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

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 2, 2 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 2, 2 > (meshPtr, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), M_comm) );

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
    // and a variable used to store the result of the integration
    // ---------------------------------------------------------------

    MatrixSmall<2, 2> identityMatrix;
    identityMatrix (0, 0) = 1;
    identityMatrix (0, 1) = 0;
    identityMatrix (1, 0) = 0;
    identityMatrix (1, 1) = 1;

    Real ETintegral (0);

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    // ---------------------------------------------------------------
    // We integrate on the domain the gradient trace (1+2 in this case)
    //
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

    M_comm->Barrier();
    M_comm->SumAll (&ETintegral, &globalIntegral, 1);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    // ---------------------------------------------------------------
    // We now compute the error as the difference between the integral
    // computed and the exact value expected ( (1+2)*4 = 12 )
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Computing the error ... " << std::flush;
    }

    Real error = std::abs (globalIntegral - 12.0);

    if (verbose)
    {
        std::cout << "Error: " << error << std::endl;
    }

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    return error;

} // run
