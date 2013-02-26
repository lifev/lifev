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
    @brief Tutorial for the assembly of the vectorial laplacian.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-10-2010

    In this tutorial, we revisit the first tutorial and explain how
    to transform it so that the vectorial laplacian is assembled. This
    gives an example of how to deal with vectorial unknowns.

    Tutorials that should be read before: 1
 */

// ---------------------------------------------------------------
// We use the same structure as the tutorial number 2, so that we
// can compare the results of the ETA framework with the classical
// assembly.
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
// We work in the LifeV namespace and define the mesh, matrix and
// vector types that we will need several times.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;


// ---------------------------------------------------------------
// As usual, we start with the definition of the MPI communicator
// and the boolean for the outputs.
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

    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // To assemble a vectorial problem, the main difference is in the
    // finite element space. To have a vectorial problem, we have to
    // indicate that the dimension of the field is not 1, but 3 (at
    // least in this example, but it could be 2 for a 2D problem).
    //
    // The field dimension appears as argument of the FESpace
    // constructor and as template argument for the ETFESpace (it is
    // the second 3, the first one indicates the dimension of the
    // space).
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the FESpace ... " << std::flush;
    }

    std::string uOrder ("P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPart, uOrder, 3, Comm) );

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

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPart, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << ETuSpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // We build now two matrices, one for each assembly procedure in
    // order to compare them in the end. There is no difference in the
    // construction of the matrices with respect to the scalar case,
    // as the map of the finite element spaces is already aware of the
    // field dimension.
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
    // With the classical assembly, nothing changes with respect to
    // a scalar laplacian.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Classical assembly ... " << std::flush;
    }

    ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;

    adrAssembler.setup (uSpace, uSpace);

    adrAssembler.addDiffusion (systemMatrix, 1.0);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // With the ETA framework, there is also no change with respect
    // to the scalar case.
    // ---------------------------------------------------------------

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

                  )
                >> ETsystemMatrix;
    }

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    // ---------------------------------------------------------------
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


