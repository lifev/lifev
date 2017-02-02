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
    @brief Tutorial about the quadrature rule in the ETA.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 28-06-2012

    In this tutorial, we show how to better select the
    quadrature rule, via the dedicated class called
    QuadratureRuleProvider. This is the right way to
    access quadrature rules.

    Tutorial 1 should be read before, to understand the basic
    concepts.
 */

// ---------------------------------------------------------------
// We include here exactly the same files that in tutorial 1,
// excepted a new one: the file containing the
// QuadratureRuleProvider
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

#include <lifev/core/fem/QuadratureRuleProvider.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#include <boost/shared_ptr.hpp>


// ---------------------------------------------------------------
// The beginning of the tutorial is exactly the same as the
// tutorial 1, so we refer to that tutorial for more comments.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;


int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);


    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }

    const UInt Nelements (10);

    std::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);

    MeshPartitioner< mesh_Type >  meshPart (fullMeshPtr, Comm);

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    if (verbose)
    {
        std::cout << " -- Building ETFESpaces ... " << std::flush;
    }

    std::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );

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
        std::cout << " -- Defining the matrix ... " << std::flush;
    }

    std::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uSpace->map() ) );

    *systemMatrix *= 0.0;

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    if (verbose)
    {
        std::cout << " -- Assembling the Laplace matrix ... " << std::flush;
    }


    {
        using namespace ExpressionAssembly;

        // ---------------------------------------------------------------
        // Now comes the interesting part. During the assembly, we need
        // to specify which quadrature rule we want to use for the
        // integration. In the tutorial 1, we used "quadRuleTetra4pt"
        // as quadrature rule.
        //
        // This is not very convenient, since
        // - One must know the quadrature rules available and their names
        // - No information about the degrees of the quadrature is given.
        //
        // The QuadratureRuleProvider class was designed exactly to
        // overcome these issues. The most useful method of that class
        // is the provideExactness: this method returns, for a given
        // element shape, a quadrature rule with the required exactness.
        // That is what we use here.
        //
        // The required exactness here is 0, since the function to be
        // integrated is constant in each element. To compute
        // the mass matrix, we would have chosen 2 (the function to
        // integrate would be quadratic in each element).
        // ---------------------------------------------------------------

        integrate (  elements (uSpace->mesh() ),
                     QuadratureRuleProvider::provideExactness (TETRA, 0),
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
                  )
                >> systemMatrix;
    }

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // Remark that you can change the behaviour of the
    // QuadratureRuleProvider with the setters that it has, e.g., to
    // obtain a warning when the precise exactness was not found, we
    // use the code
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Changing the behaviour of the QRProvider ... " << std::flush;
    }

    if (verbose)
    {
        QuadratureRuleProvider::setBehaviorNoPreciseExactness (QuadratureRuleProvider::WarningAndReturnSup );
    }

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We finish the test as in tutorial 1.
    // ---------------------------------------------------------------


    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }

    systemMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    Real matrixNorm ( systemMatrix->normInf() );

    if (verbose)
    {
        std::cout << " Matrix norm : " << matrixNorm << std::endl;
    }


#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    Real matrixNormDiff (std::abs (matrixNorm - 3.2) );

    if (verbose)
    {
        std::cout << " Error : " << matrixNormDiff << std::endl;
    }

    Real testTolerance (1e-10);

    if ( matrixNormDiff < testTolerance )
    {
        return ( EXIT_SUCCESS );
    }
    return ( EXIT_FAILURE );

}


