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
   @file ETA_ADR2DTest.cpp
   @author L. Pasquale <luca.pasquale@mail.polimi.it>
   @date 2012-11-20
 */

// ===================================================
//! Includes
// ===================================================

#include "ETA_ADR1DTest.hpp"

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

typedef RegionMesh<LinearLine> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;

// ---------------------------------------------------------------
// We define then a function, which represents a velocity field,
// of magnitude 1 in the y direction (supposing an (x,y)
// reference frame).
// ---------------------------------------------------------------


// ===================================================
//!                   Functions
// ===================================================

Real betaFct ( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& /* i */ )
{
    return 1;
}

Real fRhs ( const Real& /* t */, const Real& /*x*/, const Real& /*y*/, const Real& /* z */ , const ID& /* i */ )
{
    return  2;
}


// ===================================================
//!                  Constructors
// ===================================================

ETA_ADR1DTest::ETA_ADR1DTest ()
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
ETA_ADR1DTest::run()
{
    bool verbose (M_comm->MyPID() == 0);
    // ---------------------------------------------------------------
    // We define the mesh and parition it. We use the domain
    // (-1,1) and a structured mesh.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building the mesh ... " << std::flush;
    }

    const UInt Nelements (10);

    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh1D ( *fullMeshPtr, 0, Nelements, false,
                    1.0,
                    -1.0);

    // Mesh partitioning doesn't work yet  with 1D meshes
    //MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_comm);
    //boost::shared_ptr< mesh_Type > meshPtr (meshPart.meshPartition() );
    //fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We start by defining the finite element spaces. We use the type
    // FESpace for the classical way and copy it in an ETFESpace for
    // the ET assembly.
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
    ( new FESpace< mesh_Type, MapEpetra > (fullMeshPtr, uOrder, 1, M_comm) );

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > betaSpace
    ( new FESpace< mesh_Type, MapEpetra > (fullMeshPtr, bOrder, 1, M_comm) );

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

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 1, 1 > > ETuSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 1, 1 > (fullMeshPtr, & (uSpace->refFE() ), & (uSpace->fe().geoMap() ), M_comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 1, 1 > > ETbetaSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 1, 1 > (fullMeshPtr, & (betaSpace->refFE() ), & (betaSpace->fe().geoMap() ), M_comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << ETuSpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // We then interpolate the advection function.
    // This is can only be performed with the classical FESpace.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Interpolating the advection field ... " << std::flush;
    }

    vector_Type beta (betaSpace->map(), Repeated);
    betaSpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (betaFct), beta, 0.0);

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
    // Definition of the RHS
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Defining and interpolating the RHS ... " << std::flush;
    }

    vector_Type rhs (uSpace->map(), Repeated);
    rhs *= 0.0;
    vector_Type ETrhs (uSpace->map(), Repeated);
    ETrhs *= 0.0;

    vector_Type fInterpolated (uSpace->map(), Repeated);
    fInterpolated *= 0.0;
    uSpace->interpolate ( static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (fRhs) , fInterpolated, 0.0 );

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
    //
    // We also assemble the RHS
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

    adrAssembler.addMassRhs (rhs, fInterpolated);

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
    // Remark also that the Real constants and VectorSmall constants
    // can be used directly in the expressions to integrate.
    // ---------------------------------------------------------------

    LifeChrono ETChrono;
    ETChrono.start();

    if (verbose)
    {
        std::cout << " -- ET assembly ... " << std::flush;
    }

    // Workaround for 1D case where grad(phi_j) is a VectorSmall<1>
    // and not a double: we use a dot product with a unitary VectorSmall<1>
    // to obtain its value as adouble
    VectorSmall<1> oneVector;
    oneVector[0] = 1.0;


    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr(),
                    ETuSpace,
                    ETuSpace,

                    dot ( grad (phi_i) , grad (phi_j) )
                    + dot ( grad (phi_j) , value (oneVector) ) * value (ETbetaSpace, beta) *phi_i
                    + 2.0 * phi_i * phi_j

                  )
                >> ETsystemMatrix;

        integrate ( elements (ETuSpace->mesh() ),
                    uSpace->qr(),
                    ETuSpace,

                    value (ETuSpace, fInterpolated) *phi_i

                  )
                >> ETrhs;

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
    // We finally need to check that both yield the same matrix. In
    // that aim, we need to finalize both matrices.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Closing the matrices and vectors... " << std::flush;
    }

    systemMatrix->globalAssemble();
    ETsystemMatrix->globalAssemble();

    rhs.globalAssemble();
    ETrhs.globalAssemble();

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


    // ---------------------------------------------------------------
    // Finally, we check that the two right hand sides are similar
    // with the two assembly procedures. We also check that the two
    // integrals correspond to their exact values.
    // ---------------------------------------------------------------

    vector_Type checkRhs (ETuSpace->map(), Repeated);
    checkRhs = 0.0;

    checkRhs += rhs;
    checkRhs -= ETrhs;

    checkRhs.globalAssemble();

    vector_Type checkRhsUnique (checkRhs, Unique);

    Real errorRhsNorm ( checkRhsUnique.normInf() );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " Matrix Error : " << errorNorm << std::endl;
    }
    if (verbose)
    {
        std::cout << " Rhs error : " << errorRhsNorm << std::endl;
    }


    // ---------------------------------------------------------------
    // Finally we return the maximum of the two errors
    // ---------------------------------------------------------------


    return std::max (errorNorm, errorRhsNorm);

} // run
