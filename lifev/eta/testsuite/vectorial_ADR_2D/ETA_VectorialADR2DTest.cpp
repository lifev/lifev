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
   @file ETA_VectorialADR2DTest.cpp
   @author L. Pasquale <luca.pasquale@mail.polimi.it>
   @date 2012-11-20
 */

// ===================================================
//! Includes
// ===================================================

#include "ETA_VectorialADR2DTest.hpp"

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

// ---------------------------------------------------------------
// We define then a function, which represents a velocity field,
// of magnitude 1 in the y direction (supposing an (x,y)
// reference frame).
// ---------------------------------------------------------------


// ===================================================
//!                   Functions
// ===================================================

Real betaFct( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& i )
{
    if (i == 1)
    {
        return 1;
    }
    return 0;
}

// ===================================================
//!                  Constructors
// ===================================================

ETA_VectorialADR2DTest::ETA_VectorialADR2DTest ()
{

#ifdef EPETRA_MPI
    M_comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
#else
    M_comm.reset( new Epetra_SerialComm() );
#endif

}

// ===================================================
//!                      Methods
// ===================================================

Real
ETA_VectorialADR2DTest::run()
{
    bool verbose(M_comm->MyPID()==0);
// ---------------------------------------------------------------
// We define the mesh and parition it. We use the domain
// (-1,1)x(-1,1) and a structured mesh.
// ---------------------------------------------------------------

    if (verbose) std::cout << " -- Building and partitioning the mesh ... " << std::flush;

    const UInt Nelements(10);

    boost::shared_ptr< mesh_Type > fullMeshPtr(new mesh_Type);

    regularMesh2D( *fullMeshPtr, 0, Nelements, Nelements, false,
                   2.0,   2.0,
                   -0.0,  -0.0);

    MeshPartitioner< mesh_Type >   meshPart(fullMeshPtr, M_comm);
    boost::shared_ptr< mesh_Type > meshPtr (meshPart.meshPartition());

    fullMeshPtr.reset();

    if (verbose) std::cout << " done ! " << std::endl;


// ---------------------------------------------------------------
// We start by defining the finite element spaces. We use the type
// FESpace for the classical way and copy it in an ETFESpace for
// the ET assembly.
//
// In the end, both solution spaces display their respective
// number of degrees of freedom, which must be the same.
// ---------------------------------------------------------------

    if (verbose) std::cout << " -- Building FESpaces ... " << std::flush;

    std::string uOrder("P1");
    std::string bOrder("P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uSpace
        ( new FESpace< mesh_Type, MapEpetra >(meshPtr,uOrder, 2, M_comm));

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > betaSpace
        ( new FESpace< mesh_Type, MapEpetra >(meshPtr,bOrder, 2, M_comm));

    if (verbose) std::cout << " done ! " << std::endl;
    if (verbose) std::cout << " ---> Dofs: " << uSpace->dof().numTotalDof() << std::endl;

    if (verbose) std::cout << " -- Building ETFESpaces ... " << std::flush;

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 2, 2 > > ETuSpace
        ( new ETFESpace< mesh_Type, MapEpetra, 2, 2 >(meshPart,&(uSpace->refFE()),&(uSpace->fe().geoMap()), M_comm));

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 2, 2 > > ETbetaSpace
        ( new ETFESpace< mesh_Type, MapEpetra, 2, 2 >(meshPart,&(betaSpace->refFE()),&(betaSpace->fe().geoMap()), M_comm));

    if (verbose) std::cout << " done ! " << std::endl;
    if (verbose) std::cout << " ---> Dofs: " << ETuSpace->dof().numTotalDof() << std::endl;


// ---------------------------------------------------------------
// We interpolate then the advection function of the mesh at hand.
// This is performed with the classical FESpace only.
// ---------------------------------------------------------------

    if (verbose) std::cout << " -- Interpolating the advection field ... " << std::flush;

    vector_Type beta(betaSpace->map(),Repeated);
    betaSpace->interpolate(static_cast<FESpace< mesh_Type, MapEpetra >::function_Type>(betaFct),beta,0.0);

    if (verbose) std::cout << " done! " << std::endl;


// ---------------------------------------------------------------
// We build now two matrices, one for each assembly procedure in
// order to compare them in the end. There is no difference in the
// construction of the matrices.
// ---------------------------------------------------------------

    if (verbose) std::cout << " -- Defining the matrices ... " << std::flush;

    boost::shared_ptr<matrix_Type> systemMatrix(new matrix_Type( uSpace->map() ));
    *systemMatrix *=0.0;

    boost::shared_ptr<matrix_Type> ETsystemMatrix(new matrix_Type( ETuSpace->map() ));
    *ETsystemMatrix *=0.0;

    if (verbose) std::cout << " done! " << std::endl;


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

    if (verbose) std::cout << " -- Classical assembly ... " << std::flush;

    ADRAssembler<mesh_Type,matrix_Type,vector_Type> adrAssembler;

    adrAssembler.setup(uSpace,betaSpace);

    adrAssembler.addDiffusion(systemMatrix,1.0);

    adrAssembler.addAdvection(systemMatrix,beta);

    adrAssembler.addMass(systemMatrix,2.0);

    StdChrono.stop();

    if (verbose) std::cout << " done ! " << std::endl;
    if (verbose) std::cout << " Time : " << StdChrono.diff() << std::endl;


// ---------------------------------------------------------------
// We perform now the same assembly with the ET assembly and still
// monitor the timings required.
// 
// Remark also that the Real constants and VectorSmall constants
// can be used directly in the expressions to integrate.
// ---------------------------------------------------------------

    LifeChrono ETChrono;
    ETChrono.start();

    if (verbose) std::cout << " -- ET assembly ... " << std::flush;


    {
        using namespace ExpressionAssembly;

        integrate( elements(ETuSpace->mesh()),
                   uSpace->qr(),
                   ETuSpace,
                   ETuSpace,

                   dot( grad(phi_i) , grad(phi_j) )
                   +dot( grad(phi_j) * value(ETbetaSpace,beta), phi_i)
                   + 2.0*dot(phi_i,phi_j)
                   
                   )
            >> ETsystemMatrix;

    }

    ETChrono.stop();

    if (verbose) std::cout << " done! " << std::endl;
    if (verbose) std::cout << " Time : " << ETChrono.diff() << std::endl;


// ---------------------------------------------------------------
// We finally need to check that both yield the same matrix. In
// that aim, we need to finalize both matrices.
// ---------------------------------------------------------------

    if (verbose) std::cout << " -- Closing the matrices ... " << std::flush;

    systemMatrix->globalAssemble();
    ETsystemMatrix->globalAssemble();
    
    if (verbose) std::cout << " done ! " << std::endl;

// ---------------------------------------------------------------
// We compute now the matrix of the difference and finally the 
// norm of the difference. This should be very low if the two
// matrices are identical.
// ---------------------------------------------------------------

    if (verbose) std::cout << " -- Computing the error ... " << std::flush;

    boost::shared_ptr<matrix_Type> checkMatrix(new matrix_Type( ETuSpace->map()));
    *checkMatrix *=0.0;

    *checkMatrix += *systemMatrix;
    *checkMatrix += (*ETsystemMatrix)*(-1);

    checkMatrix->globalAssemble();

    Real errorNorm( checkMatrix->normInf() );

    if (verbose) std::cout << " done ! " << std::endl;



// ---------------------------------------------------------------
// We finally display the error norm and compare with the
// tolerance of the test.
// ---------------------------------------------------------------

    if (verbose) std::cout << " Matrix Error : " << errorNorm << std::endl;

    return errorNorm;

} // run
