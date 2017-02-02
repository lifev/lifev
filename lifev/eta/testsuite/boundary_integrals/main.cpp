/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Claudia Colciago <claudia.colciago@epfl.ch>
       Date: 2013-07-29

  Copyright (C) 2005 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file main.cpp
   \author Claudia Colciago<claudia.colciago@epfl.ch>
   \date 2013-07-29

 */


// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>


#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/core/array/MatrixBlockMonolithicEpetra.hpp>
#include <lifev/core/array/VectorBlockMonolithicEpetra.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/DOFInterface3Dto3D.hpp>
#include <lifev/core/fem/GradientRecovery.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/level_set/fem/LevelSetQRAdapter.hpp>
#include <iostream>

using namespace LifeV;

const Real pi = 3.141592653589793;

// Coefficients
Real beta ( 1.0 );
Real alpha ( 1.0 );
Real nu ( 1.0 );
UInt wall ( 30 );

Real uExact ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return std::sin ( pi * y ) * std::cos ( pi * x ) * std::exp ( z ) ;
}

Real laplacianExact ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return    ( - pi * pi * std::sin (pi * y) * std::cos (pi * x ) * std::exp ( z )
                - pi * pi * std::sin (pi * y) * std::cos (pi * x ) * std::exp ( z )
                + std::sin (pi * y) * cos (pi * x ) * exp ( z ) );
}

Real gRobinRhs ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    MatrixSmall<3, 3> hessian;

    hessian[0][0] = - pi * pi * std::sin ( pi * y ) * std::cos ( pi * x ) * std::exp ( z );

    hessian[0][1] = - pi * pi * std::sin (pi * x ) * std::cos (pi * y) * std::exp ( z );
    hessian[1][0] = hessian[0][1];

    hessian[0][2] = - pi * std::sin ( pi * y ) * std::sin ( pi * x ) * std::exp ( z );
    hessian[2][0] = hessian[0][2];

    hessian[1][1] = - pi * pi * std::sin ( pi * y ) * std::cos ( pi * x ) * std::exp ( z );
    hessian[2][2] =  std::sin ( pi * y ) * cos ( pi * x ) * exp ( z );

    hessian[1][2] =  pi * std::cos ( pi * y ) * std::cos ( pi * x ) * std::exp ( z );
    hessian[2][1] =  hessian[1][2];

    VectorSmall<3> gradient;
    gradient[0] = - pi * std::sin ( pi * y ) * std::sin ( pi * x ) * std::exp ( z );
    gradient[1] = pi * std::cos ( pi * y ) * std::cos ( pi * x ) * std::exp ( z );
    gradient[2] = std::sin ( pi * y ) * std::cos ( pi * x ) * std::exp ( z );

    VectorSmall<3> normal;

    Real traceHessian = hessian[0][0] + hessian[1][1] + hessian[2][2];

    normal[0] = 0;
    normal[1] = 0;
    normal[2] = 1;

    normal.normalize();

    return ( nu * gradient.dot ( normal ) + alpha * uExact ( 0, x, y, z, 0 )
             - beta * ( traceHessian - ( hessian * normal ).dot ( normal ) ) );
}

/* LaplacianRhs */
class LaplacianExactFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return laplacianExact ( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    LaplacianExactFunctor() {}
    LaplacianExactFunctor (const LaplacianExactFunctor&) {}
    ~LaplacianExactFunctor() {}
};

class uExactFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return uExact ( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    uExactFunctor() {}
    uExactFunctor (const uExactFunctor&) {}
    ~uExactFunctor() {}
};

class gradExactFunctor
{
public:
    typedef VectorSmall<3> return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        VectorSmall<3> gradient;
        Real x = spaceCoordinates[0];
        Real y = spaceCoordinates[1];
        Real z = spaceCoordinates[2];
        gradient[0] = - pi * std::sin ( pi * y ) * std::sin ( pi * x ) * std::exp ( z );
        gradient[1] = pi * std::cos ( pi * y ) * std::cos ( pi * x ) * std::exp ( z );
        gradient[2] = std::sin ( pi * y ) * std::cos ( pi * x ) * std::exp ( z );

        return gradient;
    }

    gradExactFunctor() {}
    gradExactFunctor (const gradExactFunctor&) {}
    ~gradExactFunctor() {}
};

class GRobinRhsFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return gRobinRhs ( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    GRobinRhsFunctor() {}
    GRobinRhsFunctor (const GRobinRhsFunctor&) {}
    ~GRobinRhsFunctor() {}
};

/* Register the preconditioner */
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}

/* Some typedef */
typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef std::shared_ptr<vector_Type> vectorPtr_Type;
typedef std::shared_ptr<matrix_Type> matrixPtr_Type;

typedef LifeV::Preconditioner             basePrec_Type;
typedef std::shared_ptr<basePrec_Type>  basePrecPtr_Type;
typedef LifeV::PreconditionerIfpack       prec_Type;
typedef std::shared_ptr<prec_Type>      precPtr_Type;

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    // a flag to see who's the leader for output purposes
    bool verbose = (Comm->MyPID() == 0);

    // Open and read the data file
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    // Load the mesh
    MeshData dataMesh;
    dataMesh.setup (dataFile, "mesh");
    std::shared_ptr < mesh_Type > fullMeshPtr (new mesh_Type);
    readMesh (*fullMeshPtr, dataMesh);

    // Partition the mesh
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
    std::shared_ptr < mesh_Type > localMeshPtr (new mesh_Type);
    localMeshPtr = meshPart.meshPartition();

    // Free the global mesh
    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " Building FESpaces  " << std::endl;
    }

    std::string uOrder ("P1");

    std::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace ( new FESpace< mesh_Type, MapEpetra > (localMeshPtr, uOrder, 1, Comm) );

    if (verbose)
    {
        std::cout << std::endl << " ### Dof Summary ###: " <<  std::endl;
    }
    if (verbose)
    {
        std::cout << " u  : " << uFESpace->map().map (Unique)->NumGlobalElements() << std::endl;
    }

    if (verbose)
    {
        std::cout << " Building EA FESpaces  " << std::endl;
    }

    std::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETuFESpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (localMeshPtr, & (uFESpace->refFE() ), Comm) );

    vectorPtr_Type uSolution ( new vector_Type ( ETuFESpace->map() , Unique) );

    if (verbose)
    {
        std::cout << " Building the solvers " << std::endl;
    }


    if ( verbose )
    {
        std::cout << "Setting up LinearSolver (AztecOO)... " << std::flush;
    }

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

    LinearSolver linearSolver;
    linearSolver.setCommunicator ( Comm );
    linearSolver.setParameters ( *belosList );
    linearSolver.setPreconditioner ( precPtr );
    if ( verbose )
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << " Building the exporter " << std::endl;
    }

    std::string const exporterFileName    =  dataFile ( "exporter/filename", "cube");
    ExporterHDF5<mesh_Type> exporter ( dataFile, meshPart.meshPartition(), exporterFileName, Comm->MyPID() );
    exporter.setMultimesh (false);

    std::shared_ptr<vector_Type> uExported ( new vector_Type (ETuFESpace->map(), exporter.mapType() ) );
    std::shared_ptr<vector_Type> errorExported ( new vector_Type (ETuFESpace->map(), exporter.mapType() ) );

    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "u", uFESpace, uExported, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "err", uFESpace, errorExported, UInt (0) );

    LifeChrono ChronoItem;
    ChronoItem.start();
    if (verbose)
    {
        std::cout << " Assembling the matrix ... " << std::flush;
    }

    QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );

    matrixPtr_Type systemMatrix (new matrix_Type ( ETuFESpace->map() ) );
    *systemMatrix *= 0.0;

    {
        using namespace ExpressionAssembly;

        integrate (
            elements (ETuFESpace->mesh() ), // Mesh

            uFESpace->qr(), // QR

            ETuFESpace,
            ETuFESpace,

            value (nu) * dot (grad (phi_i) , grad (phi_j) )

        )
                >> *systemMatrix;

        integrate ( boundary (ETuFESpace->mesh(), wall),
                    myBDQR,

                    ETuFESpace,
                    ETuFESpace,

                    value ( alpha ) * phi_j * phi_i

                    + value ( beta ) * dot ( ( grad (phi_j) - dot ( grad (phi_j) , Nface ) * Nface ) ,
                                             ( grad (phi_i) - dot ( grad (phi_i) , Nface ) * Nface ) )

                  )
                >> *systemMatrix;
    }

    ChronoItem.stop();
    if (verbose)
    {
        std::cout << ChronoItem.diff() << " s" << std::endl;
    }

    if (verbose)
    {
        std::cout << " Assembling the rhs ... " << std::flush;
    }

    ChronoItem.start();

    std::shared_ptr<LaplacianExactFunctor> laplacianFctRhs ( new LaplacianExactFunctor );
    std::shared_ptr<GRobinRhsFunctor> gRobinFctRhs ( new GRobinRhsFunctor );

    vector_Type uRhs ( ETuFESpace->map() , Repeated );
    uRhs *= 0.0;

    {
        using namespace ExpressionAssembly;

        integrate ( elements (ETuFESpace->mesh() ), // Mesh

                    uFESpace->qr(), // QR

                    ETuFESpace,

                    value (-nu) * eval ( laplacianFctRhs, X ) * phi_i

                  )
                >> uRhs;


        integrate ( boundary (ETuFESpace->mesh(), wall),
                    myBDQR,

                    ETuFESpace,

                    eval ( gRobinFctRhs, X ) * phi_i

                  )
                >> uRhs;
    }

    ChronoItem.stop();
    if (verbose)
    {
        std::cout << ChronoItem.diff() << " s" << std::endl;
    }

    vectorPtr_Type uRhsUnique ( new vector_Type ( uRhs, Unique ) );

    systemMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << "[Navier-Stokes] Applying Dirichlet boundary conditions ... " << std::flush;
    }

    ChronoItem.start();
    BCHandler bcHandler;
    BCFunctionBase uexBCFct ( uExact );

    bcHandler.addBC ("Wall2", 31, Essential, Full, uexBCFct, 1);
    bcHandler.addBC ("Corners", 32, EssentialEdges, Full, uexBCFct, 1);

    bcHandler.bcUpdate ( *meshPart.meshPartition(), uFESpace->feBd(), uFESpace->dof() );
    bcManage (*systemMatrix, *uRhsUnique,
              *uFESpace->mesh(), uFESpace->dof(),
              bcHandler, uFESpace->feBd(), 1.0, Real (0.0) );

    ChronoItem.stop();

    if (verbose)
    {
        std::cout << ChronoItem.diff() << " s" << std::endl;
    }

    if (verbose)
    {
        std::cout << " Solving the system " << std::endl;
    }


    linearSolver.setOperator ( systemMatrix );
    linearSolver.setRightHandSide ( uRhsUnique );
    linearSolver.solve ( uSolution );

    *uExported = *uSolution;
    vector_Type errorVector ( ETuFESpace->map() , Unique , Zero);
    vector_Type uexVector ( uFESpace->map(), Unique);
    uFESpace->interpolate ( uExact, uexVector, 0);
    errorVector = uexVector  -  *uSolution;
    *errorExported = errorVector;
    exporter.postProcess ( 1.0 );

    if (verbose)
    {
        std::cout << " Evaluate Errors " << std::endl;
    }

    Real errorL2SquaredLocal ( 0.0 );
    Real errorH1SquaredLocal ( 0.0 );
    Real errorL2Squared ( 0.0 );
    Real errorH1Squared ( 0.0 );
    Real errorH1BoundarySquared ( 0.0 );

    vector_Type errorH1BoundaryVector ( ETuFESpace->map(), Repeated );
    vector_Type errorH1BoundaryVectorUnique ( ETuFESpace->map() );

    std::shared_ptr<uExactFunctor> uExactFct ( new uExactFunctor );
    std::shared_ptr<gradExactFunctor> gradExactFct ( new gradExactFunctor );

    {
        using namespace ExpressionAssembly;
        integrate (
            elements (ETuFESpace->mesh() ), // Mesh

            uFESpace->qr(), // QR

            dot (  ( eval ( gradExactFct, X ) - grad ( ETuFESpace , *uSolution ) ) ,
                   ( eval ( gradExactFct, X ) - grad ( ETuFESpace , *uSolution ) ) ) /** phi_i*/

        )
                >> errorH1SquaredLocal;


        integrate (
            elements (ETuFESpace->mesh() ), // Mesh

            uFESpace->qr(), // QR

            ( eval ( uExactFct, X ) - value ( ETuFESpace , *uSolution ) ) * (  eval ( uExactFct, X ) - value ( ETuFESpace , *uSolution ) )

        )
                >> errorL2SquaredLocal;

        integrate ( boundary (ETuFESpace->mesh(), wall),
                    myBDQR,

                    ETuFESpace,

                    dot ( ( eval ( gradExactFct, X ) - dot ( eval ( gradExactFct, X ) , Nface ) * Nface )
                          -  ( grad ( ETuFESpace , *uSolution ) - dot ( grad ( ETuFESpace , *uSolution ) , Nface ) * Nface ) ,
                          ( eval ( gradExactFct, X ) - dot ( eval ( gradExactFct, X ) , Nface ) * Nface )
                          -  ( grad ( ETuFESpace , *uSolution ) - dot ( grad ( ETuFESpace , *uSolution ) , Nface ) * Nface )
                        ) * phi_i  +
                    ( eval ( uExactFct, X ) - value ( ETuFESpace , *uSolution ) )
                    * (  eval ( uExactFct, X ) - value ( ETuFESpace , *uSolution ) ) * phi_i

                  )

                >> errorH1BoundaryVector;

    }

    Comm->Barrier();
    Comm->SumAll (&errorH1SquaredLocal, &errorH1Squared, 1);
    Comm->SumAll (&errorL2SquaredLocal, &errorL2Squared, 1);

    vector_Type oneVector ( ETuFESpace->map(), Unique );
    errorH1BoundaryVectorUnique = errorH1BoundaryVector;
    oneVector *= 0;
    oneVector += 1;
    errorH1BoundarySquared = errorH1BoundaryVectorUnique.dot ( oneVector );

    if (verbose)
    {
        std::cout << " l2 error norm " <<  std::sqrt ( errorL2Squared )  << std::endl;

        std::cout << " H1 error norm " <<  sqrt ( errorH1Squared ) << std::endl;

        std::cout << " H1 Gamma error norm " << std::sqrt ( errorH1BoundarySquared ) << std::endl;
    }

    exporter.closeFile();

    Real tolerance (1e-5);
    bool success ( false );

    if ( (  abs ( sqrt (errorL2Squared) - 0.0768669 ) < tolerance )
            && (  abs ( sqrt (errorH1Squared) - 1.76249  ) < tolerance )
            && ( abs ( sqrt (errorH1BoundarySquared) - 2.35064 ) < tolerance )
       )
    {
        success = true ;
    }

    if (!success)
    {
        if (verbose)
        {
            std::cout << "End Result: TEST NOT PASSED" << std::endl;
        }
    }
    else
    {
        if (verbose)
        {
            std::cout << "End Result: TEST PASSED" << std::endl;
        }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( !success )
    {
        return ( EXIT_FAILURE );
    }
    else
    {
        return ( EXIT_SUCCESS );
    }


}


