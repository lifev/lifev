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
    @brief Laplacian problem

    @author Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
    @date 19-05-2015

    This code aims at showing how to assemble a stiffness matrix as sum of
    matrices corresponding to the restriction of the bilinear form to subdomains.
    This is done by a extension of the method integrate
    that computes the stiffness matrix only over a subset of the volumes identified by a
    specific physical flag. The volume subdivision is done through the class
    lifev/core/mesh/MeshVolumeSubdivision

    The macro TESTMESHSUB allows to check the correctness of the local matrices by comparing
    them to the ones computed by the usual integrate (simply by defining for every subdomain a
    diffusion coefficient that is an indicator function on it). The matrices computed
    in this way are saved and can be checked. Then the full matrix is assembled using the standard
    integrate and both results are exported in paraview.
    The check is coded for the mesh 4cube1.mesh which correspond to have a [0,1]^3 cube divided in
    4 regions. The diffusion coefficients are read from datafile.
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChronoManager.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/BCManage.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/BuildGraph.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <Epetra_FECrsGraph.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/filter/ExporterHDF5.hpp>

#include <boost/shared_ptr.hpp>

#include <lifev/eta/examples/laplacian/laplacianFunctor.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/ElementShapes.hpp>
#include <lifev/core/mesh/MeshEntityContainer.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshVolumeSubdivision.hpp>

using namespace LifeV;

// Dirichlet BC functions
Real zeroFunction (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.;
}

Real sourceFunction (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    return 100.;
}

const int BACK   = 101;
const int FRONT  = 102;
const int LEFT   = 103;
const int RIGHT  = 104;
const int BOTTOM = 105;
const int TOP    = 106;

#define TESTMESHSUB


Real diffusion4Function1 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    if( z<0.5 && y < 0.5 ) return 1.0;

    return 0.;
}
Real diffusion4Function2 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    if( z>0.5 && y < 0.5 ) return 1.0;

    return 0.;
}

Real diffusion4Function3 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    if( z<0.5 && y > 0.5 ) return 1.0;

    return 0.;
}
Real diffusion4Function4 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    if( z>0.5 && y > 0.5 ) return 1.0;

    return 0.;
}

Real mu1 = 1.0;
Real mu2 = 1.0;
Real mu3 = 1.0;
Real mu4 = 1.0;

Real diffusion4FunctionTotal (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    if( z<0.5 && y < 0.5 ) return mu1;
    if( z>0.5 && y < 0.5 ) return mu2;
    if( z<0.5 && y > 0.5 ) return mu3;

    return mu4;
}


int main ( int argc, char** argv )
{

    // MPI initialization
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
    boost::shared_ptr< Epetra_Comm > Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr< Epetra_Comm > Comm ( new Epetra_SerialComm );
#endif

    // Reading parameters through GetPot
    GetPot command_line ( argc, argv );
    const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
    GetPot dataFile ( dataFileName );

    const bool verbose ( Comm->MyPID() == 0 );

    if ( verbose )
    {
        std::cout << "-- Building and partitioning the mesh... " << std::flush;
    }


    // +-----------------------------------------------+
    // |                Building mesh                  |
    // +-----------------------------------------------+

    typedef RegionMesh< LinearTetra >                                       mesh_Type;

    boost::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );

    MeshData meshData;
    meshData.setup ( dataFile, "mesh");
    readMesh (*fullMeshPtr, meshData);

    boost::shared_ptr< mesh_Type > localMeshPtr;

    MeshPartitioner< mesh_Type > meshPart;

    meshPart.doPartition ( fullMeshPtr, Comm );
    localMeshPtr = meshPart.meshPartition();

    // Clearing global mesh
    fullMeshPtr.reset();

    if ( verbose )
    {
        std::cout << " done" << std::endl;
    }

    // +-----------------------------------------------+
    // |         Building FE space and matrix          |
    // +-----------------------------------------------+

    typedef FESpace< mesh_Type, MapEpetra >                                 uSpaceStd_Type;
    typedef boost::shared_ptr< uSpaceStd_Type >                             uSpaceStdPtr_Type;

    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >                         uSpaceETA_Type;
    typedef boost::shared_ptr< uSpaceETA_Type >                             uSpaceETAPtr_Type;

    typedef FESpace<mesh_Type, MapEpetra>::function_Type                    function_Type;

    if ( verbose )
    {
        std::cout << "-- Building finite elements spaces... " << std::flush;
    }

    // Defining finite elements standard and ET spaces
    uSpaceStdPtr_Type uFESpace ( new uSpaceStd_Type ( localMeshPtr, dataFile( "finite_element/degree", "P1" ), 1, Comm ) );
    uSpaceETAPtr_Type ETuFESpace ( new uSpaceETA_Type ( localMeshPtr, & ( uFESpace->refFE() ), & ( uFESpace->fe().geoMap() ), Comm ) );

    if ( verbose )
    {
        std::cout << " done. " << " (Dofs: " << uFESpace->dof().numTotalDof() << ")" << std::endl;
    }

    typedef MatrixEpetra< Real >                                            matrix_Type;
    typedef boost::shared_ptr< MatrixEpetra< Real > >                       matrixPtr_Type;
    typedef VectorEpetra                                                    vector_Type;
    typedef boost::shared_ptr<VectorEpetra>                                 vectorPtr_Type;

    if ( verbose )
    {
        std::cout << " done " << std::endl;
    }

    // Clearing problem's matrix
    if ( verbose )
    {
        std::cout << "-- Subdividing the mesh... " << std::flush;
    }

    int numSubregions = dataFile( "physicalParameters/numSubregions", 1 );

    typedef LifeV::MeshVolumeSubdivision<RegionMesh<LinearTetra> >          meshSub_Type;
    typedef boost::shared_ptr<meshSub_Type>                                 meshSubPtr_Type;

    // suppose we are running the test with 4cube1.mesh,
    // so we have 4 subregions identified by the flags 1001,1002,1003,1004
    Epetra_IntSerialDenseVector regions(numSubregions);

    regions(0) = 1001;
    regions(1) = 1002;
    regions(2) = 1003;
    regions(3) = 1004;

    meshSubPtr_Type meshSub;
    meshSub.reset( new meshSub_Type( Comm, localMeshPtr, regions, numSubregions ) );

    meshSub->makeSubDivision();

//    meshSub->printNumElementPerFlag();

    if ( verbose )
    {
        std::cout << "-- Assembling the Laplace submatrices... " << std::flush;
    }

    matrixPtr_Type systemMatrixMeshSub1( new matrix_Type( ETuFESpace->map(), 100 ) );
    matrixPtr_Type systemMatrixMeshSub2( new matrix_Type( ETuFESpace->map(), 100 ) );
    matrixPtr_Type systemMatrixMeshSub3( new matrix_Type( ETuFESpace->map(), 100 ) );
    matrixPtr_Type systemMatrixMeshSub4( new matrix_Type( ETuFESpace->map(), 100 ) );

    {
        using namespace ExpressionAssembly;

        if ( verbose )
        {
            std::cout << "Building with MeshSub" << std::endl;
        }

        integrate ( elements ( localMeshPtr, meshSub->getFlag(0),
                               meshSub->getNumElements( 0 ), meshSub->getSubmesh( 0 ),
                               true ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                value(1.0) * dot( grad( phi_j ), grad( phi_i ) )
                )
                >> systemMatrixMeshSub1;

        integrate ( elements ( localMeshPtr, meshSub->getFlag(1),
                               meshSub->getNumElements( 1 ), meshSub->getSubmesh( 1 ),
                               true ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                value(1.0) * dot( grad( phi_j ), grad( phi_i ) )
                )
                >> systemMatrixMeshSub2;

        integrate ( elements ( localMeshPtr, meshSub->getFlag(2),
                               meshSub->getNumElements( 2 ), meshSub->getSubmesh( 2 ),
                               true ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                value(1.0) * dot( grad( phi_j ), grad( phi_i ) )
                )
                >> systemMatrixMeshSub3;

        integrate ( elements ( localMeshPtr, meshSub->getFlag(3),
                               meshSub->getNumElements( 3 ), meshSub->getSubmesh( 3 ),
                               true ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                value(1.0) * dot( grad( phi_j ), grad( phi_i ) )
                )
                >> systemMatrixMeshSub4;

    }

    systemMatrixMeshSub1->globalAssemble();
    systemMatrixMeshSub2->globalAssemble();
    systemMatrixMeshSub3->globalAssemble();
    systemMatrixMeshSub4->globalAssemble();

    // diffusion parameters
    mu1 = dataFile( "physicalParameters/mu1", 1 );
    mu2 = dataFile( "physicalParameters/mu2", 1 );
    mu3 = dataFile( "physicalParameters/mu3", 1 );
    mu4 = dataFile( "physicalParameters/mu4", 1 );

    matrixPtr_Type matrixMeshSubFull( new matrix_Type( ETuFESpace->map(), 100 ) );

    *matrixMeshSubFull += ( *systemMatrixMeshSub1 ) * mu1;
    *matrixMeshSubFull += ( *systemMatrixMeshSub2 ) * mu2;
    *matrixMeshSubFull += ( *systemMatrixMeshSub3 ) * mu3;
    *matrixMeshSubFull += ( *systemMatrixMeshSub4 ) * mu4;

    if (verbose)
    {
        std::cout << " done" << std::endl;
    }

    // +-----------------------------------------------+
    // |     Initializing vectors and exporter         |
    // +-----------------------------------------------+

    vectorPtr_Type rhsLap;
    vectorPtr_Type solutionLap;

    rhsLap.reset ( new vector_Type ( uFESpace->map(), Unique ) );
    solutionLap.reset ( new vector_Type ( uFESpace->map(), Unique ) );

    rhsLap->zero();
    solutionLap->zero();

    boost::shared_ptr<laplacianFunctor< Real > >  laplacianSourceFunctor ( new laplacianFunctor< Real >( sourceFunction ) );

    {
        using namespace ExpressionAssembly;

        integrate (
                    elements ( localMeshPtr ),
                    uFESpace->qr(),
                    ETuFESpace,
                    eval(laplacianSourceFunctor, X) * phi_i
                )
                >> rhsLap;

    }

    // Setting exporter
    ExporterHDF5< mesh_Type > exporter ( dataFile, "exporter" );
    exporter.setMeshProcId( localMeshPtr, Comm->MyPID() );
    exporter.setPrefix( "mesh_volume_subdivision_laplace" );
    exporter.setPostDir( "./" );

    exporter.addVariable ( ExporterData< mesh_Type >::ScalarField, "temperature", uFESpace, solutionLap, UInt ( 0 ) );

    // +-----------------------------------------------+
    // |                  Setting BCs                  |
    // +-----------------------------------------------+

    if (verbose)
    {
        std::cout << "-- Setting boundary conditions... " << std::flush;
    }

    BCHandler bcHandler;

    BCFunctionBase ZeroBC ( zeroFunction );

    bcHandler.addBC( "Back",   BACK,   Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Left",   LEFT,   Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Top",    TOP,    Essential, Full, ZeroBC, 1 );


    bcHandler.addBC( "Front",  FRONT,  Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Right",  RIGHT,  Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Bottom", BOTTOM, Essential, Full, ZeroBC, 1 );

    bcHandler.bcUpdate( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

    matrixMeshSubFull->globalAssemble();
    rhsLap->globalAssemble();

    bcManage ( *matrixMeshSubFull, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );

    if (verbose)
    {
        std::cout << " done " << std::endl;
    }

    // +-----------------------------------------------+
    // |       Setting solver and preconditioner       |
    // +-----------------------------------------------+

    typedef LinearSolver::SolverType                                        solver_Type;

    typedef LifeV::Preconditioner                                           basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>                                basePrecPtr_Type;
    typedef PreconditionerIfpack                                            prec_Type;
    typedef boost::shared_ptr<prec_Type>                                    precPtr_Type;


    if (verbose)
    {
        std::cout << "-- Setting up the solver ... " << std::flush;
    }

    LinearSolver linearSolver ( Comm );
    linearSolver.setOperator ( matrixMeshSubFull );

    Teuchos::RCP< Teuchos::ParameterList > aztecList = Teuchos::rcp ( new Teuchos::ParameterList );
    aztecList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

    linearSolver.setParameters ( *aztecList );

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    linearSolver.setPreconditioner ( precPtr );


    if (verbose)
    {
        std::cout << " done" << std::endl;
    }


    // +-----------------------------------------------+
    // |              Solving problem                  |
    // +-----------------------------------------------+

    if (verbose)
    {
        std::cout << "-- Solving problem ... " << std::flush;
    }

    linearSolver.setRightHandSide( rhsLap );
    linearSolver.solve( solutionLap );

    exporter.postProcess( 0.0 );

#ifdef TESTMESHSUB

    matrixPtr_Type standardSystemMatrix1( new matrix_Type( ETuFESpace->map(), 100 ) );
    matrixPtr_Type standardSystemMatrix2( new matrix_Type( ETuFESpace->map(), 100 ) );
    matrixPtr_Type standardSystemMatrix3( new matrix_Type( ETuFESpace->map(), 100 ) );
    matrixPtr_Type standardSystemMatrix4( new matrix_Type( ETuFESpace->map(), 100 ) );
    matrixPtr_Type systemMatrixTotal( new matrix_Type( ETuFESpace->map(), 100 ) );

    boost::shared_ptr<laplacianFunctor< Real > >  laplacianDiffusionFunctor1 ( new laplacianFunctor< Real >( diffusion4Function1 ) );
    boost::shared_ptr<laplacianFunctor< Real > >  laplacianDiffusionFunctor2 ( new laplacianFunctor< Real >( diffusion4Function2 ) );
    boost::shared_ptr<laplacianFunctor< Real > >  laplacianDiffusionFunctor3 ( new laplacianFunctor< Real >( diffusion4Function3 ) );
    boost::shared_ptr<laplacianFunctor< Real > >  laplacianDiffusionFunctor4 ( new laplacianFunctor< Real >( diffusion4Function4 ) );
    boost::shared_ptr<laplacianFunctor< Real > >  laplacianDiffusionFunctorTotal ( new laplacianFunctor< Real >( diffusion4FunctionTotal ) );

    {
        using namespace ExpressionAssembly;
        integrate (
                elements ( localMeshPtr ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                eval(laplacianDiffusionFunctor1, X) * dot ( grad ( phi_i ) , grad ( phi_j ) )
        )
                >> standardSystemMatrix1;

        integrate (
                elements ( localMeshPtr ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                eval(laplacianDiffusionFunctor2, X) * dot ( grad ( phi_i ) , grad ( phi_j ) )
        )
                >> standardSystemMatrix2;

        integrate (
                elements ( localMeshPtr ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                eval(laplacianDiffusionFunctor3, X) * dot ( grad ( phi_i ) , grad ( phi_j ) )
        )
                >> standardSystemMatrix3;

        integrate (
                elements ( localMeshPtr ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                eval(laplacianDiffusionFunctor4, X) * dot ( grad ( phi_i ) , grad ( phi_j ) )
        )
            >> standardSystemMatrix4;

        integrate (
                elements ( localMeshPtr ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                eval(laplacianDiffusionFunctorTotal, X) * dot ( grad ( phi_i ) , grad ( phi_j ) )
        )
                >> systemMatrixTotal;

    }

    std::stringstream convertNumProc, convertProc;

    convertNumProc << Comm->NumProc();
    convertProc << Comm->MyPID();

    bcManage ( *standardSystemMatrix1, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    bcManage ( *standardSystemMatrix2, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    bcManage ( *standardSystemMatrix3, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    bcManage ( *standardSystemMatrix4, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    bcManage ( *systemMatrixTotal, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );

    bcManage ( *systemMatrixMeshSub1, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    bcManage ( *systemMatrixMeshSub2, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    bcManage ( *systemMatrixMeshSub3, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
    bcManage ( *systemMatrixMeshSub4, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );

    standardSystemMatrix1->globalAssemble();
    standardSystemMatrix2->globalAssemble();
    standardSystemMatrix3->globalAssemble();
    standardSystemMatrix4->globalAssemble();
    systemMatrixTotal->globalAssemble();

    /*
    standardSystemMatrix1->spy("matrixClassic1" + convertNumProc.str()  );
    standardSystemMatrix2->spy("matrixClassic2" + convertNumProc.str() );
    standardSystemMatrix3->spy("matrixClassic3" + convertNumProc.str()  );
    standardSystemMatrix4->spy("matrixClassic4" + convertNumProc.str() );
    systemMatrixMeshSub1->spy("matrixMeshSub1" + convertNumProc.str()) ;
    systemMatrixMeshSub2->spy("matrixMeshSub2" + convertNumProc.str() );
    systemMatrixMeshSub3->spy("matrixMeshSub3" + convertNumProc.str()) ;
    systemMatrixMeshSub4->spy("matrixMeshSub4" + convertNumProc.str() );

    systemMatrixTotal->spy("matrixClassicFull" + convertNumProc.str() );
    matrixMeshSubFull->spy("matrixMeshSubFull" + convertNumProc.str() );
	*/

    solutionLap->zero();

    linearSolver.setOperator ( systemMatrixTotal );
    linearSolver.setPreconditioner ( precPtr );
    linearSolver.solve( solutionLap );

    exporter.postProcess( 1.0 );

#endif





    exporter.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}
