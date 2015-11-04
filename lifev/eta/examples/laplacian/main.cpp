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

    @author Antonello Gerbi <antonello.gerbi@epfl.ch>
    @maintainer Niccolo' Dal Santo <niccolo.dalsanto@epfl.ch>
    @date 10-11-2014

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

using namespace LifeV;

// Dirichlet BC functions
Real nonZeroFunction (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 100.;
}

Real zeroFunction (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.;
}

Real sourceFunction (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    return 3 * M_PI * M_PI * sin( M_PI * x ) * sin( M_PI * y ) * sin( M_PI * z );
}

Real uExactFunction (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    return sin( M_PI * y ) * sin( M_PI * z ) * sin ( M_PI * x );
}

VectorSmall< 3 > uGradExactFunction (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& /*i*/)
{
    VectorSmall< 3 > v;

    v[0] = M_PI * cos( M_PI * x ) * sin( M_PI * y ) * sin( M_PI * z );
    v[1] = M_PI * sin( M_PI * x ) * cos( M_PI * y ) * sin( M_PI * z );
    v[2] = M_PI * sin( M_PI * x ) * sin( M_PI * y ) * cos( M_PI * z );

    return v;
}

// Cube's walls identifiers
const int BACK   = 1;
const int FRONT  = 2;
const int LEFT   = 3;
const int RIGHT  = 4;
const int BOTTOM = 5;
const int TOP    = 6;


int main ( int argc, char** argv )
{

    // MPI initialization
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
    std::shared_ptr< Epetra_Comm > Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    std::shared_ptr< Epetra_Comm > Comm ( new Epetra_SerialComm );
#endif

    // Initializing chronometers
    LifeChrono globalChrono;
    LifeChrono initChrono;

    globalChrono.start();


    // Reading parameters through GetPot
    GetPot command_line ( argc, argv );
    const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
    GetPot dataFile ( dataFileName );

    const bool verbose ( Comm->MyPID() == 0 );

    if ( verbose )
    {
        std::cout << "-- Building and partitioning the mesh... " << std::flush;
    }

    initChrono.start();



    // +-----------------------------------------------+
    // |                Building mesh                  |
    // +-----------------------------------------------+

    typedef RegionMesh< LinearTetra >                                       mesh_Type;

    std::shared_ptr< mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );


    // Building structured mesh (in this case a cube)
    regularMesh3D ( *fullMeshPtr, 0,
                    dataFile( "mesh/nx", 15 ), dataFile( "mesh/ny", 15 ), dataFile( "mesh/nz", 15 ),
                    dataFile ( "mesh/verbose", false ),
                    2.0, 2.0, 2.0, -1.0, -1.0, -1.0 );


    // Partitioning mesh, possibly with overlap
    const UInt overlap ( dataFile( "mesh/overlap", 0 ) );

    std::shared_ptr< mesh_Type > localMeshPtr;

    MeshPartitioner< mesh_Type > meshPart;

    if ( overlap )
    {
        meshPart.setPartitionOverlap ( overlap );
    }

    meshPart.doPartition ( fullMeshPtr, Comm );
    localMeshPtr = meshPart.meshPartition();

    // Clearing global mesh
    fullMeshPtr.reset();

    initChrono.stop();

    if ( verbose )
    {
        std::cout << " done in " << initChrono.diff() << " s!" << std::endl;
    }

    initChrono.reset();
    initChrono.start();



    // +-----------------------------------------------+
    // |         Building FE space and matrix          |
    // +-----------------------------------------------+

    typedef FESpace< mesh_Type, MapEpetra >                                 uSpaceStd_Type;
    typedef std::shared_ptr< uSpaceStd_Type >                             uSpaceStdPtr_Type;

    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >                         uSpaceETA_Type;
    typedef std::shared_ptr< uSpaceETA_Type >                             uSpaceETAPtr_Type;

    typedef FESpace<mesh_Type, MapEpetra>::function_Type                    function_Type;

    if ( verbose )
    {
        std::cout << "-- Building finite elements spaces... " << std::flush;
    }

    // Defining finite elements standard and ET spaces
    uSpaceStdPtr_Type uFESpace ( new uSpaceStd_Type ( localMeshPtr, dataFile( "finite_element/degree", "P1" ), 1, Comm ) );
    uSpaceETAPtr_Type ETuFESpace ( new uSpaceETA_Type ( localMeshPtr, & ( uFESpace->refFE() ), & ( uFESpace->fe().geoMap() ), Comm ) );

    initChrono.stop();

    if ( verbose )
    {
        std::cout << " done in " << initChrono.diff() << " s! (Dofs: " << uFESpace->dof().numTotalDof() << ")" << std::endl;
    }

    typedef Epetra_FECrsGraph                                               graph_Type;
    typedef std::shared_ptr<Epetra_FECrsGraph>                            graphPtr_Type;

    typedef MatrixEpetra< Real >                                            matrix_Type;
    typedef std::shared_ptr< MatrixEpetra< Real > >                       matrixPtr_Type;
    typedef VectorEpetra                                                    vector_Type;
    typedef std::shared_ptr<VectorEpetra>                                 vectorPtr_Type;

    // Declaring the problem's graph and matrix
    graphPtr_Type systemGraph;
    matrixPtr_Type systemMatrix;

    initChrono.reset();
    initChrono.start();

    if ( overlap )
    {
        systemGraph.reset ( new graph_Type ( Copy, * ( uFESpace->map().map( Unique ) ), 50, true ) );
    }
    else
    {
        systemGraph.reset ( new graph_Type ( Copy, * ( uFESpace->map().map( Unique ) ), 50 ) );
    }


    if ( verbose )
    {
        std::cout << "-- Assembling matrix graph... " << std::flush;
    }

    {

      using namespace ExpressionAssembly;

      buildGraph (
            elements ( localMeshPtr ),
            uFESpace->qr(),
            ETuFESpace,
            ETuFESpace,
            dot ( grad ( phi_i ) , grad ( phi_j ) )
            )
            >> systemGraph;

    }

    systemGraph->GlobalAssemble();


    if ( overlap )
    {

      systemMatrix.reset( new matrix_Type ( ETuFESpace->map(), *systemGraph, true ) );

    }
    else
    {

      systemMatrix.reset( new matrix_Type ( ETuFESpace->map(), *systemGraph ) );

    }

    initChrono.stop();

    if ( verbose )
    {
        std::cout << " done in " << initChrono.diff() << " s!" << std::endl;
    }

    // Clearing problem's matrix
    systemMatrix->zero();

    initChrono.reset();
    initChrono.start();

    if ( verbose )
    {
        std::cout << "-- Assembling the Laplace matrix... " << std::flush;
    }

    {
        using namespace ExpressionAssembly;

        integrate (
                elements ( localMeshPtr ),
                uFESpace->qr(),
                ETuFESpace,
                ETuFESpace,
                dot ( grad ( phi_i ) , grad ( phi_j ) )
        )
                >> systemMatrix;
    }

    initChrono.stop();

    if (verbose)
    {
        std::cout << " done in " << initChrono.diff() << " s!" << std::endl;
    }

    // +-----------------------------------------------+
    // |     Initializing vectors and exporter         |
    // +-----------------------------------------------+

    vectorPtr_Type rhsLap;
    vectorPtr_Type solutionLap;

    if ( overlap )
    {

      rhsLap.reset ( new vector_Type ( uFESpace->map(), Unique, Zero ) );
      solutionLap.reset ( new vector_Type ( uFESpace->map(), Unique, Zero ) );

    }
    else
    {

      rhsLap.reset ( new vector_Type ( uFESpace->map(), Unique ) );
      solutionLap.reset ( new vector_Type ( uFESpace->map(), Unique ) );

    }

    rhsLap->zero();
    solutionLap->zero();


    std::shared_ptr<laplacianFunctor< Real > >  laplacianSourceFunctor ( new laplacianFunctor< Real >( sourceFunction ) );

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
    exporter.setPrefix( "laplace" );
    exporter.setPostDir( "./" );

    exporter.addVariable ( ExporterData< mesh_Type >::ScalarField, "temperature", uFESpace, solutionLap, UInt ( 0 ) );


    // +-----------------------------------------------+
    // |                  Setting BCs                  |
    // +-----------------------------------------------+

    if (verbose)
    {
        std::cout << "-- Setting boundary conditions... " << std::flush;
    }
    systemMatrix->globalAssemble();
    rhsLap->globalAssemble();

    initChrono.reset();
    initChrono.start();

    BCHandler bcHandler;

    BCFunctionBase ZeroBC ( zeroFunction );
    BCFunctionBase OneBC ( nonZeroFunction );

    bcHandler.addBC( "Back",   BACK,   Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Left",   LEFT,   Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Top",    TOP,    Essential, Full, ZeroBC, 1 );

    bcHandler.addBC( "Front",  FRONT,  Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Right",  RIGHT,  Essential, Full, ZeroBC, 1 );
    bcHandler.addBC( "Bottom", BOTTOM, Essential, Full, ZeroBC, 1 );

    bcHandler.bcUpdate( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );


    bcManage ( *systemMatrix, *rhsLap, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );

    initChrono.stop();

    if (verbose)
    {
        std::cout << " done in " << initChrono.diff() << " s!" << std::endl;
    }

    // +-----------------------------------------------+
    // |       Setting solver and preconditioner       |
    // +-----------------------------------------------+

    typedef LinearSolver::SolverType                                        solver_Type;

    typedef LifeV::Preconditioner                                           basePrec_Type;
    typedef std::shared_ptr<basePrec_Type>                                basePrecPtr_Type;
    typedef PreconditionerIfpack                                            prec_Type;
    typedef std::shared_ptr<prec_Type>                                    precPtr_Type;


    if (verbose)
    {
        std::cout << "-- Setting up the solver ... " << std::flush;
    }

    initChrono.reset();
    initChrono.start();

    LinearSolver linearSolver ( Comm );
    linearSolver.setOperator ( systemMatrix );

    Teuchos::RCP< Teuchos::ParameterList > aztecList = Teuchos::rcp ( new Teuchos::ParameterList );
    aztecList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

    linearSolver.setParameters ( *aztecList );

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    linearSolver.setPreconditioner ( precPtr );

    initChrono.stop();

    if (verbose)
    {
        std::cout << " done in " << initChrono.diff() << " s!" << std::endl;
    }


    // +-----------------------------------------------+
    // |              Solving problem                  |
    // +-----------------------------------------------+

    if (verbose)
    {
        std::cout << "-- Solving problem ... " << std::flush;
    }

    initChrono.reset();
    initChrono.start();

    linearSolver.setRightHandSide( rhsLap );
    linearSolver.solve( solutionLap );

    initChrono.stop();


    function_Type uEx(uExactFunction);

    // Save exact solution
    vectorPtr_Type uExPtrInterpolated( new vector_Type (uFESpace->map(), Repeated) );

    uExPtrInterpolated->zero();

    uFESpace->interpolate (uEx, *uExPtrInterpolated, 0.0);

    exporter.addVariable ( ExporterData< mesh_Type >::ScalarField, "realTemperature", uFESpace, uExPtrInterpolated, UInt ( 0 ) );

    // using general functors

    Real L2ErrorLap = 0.0;
    Real TotL2ErrorLap = 0.0;

    Real H1SeminormLap = 0.0;
    Real TotH1SeminormLap = 0.0;

    std::shared_ptr<laplacianFunctor< Real > >  laplacianExactFunctor ( new laplacianFunctor< Real >( uExactFunction ) );
    std::shared_ptr<laplacianFunctor< VectorSmall<3> > >  laplacianExactGradientFunctor ( new laplacianFunctor< VectorSmall<3> >( uGradExactFunction ) );

    {
        using namespace ExpressionAssembly;

        integrate (
                    elements ( localMeshPtr ),
                    uFESpace->qr(),
                    ( eval(laplacianExactFunctor, X) - value (ETuFESpace, *solutionLap) )
                    * ( eval(laplacianExactFunctor, X) - value (ETuFESpace, *solutionLap) )
                ) >> L2ErrorLap;

    }

    {
        using namespace ExpressionAssembly;

        integrate (
                    elements ( localMeshPtr ),
                    uFESpace->qr(),
                    dot( eval(laplacianExactGradientFunctor, X) - grad( ETuFESpace, *solutionLap ),
                         eval(laplacianExactGradientFunctor, X) - grad( ETuFESpace, *solutionLap ) )
                ) >> H1SeminormLap;

    }

    Comm->Barrier();
    Comm->SumAll (&L2ErrorLap, &TotL2ErrorLap, 1);
    Comm->SumAll (&H1SeminormLap, &TotH1SeminormLap, 1);

    if (verbose)
    {
        std::cout << "TotError General in L2 norm is " << sqrt( TotL2ErrorLap ) << std::endl;
        std::cout << "TotError General in H1 norm is " << sqrt( TotL2ErrorLap + TotH1SeminormLap ) << std::endl;

    }

    if (verbose)
    {
        std::cout << " done in " << initChrono.diff() << " s!" << std::endl;
    }

    exporter.postProcess( 0 );
    exporter.closeFile();

    globalChrono.stop();

    if (verbose)
    {
        std::cout << std::endl << "Problem solved in  " << globalChrono.diff() << " s!" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}
