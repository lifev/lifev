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
    @brief

    @author Claudia Colciago <claudia.colciago@epfl.ch>
    @date 08-11-2012
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

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/fem/BCManage.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEnsight.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>

#include "ETRobinMembraneSolver.hpp"
#include "ud_functions.hpp"

#include <iostream>

#define FLUX 1

#define OUTLET 3
#define INLET 2
#define WALL 1
#define RING 20
#define RING2 30
#define OLDRING 40
#define OLDARTERY 4

using namespace LifeV;


const Real PI=3.141592653589793;


struct ETRobinMembraneSolver::Private
{
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    std::string    data_file_name;

    // Data for the fluid
    Real         Re;
    Real         nu;         /**< viscosity (in m^2/s) */
    Real         H;          /**< height and width of the domain (in m) */
    Real         D;          /**< diameter of the cylinder (in m) */
    Real         density;
    Real         R;          //radius


    //Data for the membrnae
    Real rhos;
    Real Hs;
    Real ni;
    Real E;


    //Data of the specific problem
    std::string initial_sol;
    Real numLM;
    boost::shared_ptr<Epetra_Comm>   comm;

};

ETRobinMembraneSolver::ETRobinMembraneSolver( int argc, char** argv )
    :
    M_d( new Private )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    M_d->data_file_name  = data_file_name;

    M_d->Re              = dataFile( "fluid/physics/Re", 1. );
    M_d->nu              = dataFile( "fluid/physics/viscosity", 1. );
    M_d->H               = dataFile( "fluid/physics/H", 20. );
    M_d->D               = dataFile( "fluid/physics/D", 1. );
    M_d->R               = M_d->D/2;
    M_d->density         = dataFile( "fluid/physics/density", 1. );

    M_d->rhos            = dataFile( "membrane/physics/density_sol", 1. );
    M_d->E               = dataFile( "membrane/physics/young", 1. );
    M_d->ni              = dataFile( "membrane/physics/poisson", 1. );
    M_d->Hs              = dataFile( "membrane/physics/wall_thickness", 1. );


    M_d->initial_sol = (std::string) dataFile( "fluid/problem/initial_sol", "none");
    M_d->numLM       = dataFile( "fluid/problem/numLM"      , 1     );

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    int ntasks = 0;
    M_d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    if (!M_d->comm->MyPID())
	{
	    std::cout << "My PID            = " << M_d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
	    std::cout << "fluid density     = " << M_d->density       << std::endl
		      << "nu                = " << M_d->nu            << std::endl
		      << "poisson           = " << M_d->ni            << std::endl
		      << "initial solution  = " << M_d->initial_sol   << std::endl
		      << "Young modulus     = " << M_d->E             << std::endl
		      << "Wall thickness    = " << M_d->Hs            << std::endl;

	}

#else
    M_d->comm.reset( new Epetra_SerialComm() );
#endif

}

void
ETRobinMembraneSolver::run()

{

    bool verbose = (M_d->comm->MyPID() == 0);

    //------------Creating file to store data------------------------------

    GetPot dataFile( M_d->data_file_name );

    boost::shared_ptr<OseenData> oseenData(new OseenData());
    oseenData->setup( dataFile );

    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");

    Real LameI = (M_d->Hs * M_d->E * M_d->ni )/( ( 1 - M_d->ni*M_d->ni ) );
    Real LameII = M_d->Hs * M_d->E / ( 2 * ( 1 + M_d->ni ) );

    if (verbose) std::cout<<" LameI : "<<LameI<<std::endl;
    if (verbose) std::cout<<" LameII : "<<LameII<<std::endl;

    //-----------------------------------------------------------------------

    //----------Creating Mesh and FESpaces-----------------------------------

    boost::shared_ptr< mesh_type > fullMeshPtr (new mesh_type);
    readMesh(*fullMeshPtr, meshData);

    MeshPartitioner< mesh_type >   meshPart(fullMeshPtr, M_d->comm);

    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");
    std::string mOrder =  dataFile( "mesh_motion/space_discretization/m_order", "P1");

    if (verbose)
        std::cout << "Building the FE spaces ... " << std::flush;

    M_uFESpace.reset( new FESpace< mesh_type, MapEpetra >( meshPart, uOrder, 3,  M_d->comm ) );
    M_pFESpace.reset( new FESpace< mesh_type, MapEpetra >( meshPart, pOrder, 1,  M_d->comm ) );

    M_ETuFESpace.reset( new ETFESpace< mesh_type, MapEpetra, 3, 3 >( meshPart, &(M_uFESpace->refFE()), M_d->comm ) );
    M_ETpFESpace.reset( new ETFESpace< mesh_type, MapEpetra, 3, 1 >( meshPart, &(M_pFESpace->refFE()), M_d->comm ) );

    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalVelDof   = M_uFESpace->map().map( Unique )->NumGlobalElements();
    UInt totalPressDof = M_pFESpace->map().map( Unique )->NumGlobalElements();

    if (verbose) std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;
    if (verbose) std::cout << "Total FLux DOF = " << M_d->numLM << std::endl;

    MapEpetra fluidMap( M_uFESpace->map() );
    MapEpetra fluxMap ( M_d->numLM, M_d->comm );

#ifdef FLUX
   MapEpetra fullMap ( M_uFESpace->map() + M_pFESpace->map() + fluxMap );
#else
   MapEpetra fullMap ( M_uFESpace->map() + M_pFESpace->map() );
#endif

   vectorPtr_type checkVector(new vector_type(M_uFESpace->map(), Unique));

    DOFInterface3Dto3D interfaceDOF(M_uFESpace->refFE(), M_uFESpace->dof());
    interfaceDOF.update(*M_uFESpace->mesh(), WALL, *M_uFESpace->mesh(), WALL, 0);
    createInterfaceMap( checkVector ,  meshPart.meshPartition() ,  M_uFESpace->dof());

    //------------------FLUID PROBLEM--------------------------------

    if (verbose) std::cout << "Calling the solver constructors ... ";

    SolverAztecOO NSSolver;
    NSSolver.setCommunicator(M_d->comm);
    NSSolver.setDataFromGetPot(dataFile,"solver");
    NSSolver.setupPreconditioner(dataFile,"prec");

    if (verbose) std::cout << "done." << std::endl;

    //---------------------------------------------------------

    //--------------------Creating BDF objects---------------------------

    if (verbose) std::cout << "Calling the TimeAdvance constructors ... ";

    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();

    TimeAdvanceBDFNavierStokes<vector_type> fluidTimeAdvance;
    fluidTimeAdvance.setup(oseenData->dataTimeAdvance()->orderBDF());

    TimeAdvanceBDF<vector_type> dispTimeAdvance;
    dispTimeAdvance.setup(oseenData->dataTimeAdvance()->orderBDF());

    //TimeAdvance for the fluid
    fluidTimeAdvance.bdfVelocity().setTimeStep(oseenData->dataTime()->timeStep());
    if(verbose)
	fluidTimeAdvance.showMe();

    //TimeAdvance for the displacement
    dispTimeAdvance.setTimeStep(oseenData->dataTime()->timeStep());
    if(verbose)
	dispTimeAdvance.showMe();

    Real alpha = fluidTimeAdvance.bdfVelocity().coefficientFirstDerivative(t0);

    if(verbose) std::cout<<" done. \n";

    //---------------------------------------------------------

    //--------------Boundary Conditions-------------------------

    if (verbose) std::cout << "Calling the BCHandler constructor ... ";

    //creating BCHandler objects
    BCHandler bcHFluid;
    BCFunctionBase flow_in (linearVelInletCylinder);
    //BCFunctionBase flux_in (linearPontdist);
    BCFunctionBase flux_in (linearInletCylinder);
    //BCFunctionBase flow_out (linearPopliteal);
    BCFunctionBase uZero( fZero );

#ifdef FLUX
    bcHFluid.addBC("InFlow" ,   INLET,    Flux,                Normal, flux_in    );
    bcHFluid.addBC("InFlow_2" , OUTLET,   Natural,             Full,   uZero,   3 );
    bcHFluid.addBC("Old" ,   OLDARTERY,    Essential,           Full,   uZero, 3 );
    //In order to have a well defined problem you need to set conditions on boundaries of boundaries
    bcHFluid.addBC("Ring" ,     RING,     EssentialEdges,   Full,   uZero,   3 );
    bcHFluid.addBC("Ring" ,     OLDRING,     EssentialEdges,   Full,   uZero,   3 );
    bcHFluid.addBC("Ring3" ,    RING2,    EssentialEdges,   Full,   uZero,   3 );
    bcHFluid.setOffset("InFlow", totalVelDof + totalPressDof);

#else
    bcHFluid.addBC("InFlow" ,   INLET,    Essential,           Full,   flow_in, 3 );
    bcHFluid.addBC("InFlow_2" , OUTLET,   Natural,             Full,   uZero,   3 );
    //In order to have a well defined problem you need to set conditions on boundaries of boundaries
    bcHFluid.addBC("Ring" ,     RING,     EssentialVertices,   Full,   uZero,   3 );
    bcHFluid.addBC("Ring3" ,    RING2,    EssentialVertices,   Full,   uZero,   3 );
#endif

    if(verbose)
	std::cout<<" BC done \n";

    //----------------------------------------------------------------------------

#ifdef FLUX
     vector_block_type NSSolution(M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap, Unique);
#else
    vector_block_type NSSolution(M_ETuFESpace->map() | M_ETpFESpace->map(), Unique);
#endif

    vector_type velocitySolution(M_ETuFESpace->map(),Repeated);
    vector_type dispSolution(M_ETuFESpace->map(),Repeated);

    NSSolution *= 0.0;

    //-------------Creating Exporter Object------------------------------------------------

    if (verbose) std::cout << "Calling the Exporter constructor ... ";

    std::string const exporterFileName    =  dataFile( "exporter/filename", "cylinder");
    UInt  exportEach = dataFile("exporter/each",1);

    ExporterHDF5< mesh_type > exporter( dataFile, meshPart.meshPartition(), exporterFileName, M_d->comm->MyPID());

    vectorPtr_type velAndPressureExporter ( new vector_type(NSSolution,Repeated ));
    vectorPtr_type dispExporter (new vector_type(dispSolution,Repeated ));

    exporter.addVariable( ExporterData< mesh_type >::VectorField, "f-velocity",
			  M_uFESpace, velAndPressureExporter, UInt(0) );

    exporter.addVariable( ExporterData< mesh_type >::ScalarField, "f-pressure",
			  M_pFESpace, velAndPressureExporter, UInt(3*M_uFESpace->dof().numTotalDof() ) );

    exporter.addVariable( ExporterData< mesh_type >::VectorField, "s-displacement",
			  M_uFESpace, dispExporter, UInt(0) ) ;
    if(verbose)
	std::cout<<" done. \n";

    //-----------------------------------------------------------------------------

    //--------------------------Initialization--------------------------------------



    if (M_d->initial_sol == "restart")
	{
	    ASSERT(oseenData->dataTimeAdvance()->orderBDF() == 1, "Restart works only for BDF1");

	    if (verbose) std::cout << std::endl;
	    if (verbose) std::cout << "Restoring the previous solution ... " << std::endl;

	    std::string filename = dataFile("importer/filename", "cylinder");

	    LifeV::ExporterHDF5<mesh_type> importer( dataFile, filename);
	    importer.setMeshProcId(M_uFESpace->mesh(), M_d->comm->MyPID());

	    importer.addVariable( ExporterData<mesh_type>::VectorField,
				  "f-velocity",
				  M_uFESpace,
				  velAndPressureExporter,
				  UInt ( 0 ) );

	    importer.addVariable( ExporterData<mesh_type>::ScalarField,
				  "f-pressure",
				  M_pFESpace,
				  velAndPressureExporter,
				  3*M_uFESpace->dof().numTotalDof() );

	    importer.addVariable( ExporterData<mesh_type>::VectorField,
				  "s-displacement",
				  M_uFESpace,
				  dispExporter,
				  UInt( 0 ) );

	    exporter.setTimeIndex(importer.importFromTime(t0));

	    Real norm = velAndPressureExporter->norm2();
	    if (verbose)
		std::cout << "   f- restart solution norm = " << norm << std::endl;


	}

    fluidTimeAdvance.bdfVelocity().setInitialCondition( *(velAndPressureExporter) );
    dispTimeAdvance.setInitialCondition( dispSolution );

    velocitySolution = *velAndPressureExporter;
    dispSolution = *dispExporter;

    fluidTimeAdvance.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
    dispTimeAdvance.updateRHSContribution(oseenData->dataTime()->timeStep());

    if(verbose)
	std::cout<<" done \n";

    exporter.postProcess( t0 );

    //-------------------------------------------------------------------------------------

#ifdef FLUX
    boost::shared_ptr<matrix_block_type> NSMatrixConstant(new matrix_block_type( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap ) );
    *NSMatrixConstant *= 0.0;

#else
    boost::shared_ptr<matrix_block_type> NSMatrixConstant(new matrix_block_type( M_ETuFESpace->map() | M_ETpFESpace->map() ) );
    *NSMatrixConstant *= 0.0;
#endif


    //----------------------Temporal Loop----------------------------------

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << " ### Simulation times ### " << std::endl;
    if (verbose) std::cout << " From " << t0 << " to " << tFinal << std::endl;
    if (verbose) std::cout << " Time step: " << dt << std::endl;
    if (verbose) std::cout << std::endl;

    Real currentTime(t0);
    UInt niter(0);

    *NSMatrixConstant *=0.0;

    QuadratureBoundary myBDQR(buildTetraBDQR(quadRuleTria4pt));

    {

	using namespace ExpressionAssembly;

	integrate(
		  elements(M_ETuFESpace->mesh()), // Mesh

		  M_uFESpace->qr(), // QR

		  M_ETuFESpace,
		  M_ETuFESpace,

		  // Intertial Term
		  M_d->density
		  * ( value(alpha/dt) * dot(phi_i,phi_j) )

		  // Viscous Term
		  + M_d->nu * dot(grad(phi_i) , grad(phi_j))

		  )
	    >> NSMatrixConstant->block(0,0);

	integrate(
		  elements(M_ETuFESpace->mesh()), // Mesh

		  M_uFESpace->qr(), // QR

		  M_ETuFESpace,
		  M_ETpFESpace,

		  value(-1.0)*phi_j*div(phi_i)


		  )
	    >> NSMatrixConstant->block(0,1);

	integrate(
		  elements(M_ETuFESpace->mesh()), // Mesh

		  M_uFESpace->qr(), // QR

		  M_ETpFESpace,
		  M_ETuFESpace,

		  value(-1.0)*phi_i*div(phi_j)

		  )
	    >> NSMatrixConstant->block(1,0);


	integrate( boundary(M_ETuFESpace->mesh(),WALL),
		   myBDQR,

		   M_ETuFESpace,
		   M_ETuFESpace,

		   //Boundary Mass
		   value(M_d->rhos * M_d->Hs * alpha / dt) * dot(phi_j,phi_i)

		   )
	    >> NSMatrixConstant->block(0,0);

    }


    MatrixSmall<3,3> Eye;
    Eye *= 0.0;
    Eye[0][0]=1;
    Eye[1][1]=1;
    Eye[2][2]=1;

    {
	using namespace ::LifeV::ExpressionAssembly;


	integrate( boundary(M_ETuFESpace->mesh(),WALL),
		   myBDQR,

		   M_ETuFESpace,
		   M_ETuFESpace,

		   //Boundary Stiffness
		   ( dt / alpha ) *
		   2  * LameII  *
		   0.5 * dot( ( grad(phi_j) + (-1) * grad(phi_j) * outerProduct( Nface, Nface ) )
			      + transpose(grad(phi_j) + (-1) * grad(phi_j) * outerProduct( Nface, Nface )),
			      ( grad(phi_i) + ( (-1) * grad(phi_i) * outerProduct( Nface, Nface ) ) ) ) +

		   ( dt / alpha ) *
		   LameI  * dot( value( Eye ) , ( grad(phi_j) + (-1)*  grad(phi_j) * outerProduct( Nface, Nface ) ) )
		   * dot( value( Eye ) ,  ( grad(phi_i) + (-1)* grad(phi_i) * outerProduct( Nface, Nface ) )  )

		   )
	    >> NSMatrixConstant->block(0,0);

    }

	NSMatrixConstant->globalAssemble();



	while ( currentTime < tFinal)
	    {
		LifeChrono ChronoIteration;
		ChronoIteration.start();

		currentTime += dt;
		niter += 1;

		if (verbose) std::cout << std::endl;
		if (verbose) std::cout << "----------------------------" << std::endl;
		if (verbose) std::cout << " Time : " << currentTime << std::endl;
		if (verbose) std::cout << " Iter : " << niter << std::endl;
		if (verbose) std::cout << "----------------------------" << std::endl;
		if (verbose) std::cout << std::endl;

		vector_type velocityExtrapolated(velocitySolution,Repeated);
		fluidTimeAdvance.bdfVelocity().extrapolation( velocityExtrapolated );
		vector_type velocityBdfRHS(velocitySolution,Repeated);
		velocityBdfRHS = fluidTimeAdvance.bdfVelocity().rhsContributionFirstDerivative();

		vector_type dispExtrapolated(dispSolution,Repeated);
		dispTimeAdvance.extrapolation( dispExtrapolated );
		vector_type dispBdfRHS(fluidMap,Repeated);
		dispBdfRHS = dispTimeAdvance.rhsContributionFirstDerivative();

#ifdef FLUX
		boost::shared_ptr<matrix_block_type> NSMatrix(new matrix_block_type( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap ) );
		*NSMatrix *= 0.0;

#else

		boost::shared_ptr<matrix_block_type> NSMatrix(new matrix_block_type( M_ETuFESpace->map() | M_ETpFESpace->map() ) );
		*NSMatrix *= 0.0;

#endif

	    if (verbose) std::cout << "[Navier-Stokes] Assembling the matrix ... " << std::flush;

	    LifeChrono ChronoItem;
	    ChronoItem.start();

	    {

		using namespace ::LifeV::ExpressionAssembly;

		integrate(
			  elements(M_ETuFESpace->mesh()), // Mesh
			  M_uFESpace->qr(), // QR

			  M_ETuFESpace,
			  M_ETuFESpace,

			  // Advection Term
			  M_d->density
			  * ( dot(grad(phi_j) * value(M_ETuFESpace,velocityExtrapolated) , phi_i) )

			  )
		    >> NSMatrix->block(0,0);
	    }


	    ChronoItem.stop();
	    if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


	    if (verbose) std::cout << "[Navier-Stokes] Adding constant parts ... " << std::flush;
	    ChronoItem.start();

	    *NSMatrix += *NSMatrixConstant;

	    ChronoItem.stop();
	    if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


	    if (verbose) std::cout << "[Navier-Stokes] Assembling the rhs ... " << std::flush;
	    ChronoItem.start();

#ifdef FLUX
	    vector_block_type NSRhs( M_ETuFESpace->map() | M_ETpFESpace->map() | fluxMap, Repeated );
#else
	    vector_block_type NSRhs( M_ETuFESpace->map() | M_ETpFESpace->map(), Repeated );
#endif
	    NSRhs *= 0.0;

	    {
		using namespace ExpressionAssembly;

		integrate(
			  elements(M_ETuFESpace->mesh()), // Mesh

			  M_uFESpace->qr(), // QR

			  M_ETuFESpace,

			  //Inertial Term
			  M_d->density * dot(value(M_ETuFESpace,velocityBdfRHS), phi_i )

			  )
		    >> NSRhs.block(0);

	    }

	    ChronoItem.stop();
	    if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;

	    if (verbose) std::cout << "[Navier-Stokes] Boundary Intergrals in the rhs ... " << std::flush;
	    ChronoItem.start();

	    {
		using namespace ExpressionAssembly;

		integrate( boundary(M_ETuFESpace->mesh(),WALL),
			   myBDQR,

			   M_ETuFESpace,

			   //BOUNDARY STIFFNESS
			   value(-1)*( dt / alpha ) *
			   2  * LameII  *
			   0.5 * dot( ( grad(M_ETuFESpace, dispBdfRHS) + (-1) * grad(M_ETuFESpace, dispBdfRHS) * outerProduct( Nface, Nface ) )
				      + transpose(grad(M_ETuFESpace, dispBdfRHS) + (-1) * grad(M_ETuFESpace, dispBdfRHS) * outerProduct( Nface, Nface )),
				      ( grad(phi_i) + ( (-1) * grad(phi_i) * outerProduct( Nface, Nface ) ) ) ) +

			   value(-1)*( dt / alpha ) *
			   LameI  * dot( value( Eye ) , ( grad(M_ETuFESpace, dispBdfRHS) + (-1)*  grad(M_ETuFESpace, dispBdfRHS) * outerProduct( Nface, Nface ) ) )
			   * dot( value( Eye ) ,  ( grad(phi_i) + (-1)* grad(phi_i) * outerProduct( Nface, Nface ) )  )
			   ) >> NSRhs.block(0);

		integrate(
			  boundary(M_ETuFESpace->mesh(),WALL), // Mesh

			  myBDQR, // QR

			  M_ETuFESpace,

			  //BOUNDARY MASS
			  ( M_d->Hs * M_d->rhos * alpha) * dot(value(M_ETuFESpace,velocityBdfRHS), phi_i )

			  )
		    >> NSRhs.block(0);

	    }

	    ChronoItem.stop();
	    if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;

	    if (verbose) std::cout << "[Navier-Stokes] Closing the matrix and the rhs ... " << std::flush;

	    ChronoItem.start();

	    NSMatrix->globalAssemble();

	    NSRhs.globalAssemble();
	    vector_block_type NSRhsUnique( NSRhs, Unique );

	    ChronoItem.stop();
	    if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


	    if (verbose) std::cout << "[Navier-Stokes] Applying boundary conditions ... " << std::flush;

	    bcHFluid.bcUpdate( *meshPart.meshPartition(), M_uFESpace->feBd(), M_uFESpace->dof() );
	    bcManage(*NSMatrix, NSRhsUnique,
		     *M_uFESpace->mesh(), M_uFESpace->dof(),
		     bcHFluid, M_uFESpace->feBd(), 1.0, currentTime);

	    ChronoItem.stop();
	    if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;

	    if (verbose) std::cout << "[Navier-Stokes] Solving the system " << std::endl;

	    NSSolver.setMatrix(*NSMatrix);

	    boost::shared_ptr<matrix_type> NSMatrixNoBlock(new matrix_type( NSMatrix->matrixPtr() ));

	    NSSolver.solveSystem(NSRhsUnique,NSSolution,NSMatrixNoBlock);

	    NSSolution.spy("solution");

	    if (verbose) std::cout << "[Navier-Stokes] Time advancing ... " << std::flush;

	    ChronoItem.start();
	    fluidTimeAdvance.bdfVelocity().shiftRight( NSSolution );

	    velocitySolution.subset(NSSolution);

	    vector_type dispInterface(M_interfaceMap, Repeated);

	    dispInterface = dt/alpha * ( velocitySolution + dispBdfRHS );
	    dispInterface.spy("dispInterface");

	    dispSolution *= 0.0;


	    dispSolution.subset(dispInterface, *M_interfaceMap, 0, 0);

	    dispTimeAdvance.shiftRight( dispSolution );

	    fluidTimeAdvance.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
	    dispTimeAdvance.updateRHSContribution(oseenData->dataTime()->timeStep());

	    ChronoItem.stop();
	    if (verbose) std::cout << ChronoItem.diff() << " s" << std::endl;


	    if (verbose) std::cout << std::endl;

	    if (verbose) std::cout << " Exporting " << std::endl;


	    *velAndPressureExporter = NSSolution;
	    *dispExporter = dispSolution;
	    //*dispExporter += 100;

	    dispExporter->spy("dispExporter");

	    if (niter%exportEach == 0)
		{
		    exporter.postProcess(currentTime);
		}

	    ChronoIteration.stop();
	    if (verbose) std::cout << std::endl << " Total iteration time : " << ChronoIteration.diff() << " s" << std::endl;

	} // end time loop


    exporter.closeFile();

}


// void ETRobinMembraneSolver::createInterfaceMap( std::map<ID, ID> const& locDofMap, const DOF& dof )
// {
//     Displayer disp(M_d->comm);
//     disp.leaderPrint("Building the Interface Map ...             ");

//     std::vector<int> dofInterfaceFluid;

//     typedef std::map<ID, ID>::const_iterator iterator_Type;

//     //std::map<ID, ID> const& locDofMap = M_dofStructureToHarmonicExtension->locDofMap();
//     dofInterfaceFluid.reserve( locDofMap.size() );

//     for (UInt dim = 0; dim < nDimensions; ++dim)
// 	for ( iterator_Type i = locDofMap.begin(); i != locDofMap.end(); ++i )
// 	    dofInterfaceFluid.push_back(i->second + dim * dof.numTotalDof()); // in solid numerotation

//     int* pointerToDofs(0);
//     if (dofInterfaceFluid.size() > 0)
// 	pointerToDofs = &dofInterfaceFluid[0];

//     M_interfaceMap.reset( new MapEpetra( -1,
// 					 static_cast<int>(dofInterfaceFluid.size()),
// 					 pointerToDofs,
// 					 M_d->comm ));
//     disp.leaderPrint("done\n");
//     M_d->comm->Barrier();

// }

void ETRobinMembraneSolver::createInterfaceMap( vectorPtr_type checkVector, meshPtr_type& mesh, const DOF& dof)
{
    std::set<ID> GID_nodes;
    for (UInt dim = 0; dim < nDimensions; ++dim)
	{

	    for ( ID iBoundaryElement = 0 ; iBoundaryElement < mesh->numBoundaryFacets(); ++iBoundaryElement )
		{
		    if( mesh->boundaryFacet( iBoundaryElement ).markerID() == WALL )
			{

			    std::vector<ID> localToGlobalMapOnBElem = dof.localToGlobalMapOnBdFacet(iBoundaryElement);

			    for (ID lDof = 0; lDof< localToGlobalMapOnBElem.size(); lDof++)
				{
				    ID gDof = dof.localToGlobalMapByBdFacet( iBoundaryElement, lDof);
				    GID_nodes.insert(gDof +  dim * dof.numTotalDof());
				}
			}
		}

	    for ( UInt i = 0; i < mesh->numVertices(); i++ )
		{
		    if( mesh->point(i).markerID() == WALL )
			{
			    GID_nodes.insert(mesh->point(i).id() +  dim * dof.numTotalDof());
			}
		}
	}

    int LocalNodesNumber = GID_nodes.size();
    int TotalNodesNumber = 0;

    MPI_Allreduce(&LocalNodesNumber, &TotalNodesNumber, 1,  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int * GlobalID = new int[LocalNodesNumber];
    int k = 0;
    for(std::set<ID>::iterator it = GID_nodes.begin(); it != GID_nodes.end(); ++it)
	{
	    GlobalID[k] = *it;
	    ++k;
	}
    M_interfaceMap.reset( new MapEpetra(-1, LocalNodesNumber, GlobalID, M_d->comm));

}
