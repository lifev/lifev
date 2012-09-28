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

    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @author Claudia Colciago <claudia.colciago@epfl.ch>
    @date 08-03-2011
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
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDFNavierStokes.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEnsight.hpp>


#include "OseenSolverBoundaryDerivative.hpp"
#include "ud_functions.hpp"
#include "RobinMembraneSolver.hpp"

#include <iostream>

#define OUTLET 3
#define INLET 2
#define WALL 1
//#define OUTERWALL 10
#define OUTERWALL 9 //(for Philippe's aorta)
#define SOLIDINTERFACE 1
//#define RING2 22
#define RING 2
#define RING3 3
#define RING4 4
#define RING5 5
#define RING6 6
#define RING7 7
//#define RING8 8
#define EDGE1 20
#define EDGE2 30
#define EDGE3 40
#define EDGE4 50
#define EDGE5 60
#define EDGE6 70

using namespace LifeV;



/*const int INLET       = 2;
const int WALL        = 1;
const int OUTLET      = 3;
const int OLD_ARTERY  = 4;
const int RINGIN      = 20;
const int RINGOUT     = 30;
const int OLD_RING    = 40;
*/


const Real PI=3.141592653589793;


struct RobinMembraneSolver::Private
{
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_Type;

    double Re;

    std::string data_file_name;

    double      nu;  /**< viscosity (in m^2/s) */
    double      H;   /**< height and width of the domain (in m) */
    double      D;   /**< diameter of the cylinder (in m) */
    bool        centered; /**< true if the cylinder is at the origin */
    double density;  
    double R;//radius
    //membrane
    double rhos;
    double Hs;
    double ni;
    double E;

    std::string initial_sol;

    boost::shared_ptr<Epetra_Comm>   comm;

    /**
     * u2d flat 2D velocity profile.
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real u2d( const Real& t,
              const Real& /*x*/,
              const Real& /*y*/,
              const Real& /*z*/,
              const ID&   id ) const
    {

      /*switch (id)
        {
        case 0: // x component
            return 0.0;
            break;
	    case 2: // z component*/
      if ( t <= 0.0025 )
	      //return 6.6e3*(1-cos(3.14159*t/0.0025));*/
	      //return -814*(sin(2*PI*t/0.005)*sin(2*PI*t/0.005));
	    return  -5e4*(sin(2*PI*t/0.005)*sin(2*PI*t/0.005));
            //      return 0.01;
	    /* return 0.0;
            break;
        case 1: // y component
            return 0.0;
            //      return 1.3332e4;
            //    else
            //      return 0.0;
            break;
	    }*/
        return 0;
    }

    fct_Type getU_2d()
    {
        fct_Type f;
        f = boost::bind(&RobinMembraneSolver::Private::u2d, this, _1, _2, _3, _4, _5);
        return f;
    }


};

RobinMembraneSolver::RobinMembraneSolver( int argc,
                    char** argv )
        :
        M_d( new Private )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    M_d->data_file_name = data_file_name;

    M_d->Re          = dataFile( "fluid/problem/Re", 1. );
    M_d->ni          = dataFile( "fluid/physics/poisson", 1. );
    M_d->rhos        = dataFile( "fluid/physics/density_sol", 1. );
    M_d->E           = dataFile( "fluid/physics/young", 1. );
    M_d->nu          = dataFile( "fluid/physics/viscosity", 1. );
    M_d->Hs          = dataFile( "fluid/physics/wall_thickness", 1. );
    M_d->H           = dataFile( "fluid/problem/H", 20. );
    M_d->D           =               dataFile( "fluid/problem/D", 1. );
    M_d->R = M_d->D/2;
    M_d->centered    = (bool)        dataFile( "fluid/problem/centered", 0 );
    M_d->density = dataFile( "fluid/physics/density", 1. );
    M_d->initial_sol = (std::string) dataFile( "fluid/problem/initial_sol", "none");
    std::cout << M_d->initial_sol << std::endl;


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    int ntasks = 0;
    M_d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    if (!M_d->comm->MyPID())
    {
        std::cout << "My PID = " << M_d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "fluid density = " << M_d->density << std::endl
                  << "nu = " << M_d->nu << std::endl
	          << "poisson = " << M_d->ni << std::endl
                  << "initial solution  = " << M_d->initial_sol  << std::endl
                  << "Young modulus  = " << M_d->E  << std::endl
                  << "Wall thickness  = " << M_d->Hs  << std::endl;

    }

#else
    M_d->comm.reset( new Epetra_SerialComm() );
#endif

}

void
RobinMembraneSolver::run()

{

    bool verbose = (M_d->comm->MyPID() == 0);

    //------------Creating file to store data------------------------------

    typedef OseenSolver< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    GetPot dataFile( M_d->data_file_name );

    boost::shared_ptr<OseenData> oseenData(new OseenData());
    oseenData->setup( dataFile );

    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");

    //-----------------------------------------------------------------------

    //----------Creating Mesh and FESpaces-----------------------------------

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra>);
    readMesh(*fullMeshPtr, meshData);

    MeshPartitioner< RegionMesh<LinearTetra> >   meshPart(fullMeshPtr, M_d->comm);

    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");
    std::string mOrder =  dataFile( "mesh_motion/space_discretization/m_order", "P1");
    if (verbose)
        std::cout << "Building the FE spaces ... " << std::flush;

    M_uFESpace.reset(new FESpace<mesh_Type, MapEpetra>(meshPart,  uOrder,3,  M_d->comm));
    M_pFESpace.reset( new FESpace<mesh_Type, MapEpetra>(meshPart, pOrder,1,  M_d->comm));
    M_mFESpace.reset(new FESpace<mesh_Type, MapEpetra>(meshPart,  mOrder,3,  M_d->comm));

    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalVelDof   = M_uFESpace->map().map(Unique)->NumGlobalElements();
    UInt totalPressDof = M_pFESpace->map().map(Unique)->NumGlobalElements();

    if (verbose) std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;
   
    
    //------------------FLUID PROBLEM--------------------------------

    if (verbose) std::cout << "Calling the fluid constructor ... ";
    int numLM = 1;

    M_fluid.reset(new solver_Type( oseenData,
                                              *M_uFESpace,
                                              *M_pFESpace,
				   M_d->comm,  WALL, numLM ) );
    MapEpetra fluidMap(M_fluid->getMap());

    if (verbose) std::cout << "ok." << std::endl;

    Real LameI = (M_d->Hs * M_d->E * M_d->ni )/( ( 1 - M_d->ni*M_d->ni ) );
    //Real LameI = (  M_d->Hs * M_d->E * M_d->ni )/( ( 1 - 2 * M_d->ni ) * ( 1 + M_d->ni ) );
    Real LameII = M_d->Hs * M_d->E / ( 2 * ( 1 + M_d->ni ) );

    std::cout<<" LameI : "<<LameI<<std::endl;
    std::cout<<" LameII : "<<LameII<<std::endl;

    M_fluid->setUp(dataFile);
    M_fluid->buildSystem( 2 * LameII * oseenData->dataTime()->timeStep(), LameI * oseenData->dataTime()->timeStep(), M_d->rhos* M_d->Hs );

    //---------------------------------------------------------

    //--------------------Creating BDF objects---------------------------

    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();

    TimeAdvanceBDFNavierStokes<vector_Type> fluidTimeAdvance;
    fluidTimeAdvance.setup(oseenData->dataTime()->orderBDF());

    TimeAdvanceBDF<vector_Type> meshTimeAdvance;
    meshTimeAdvance.setup(oseenData->dataTime()->orderBDF());

    TimeAdvanceBDF<vector_Type> dispTimeAdvance;
    dispTimeAdvance.setup(oseenData->dataTime()->orderBDF(),UInt(2));
  
    //Dt fluido
    fluidTimeAdvance.bdfVelocity().setTimeStep(oseenData->dataTime()->timeStep());
    //fluidTimeAdvance.showMe();

    //Dt mesh
    meshTimeAdvance.setTimeStep(oseenData->dataTime()->timeStep());
    //meshTimeAdvance.showMe();

    //Dt mesh
    dispTimeAdvance.setTimeStep(oseenData->dataTime()->timeStep());
    dispTimeAdvance.showMe();
    if(verbose)
	  std::cout<<" done \n"; 


    //---------------------------------------------------------


    //--------------Boundary Conditions-------------------------

    //creating BCHandler objects
    BCHandler bcH_fluid;
    BCHandler bcH_stokes;
    BCFunctionBase flow_in (linearCylinderFlux);
    BCFunctionBase vel_in (aortaVelocityIn);
    BCFunctionBase flow_3 (linearFluxThoracicAorta);
    BCFunctionBase flow_4 (linearFluxLeftSubclavian);
    BCFunctionBase flow_5 (linearFluxRightCarotid);
    BCFunctionBase flow_6 (linearFluxRightVertebral);
    BCFunctionBase flow_7 (linearFluxRightSubclavian);
    BCFunctionBase flow_8 (linearFluxCommonCarotid);
    BCFunctionBase flow_9 (linearFluxLeftVertebral);

    BCFunctionBase uZero( fZero );


    std::vector <ID> compXY(2),compZ(1);
    compXY[0]=0;
    compXY[1]=1;
    compZ[0]=2;
        
    //Set Robin Boundary Condition
    vector_Type robinCoeff( fluidMap, Repeated );
    vector_Type robinExt(fluidMap, Repeated);
    M_uFESpace->interpolate( fZero , robinExt, 0);
    robinCoeff *= 0.0;
    robinCoeff +=  ( M_d->rhos* M_d->Hs * fluidTimeAdvance.bdfVelocity().coefficientFirstDerivative( 0 ) / dt + robinExt * (0.1 + dt));

    vector_Type robinRHS(fluidMap, Repeated, Zero  );

    BCVector bcRobin(robinRHS, M_uFESpace->dof().numTotalDof(), 0);
    bcRobin.setRobinCoeffVector( robinCoeff );
    

    Real  M_resistanceCoef;
    vector_Type M_vecResistance(fluidMap, Repeated);
    Real betaCurvCoeff=((M_d->Hs*M_d->E)/(1-M_d->ni*M_d->ni))*((1/(M_d->R*M_d->R)));
    Real area0 = M_fluid->area(3);

    M_vecResistance *=0;

    //resistance coeff
    M_vecResistance += - pow(sqrt(M_d->rhos)/(2.*sqrt(2))*(M_fluid->flux(3))/area0 + sqrt(betaCurvCoeff*sqrt(area0)),2)
            + betaCurvCoeff*sqrt(area0); 

    BCVector M_bcvResistance(M_vecResistance, M_uFESpace->dof().numTotalDof(), 1);

    bcH_fluid.addBC("InFlow" , INLET,  Flux, /*Full*/Normal, flow_in );
    bcH_fluid.addBC("InFlow_2" , OUTLET,  Natural, Full, uZero, 3);
    /*bcH_fluid.addBC("InFlow_3" , 9,  Flux, Normal, flow_4);
    bcH_fluid.addBC("InFlow_4" , 5,  Flux, Normal, flow_5);
    bcH_fluid.addBC("InFlow_5" , 6,  Flux, Normal, flow_6);
    bcH_fluid.addBC("InFlow_6" , 7,  Flux, Normal, flow_7);
    bcH_fluid.addBC("InFlow_7" , 4,  Flux, Normal, flow_8);
    bcH_fluid.addBC("InFlow_8" , 8,  Flux, Normal, flow_9);*/
    bcH_fluid.addBC("Ring2" , 20,  EssentialVertices, Full, uZero,3);
    bcH_fluid.addBC( "Wall",     WALL,     Robin,   Full,     bcRobin, 3 );  
    bcH_fluid.setOffset("InFlow", totalVelDof + totalPressDof);
    /*bcH_fluid.setOffset("InFlow_2", totalVelDof + totalPressDof +1);
    bcH_fluid.setOffset("InFlow_3", totalVelDof + totalPressDof + 2);
    bcH_fluid.setOffset("InFlow_4", totalVelDof + totalPressDof + 3);
    bcH_fluid.setOffset("InFlow_5", totalVelDof + totalPressDof + 4);
    bcH_fluid.setOffset("InFlow_6", totalVelDof + totalPressDof + 5);
    bcH_fluid.setOffset("InFlow_7", totalVelDof + totalPressDof + 6);
    bcH_fluid.setOffset("InFlow_8", totalVelDof + totalPressDof + 7);
    */       


    bcH_stokes.addBC("InFlow" , INLET,  Flux, /*Full*/Normal, flow_in);
    bcH_stokes.addBC("InFlow_2" , OUTLET,  Flux/*Essential*/, Normal, flow_3);
    bcH_stokes.addBC("InFlow_3" , 9,  Flux/*Essential*/, Normal, flow_4);
    bcH_stokes.addBC("InFlow_4" , 5,  Flux/*Essential*/, Normal, flow_5);
    bcH_stokes.addBC("InFlow_5" , 6,  Flux/*Essential*/, Normal, flow_6);
    bcH_stokes.addBC("InFlow_6" , 7,  Flux/*Essential*/, Normal, flow_7);
    bcH_stokes.addBC("InFlow_7" , 4,  Flux/*Essential*/, Normal, flow_8);
    bcH_stokes.addBC("InFlow_8" , 8,  Flux/*Essential*/, Normal, flow_9);
    bcH_stokes.addBC("Ring2" , 30,  EssentialVertices, Full, uZero,3);
    bcH_stokes.addBC( "Wall",     WALL,     Essential,   Full,     uZero, 3 );  
    bcH_stokes.setOffset("InFlow", totalVelDof + totalPressDof);
    bcH_stokes.setOffset("InFlow_2", totalVelDof + totalPressDof + 1);
    bcH_stokes.setOffset("InFlow_3", totalVelDof + totalPressDof + 2);
    bcH_stokes.setOffset("InFlow_4", totalVelDof + totalPressDof + 3);
    bcH_stokes.setOffset("InFlow_5", totalVelDof + totalPressDof + 4);
    bcH_stokes.setOffset("InFlow_6", totalVelDof + totalPressDof + 5);
    bcH_stokes.setOffset("InFlow_7", totalVelDof + totalPressDof + 6);
    bcH_stokes.setOffset("InFlow_8", totalVelDof + totalPressDof + 7);
  
    if(verbose)
	  std::cout<<" BC done \n"; 

    //----------------------------------------------------------------------------


    //-------------Creating Exporter Object------------------------------------------------

    std::string const exporterFileName    =  dataFile( "exporter/filename", "cylinder");

    ExporterHDF5<RegionMesh<LinearTetra> > exporter( dataFile, meshPart.meshPartition(), exporterFileName, M_d->comm->MyPID());

    vectorPtr_Type velAndPressure ( new vector_Type( fluidMap, exporter.mapType(), Zero ));
    vectorPtr_Type eta (new vector_Type(M_uFESpace->map(), exporter.mapType(), Zero ));
    vectorPtr_Type vectorUseful (new vector_Type(M_uFESpace->map(), exporter.mapType(), Zero ));


    exporter.addVariable( ExporterData< mesh_Type >::VectorField, "f-velocity",
			  M_uFESpace, velAndPressure, UInt(0) );

    exporter.addVariable( ExporterData< mesh_Type >::ScalarField, "f-pressure",
			  M_pFESpace, velAndPressure, UInt(3*M_uFESpace->dof().numTotalDof() ) );  

    exporter.addVariable( ExporterData< mesh_Type >::VectorField, "s-displacement",
			  M_uFESpace, eta, UInt(0) ) ;
    if(verbose)
	  std::cout<<" Exporter done \n"; 

    //-----------------------------------------------------------------------------

    
    //--------------------------Initialization--------------------------------------

    vector_Type beta( fluidMap );    
    vector_Type rhs ( fluidMap );
    vector_Type etaNew( fluidMap );
    vector_Type etaOld( fluidMap );

    etaOld *= 0;
    etaNew *= 0;

    if (M_d->initial_sol == "stokes")
    {
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

        oseenData->dataTime()->setTime(t0);

        beta *= 0.;
        rhs  *= 0.;

        M_fluid->updateSystem(0, beta, rhs );
        M_fluid->iterate( bcH_stokes );

        *velAndPressure = *M_fluid->solution();
        M_fluid->resetPreconditioner();

	etaOld= etaNew;	
	etaNew = (dt*(*(M_fluid->solution())) + etaOld);
	*eta = etaNew;

    }

    if (M_d->initial_sol == "restart")
    {
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "Restoring the previous solution ... " << std::endl;
       
	
	std::string start    = dataFile( "fluid/importer/start", "00000");
        std::string filename = dataFile("importer/filename", "cylinder");

        LifeV::ExporterHDF5<mesh_Type> importer( dataFile, filename);
        importer.setMeshProcId(M_uFESpace->mesh(), M_d->comm->MyPID());

        importer.addVariable( ExporterData<mesh_Type>::VectorField,
                              "f-velocity",
                              M_uFESpace,
                              velAndPressure,
                              UInt ( 0 ) );

        importer.addVariable( ExporterData<mesh_Type>::ScalarField,
                              "f-pressure",
                              M_pFESpace,
                              velAndPressure,
                              3*M_uFESpace->dof().numTotalDof() );

        importer.addVariable( ExporterData<mesh_Type>::VectorField,
                              "s-displacement",
                              M_uFESpace,
                              eta,
                              UInt( 0 ) );

        exporter.setTimeIndex(importer.importFromTime(t0));

        double norm = velAndPressure->norm2();
        if (verbose)
            std::cout << "   f- restart solution norm = " << norm << std::endl;
        M_fluid->initialize(*velAndPressure);
	etaNew = *eta;


    }

    fluidTimeAdvance.bdfVelocity().setInitialCondition( *(velAndPressure) );
    //*velAndPressure = *(M_fluid->solution());
 
    std::vector <vector_Type> initialCondition;
    for ( int previousPass=0; previousPass <(int) oseenData->dataTime()->orderBDF() ; previousPass++)
	{	    
	    initialCondition.push_back(etaNew);
	}
    if(verbose)
	  std::cout<<" Initializing displacement BDF \n";     
    dispTimeAdvance.setInitialCondition( initialCondition );
    if(verbose)
	  std::cout<<" done \n"; 
    //exportGID(*vectorUseful); 

    exporter.postProcess( t0 );
    
    //-------------------------------------------------------------------------------------


    // Temporal loop

    LifeChrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        oseenData->dataTime()->setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time "<< oseenData->dataTime()->time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        //--------------Solving Fluid Problem----------------------

        //update Boundary Conditions
        if(verbose)
	  std::cout<<" Update and Solve Fluid Problem \n";
	
        fluidTimeAdvance.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
	dispTimeAdvance.updateRHSContribution(oseenData->dataTime()->timeStep());

	robinRHS = ( ( M_d->rhos * M_d->Hs ) ) * fluidTimeAdvance.bdfVelocity().rhsContributionFirstDerivative()
	                -robinExt*dt*dispTimeAdvance.rhsContributionFirstDerivative() ;
	
	M_vecResistance *= 0.0;
	
	M_vecResistance +=  - pow(sqrt(M_d->rhos)/(2.*sqrt(2))*M_fluid->flux(3)/area0 + sqrt(betaCurvCoeff*sqrt(area0)),2)
            + betaCurvCoeff*sqrt(area0);

        if(verbose)
	  std::cout<<"Update Fluid Convective Term\n ";

	beta*=0.0;

        beta += fluidTimeAdvance.bdfVelocity().extrapolation();

	rhs  = M_fluid->matrixMass()*fluidTimeAdvance.bdfVelocity().rhsContributionFirstDerivative();

        Real alpha = fluidTimeAdvance.bdfVelocity().coefficientFirstDerivative( 0 ) / oseenData->dataTime()->timeStep();
        
        M_fluid->updateSystem( alpha, beta, rhs, -1.0*dispTimeAdvance.rhsContributionFirstDerivative() );

        M_fluid->iterate( bcH_fluid );
        
        fluidTimeAdvance.bdfVelocity().shiftRight( *(M_fluid->solution()) );

        *velAndPressure = *(M_fluid->solution());

	//---------------------------------------------------------------------------

	//------------------Updating solid displacement-----------------------------

	etaOld = etaNew;

	etaNew = (dt*(*M_fluid->solution()) + etaOld);

	dispTimeAdvance.shiftRight( etaNew );

	*eta = etaNew;

	//------------------------------------------------------------------------------------

        exporter.postProcess( time );

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//********************************************************************************************

void RobinMembraneSolver::exportGID( vector_Type& vectorGID )
{
    
    Int numMyEntries = vectorGID.epetraVector().MyLength();
    const Int* gids  = vectorGID.blockMap().MyGlobalElements();

    for ( Int i = 0; i < numMyEntries; ++i )
    {
        vectorGID(gids[i]) = gids[i];
    }

}

void RobinMembraneSolver::initialize( GetPot const& data_file, vectorPtr_Type velAndPressure, vectorPtr_Type eta)
{

    Real Tstart=data_file( "fluid/time_discretization/initialtime", 0.);

    std::cout<<Tstart<<std::endl;

    std::string const importerType =  data_file( "importer/type", "hdf5");
    std::string const fluidName    =  data_file( "importer/filename", "testRestart");

    std::cout<<"Il nome del file che voglio leggere: "<<fluidName<<std::endl;

    ExporterHDF5<RegionMesh<LinearTetra> > importer( data_file, fluidName);

    vectorPtr_Type  velAndPressImport( new vector_Type( M_fluid->getMap(), importer.mapType(), Zero ) );
    vectorPtr_Type  etaImport( new vector_Type( M_uFESpace->map(), importer.mapType(), Zero ) );

    importer.setMeshProcId(M_uFESpace->mesh(), M_uFESpace->map().comm().MyPID());
    
    importer.addVariable( ExporterData<mesh_Type>::VectorField, "f-velocity",
                                  M_uFESpace, velAndPressImport, UInt(0) );

    importer.addVariable( ExporterData<mesh_Type>::ScalarField, "f-pressure",
                                  M_pFESpace, velAndPressImport,
                                  UInt(3*M_uFESpace->dof().numTotalDof()) );

    importer.addVariable( ExporterData<mesh_Type>::VectorField, "s-displacement",
                                  M_mFESpace, etaImport, UInt(0) );



    boost::shared_ptr<LifeV::VectorEpetra> initSol(new LifeV::VectorEpetra(*velAndPressImport));
    boost::shared_ptr<LifeV::VectorEpetra> UniqueV;
    
    importer.importFromTime(Tstart);
    
    UniqueV.reset( new vector_Type(*velAndPressImport, Unique, Zero));
    *initSol=*UniqueV;

    *velAndPressure=*velAndPressImport;
    *eta=*etaImport;

    M_fluid->initialize(*velAndPressImport);
    
}

