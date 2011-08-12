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


// ===================================================
//! Includes
// ===================================================
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

#include <life/lifearray/MapEpetra.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>

#include <life/lifesolver/VenantKirchhoffElasticData.hpp>
#include <life/lifesolver/StructuralMaterial.hpp>
#include <life/lifesolver/StructuralSolver.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialLinear.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialNonLinear.hpp>

#include <life/lifefem/TimeData.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifefem/TimeAdvance.hpp>
#include <life/lifefem/TimeAdvanceNewmark.hpp>
#include <life/lifefem/TimeAdvanceBDF.hpp>

#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>

#include <iostream>

#include "rotation.hpp"


// ===================================================
//! Namespaces & Define
// ===================================================
using namespace LifeV;

const int TOP    = 3;
const int BOTTOM = 2;
const int WALL   = 1;
const int RING  = 20;

// ===================================================
//! Private members
// ===================================================
struct problem::Private
{
    Private():
            rho (1), young(), poisson()
    {}


    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;
    Real rho, young, poisson;
    std::string    data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

    Real alpha ;
    Real theta , dtheta, dtheta2, ddtheta;
    Real c1, c2, c3, dc1, dc2, dc3, ddc1, ddc2, ddc3;

    Real displacementExact(const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
  {
  switch(i) {
    case 0:
      return  X * ( cos(theta) - 1 ) - Y * sin( theta ) + c1;
        break;
    case 1:
      return  X * sin( theta ) + Y * ( cos(theta) - 1 )+c2;
        break;
    case 2:
        return c3;
        break;
    default:
        ERROR_MSG("This entrie is not allowed");
        break;
    }
}


Real velocityExact(const Real& t, const Real& X, const Real& Y, const Real& Z, const ID& i)
{
  switch(i) {
    case 0:
      return - dtheta*( X*sin(theta) + Y*cos(theta) ) + dc1;
        break;
    case 1:
      return dtheta *( X*cos(theta) - Y*sin(theta) ) + dc2;
        break;
    case 2:
        return dc3;
        break;
    default:
        ERROR_MSG("This entrie is not allowed");
        break;
    }
}


Real accelerateExact(const Real& t, const Real& X, const Real& Y, const Real& Z, const ID& i)
{
  switch(i) {
    case 0:
      return - ddtheta*( X*sin(theta)+Y*cos(theta) ) + ddc1
	     - dtheta2 *( X*cos(theta) - Y*sin(theta) )	;
        break;
    case 1:
      return ddtheta * ( X*cos(theta) - Y*sin(theta) ) + ddc2
	     - dtheta2 *( X*sin(theta) + Y*cos(theta) )	;
        break;
    case 2:
        return dc3;
        break;
    default:
        ERROR_MSG("This entrie is not allowed");
        break;
    }
}


Real sourceTerm(const Real& t, const Real& X, const Real& Y, const Real& Z, const ID& i)
{
    switch(i) {
  case 0:
    return -rho*(ddtheta * ( X*sin(theta) + Y*cos(theta) )+ dtheta2 * ( X*cos(theta) - Y*sin(theta))-ddc1);
    break;
  case 1:
    return rho*( ddtheta* (X*cos(theta)-Y*sin(theta) ) - dtheta2 * ( X*sin(theta)+Y*cos(theta))+ddc2);
    break;
  case 2:
    return rho*ddc3;
    break;
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  } 
}

    fct_type getDisplacementExact()
    {
        fct_type f;
        f = boost::bind(&problem::Private::displacementExact, this, _1, _2, _3, _4, _5);
        return f;
    }

    fct_type getVelocityExact()
    {
        fct_type f;
        f = boost::bind(&problem::Private::velocityExact, this, _1, _2, _3, _4, _5);
        return f;
    }

    fct_type getAccelerateExact()
    {
        fct_type f;
        f = boost::bind(&problem::Private::accelerateExact, this,  _1, _2, _3, _4, _5);
        return f;
    }

    fct_type getSourceTerm()
    {
        fct_type f;
        f = boost::bind(&problem::Private::sourceTerm, this, _1, _2, _3, _4, _5);
        return f;
    }

};





// ===================================================
//! Constructors
// ===================================================
problem::problem( int          argc,
                  char**                argv,
                  boost::shared_ptr<Epetra_Comm>        structComm):
        members( new Private() )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    members->data_file_name = data_file_name;

    members->rho = dataFile( "solid/physics/density", 1. );
    members->young = dataFile( "solid/physics/young", 1. );
    members->poisson  = dataFile( "solid/physics/poisson", 1. );

    std::cout << "density = " << members->rho << std::endl
              << "young   = " << members->young << std::endl
              << "poisson = " << members->poisson << std::endl;

    members->comm = structComm;
    int ntasks = members->comm->NumProc();

    if (!members->comm->MyPID()) std::cout << "My PID = " << members->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
}

// ===================================================
//! Methods
// ===================================================
void
problem::run()
{
    typedef RegionMesh3D<LinearTetra>                                   mesh_Type;
    typedef StructuralSolver< mesh_Type >::vector_Type vector_type;
    typedef boost::shared_ptr<vector_type>                              vector_ptrtype;

    typedef boost::shared_ptr< TimeAdvance< vector_type > >             TimeAdvance_type;
    typedef FESpace< mesh_Type, MapEpetra >                             FESpace_type;
    typedef  boost::shared_ptr<FESpace_type>                            FESpace_ptrtype;

    bool verbose = (members->comm->MyPID() == 0);

    //
    // dataElasticStructure
    //

    GetPot dataFile( members->data_file_name.c_str() );
    boost::shared_ptr<VenantKirchhoffElasticData> dataProblem(new VenantKirchhoffElasticData( ));
    dataProblem->setup(dataFile, "solid");

    dataProblem->showMe();

    MeshData             meshData;
    meshData.setup(dataFile, "solid/space_discretization");

    meshData.showMe();

    boost::shared_ptr<mesh_Type > fullMeshPtr(new mesh_Type);
    readMesh(*fullMeshPtr, meshData);

    MeshPartitioner< mesh_Type > meshPart( fullMeshPtr, members->comm );

    //
    // The Problem Solver
    //

    // Scalar Solution Space:

    std::string Order =  dataFile( "solid/space_discretization/order", "P1");
    const ReferenceFE*    refFE(0);
    const QuadratureRule* qR(0);
    const QuadratureRule* bdQr(0);

    if ( Order.compare("P1") == 0 )
    {
        if (verbose) std::cout << "  Space order : P1" << std::flush;

        refFE = &feTetraP1;
        qR    = &quadRuleTetra15pt; // DoE 5
        bdQr  = &quadRuleTria4pt;   // DoE 2

    }
    else if ( Order.compare("P2") == 0 )
    {
        if (verbose) std::cout << " Space order : P2";

        refFE = &feTetraP2;
        qR    = &quadRuleTetra64pt; // DoE 6
        bdQr = &quadRuleTria4pt;   // DoE 2

    }
    if (verbose) std::cout << std::endl;


    // finite element space of the solution

    std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");

    FESpace_ptrtype feSpace( new FESpace_type(meshPart,dOrder,1,members->comm) );

    // instantiation of the VenantKirchhoffViscoelasticSolver class

    StructuralSolver< mesh_Type > problem;

    problem.setup(dataProblem,  feSpace,
                  members->comm);

    problem.setDataFromGetPot(dataFile);

    // the boundary conditions
    boost::shared_ptr<BCHandler> BCh( new BCHandler() );
 
    //
    // Boundary conditions for the displacement
    //
    BCFunctionBase sol(members->getDisplacementExact());

    BCh->addBC("EdgesIn",      2,  Essential, Full, sol,  3);
    BCh->addBC("EdgesIn",      3,  Essential, Full, sol,  3);
    BCh->addBC("EdgesIn",      20, Essential, Full, sol,  3);
    BCh->addBC("EdgesIn",      30, Essential, Full, sol,  3);
    BCh->addBC("EXternalWall", 10, Essential, Full, sol,  3);
    BCh->addBC("InternalWall", 1,  Essential, Full, sol,  3);

    std::ofstream out_norm;
    if (verbose)
    {
        out_norm.open("norm.txt");
        out_norm << "  time   "

        <<"  displacement L2_Error    "
        <<"  velocityL2_Error    "
        <<"  accelerateL2_Error  "
	<<"  \n";


        out_norm.close();
    }

    LifeChrono chrono;

    std::string TimeAdvanceMethod =  dataFile( "solid/time_discretization/method", "Newmark");

    TimeAdvance_type  timeAdvance( TimeAdvanceFactory::instance().createObject( TimeAdvanceMethod ) );

    UInt OrderDev = 2;

std::cout<<"dataProblem->dataTime()->orderBDF() "<<dataProblem->dataTime()->orderBDF()<<"\n";

    //! initialization of parameters of time Advance method:
    if (TimeAdvanceMethod =="Newmark")
        timeAdvance->setup( dataProblem->dataTime()->coefficientsNewmark() , OrderDev);

    if (TimeAdvanceMethod =="BDF")
        timeAdvance->setup(dataProblem->dataTime()->orderBDF() , OrderDev);

    timeAdvance->setTimeStep(dataProblem->dataTime()->timeStep());
    timeAdvance->showMe();


     const Real Pi = 3.14159265358979323846264338328;

   
     /*
     // Traslate function is a Sino:

     dtheta = 0, dtheta2 = 0 , ddtheta = 0, theta= 0;
     c1=sin(alpha* t), c2=cos( alpha * t ), c3=0.0,
     dc1=alpha * cos(alpha*t ), dc2= -alpha * sin(alpha*t),
     dc3=0.0, ddc1 = -alpha * alpha * sin(alpha*t),
     ddc2= -alpha * alpha * cos(alpha*t), ddc3=.0;
     */

   // Rotation 1:
   members->alpha = 10;
   members->theta = members->alpha*Pi, members->dtheta= members->alpha * Pi, members->ddtheta=0, members->dtheta2=members->dtheta*members->dtheta;


   // Rotation 2:
   /*
   Real theta  = -1./100*cos(100*Pi*t),  
            dtheta =Pi *sin(100*Pi*t),
           ddtheta =100* Pi*Pi*cos(100*Pi*t),
           dtheta2 = dtheta*dtheta;
    */

   //Rotation 3:
   /*
   Real theta  = 1./5.*(1-cos(50*Pi*t)),  
      dtheta = 10*Pi *sin(50*Pi* t),
      ddtheta =500 * Pi*Pi*cos(50*Pi*t),
      dtheta2 = dtheta*dtheta;
    */

  members->c1 = 0; 
  members->c2 = 0;
  members->c3 = 0;
  members->dc1 = 0; 
  members->dc2 = 0;  
  members->dc3 = 0; 
  members->ddc1 = 0;
  members->ddc2 = 0;
  members->ddc3 = 0;

  //Linear Rotation :

  //  Real c1=1*t, c2=1*t, c3=0, dc1=1, dc2=1, dc3=0, ddc1 =0, ddc2=0, ddc3=0;

  // Quadratic Rotation:

  //Real c1=t*t, c2=t*t, c3=0, dc1=2*t, dc2=2*t, dc3=0, ddc1 =2, ddc2=2, ddc3=0;




    dataProblem->setup(dataFile, "solid");

    chrono.start();

    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative( 0 ) / ( dataProblem->dataTime()->timeStep()*dataProblem->dataTime()->timeStep());

    problem.buildSystem(timeAdvanceCoefficient);


    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose ) std::cout << "ok." << std::endl;

    MapEpetra structMap(feSpace->refFE(), meshPart, members->comm);

    MapEpetra fullMap;

    for (UInt ii = 0; ii < nDimensions; ++ii)
    {
        fullMap += structMap;
    }

   // computing the rhs
    vector_type rhs ( fullMap, Unique );


    // postProcess
    boost::shared_ptr< Exporter<mesh_Type > > exporter;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
        exporter.reset( new ExporterHDF5<mesh_Type > ( dataFile, "problem" ) );
    else
#endif
    {
        if (exporterType.compare("none") == 0)
            exporter.reset( new ExporterEmpty<mesh_Type > ( dataFile, meshPart.meshPartition(), "problem", members ->comm->MyPID()) );
        else
            exporter.reset( new ExporterEnsight<mesh_Type > ( dataFile, meshPart.meshPartition(), "problem",   members->comm->MyPID()) );
    }

    exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( meshPart.meshPartition(),  members->comm->MyPID() );

    vector_ptrtype disp ( new vector_type(problem.displacement(), exporter->mapType() ) );
    vector_ptrtype vel  ( new vector_type(problem.displacement(),  exporter->mapType() ) );
    vector_ptrtype acc ( new vector_type(problem.displacement(),  exporter->mapType() ) );
    
    exporter->addVariable( ExporterData<mesh_Type>::VectorField, "displacement", feSpace,  disp, UInt(0) );
	
    exporter->addVariable( ExporterData<mesh_Type>::VectorField, "velocity", feSpace, vel, UInt(0) );

    exporter->addVariable( ExporterData<mesh_Type>::VectorField, "acceleration", feSpace, acc,  UInt(0) );

    exporter->postProcess( 0 );

    //evaluate disp and vel as interpolate the bcFunction d0 and v0
    feSpace->interpolate(members->getDisplacementExact(), *disp, 0.0);
    feSpace->interpolate(members->getVelocityExact(), *vel , 0.0);
    feSpace->interpolate(members->getAccelerateExact(), *acc, 0.0);

    //evaluate disp and vel as interpolate the bcFunction d0 and v0 and w0

    std::vector<vector_type> uv0;

    Real dt = dataProblem->dataTime()->timeStep();
    Real T  = dataProblem->dataTime()->endTime();

    if (TimeAdvanceMethod =="Newmark")
    {
        uv0.push_back(*disp);
        uv0.push_back(*vel);
        uv0.push_back(*acc);
    }
    if (TimeAdvanceMethod =="BDF")
    {
        for ( UInt previousPass=0; previousPass < dataProblem->orderBDF() ; previousPass++)
        {
            Real previousTimeStep = -previousPass*dt;
            feSpace->interpolate(members->getDisplacementExact(), *disp, previousTimeStep );
            uv0.push_back(*disp);
        }
    }

    timeAdvance->setInitialCondition(uv0);

    timeAdvance-> setTimeStep(dataProblem->dataTime()->timeStep());

    timeAdvance->updateRHSContribution(dataProblem->dataTime()->timeStep());

    vector_type uComputed(fullMap, Repeated );
    vector_type uExa(fullMap, Repeated );
    vector_type vExa(fullMap, Repeated );
    vector_type wExa(fullMap, Repeated );

    *disp = timeAdvance->solution();
    *vel = timeAdvance->velocity();
    *acc = timeAdvance->accelerate();

    exporter->postProcess( 0 );

    for (Real time = dt; time <= T; time += dt)
    {
        dataProblem->setTime(time);

   // Rotation 1:
   members->alpha = 10;
   members->theta = members->alpha*Pi*time, members->dtheta= members->alpha * Pi,   members->ddtheta=0, members->dtheta2=members->dtheta*members->dtheta;



        if (verbose)
        {
            std::cout << std::endl;
            std::cout << " P - Now we are at time " << dataProblem->dataTime()->time() << " s." << std::endl;
        }

        rhs *=0;

        timeAdvance->updateRHSContribution( dt );

        //evaluate rhs

        feSpace->l2ScalarProduct(members->getSourceTerm(), rhs, time);
        rhs += *problem.Mass() *timeAdvance->rhsContributionSecondDerivative();

        //update system
        problem.updateRightHandSide(rhs );

        //solver system
        problem.iterate( BCh );    // Computes the matrices and solves the system

        //update unknowns of timeAdvance
        timeAdvance->shiftRight(problem.displacement());

        //evaluate uexact solution
      //  feSpace->interpolate(members->getDisplacementExact(), Exact , time);
      //  feSpace->interpolate(members->getVelocityExact(), vExact , time);
      //  feSpace->interpolate(members->getAccelerateExact(), wExact , time);

        *disp =  timeAdvance->solution();
        *vel = timeAdvance->velocity();
        *acc = timeAdvance->accelerate();

        //postProcess
        exporter->postProcess( time );

        // Error L2

        Real displacementL2_Error, velocityL2_Error, accelerateL2_Error;

        displacementL2_Error = feSpace->l2Error(members->getDisplacementExact(),*disp , time);
        velocityL2_Error = feSpace->l2Error(members->getVelocityExact(), *vel , time );
        accelerateL2_Error = feSpace->l2Error(members->getAccelerateExact(), *acc, time );


        //save the norm
        out_norm.open("norm.txt", std::ios::app);
        out_norm << time  << "   "
        << displacementL2_Error << "   "
        << velocityL2_Error << "   "
        << accelerateL2_Error << " \n";
      
        out_norm.close();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    chrono.stop();
    if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
}
