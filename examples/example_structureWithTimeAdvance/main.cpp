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
   \file main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-16
 */
#ifdef TWODIM
#error test_structure cannot be compiled in 2D
#endif

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

#include <life/lifecore/LifeV.hpp>
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>


//Include fils which were in the structure.cpp file
#include <life/lifearray/MapEpetra.hpp>

#include <life/lifefem/TimeAdvance.hpp>
#include <life/lifefem/TimeAdvanceNewmark.hpp>
#include <life/lifefem/TimeAdvanceBDF.hpp>

#include <life/lifemesh/MeshData.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>

#include <life/lifesolver/VenantKirchhoffElasticData.hpp>

#include <life/lifesolver/StructuralMaterial.hpp>
#include <life/lifesolver/StructuralSolver.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialLinear.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialNonLinear.hpp>

#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>

#include <iostream>


using namespace LifeV;

int returnValue = EXIT_SUCCESS; // For the final check
enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}

std::set<UInt> parseList( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    size_t commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );
        setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
        stringList = stringList.substr( commaPos+1 );
    }
    setList.insert( atoi( stringList.c_str() ) );
    return setList;
}

class Structure
{
public:
  //@}
  /** @name Constructors, destructor
   */
  //@{

    Structure( int                                   argc,
               char**                                argv,
               boost::shared_ptr<Epetra_Comm>        structComm );

    ~Structure()
    {}
    //@}

    //@{
    void run()
    {
        run3d();
    }
    void CheckResults(const Real& dispNorm, const Real& time);
    void resultChanged(Real time);

    //@}

protected:

private:

    /**
     * run the 2D driven cylinder simulation
     */
    void run2d();

    /**
     * run the 3D driven cylinder simulation
     */
    void run3d();

private:
    struct Private;
    boost::shared_ptr<Private> parameters;
};

struct Structure::Private
{
    Private() :
            rho(1), poisson(1), young(1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;
    double rho, poisson, young;
    Real alpha ;
    Real theta , dtheta, dtheta2, ddtheta;
    Real c1, c2, c3, dc1, dc2, dc3, ddc1, ddc2, ddc3;


    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

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
        f = boost::bind(&Structure::Private::displacementExact, this, _1, _2, _3, _4, _5);
        return f;
    }

    fct_type getVelocityExact()
    {
        fct_type f;
        f = boost::bind(&Structure::Private::velocityExact, this, _1, _2, _3, _4, _5);
        return f;
    }

    fct_type getAccelerateExact()
    {
        fct_type f;
        f = boost::bind(&Structure::Private::accelerateExact, this,  _1, _2, _3, _4, _5);
        return f;
    }

    fct_type getSourceTerm()
    {
        fct_type f;
        f = boost::bind(&Structure::Private::sourceTerm, this, _1, _2, _3, _4, _5);
        return f;
    }

};


Structure::Structure( int                                   argc,
                      char**                                argv,
                      boost::shared_ptr<Epetra_Comm>        structComm):
        parameters( new Private() )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->rho = dataFile( "solid/physics/density", 1. );
    parameters->young = dataFile( "solid/physics/young", 1. );
    parameters->poisson  = dataFile( "solid/physics/poisson", 1. );

    std::cout << "density = " << parameters->rho << std::endl
              << "young   = " << parameters->young << std::endl
              << "poisson = " << parameters->poisson << std::endl;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID()) std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;

}

void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}

void
Structure::run3d()
{
    typedef StructuralSolver< RegionMesh3D<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    typedef boost::shared_ptr< TimeAdvance< vector_Type > >       timeAdvance_type;

    bool verbose = (parameters->comm->MyPID() == 0);
    // Number of boundary conditions for the velocity and mesh motion
    //
    boost::shared_ptr<BCHandler> BCh( new BCHandler() );

    //
    // dataElasticStructure
    //

    GetPot dataFile( parameters->data_file_name.c_str() );

    boost::shared_ptr<VenantKirchhoffElasticData> dataStructure(new VenantKirchhoffElasticData( ));
    dataStructure->setup(dataFile);

    MeshData             meshData;
    meshData.setup(dataFile, "solid/space_discretization");

    boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);
    readMesh(*fullMeshPtr, meshData);


    MeshPartitioner< RegionMesh3D<LinearTetra> > meshPart( fullMeshPtr, parameters->comm );

    std::string dOrder =  dataFile( "solid/space_discretization/order", "P1");

    typedef FESpace< RegionMesh3D<LinearTetra>, MapEpetra > solidFESpace_type;
    typedef boost::shared_ptr<solidFESpace_type> solidFESpace_ptrtype;
    solidFESpace_ptrtype dFESpace( new solidFESpace_type(meshPart,dOrder,3,parameters->comm) );
    if (verbose) std::cout << std::endl;

    MapEpetra structMap(dFESpace->refFE(), meshPart, parameters->comm);

    MapEpetra fullMap;



  /*
     // Traslate function is a Sino:

     dtheta = 0, dtheta2 = 0 , ddtheta = 0, theta= 0;
     c1=sin(alpha* t), c2=cos( alpha * t ), c3=0.0,
     dc1=alpha * cos(alpha*t ), dc2= -alpha * sin(alpha*t),
     dc3=0.0, ddc1 = -alpha * alpha * sin(alpha*t),
     ddc2= -alpha * alpha * cos(alpha*t), ddc3=.0;
     */
     const Real Pi = 3.14159265358979323846264338328;
   // Rotation 1:
   parameters->alpha = 10;
   parameters->theta = parameters->alpha*Pi, parameters->dtheta= parameters->alpha * Pi, parameters->ddtheta=0, parameters->dtheta2=parameters->dtheta*parameters->dtheta;


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

  parameters->c1 = 10;
  parameters->c2 = 0;
  parameters->c3 = 0;
  parameters->dc1 = 0;
  parameters->dc2 = 0;
  parameters->dc3 = 0;
  parameters->ddc1 = 0;
  parameters->ddc2 = 0;
  parameters->ddc3 = 0;


  //Linear Rotation :

  //  Real c1=1*t, c2=1*t, c3=0, dc1=1, dc2=1, dc3=0, ddc1 =0, ddc2=0, ddc3=0;

  // Quadratic Rotation:

  //Real c1=t*t, c2=t*t, c3=0, dc1=2*t, dc2=2*t, dc3=0, ddc1 =2, ddc2=2, ddc3=0;



   std::string timeAdvanceMethod =  dataFile( "solid/time_discretization/method", "Newmark");

    timeAdvance_type  timeAdvance( TimeAdvanceFactory::instance().createObject( timeAdvanceMethod ) );

    UInt OrderDev = 2;

    //! initialization of parameters of time Advance method:
    if (timeAdvanceMethod =="Newmark")
        timeAdvance->setup( dataStructure->dataTime()->coefficientsNewmark() , OrderDev);

    if (timeAdvanceMethod =="BDF")
        timeAdvance->setup(dataStructure->dataTime()->orderBDF() , OrderDev);

    timeAdvance->setTimeStep(dataStructure->dataTime()->timeStep());
    timeAdvance->showMe();

    for (UInt ii = 0; ii < nDimensions; ++ii)
    {
        fullMap += structMap;
    }


    // BC for cyl1x02_1796_edge.mesh
    vector <ID> compx(1), compy(1), compz(1);
    compx[0]=0; compy[0]=1, compz[0]=2;

    BCFunctionBase sol(parameters->getDisplacementExact());

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


    StructuralSolver< RegionMesh3D<LinearTetra> > solid;
    solid.setup(dataStructure,
                dFESpace,
		BCh,
                parameters->comm);

    solid.setDataFromGetPot(dataFile);

    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative( 0 ) / (dataStructure->dataTime()->timeStep()*dataStructure->dataTime()->timeStep());
    solid.buildSystem(timeAdvanceCoefficient);

    //
    // Temporal data and initial conditions
    //

    Real dt = dataStructure->dataTime()->timeStep();
    Real T  = dataStructure->dataTime()->endTime();

    vectorPtr_Type rhs(new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type disp(new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type vel(new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type acc(new vector_Type(solid.displacement(), Unique));

    dFESpace->interpolate(parameters->getDisplacementExact(), *disp, 0.0);
    dFESpace->interpolate(parameters->getVelocityExact(), *vel , 0.0);
    dFESpace->interpolate(parameters->getVelocityExact(), *acc , 0.0);

    if (verbose) std::cout << "S- initialization ... ";

    //solid.initialize(d0,w0,a0); // displacement, velocity, acceleration

    std::vector<vectorPtr_Type> uv0;


    if (timeAdvanceMethod =="Newmark")
      {
        uv0.push_back(disp);
        uv0.push_back(vel);
        uv0.push_back(acc);
      }
    if (timeAdvanceMethod =="BDF")
    {
        for ( UInt previousPass=0; previousPass < dataStructure->dataTime()->orderBDF() ; previousPass++)
        {
	  Real previousTimeStep = -previousPass*dt;
	  std::cout<<"BDF " <<previousTimeStep<<"\n";
	  dFESpace->interpolate(parameters->getDisplacementExact(), *disp, previousTimeStep );
	  uv0.push_back(disp);
        }
    }

    timeAdvance->setInitialCondition(uv0);

    timeAdvance->setTimeStep(dataStructure->dataTime()->timeStep());

    timeAdvance->updateRHSContribution(dataStructure->dataTime()->timeStep());

    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose ) std::cout << "ok." << std::endl;
    //if (parameters->comm->NumProc() == 1 )  solid.postProcess();


    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
      exporter.reset( new ExporterHDF5<RegionMesh3D<LinearTetra> > ( dataFile, "structure" ) );
    }
    else
#endif
      {
        if (exporterType.compare("none") == 0)
	  {
	    exporter.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID()) );
	  }
        else
	  {
	    exporter.reset( new ExporterEnsight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID()) );
	  }
      }

	exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
	exporter->setMeshProcId( meshPart.meshPartition(), parameters->comm->MyPID() );

	//vectorPtr_Type solidDisp ( new vector_Type(solid.getDisplacement(), exporter->mapType() ) );
	//vectorPtr_Type solidVel  ( new vector_Type(solid.getVelocity(),  exporter->mapType() ) );

	vectorPtr_Type solidDisp ( new vector_Type(solid.displacement(), exporter->mapType() ) );
	vectorPtr_Type solidVel  ( new vector_Type(solid.displacement(),  exporter->mapType() ) );
	vectorPtr_Type solidAcc  ( new vector_Type(solid.displacement(),  exporter->mapType() ) );

	exporter->addVariable( ExporterData<RegionMesh3D<LinearTetra> >::VectorField, "displacement",
                           dFESpace, solidDisp, UInt(0) );

	exporter->addVariable( ExporterData<RegionMesh3D<LinearTetra> >::VectorField, "velocity",
                           dFESpace, solidVel, UInt(0) );

	exporter->addVariable( ExporterData<RegionMesh3D<LinearTetra> >::VectorField, "acceleration",
                           dFESpace, solidAcc, UInt(0) );



	exporter->postProcess( 0 );
	solid.updateSystem();

	//
	// Temporal loop
	//

	for (Real time = dt; time <= T; time += dt)
	  {
	    dataStructure->dataTime()->setTime(time);

	    if (verbose)
	      {
		std::cout << std::endl;
		std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
	      }

	    *rhs *=0;
	    timeAdvance->updateRHSContribution( dt );
	    dFESpace->l2ScalarProduct(parameters->getSourceTerm(), *rhs, time);
	    *rhs += *solid.Mass() *timeAdvance->rhsContributionSecondDerivative()/timeAdvanceCoefficient;
	    solid.updateRightHandSide( *rhs );

	    solid.iterate( BCh );                  // Computes the matrices and solves the system
	    //if (parameters->comm->NumProc() == 1 )  solid.postProcess(); // Post-presssing

	    timeAdvance->shiftRight(solid.displacement());

	    *solidDisp = solid.displacement();
	    *solidVel  = timeAdvance->velocity();
	    *solidAcc  = timeAdvance->accelerate();

	    //if (parameters->comm->NumProc() == 1 )  solid.postProcess(); // Post-presssing

            //this->CheckResults(solid.displacement().norm2(),time);
	    exporter->postProcess( time );
        // Error L2

        Real displacementL2_Error, velocityL2_Error, accelerateL2_Error;

        displacementL2_Error = dFESpace->l2Error(parameters->getDisplacementExact(),*disp , time);
        velocityL2_Error = dFESpace->l2Error(parameters->getVelocityExact(), *vel , time );
        accelerateL2_Error = dFESpace->l2Error(parameters->getAccelerateExact(), *acc, time );


        //save the norm
        out_norm.open("norm.txt", std::ios::app);
        out_norm << time  << "   "
        << displacementL2_Error << "   "
        << velocityL2_Error << "   "
        << accelerateL2_Error << " \n";

        out_norm.close();
        MPI_Barrier(MPI_COMM_WORLD);

	  }

}

void Structure::CheckResults(const Real& dispNorm,const Real& time)
{
    if ( time == 0.001  && std::fabs(dispNorm-1.18594)>1e-4 )
        this->resultChanged(time);
    if ( time == 0.002  && std::fabs(dispNorm-1.10232)>1e-4 )
        this->resultChanged(time);
    if ( time == 0.003  && std::fabs(dispNorm-0.808509)>1e-4 )
        this->resultChanged(time);
}

void Structure::resultChanged(Real time)
{
  std::cout << "Some modifications led to changes in the l2 norm of the solution at time " << time << std::endl;
  returnValue = EXIT_FAILURE;
}



int
main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);

    boost::shared_ptr<Epetra_MpiComm> Comm(new Epetra_MpiComm( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
        cout << "% using MPI" << endl;
#else
    boost::shared_ptr<Epetra_SerialComm> Comm( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif

    Structure structure( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}

