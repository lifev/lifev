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
#include <life/lifesolver/NeoHookeanMaterialNonLinear.hpp>
#include <life/lifesolver/ExponentialMaterialNonLinear.hpp>

#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>

#include <iostream>


using namespace LifeV;

int returnValue = EXIT_SUCCESS;
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
            rho(1), poisson(1), young(1), bulk(1), alpha(1), gamma(1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;
    double rho, poisson, young, bulk, alpha, gamma;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

static Real bcZero(const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
{
    return  0.;
}

static Real bcNonZero(const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
{
    return  300000.;
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

    parameters->rho     = dataFile( "solid/physics/density", 1. );
    parameters->young   = dataFile( "solid/physics/young",   1. );
    parameters->poisson = dataFile( "solid/physics/poisson", 1. );
    parameters->bulk    = dataFile( "solid/physics/bulk",    1. );
    parameters->alpha   = dataFile( "solid/physics/alpha",   1. );
    parameters->gamma   = dataFile( "solid/physics/gamma",   1. );

    std::cout << "density = " << parameters->rho     << std::endl
              << "young   = " << parameters->young   << std::endl
              << "poisson = " << parameters->poisson << std::endl
              << "bulk    = " << parameters->bulk    << std::endl
              << "alpha   = " << parameters->alpha   << std::endl
              << "gamma   = " << parameters->gamma   << std::endl;

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

    //! Number of boundary conditions for the velocity and mesh motion
    boost::shared_ptr<BCHandler> BCh( new BCHandler() );

    //! dataElasticStructure
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


    //! =================================================================================
    //! BC for StructuredCube4.mesh
    //! =================================================================================
    vector <ID> compx(1), compy(1), compz(1);
    compx[0]=0; compy[0]=1, compz[0]=2;

    BCFunctionBase zero(Private::bcZero);
    BCFunctionBase nonZero(Private::bcNonZero);

    BCh->addBC("EdgesIn",      2,  Natural,           Component, nonZero, compx);
    BCh->addBC("EdgesIn",      4,  EssentialVertices, Component, zero,    compx);
    //! =================================================================================


    //! 1. Constructor of the structuralSolver
    StructuralSolver< RegionMesh3D<LinearTetra> > solid;

    //! 2. Setup of the structuralSolver
    solid.setup(dataStructure,
                dFESpace,
		BCh,
                parameters->comm);

    //! 3. Setting data from getPot
    solid.setDataFromGetPot(dataFile);

    //! 4. Building system using TimeAdvance class
    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative( 0 ) / (dataStructure->dataTime()->timeStep()*dataStructure->dataTime()->timeStep());
    solid.buildSystem(timeAdvanceCoefficient);



    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    //! 5. Initial data
    Real dt = dataStructure->dataTime()->timeStep();
    Real T  = dataStructure->dataTime()->endTime();

    vectorPtr_Type rhs (new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type disp(new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type vel (new vector_Type(solid.displacement(), Unique));
    vectorPtr_Type acc (new vector_Type(solid.displacement(), Unique));

    if (verbose) std::cout << "S- initialization ... ";

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
	  uv0.push_back(disp);
        }
    }

    timeAdvance->setInitialCondition(uv0);

    timeAdvance->setTimeStep(dataStructure->dataTime()->timeStep());

    timeAdvance->updateRHSContribution(dataStructure->dataTime()->timeStep());

    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose ) std::cout << "ok." << std::endl;

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

    exporter->setPostDir( "./" );
    exporter->setMeshProcId( meshPart.meshPartition(), parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type(solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidVel  ( new vector_Type(solid.displacement(),  exporter->mapType() ) );
    vectorPtr_Type solidAcc  ( new vector_Type(solid.displacement(),  exporter->mapType() ) );

    exporter->addVariable( ExporterData<RegionMesh3D<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt(0) );
    exporter->addVariable( ExporterData<RegionMesh3D<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt(0) );
    exporter->addVariable( ExporterData<RegionMesh3D<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt(0) );

    exporter->postProcess( 0 );
    //! =================================================================================   



    //! =============================================================================
    //! Temporal loop
    //! =============================================================================
    for (Real time = dt; time <= T; time += dt)
    {
	dataStructure->dataTime()->setTime(time);

	if (verbose)
        {
		std::cout << std::endl;
		std::cout << "S- Now we are at time " << dataStructure->dataTime()->time() << " s." << std::endl;
	}

        //! 6. Updating right-hand side
	*rhs *=0;
	timeAdvance->updateRHSContribution( dt );
	*rhs += *solid.Mass() *timeAdvance->rhsContributionSecondDerivative()/timeAdvanceCoefficient;
	solid.updateRightHandSide( *rhs );

        //! 7. Iterate --> Calling Newton
	solid.iterate( BCh );

        timeAdvance->shiftRight( solid.displacement() );

        *solidDisp = solid.displacement();
	*solidVel  = timeAdvance->velocity();
	*solidAcc  = timeAdvance->accelerate();

	exporter->postProcess( time );

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


