#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>

#include <ChorinTemam.hpp>

#include <ct.hpp>
#include <iostream>

using namespace LifeV;

/* References for tube3D_xxxM.mesh */
const int INLET        = 2;
const int OUTLET       = 3;
const int WALL         = 1;
/* Dimensions for tube3D_xxxM.mesh */
const Real RADIUS      = 4.9;
const Real HEIGHT      = 100.0;
/* References for cylinder */
const int CYL_INLET    = 40;
const int CYL_WALL     = 60;
const int CYL_SLIPWALL = 61;
const int CYL_OUTLET   = 50;
const int CYL_CYLINDER = 70;


/* Define for different meshes, bcs */
#undef __CT_VELOCITY_TUBE	/* for velocity excitation on tube */
#undef __CT_PRESSURE_TUBE	/* for pressure excitation on tube */
#undef __CT_CYLINDER_CASE	/* for cylinder to compare w/ Oseen */
#define __CT_CYLINDER_CASE 23


/*
 * The CT::Private struct contains mainly the functions that will be used
 * for computing boundary conditions 
 * We follow the encapsulating mechanism of test_cylinder.
 */

struct CT::Private
{
    Private() 
	{}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    std::string data_file_name;

    Epetra_Comm*   comm;

    Real u3DIn( const Real& t, 
		const Real& x, 
		const Real& y, 
		const Real& z, 
		const ID& id ) const
    {
            if ( id == 3 ) {
                  return 100.0 * t + 0.1;
		}
            else {
		  return 0.0;
            }
    }

    fct_type get_u3DIn()
        {
            fct_type f;
            f = boost::bind(&CT::Private::u3DIn, this, _1, _2, _3, _4, _5);
            return f;
        }

    Real u3DZero( const Real& t,
		  const Real& x,
		  const Real& y,
		  const Real& z,
		  const ID&   id ) const
    {
                 return 0.0;
    }

    fct_type get_u3DZero()
        {
            fct_type f;
            f = boost::bind(&CT::Private::u3DZero, this, _1, _2, _3, _4, _5);
            return f;
        }

    Real p3DIn( const Real& t,
		const Real& x, 
		const Real& y,
		const Real& z,
		const ID& id ) const
    {
	return -0.1 * cos(t) + 0.1;
    }
 
    fct_type get_p3DIn()
	{
	    fct_type f;
	    f = boost::bind(&CT::Private::p3DIn, this, _1, _2, _3, _4, _5);
	    return f;
  	}

    Real u3Dcyl( const Real& t,
		 const Real& x,
		 const Real& y,
		 const Real& z,
		 const ID& id ) const
    {
        if ( id == 1 ) {
            return 1./(20.*20.)*(y + 20.)*(20. -y);
        } else {
            return 0.;
        }
    }

    fct_type get_u3Dcyl()
        {
            fct_type f;
            f = boost::bind(&CT::Private::u3Dcyl, this, _1, _2, _3, _4, _5);
            return f;
        }

};

CT::CT( int argc,
                    char** argv,
                    LifeV::AboutData const& /*ad*/,
                    LifeV::po::options_description const& /*od*/ )
    :
    d( new Private )
{
    GetPot command_line(argc, argv);
    const char* data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    d->data_file_name = data_file_name;

#ifdef EPETRA_MPI
    d->comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    std::cout << "  t-  MPI Initialization from PID = " << d->comm->MyPID() 
	<< " among " << ntasks << " running." << std::endl;
#else
    d->comm = new Epetra_SerialComm();
#endif

}

/*
 * CT::run() 
 * Inspired from test_cylinder.
 */

void
CT::run()
{

    typedef ChorinTemam< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    

    // Reading from data file
    GetPot dataFile( d->data_file_name.c_str() );

    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (d->comm->MyPID() == 0);

#ifdef __CT_VELOCITY_TUBE
    // Boundary conditions for tube
    BCHandler bcHu( 3, BCHandler::HINT_BC_NONE );
    BCHandler bcHp( 3, BCHandler::HINT_BC_NONE );
    BCFunctionBase uZero( d->get_u3DZero() );
    BCFunctionBase uIn  ( d->get_u3DIn() );

    // bc for the velocity
    bcHu.addBC( "Inlet",    INLET,    Essential, Full,      uIn,   3 );
    bcHu.addBC( "Outlet",   OUTLET,   Natural,   Full,      uZero, 3 );
    bcHu.addBC( "Wall",     WALL,     Essential, Full,      uZero, 3 );
    // bc for the pressure 
    bcHp.addBC( "Inlet",    INLET,    Natural,   Scalar,    uZero    );
    bcHp.addBC( "Outlet",   OUTLET,   Essential, Scalar,    uZero    );
    bcHp.addBC( "Wall",	    WALL,     Natural,   Scalar,    uZero    );
#endif
#ifdef __CT_PRESSURE_TUBE
    // Boundary conditions for tube
    BCHandler bcHu( 3, BCHandler::HINT_BC_NONE );
    BCHandler bcHp( 3, BCHandler::HINT_BC_NONE );
    BCFunctionBase uZero( d->get_u3DZero() );
    BCFunctionBase pIn ( d->get_p3dIn() );

    // bc for the velocity
    bcHu.addBC ("Inlet",    INLET,    Natural,   Full,      uZero, 3);
    bcHu.addBC ("Outlet",   OUTLET,   Natural,   Full,      uZero, 3);
    bcHu.addBC ("Wall",     WALL,     Essential, Full,      uZero, 3);
    // bc for the pressure
    bcHp.addBC ("Inlet",    INLET,    Essential, Scalar,    pIn);
    bcHp.addBC ("Outlet",   OUTLET,   Essential, Scalar,    uZero);
    bcHp.addBC ("Wall",     WALL,     Natural,   Scalar,    uZero);
#endif
#ifdef __CT_CYLINDER_CASE
    // Boundary conditions for cylinder (same as Oseen w/ test_cylinder)
    BCHandler bcHu( 5, BCHandler::HINT_BC_NONE);
    BCHandler bcHp( 5, BCHandler::HINT_BC_NONE);
    BCFunctionBase uIn ( d->get_u3Dcyl() );
    BCFunctionBase uZero ( d->get_u3DZero() );
    std::vector<ID> zComp(1);
    zComp[0] = 3;

    // bc for the velocity
    bcHu.addBC( "Inlet",    CYL_INLET,    Essential, Full,      uIn,   3);
    bcHu.addBC( "Outlet",   CYL_OUTLET,   Natural,   Full,      uZero, 3);
    bcHu.addBC( "Wall",     CYL_WALL,     Essential, Full,      uZero, 3);
    bcHu.addBC( "SlipWall", CYL_SLIPWALL, Essential, Component, uZero, zComp);
    bcHu.addBC( "Cylinder", CYL_CYLINDER, Essential, Full,      uZero, 3);
    // bc for the pressure
    bcHp.addBC( "Inlet",    CYL_INLET,    Natural,   Scalar,     uZero);
    bcHp.addBC( "Outlet",   CYL_OUTLET,   Essential, Scalar,     uZero);
    bcHp.addBC( "Wall",     CYL_WALL,     Natural,   Scalar,     uZero);
    bcHp.addBC( "SlipWall", CYL_SLIPWALL, Natural,   Scalar,     uZero);
    bcHp.addBC( "Cylinder", CYL_CYLINDER, Natural,   Scalar,     uZero);

#endif

    // fluid solver

    const RefFE*    refFE_vel;
    const QuadRule* qR_vel;
    const QuadRule* bdQr_vel;

    const RefFE*    refFE_press;
    const QuadRule* qR_press;
    const QuadRule* bdQr_press;

    DataNavierStokes<RegionMesh3D<LinearTetra> > dataNavierStokes( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> > meshPart(*dataNavierStokes.mesh(), *d->comm);

    // fill in the space and time discretization orders
    std::string uOrder = dataFile( "fluid/discretization/vel_order", "P1");

    if ( uOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "  t-  P2 velocity " << std::flush;
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( uOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "  t-  P1 velocity ";
            refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;  // DoE 2
            bdQr_vel  = &quadRuleTria3pt;   // DoE 2
        }
        else
            if ( uOrder.compare("P1Bubble") == 0 )
            {
                if (verbose) std::cout << "  t-  P1-bubble velocity " << std::flush;
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;  // DoE 2
                bdQr_vel  = &quadRuleTria3pt;   // DoE 2
            }

    Dof uDof(*dataNavierStokes.mesh(), *refFE_vel);

    std::string pOrder =  dataFile( "fluid/discretization/press_order", "P1");
    if ( pOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 pressure " << std::flush;
        refFE_press = &feTetraP2;
        qR_press    = &quadRuleTetra15pt; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( pOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 pressure";
            refFE_press = &feTetraP1;
            qR_press    = &quadRuleTetra4pt;  // DoE 2
            bdQr_press  = &quadRuleTria3pt;   // DoE 2
        }

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "  t-  Time discretization order : " << dataNavierStokes.order_bdf() << std::endl;

    dataNavierStokes.setMesh(meshPart.mesh());

    // building velocity and pressure FE spaces
    if (verbose)
        std::cout << "  t-  Building the velocity FE space ... " << std::flush;
    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,
                                                             *refFE_vel,
                                                             *qR_vel,
                                                             *bdQr_vel,
                                                             3,
                                                             *d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "  t-  Building the pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > pFESpace(meshPart,
                                                             *refFE_press,
                                                             *qR_press,
                                                             *bdQr_press,
                                                             1,
                                                             *d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;



    UInt totalVelDof   = uFESpace.map().getMap(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpace.map().getMap(Unique)->NumGlobalElements();


    if (verbose) std::cout << "  t-  Total Velocity Dof = " << totalVelDof << std::endl;
    if (verbose) std::cout << "  t-  Total Pressure Dof = " << totalPressDof << std::endl;

    if (verbose) std::cout << "  t-  Calling the fluid constructor ... ";

    ChorinTemam< RegionMesh3D<LinearTetra> > fluid (dataNavierStokes,
                                              uFESpace,
                                              pFESpace,
                                              bcHu,
					      bcHp,
                                              *d->comm);
    EpetraMap fullMap_u(fluid.getMap_u());
    EpetraMap fullMap_p(fluid.getMap_p());

    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "  tt- Setting up the NS data ...";
    fluid.setUp(dataFile);
    if (verbose) std::cout << "  ok." << std::endl;

    // sync
    MPI_Barrier(MPI_COMM_WORLD);

    // Initialization
    Real dt     = dataNavierStokes.timestep();
    Real t0     = dataNavierStokes.inittime();
    Real tFinal = dataNavierStokes.endtime ();


    // initialization with stokes solution

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "  tt- Computing the stokes solution ... " << std::endl << std::endl;

    dataNavierStokes.setTime(t0);

    vector_type init_u ( fullMap_u );
    vector_type init_p ( fullMap_p );

    MPI_Barrier(MPI_COMM_WORLD);

    init_u *= 0.;
    init_p *= 0.;

    fluid.initialize(init_u, init_p);

    Ensight<RegionMesh3D<LinearTetra> > ensight( dataFile, meshPart.mesh(), "tube", d->comm->MyPID());

    vector_ptrtype vel ( new vector_type(fluid.solution_u(), Repeated ) );
    vector_ptrtype press ( new vector_type(fluid.solution_p(), Repeated ) );

    ensight.addVariable( ExporterData::Vector, "velocity", vel,
                         UInt(0), uFESpace.dof().numTotalDof() );

    ensight.addVariable( ExporterData::Scalar, "pressure", press,
                         UInt(0),
                         UInt(pFESpace.dof().numTotalDof()) );
    ensight.postProcess( 0 );

    // Temporal loop

    Chrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        dataNavierStokes.setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "l-  We are now at time "<< dataNavierStokes.time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

	fluid.time_advance(time);

        fluid.iterate_u(bcHu);
	fluid.iterate_p(bcHp);

        *vel = fluid.solution_u();
        *press = fluid.solution_p();
        ensight.postProcess( time );

        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "\n l-  Total iteration time : " << chrono.diff() << " s." << std::endl;
    }

}
