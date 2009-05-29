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
#include <life/lifefem/bdf_template.hpp>
#include <life/lifefilters/ensight.hpp>

#include <ChorinTemam.hpp>

#include <ct.hpp>
#include <iostream>
#include <string>

// Include user specific test study
#include <ctUserCase.hpp>

using namespace LifeV;

CT::CT( int argc,
        char** argv,
        LifeV::AboutData const& /*ad*/,
        LifeV::po::options_description const& /*od*/ )
    :
    C_case (new CTcaseUser)
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

#ifdef EPETRA_MPI
    M_comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    std::cout << "  t-  MPI Initialization from PID = " << M_comm->MyPID()
	<< " among " << ntasks << " running." << std::endl;
#else
    M_comm = new Epetra_SerialComm();
#endif

    C_case->set_base_data(dataFile, M_comm);
    C_case->set_user_data();
    C_case->create_bcs();
    C_case->set_bcs();
}

/*
 * CT::run()
 */

void
CT::run()
{

    typedef ChorinTemam< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    // Reading from data file
    GetPot dataFile( C_case->get_data_hdl() );

    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (M_comm->MyPID() == 0);

    // retrieve boundary conditions from the CTcase
    boost::shared_ptr<BCHandler> bcHu = C_case->get_bcHu();
    boost::shared_ptr<BCHandler> bcHp = C_case->get_bcHp();


    // fluid solver

    const RefFE*    refFE_vel;
    const QuadRule* qR_vel;
    const QuadRule* bdQr_vel;

    const RefFE*    refFE_press;
    const QuadRule* qR_press;
    const QuadRule* bdQr_press;

    DataNavierStokes<RegionMesh3D<LinearTetra> > dataNavierStokes( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> > meshPart(*dataNavierStokes.mesh(), *M_comm);

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
            // qR_press    = &quadRuleTetra4pt;  // DoE 2
            qR_press    = qR_vel;    // test purpose
	    // because we need same qrule for u and p wrt coupling CT terms
	    // bdQr_press  = &quadRuleTria3pt;   // DoE 2
            bdQr_press  = bdQr_vel;	 // test purpose
        }

    int uBdfOrder = dataFile( "fluid/discretization/vel_order_bdf", 1 );
    int pBdfOrder = dataFile( "fluid/discretization/press_order_bdf", 1 );

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "  t-  Velocity time discretization order : "
		<< uBdfOrder << std::endl;
    if (verbose) std::cout << "  t-  Pressure time discretization order : "
		<< pBdfOrder << std::endl;

    dataNavierStokes.setMesh(meshPart.mesh());

    // building velocity and pressure FE spaces
    if (verbose)
        std::cout << "  t-  Building the velocity FE space ... " << std::flush;
    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,
                                                             *refFE_vel,
                                                             *qR_vel,
                                                             *bdQr_vel,
                                                             3,
                                                             *M_comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "  t-  Building the pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > pFESpace(meshPart,
                                                             *refFE_press,
                                                             *qR_press,
                                                             *bdQr_press,
                                                             1,
                                                             *M_comm);

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
                                              uBdfOrder,
                                              pBdfOrder,
                                              *bcHu,
                      		              *bcHp,
                                              *M_comm);

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

    Ensight<RegionMesh3D<LinearTetra> > ensight( dataFile, meshPart.mesh(), "testCT", M_comm->MyPID());

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

        fluid.iterate_u(*bcHu);
	fluid.iterate_p(*bcHp);

        *vel = fluid.solution_u();
        *press = fluid.solution_p();
        ensight.postProcess( time );

        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "\n l-  Total iteration time : " << chrono.diff() << " s." << std::endl;
    }

}

