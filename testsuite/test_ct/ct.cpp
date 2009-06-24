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

    int uBdfOrder = dataFile( "fluid/time_discretization/BDF_order_vel", 1 );
    int pBdfOrder = dataFile( "fluid/time_discretization/BDF_order_press", 1 );

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "  t-  Velocity time discretization order : "
		<< uBdfOrder << std::endl;
    if (verbose) std::cout << "  t-  Pressure time discretization order : "
		<< pBdfOrder << std::endl;

    dataNavierStokes.setMesh(meshPart.mesh());

    std::string uOrder = dataFile( "fluid/space_discretization/vel_order", "P1");
    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");

    // building velocity and pressure FE spaces
    if (verbose)
        std::cout << "  t-  Building the velocity FE space ... " << std::flush;
    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,
                                                             uOrder,
                                                             3,
                                                             *M_comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "  t-  Building the pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > pFESpace(meshPart,
                                                             pOrder,
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
    Real dt     = dataNavierStokes.getTimeStep();
    Real t0     = dataNavierStokes.getInitialTime();
    Real tFinal = dataNavierStokes.getEndTime();


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
            std::cout << "l-  We are now at time "<< dataNavierStokes.getTime() << " s. " << std::endl;
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

