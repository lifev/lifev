/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-19

  Copyright (C) 2005 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file ethiersteiman.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-19
 */

// #include <lifeconfig.h>

// #include <life/lifecore/life.hpp>
// #include <life/lifecore/GetPot.hpp>
// #include <life/lifecore/debug.hpp>

// #include <life/lifefilters/importer.hpp>

//#include "NavierStokesSolverBlockIP.hpp"

//#include "Epetra_SerialComm.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
//#include "life/lifesolver/NavierStokesSolver.hpp"
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
//#include <life/lifefilters/hdf5exporter.hpp>

#include <life/lifesolver/Oseen.hpp>

#include <iostream>

#include "ethiersteinman.hpp"
#include "ESSteady_function.hpp"
#include "ESUnsteady_function.hpp"


using namespace LifeV;


const int INLET    = 2;
const int WALL     = 1;
const int SLIPWALL = 20;
const int OUTLET   = 3;



Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}



std::set<UInt> parseList( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    std::string::size_type commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );
        setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
        stringList = stringList.substr( commaPos+1 );
    }
    setList.insert( atoi( stringList.c_str() ) );
    return setList;
}



struct Ethiersteinman::Private
{
    Private() :
        nu    (1),
        steady(0)
        {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    double         Re;

    std::string    data_file_name;

    double         nu;  /**< viscosity (in m^2/s) */
    //const double rho; /**< density is constant (in kg/m^3) */

    bool           steady;
    Epetra_Comm*   comm;
};








Ethiersteinman::Ethiersteinman( int argc,
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

    d->Re = dataFile( "fluid/problem/Re", 1. );
    d->nu = dataFile( "fluid/physics/viscosity", 1. ) /
        dataFile( "fluid/physics/density", 1. );

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    d->comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm = new Epetra_SerialComm();
#endif

    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "Re = " << d->Re << std::endl
                  << "nu = " << d->nu << std::endl;
    }
}

void
Ethiersteinman::run()
{
    typedef Oseen< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    // Reading from data file
    //
    GetPot dataFile( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);



    // Problem definition
    typedef EthierSteinmanUnsteady Problem;
    Problem::setParamsFromGetPot( dataFile );

    // Boundary conditions
    std::string dirichletList = dataFile( "fluid/problem/dirichletList", "" );
    std::set<UInt> dirichletMarkers = parseList( dirichletList );
    std::string neumannList = dataFile( "fluid/problem/neumannList", "" );
    std::set<UInt> neumannMarkers = parseList( neumannList );

    BCHandler::BCHints hint = neumannMarkers.size() != 0 ?
        BCHandler::HINT_BC_NONE : BCHandler::HINT_BC_ONLY_ESSENTIAL;
    BCHandler bcH( 0, hint );
    BCFunctionBase uWall( Problem::uexact );
    BCFunctionBase uNeumann( Problem::fNeumann );

    for (std::set<UInt>::const_iterator it = dirichletMarkers.begin();
         it != dirichletMarkers.end(); ++it)
    {
        bcH.addBC( "Wall", *it, Essential, Full, uWall, 3 );
    }
    for (std::set<UInt>::const_iterator it = neumannMarkers.begin();
         it != neumannMarkers.end(); ++it)
    {
        bcH.addBC( "Flux", *it, Natural, Full, uNeumann, 3 );
    }


    // fluid solver



    DataNavierStokes<RegionMesh3D<LinearTetra> > dataNavierStokes( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(*dataNavierStokes.mesh(), *d->comm);

    std::string uOrder =  dataFile( "fluid/discretization/vel_order", "P1");

    const RefFE*    refFE_vel(0);
    const QuadRule* qR_vel(0);
    const QuadRule* bdQr_vel(0);

    const RefFE*    refFE_press;
    const QuadRule* qR_press;
    const QuadRule* bdQr_press;

    if ( uOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 velocity " << std::flush;
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( uOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 velocity ";
            refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;  // DoE 2
            bdQr_vel  = &quadRuleTria3pt;   // DoE 2
        }
        else
            if ( uOrder.compare("P1Bubble") == 0 )
            {
                if (verbose) std::cout << "P1-bubble velocity " << std::flush;
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;  // DoE 2
                bdQr_vel  = &quadRuleTria3pt;   // DoE 2
            }

//    Dof uDof(*dataNavierStokes.mesh(), *refFE_vel);

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
    if (verbose) std::cout << "Time discretization order " << dataNavierStokes.order_bdf() << std::endl;

    dataNavierStokes.setMesh(meshPart.mesh());

    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,
                                                             *refFE_vel,
                                                             *qR_vel,
                                                             *bdQr_vel,
                                                             3,
                                                             *d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

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


    if (verbose) std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    Oseen< RegionMesh3D<LinearTetra> > fluid (dataNavierStokes,
                                              uFESpace,
                                              pFESpace,
                                              *d->comm);
    EpetraMap fullMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;

    fluid.setUp(dataFile);
    fluid.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);

    // Initialization

    Real dt     = dataNavierStokes.timestep();
    Real t0     = dataNavierStokes.inittime();
    Real tFinal = dataNavierStokes.endtime ();


    // bdf object to store the previous solutions

    BdfTNS<vector_type> bdf(dataNavierStokes.order_bdf());

    // initialization with stokes solution

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Computing the initial solution ... " << std::endl << std::endl;

    dataNavierStokes.setTime(t0);

    vector_type beta( fullMap );
    vector_type rhs ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    beta *= 0.;
    rhs  *= 0.;

    std::string const proj =  dataFile( "fluid/discretization/initialization", "proj");
    bool const L2proj( !proj.compare("proj") );
    if (L2proj)
        {
            uFESpace.L2ScalarProduct(Problem::uderexact, rhs, 0.);
            rhs *= -1.;
        }

    fluid.initialize( Problem::uexact, Problem::pexact );
    //fluid.initialize( zero_scalar, zero_scalar);

    beta = fluid.solution();

    std::cout << "norm beta " << beta.Norm2() << std::endl;
    std::cout << "norm rhs  " << rhs.Norm2() << std::endl;


    if (L2proj)
        {
            fluid.updateSystem( 0., beta, rhs );
            fluid.iterate(bcH);
        }

    //    fluid.postProcess();

    std::ofstream out_norm("norm.txt");

    out_norm << "time / u L2 error / L2 rel error   p L2 error / L2 rel error \n" << std::flush;

    vector_type vel  (uFESpace.map(), Repeated);
    vector_type press(pFESpace.map(), Repeated);

    vector_type velpressure ( fluid.solution(), Repeated );

    vel.subset(velpressure);
    press.subset(velpressure, uFESpace.dim()*uFESpace.fieldDim());

    double urelerr;
    double prelerr;

    double ul2error = uFESpace.L2Error(Problem::uexact, vel  , 0., &urelerr );
    double pl2error = pFESpace.L20Error(Problem::pexact, press, 0., &prelerr );

//         if (verbose)
//         {
    out_norm << 0. << " "
             << ul2error << " "
             << urelerr << " "
             << pl2error << " "
             << prelerr << "\n" << std::flush;



//     fluid.updateSystem(0, beta, rhs );
//     fluid.iterate();


    bdf.bdf_u().initialize_unk( fluid.solution() );

    fluid.resetPrec();

    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;

    exporter.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "ethiersteinman", d->comm->MyPID()) );
    // hdf5 exporter, still under development
    //    exporter.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "ethiersteinman", d->comm->MyPID()) );

    //    Ensight<RegionMesh3D<LinearTetra> > exporter( dataFile, meshPart.mesh(), "ethiersteinman", d->comm->MyPID());

    vector_ptrtype velAndPressure ( new vector_type(fluid.solution(), Repeated ) );

    exporter->addVariable( ExporterData::Vector, "velocity", velAndPressure,
                         UInt(0), uFESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::Scalar, "pressure", velAndPressure,
                         UInt(3*uFESpace.dof().numTotalDof()),
                         UInt(3*uFESpace.dof().numTotalDof()+pFESpace.dof().numTotalDof()) );
    exporter->postProcess( 0 );

    // Temporal loop

    Chrono chrono;
    int iter = 1;


    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        dataNavierStokes.setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time "<< dataNavierStokes.time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        double alpha = bdf.bdf_u().coeff_der( 0 ) / dataNavierStokes.timestep();

        beta = bdf.bdf_u().extrap();

        rhs  = fluid.matrMass()*bdf.bdf_u().time_der( dataNavierStokes.timestep() );
//        rhs *= alpha;
//        rhs  = bdf.bdf_u().time_der( dataNavierStokes.timestep() );

        std::cout << "alpha " << alpha << std::endl;
        std::cout << "norm beta " << beta.Norm2() << std::endl;
        std::cout << "norm rhs  " << rhs.Norm2() << std::endl;

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdf_u().shift_right( fluid.solution() );

        velpressure   = fluid.solution();

        vector_type vel  (uFESpace.map(), Repeated);
        vector_type press(pFESpace.map(), Repeated);

        vel.subset(velpressure);
        press.subset(velpressure, uFESpace.dim()*uFESpace.fieldDim());


        ul2error = uFESpace.L2Error(Problem::uexact, vel  , time, &urelerr );
        pl2error = pFESpace.L20Error(Problem::pexact, press, time, &prelerr );

//         if (verbose)
//         {
        out_norm << time << " "
                 << ul2error << " "
                 << urelerr << " "
                 << pl2error << " "
                 << prelerr << "\n" << std::flush;
//         }
//         if (((iter % save == 0) || (iter == 1 )))
//         {
        *velAndPressure = fluid.solution();
        exporter->postProcess( time );
        //            fluid.postProcess();
//         }


        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//////////////////////


