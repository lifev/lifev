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
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif
//#include "life/lifesolver/NavierStokesSolver.hpp"
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifefilters/noexport.hpp>


#include <iostream>

#include "ethiersteinman.hpp"


using namespace LifeV;



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
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    d->data_file_name = data_file_name;

    d->Re = dataFile( "fluid/problem/Re", 1. );
    d->nu = dataFile( "fluid/physics/viscosity", 1. ) /
        dataFile( "fluid/physics/density", 1. );

#ifdef EPETRA_MPI

    //    MPI_Init(&argc,&argv);

    d->comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm = new Epetra_SerialComm();
#endif

    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "Re = " << d->Re << std::endl
                  << "nu = " << d->nu << std::endl;

        out_norm.open("norm.txt");
        out_norm << "% time / u L2 error / L2 rel error   p L2 error / L2 rel error \n" << std::flush;
    }

}


void
Ethiersteinman::computeError( double     const& time,
                              fespace_type&     uFESpace,
                              fespace_type&     pFESpace,
                              fluid_type const& fluid )
{
        vector_type vel  (uFESpace.map(), Repeated);
        vector_type press(pFESpace.map(), Repeated);
        vector_type velpressure ( fluid.solution(), Repeated );

        velpressure = fluid.solution();
        vel.subset(velpressure);
        press.subset(velpressure, uFESpace.dim()*uFESpace.fieldDim());

        double urelerr;
        double prelerr;
        double ul2error;
        double pl2error;

        ul2error = uFESpace.L2Error (Problem::uexact, vel  , time, &urelerr );
        pl2error = pFESpace.L20Error(Problem::pexact, press, time, &prelerr );


        bool verbose = (d->comm->MyPID() == 0);

         if (verbose)
         {
             out_norm << time  << " "
                      << ul2error << " "
                      << urelerr << " "
                      << pl2error << " "
                      << prelerr << "\n" << std::flush;
         }
         checkResult(time, ul2error, pl2error);

}

void Ethiersteinman::checkResult( LifeV::Real const& time,
                                  double const&      ul2error,
                                  double const&      pl2error )
{

    double ul2stored(0);
    double pl2stored(0);

    if (std::abs(time) < 1e-8)
    {
        ul2stored =  0.0434792;
        pl2stored =  0.11569;
    }
    if (std::abs(time - 0.01) < 1e-8)
    {
        ul2stored = 0.075432;
        pl2stored =  1.92917;
    }
    if (std::abs(time - 0.02) < 1e-8)
    {
        ul2stored =  0.0936949;
        pl2stored =  0.706018;
    }
    if (std::abs(time - 0.03) < 1e-8)
    {
        ul2stored =  0.104034;
        pl2stored =  0.445147;
    }
    if (std::abs(time - 0.04) < 1e-8)
    {
        ul2stored =  0.110704;
        pl2stored =  0.354164;
    }

    if ( ul2stored > 0 &&
         ( std::abs(ul2stored - ul2error) > 1e-4*ul2stored || (std::abs(pl2stored - pl2error) > 1e-4*pl2stored ) ) )
         throw Ethiersteinman::RESULT_CHANGED_EXCEPTION(time);

}


void
Ethiersteinman::run()
{
    Chrono chrono;

    // Reading from data file
    //
    GetPot dataFile( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);

    chrono.start();

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

    DataNavierStokes<RegionMesh3D<LinearTetra> > dataNavierStokes;
    dataNavierStokes.setup( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(*dataNavierStokes.dataMesh()->mesh(), *d->comm);

    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Time discretization order " << dataNavierStokes.dataTime()->getBDF_order() << std::endl;

    dataNavierStokes.dataMesh()->setMesh(meshPart.mesh());

    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,uOrder,3,*d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > pFESpace(meshPart,pOrder,1,*d->comm);

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

    if (verbose)
        std::cout << d->comm->MyPID() << " Init time (partial) " << chrono.diff() << " s." << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // Initialization

    Real dt     = dataNavierStokes.dataTime()->getTimeStep();
    Real t0     = dataNavierStokes.dataTime()->getInitialTime();
    Real tFinal = dataNavierStokes.dataTime()->getEndTime ();


    // bdf object to store the previous solutions

    BdfTNS<vector_type> bdf(dataNavierStokes.dataTime()->getBDF_order());

    // initialization with exact solution: either interpolation or "L2-NS"-projection
    t0 -= dt * bdf.bdf_u().order();

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Computing the initial solution ... " << std::endl << std::endl;

    vector_type beta( fullMap );
    vector_type rhs ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    dataNavierStokes.dataTime()->setTime(t0);
    fluid.initialize( Problem::uexact, Problem::pexact );
    bdf.bdf_u().initialize_unk( fluid.solution() );

    std::string const proj =  dataFile( "fluid/space_discretization/initialization", "proj");
    bool const L2proj( !proj.compare("proj") );



    Real time = t0 + dt;
    for (  ; time <=  dataNavierStokes.dataTime()->getInitialTime() + dt/2.; time += dt)
    {

        dataNavierStokes.dataTime()->setTime(time);

        beta *= 0.;
        rhs  *= 0.;

        if (L2proj)
        {
            uFESpace.L2ScalarProduct(Problem::uderexact, rhs, time);
            rhs *= -1.;
        }

        fluid.initialize( Problem::uexact, Problem::pexact );

        beta = fluid.solution();

        fluid.getDisplayer().leaderPrint("norm beta ", beta.Norm2());
        fluid.getDisplayer().leaderPrint("\n");
        fluid.getDisplayer().leaderPrint("norm rhs  ", rhs.Norm2());
        fluid.getDisplayer().leaderPrint("\n");


        if (L2proj)
        {
            fluid.updateSystem( 0., beta, rhs );
            fluid.iterate(bcH);
        }

        computeError( time,
                      uFESpace,
                      pFESpace,
                      fluid );


         bdf.bdf_u().shift_right( fluid.solution() );

    }

    chrono.stop();
    if (verbose)
	std::cout << d->comm->MyPID() << " Total init time" << chrono.diff() << " s." << std::endl;
    // end initialization step

    fluid.resetPrec();

    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;

    vector_ptrtype velAndPressure;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
    {
        exporter.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, "ethiersteinman" ) );
        exporter->setDirectory( "./" ); // This is a test to see if M_post_dir is working
        exporter->setMeshProcId( meshPart.mesh(), d->comm->MyPID() );
    }
    else
#endif
    {
        if (exporterType.compare("none") == 0)
        {
            exporter.reset( new NoExport<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "ethiersteinman", d->comm->MyPID()) );
        } else {
            exporter.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "ethiersteinman", d->comm->MyPID()) );
        }
    }

    velAndPressure.reset( new vector_type(fluid.solution(), exporter->mapType() ) );

    exporter->addVariable( ExporterData::Vector, "velocity", velAndPressure,
                         UInt(0), uFESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::Scalar, "pressure", velAndPressure,
                         UInt(3*uFESpace.dof().numTotalDof()),
                         UInt(pFESpace.dof().numTotalDof()) );
    exporter->postProcess( 0 );

    std::cout << "uDOF: " << uFESpace.dof().numTotalDof() << std::endl;
    std::cout << "pDOF: " << pFESpace.dof().numTotalDof() << std::endl;

    // Temporal loop

    int iter = 1;

    Chrono chronoGlobal;
    chronoGlobal.start();

    for ( ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        dataNavierStokes.dataTime()->setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time "<< dataNavierStokes.dataTime()->getTime() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        double alpha = bdf.bdf_u().coeff_der( 0 ) / dataNavierStokes.dataTime()->getTimeStep();

        beta = bdf.bdf_u().extrap();

        rhs  = fluid.matrMass()*bdf.bdf_u().time_der( dataNavierStokes.dataTime()->getTimeStep() );
//        rhs *= alpha;
//        rhs  = bdf.bdf_u().time_der( dataNavierStokes.getTimeStep() );

        fluid.getDisplayer().leaderPrint("alpha ", alpha);
        fluid.getDisplayer().leaderPrint("\n");
        fluid.getDisplayer().leaderPrint("norm beta ", beta.Norm2());
        fluid.getDisplayer().leaderPrint("\n");
        fluid.getDisplayer().leaderPrint("norm rhs  ", rhs.Norm2());
        fluid.getDisplayer().leaderPrint("\n");

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdf_u().shift_right( fluid.solution() );


        computeError( time,
                      uFESpace,
                      pFESpace,
                      fluid );

        *velAndPressure = fluid.solution();
        exporter->postProcess( time );


        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

    chronoGlobal.stop();
    if (verbose) std::cout << "Total simulation time (time loop only) " << chronoGlobal.diff() << " s." << std::endl;

}


//////////////////////


