/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
             Christian Vergara <>
       Date: 2008-10-30

  Copyright (C) 2008 EPFL

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
   \file cylinder.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2008-06-13
 */

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

#include <boost/shared_ptr.hpp>
#include <vector>

#include <life/lifesolver/Oseen.hpp>

#include <cylinder.hpp>
#include <iostream>
#include <math.h>


using namespace LifeV;


//cylinder

#define TUBE20_MESH_SETTINGS

#ifdef TUBE20_MESH_SETTINGS
 const int INLET    = 2;
 const int WALL     = 1;
 const int SLIPWALL = 20;
 const int OUTLET   = 3;
#endif



Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

Real flux(const Real& t, const ID& i)
{

    // we have to impose (q*n)
    // here : n = ( 0, 0, 1)
    Real const pi(3.14159265);

    switch(i) {
    case 1:
        return 0.0;
        break;
    case 3:
        if ( t <= 1 )
            return sin(pi*t);
        return 0.0;
        break;
    case 2:
        return 0.0;
        break;
    }
    return 0;
}



struct Cylinder::Private
{
    Private() :
        nu(1), D(1), period(1),
        lambda()
        {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    std::string data_file_name;

    Real Re;
    Real nu;  /**< viscosity (in m^2/s) */
    Real D;
    Real period;
    boost::shared_ptr< std::vector<Real> > lambda;

    void setLambda(boost::shared_ptr< std::vector<Real> >& lam) { lambda = lam;}

    Epetra_Comm*   comm;
    /**
     * get the characteristic velocity
     *
     * @return the characteristic velocity
     */
    Real Ubar() const { return nu*Re/D; }


    /**
     * u3d 3D velocity profile.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real flux3d( const ID&   id, const Real t  ) const
        {
            return (Ubar() * flux( period*t, id ));
        }

    fct_type get_flux3d()
        {
            fct_type f;
            f = boost::bind(&Cylinder::Private::flux3d, this, _1, _2);
            return f;
        }

    /**
     * u3d 3D lagrabge multiplier.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real lambda3d( const Real& t,
                 const Real& /* x */,
                 const Real& /* y */,
                 const Real& /* z */,
                 const ID&   id ) const
        {
            return ( (*lambda)[0] );
        }

    fct_type get_lambda3d()
        {
            fct_type f;
            f = boost::bind(&Cylinder::Private::lambda3d, this, _1, _2, _3, _4, _5);
            return f;
        }


};

Cylinder::Cylinder( int argc,
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
    d->period  = dataFile( "fluid/problem/period", 20. );
    d->D  = dataFile( "fluid/problem/D", 1. );

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

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
                  << "nu = " << d->nu << std::endl
                  << "period  = " << d->period  << std::endl
                  << "D  = " << d->D  << std::endl;
    }
}

void
Cylinder::run()

{

    typedef Oseen< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    // Reading from data file
    //
    GetPot dataFile( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);

    // Boundary conditions
    BCHandler bcH( 4, BCHandler::HINT_BC_NONE );
    BCFunctionBase uZero( zero_scalar );
    std::vector<ID> zComp(1);
    zComp[0] = 3;

    boost::shared_ptr< std::vector<Real> > lambda;
    lambda.reset( new std::vector<Real>(1) );
    (*lambda)[1] = INLET;
    d-> setLambda(lambda);

    BCFunctionBase lambdaIn  (  d->get_lambda3d() );


#ifdef TUBE20_MESH_SETTINGS

    //cylinder
    bcH.addBC( "Inlet",    INLET,    Natural,   Full,      lambdaIn,   3 );
    bcH.addBC( "Outlet",   OUTLET,   Natural,   Full,      uZero, 3 );
    bcH.addBC( "Wall",     WALL,     Essential,   Full,    uZero, 3 );
    bcH.addBC( "Slipwall", SLIPWALL, Essential,   Full,    uZero, 3 );

#endif


    // fluid solver

    const RefFE*    refFE_vel;
    const QuadRule* qR_vel;
    const QuadRule* bdQr_vel;

    const RefFE*    refFE_press;
    const QuadRule* qR_press;
    const QuadRule* bdQr_press;

    DataNavierStokes<RegionMesh3D<LinearTetra> > dataNavierStokes( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(*dataNavierStokes.mesh(), *d->comm);

    std::string uOrder =  dataFile( "fluid/discretization/vel_order", "P1");

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

    Dof uDof(*dataNavierStokes.mesh(), *refFE_vel);

    std::string pOrder =  dataFile( "fluid/discretization/press_order", "P1");
    if ( pOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 pressure " << std::flush;
        refFE_press = &feTetraP2;
        qR_press    = qR_vel;
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( pOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 pressure";
            refFE_press = &feTetraP1;
            qR_press    = qR_vel;
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

    // Lagrange multipliers for flux imposition
    std::vector<int> lagrangeMultipliers(0);

    /*
    if  (d->comm->MyPID() == 0) // Adding lagrange multipliers in the first processor
        {
    */
            lagrangeMultipliers.resize(1); // just one flux to impose (n);
            lagrangeMultipliers[0] = 1;    // just take a numbering of the lagrange multipliers [1:n];
    /*
        }
    */

    if (verbose) std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    Oseen< RegionMesh3D<LinearTetra> > fluid (dataNavierStokes,
                                              uFESpace,
                                              pFESpace,
                                              lagrangeMultipliers,
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
    if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

    dataNavierStokes.setTime(t0);

    vector_type beta( fullMap );
    vector_type rhs ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    beta *= 0.;
    rhs  *= 0.;

    fluid.updateSystem(0, beta, rhs );
    fluid.iterate( bcH );


//    fluid.postProcess();
    bdf.bdf_u().initialize_unk( fluid.solution() );

    fluid.resetPrec();

    Ensight<RegionMesh3D<LinearTetra> > ensight( dataFile, meshPart.mesh(), "cylinder", d->comm->MyPID());

    vector_ptrtype velAndPressure ( new vector_type(fluid.solution(), Repeated ) );

    ensight.addVariable( ExporterData::Vector, "velocity", velAndPressure,
                         UInt(0), uFESpace.dof().numTotalDof() );

    ensight.addVariable( ExporterData::Scalar, "pressure", velAndPressure,
                         UInt(3*uFESpace.dof().numTotalDof()),
                         UInt(3*uFESpace.dof().numTotalDof()+pFESpace.dof().numTotalDof()) );
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
            std::cout << "We are now at time "<< dataNavierStokes.time() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

        Real alpha = bdf.bdf_u().coeff_der( 0 ) / dataNavierStokes.timestep();

        beta = bdf.bdf_u().extrap();

        rhs  = fluid.matrMass()*bdf.bdf_u().time_der( dataNavierStokes.timestep() );

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdf_u().shift_right( fluid.solution() );

//         if (((iter % save == 0) || (iter == 1 )))
//         {
        *velAndPressure = fluid.solution();
        ensight.postProcess( time );
//        fluid.postProcess();
//         }


        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//////////////////////


