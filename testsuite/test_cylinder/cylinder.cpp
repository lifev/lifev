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
   \file cylinder.cpp
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

#include <life/lifesolver/Oseen.hpp>

#include <cylinder.hpp>
#include <iostream>

using namespace LifeV;


//cylinder

// //#define CYL2D_MESH_SETTINGS
#define CYL3D_2_MESH_SETTINGS
//#define TUBE20_MESH_SETTINGS


// #ifdef CYL2D_MESH_SETTINGS // cyl2D.mesh
// const int INLET    = 1;
// const int WALL     = 1;
// const int SLIPWALL = 20;
// const int OUTLET   = 0;
// const int CYLINDER = 2;
// #endif
#ifdef CYL3D_2_MESH_SETTINGS // cyl3D-2.mesh
const int INLET    = 40;
const int WALL     = 60;
const int SLIPWALL = 61;
const int OUTLET   = 50;
const int CYLINDER = 70;
#endif
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

Real u2(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
  switch(i) {
  case 1:
    return 0.0;
    break;
  case 3:
       if ( t <= 0.003 )
          return 1.3332e4;
//      return 0.01;
      return 0.0;
      break;
  case 2:
      return 0.0;
//      return 1.3332e4;
//    else
//      return 0.0;
    break;
  }
  return 0;
}



struct Cylinder::Private
{
    Private() :
        //check(false),
        nu(1),
        //rho(1),
        H(1), D(1)
        //H(20), D(1)
        //H(0.41), D(0.1)
        {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    double Re;

    std::string data_file_name;

    double nu;  /**< viscosity (in m^2/s) */
    //const double rho; /**< density is constant (in kg/m^3) */
    double H;   /**< height and width of the domain (in m) */
    double D;   /**< diameter of the cylinder (in m) */
    bool centered; /**< true if the cylinder is at the origin */

    Epetra_Comm*   comm;
    /**
     * get the characteristic velocity
     *
     * @return the characteristic velocity
     */
    double Ubar() const { return nu*Re/D; }

    /**
     * get the magnitude of the profile velocity
     *
     *
     * @return the magnitude of the profile velocity
     */
    double Um_3d() const { return 9*Ubar()/4; }

    double Um_2d() const { return 3*Ubar()/2; }


    /**
     * u3d 3D velocity profile.
     *
     * Define the velocity profile at the inlet for the 3D cylinder
     */
    Real u3d( const Real& /* t */,
              const Real& /* x */,
              const Real& y,
              const Real& z,
              const ID&   id ) const
        {
            if ( id == 1 ) {
                if ( centered ) {
                    return Um_3d() * (H+y)*(H-y) * (H+z)*(H-z) / pow(H,4);
                } else {
                    return 16 * Um_3d() * y * z * (H-y) * (H-z) / pow(H,4);
                }
            } else {
                return 0;
            }
        }

    fct_type getU_3d()
        {
            fct_type f;
            f = boost::bind(&Cylinder::Private::u3d, this, _1, _2, _3, _4, _5);
            return f;
        }

    /**
     * u2d flat 2D velocity profile.
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real u2d( const Real& /*t*/,
              const Real& /*x*/,
              const Real& y,
              const Real& /*z*/,
              const ID&   id ) const
        {
             if ( id == 1 )
              {
                  return 1./(20.*20.)*(y + 20.)*(20. - y);
// 	         if ( centered ) {
//                      return Um_2d() * (y+H)*(H-y) / (H*H);
//                  } else {
//                      return 4 * Um_2d() * y * (H-y) / (H*H);
//                  }
             } else {
                 return 0;
             }
        }

    fct_type getU_2d()
        {
            fct_type f;
            f = boost::bind(&Cylinder::Private::u2d, this, _1, _2, _3, _4, _5);
            return f;
        }

    /**
     * one flat (1,1,1)
     *
     * Define the velocity profile at the inlet for the 2D cylinder
     */
    Real oneU( const Real& /*t*/,
               const Real& /*x*/,
               const Real& y,
               const Real& /*z*/,
               const ID&   id ) const
        {
            return 1.;
        }

    fct_type getU_one()
        {
            fct_type f;
            f = boost::bind(&Cylinder::Private::oneU, this, _1, _2, _3, _4, _5);
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
    const char* data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    d->data_file_name = data_file_name;

    d->Re = dataFile( "fluid/problem/Re", 1. );
    d->nu = dataFile( "fluid/physics/viscosity", 1. ) /
        dataFile( "fluid/physics/density", 1. );
    d->H  = 20.;//dataFile( "fluid/problem/H", 20. );
    d->D  = dataFile( "fluid/problem/D", 1. );
    d->centered = (bool)dataFile( "fluid/problem/centered", 0 );

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    d->comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm = new Epetra_SerialComm();
#endif

    if (!d->comm->MyPID()) {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
        std::cout << "Re = " << d->Re << std::endl
                  << "nu = " << d->nu << std::endl
                  << "H  = " << d->H  << std::endl
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

    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (d->comm->MyPID() == 0);

    // Boundary conditions
    BCHandler bcH( 5, BCHandler::HINT_BC_NONE );
    BCFunctionBase uZero( zero_scalar );
    std::vector<ID> zComp(1);
    zComp[0] = 3;


#ifdef CYL3D_2_MESH_SETTINGS // cyl3D-2.mesh
    BCFunctionBase uIn  (  d->getU_2d() );

    //cylinder
    bcH.addBC( "Inlet",    INLET,    Essential,   Full,      uIn,   3 );
    bcH.addBC( "Outlet",   OUTLET,   Natural,   Full,      uZero, 3 );
//     if ( WALL != INLET )
    bcH.addBC( "Wall",     WALL,     Essential, Full,      uZero, 3 );
    bcH.addBC( "Slipwall", SLIPWALL, Essential, Component, uZero, zComp );
    bcH.addBC( "Cylinder", CYLINDER, Essential, Full,      uZero, 3 );
//    bcH.addBC( "Slipwall", SLIPWALL, Essential, Full, uZero , 3 );
#endif
#ifdef TUBE20_MESH_SETTINGS
    BCFunctionBase uIn  (  d->getU_one() );
    //BCFunctionBase unormal(  d->get_normal() );

    //cylinder
    bcH.addBC( "Inlet",    INLET,    Essential,   Full,      uIn, 3 );
    bcH.addBC( "Outlet",   OUTLET,   Essential,   Full,      uIn, 3 );

    //bcH.addBC( "Wall",     WALL,     Natural,     Full,      uNormal, 3 );
    //    bcH.addBC( "Wall",     WALL,     Natural,     Full,      uZero, 3 );

    bcH.addBC( "Slipwall", SLIPWALL, Essential,   Full,      uIn, 3 );
#endif


//    bcH.showMe();

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
                                              bcH,
                                              *d->comm);
    EpetraMap fullMap(fluid.getMap());


#ifdef TUBE20_MESH_SETTINGS
    vector_type vec_lambda_aux(fullMap),	mixtevec_aux(fullMap);
    vec_lambda_aux.getEpetraVector().PutScalar(1);
    vec_lambda_aux.GlobalAssemble();
    mixtevec_aux.getEpetraVector().PutScalar(1);
    mixtevec_aux.GlobalAssemble();

    vector_type vec_lambda(fluid.getMap(), Repeated), mixtevec(fluid.getMap(), Repeated);

    vec_lambda = vec_lambda_aux;
    mixtevec = mixtevec_aux;

    // Robin BC
    BCVector robin_wall(vec_lambda, uFESpace.dof().numTotalDof());

    // Neumann BC, normal component
    //BCVector robin_wall(vec_lambda, uFESpace.dof().numTotalDof(),1);

    robin_wall.setMixteCoef(0);
    robin_wall.setMixteVec(mixtevec);

    vec_lambda_aux.spy("lambda");
    vec_lambda.spy("lambdaRep");



    bcH.addBC( "Wall",      1, Mixte, Full,  robin_wall, 3 );
    //bcH.addBC( "Wall",      1, Natural, Full,  robin_wall, 3 );
#endif




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

        double alpha = bdf.bdf_u().coeff_der( 0 ) / dataNavierStokes.timestep();

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


