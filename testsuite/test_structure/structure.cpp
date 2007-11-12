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

//#include "NavierStokesSolverIP.hpp"
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataElasticStructure.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>

#include <life/lifefilters/ensight.hpp>


#include "structure.hpp"
#include <iostream>

using namespace LifeV;



Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    if (t == 0)
    {
        switch(i) {
            case 1:
                return -z*(z - 5.)*x/50.;
                return 0;
                break;
            case 2:
                return -z*(z - 5.)*y/50.;
                return 0;
                break;
            case 3:
                return 0.0;
                break;
            default:
                ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
                break;
        }
    }
    else
    {
        return 0.0;
    }
}

Real w0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{

  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break;
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }
}



Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}

struct Structure::Private
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
    //double h;
    //double T0;
    //double T;
    //double dt;
    //TimeScheme time_order;
    double Re;
    //bool check;

    //std::string exporter_str;
    std::string data_file_name;

    double nu;  /**< viscosity (in m^2/s) */
    //const double rho; /**< density is constant (in kg/m^3) */
    double H;   /**< height and width of the domain (in m) */
    double D;   /**< diameter of the cylinder (in m) */
    bool centered; /**< true if the cylinder is at the origin */

    Epetra_Comm*     comm;
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
    Real d3d( const Real& /* t */,
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

    fct_type getD_3d()
        {
            fct_type f;
            f = boost::bind(&Structure::Private::d3d, this, _1, _2, _3, _4, _5);
            return f;
        }

};

Structure::Structure( int                                   argc,
                      char**                                argv,
                      Epetra_Comm &                         structComm,
                      LifeV::AboutData const&               ad,
                      LifeV::po::options_description const& od ):
    d( new Private )
{
    GetPot command_line(argc, argv);
    const char* data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    d->data_file_name = data_file_name;

    d->Re = dataFile( "fluid/problem/Re", 1. );
    d->nu = dataFile( "fluid/physics/viscosity", 1. ) /
        dataFile( "fluid/physics/density", 1. );
    d->H  = dataFile( "fluid/problem/H", 1. );
    d->D  = dataFile( "fluid/problem/D", 1. );
    d->centered = (bool)dataFile( "fluid/problem/centered", 0 );

    std::cout << "Re = " << d->Re << std::endl
              << "nu = " << d->nu << std::endl
              << "H  = " << d->H  << std::endl
              << "D  = " << d->D  << std::endl;

    d->comm = &structComm;
    int ntasks = d->comm->NumProc();

    if (!d->comm->MyPID()) std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;

}

void
Structure::run2d()
{
    std::cout << "2D cylinder test case is not available yet\n";
}

void
Structure::run3d()
{

    typedef VenantKirchhofSolver< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    bool verbose = (d->comm->MyPID() == 0);
    // Number of boundary conditions for the velocity and mesh motion
    //
    BCHandler BCh(2);
    BCFunctionBase dZero( zero_scalar );
    //
    // dataElasticStructure
    //

    GetPot dataFile( d->data_file_name.c_str() );

    DataElasticStructure< RegionMesh3D<LinearTetra > >    dataStructure( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >            meshPart( *dataStructure.mesh(),
                                                                    *d->comm );

//    meshPart.rebuildMesh();


//    exit(0);
    const RefFE*    refFE(0);
    const QuadRule* qR(0);
    const QuadRule* bdQr(0);

    std::string pOrder =  dataFile( "solid/discretization/order", "P1");
    if ( pOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 displacement " << std::flush;
        refFE = &feTetraP2;
        qR    = &quadRuleTetra15pt; // DoE 5
        bdQr  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( pOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 displacement";
            refFE = &feTetraP1;
            qR    = &quadRuleTetra4pt;  // DoE 2
            bdQr  = &quadRuleTria3pt;   // DoE 2
        }

    dataStructure.setMesh(meshPart.mesh());

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > dFESpace(meshPart,
                                                             *refFE,
                                                             *qR,
                                                             *bdQr,
                                                             3,
                                                             *d->comm);
    if (verbose) std::cout << std::endl;

    EpetraMap structMap(*refFE, meshPart, *d->comm);

    EpetraMap fullMap;

    for (UInt ii = 0; ii < nDimensions; ++ii)
    {
         fullMap += structMap;
    }

    VenantKirchhofSolver< RegionMesh3D<LinearTetra> > solid( dataStructure,
                                                             dFESpace,
                                                             *d->comm);
    solid.setUp(dataFile);
    solid.buildSystem();
    //
    // Boundary conditions for the displacement
    //
    BCFunctionBase fixed(dZero);

    BCh.addBC("Base2 ", 2 , Essential, Full, fixed, 3);
    BCh.addBC("Base3 ", 3 , Essential, Full, fixed, 3);
    //
    // Temporal data and initial conditions
    //
    Real dt = dataStructure.timestep();
    Real T  = dataStructure.endtime();

    EpetraVector<double> disp(solid.disp(),*fullMap.getRepeatedEpetra_Map());
    EpetraVector<double> vel (solid.vel(),*fullMap.getRepeatedEpetra_Map());

    dFESpace.interpolate(d0, disp, 0.0);
    dFESpace.interpolate(w0, vel , 0.0);

    if (verbose) std::cout << "S- initialization ... ";

//    solid.initialize(d0, w0);
    solid.initialize(disp,vel); // displacement and velocity

    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose ) std::cout << "ok." << std::endl;
    //if (d->comm->NumProc() == 1 )  solid.postProcess();


    Ensight<RegionMesh3D<LinearTetra> > ensight( dataFile, meshPart.mesh(), "structure", d->comm->MyPID());

    vector_ptrtype solidDisp ( new vector_type(solid.disp(), solid.getRepeatedEpetraMap() ) );
    vector_ptrtype solidVel  ( new vector_type(solid.vel(), solid.getRepeatedEpetraMap() ) );

    ensight.addVariable( ExporterData::Vector, "displacement", solidDisp,
                         UInt(0), dFESpace.dof().numTotalDof() );

    ensight.addVariable( ExporterData::Vector, "velocity", solidVel,
                         UInt(0), dFESpace.dof().numTotalDof() );

    ensight.postProcess( 0 );


    //
    // Temporal loop
    //

    for (Real time = dt; time <= T; time += dt)
    {
        dataStructure.setTime(time);

        if (verbose)
            {
                std::cout << std::endl;
                std::cout << "S- Now we are at time " << dataStructure.time() << " s." << std::endl;
            }

        //solid.updateSystem(dZero);    // Computes the rigth hand side
        solid.updateSystem();    // Computes the rigth hand side
        solid.iterate( BCh );                  // Computes the matrices and solves the system
        //if (d->comm->NumProc() == 1 )  solid.postProcess(); // Post-presssing

        *solidDisp = solid.disp();
        *solidVel  = solid.vel();
        ensight.postProcess( time );

    }


}


//////////////////////


