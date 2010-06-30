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
   \file poiseuille.cpp
   \author F. Nobile, M. Pozzoli, C. Vergara
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
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifesolver/Oseen.hpp>

#include <iostream>

#include "resistance.hpp"


using namespace LifeV;
	using namespace std;
  std::ofstream outfile;


//cyltetra.mesh
 const int INLET    = 1;
 const int WALL     = 2;
 const int OUTLET  =3;
 const int INEDGES = 4;
 const int OUTEDGES   = 5;

Real zero_scalar(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
  return 0.0;
}


Real p0(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
  return 10;
}



Real velocity(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    switch(i) {
    case 1:
      return 0;
      break;
    case 3:
      return   20*( 1-(x*x+y*y)/(0.5*0.5) );
      break;
    case 2:
      return 0;
      break;
    }
  return 0;
  }




struct  ResistanceProblem::Private
{
  Private() :
         nu(1)//,
         {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

   double Re;

   std::string data_file_name;

   double nu;  /**< viscosity (in cm^2/s) */
   double L;   /**< height and width of the domain (in cm) */
   double R;   /**< radius of the cylinder (in cm) */
   double resistance;
   Epetra_Comm* comm;
};



 ResistanceProblem:: ResistanceProblem( int argc,
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

     d->nu = dataFile( "fluid/physics/viscosity", 0.03 ) /
     dataFile( "fluid/physics/density", 1. );
     d->L = dataFile( "fluid/problem/L", 5. );
     d->R  = 0.5;
     Real  pi = 3.14159265358979 ;

     d->resistance = 8*d->nu*d->L/(pi*pow(d->R,4));


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

      std::cout <<"\n######################\n"
	        << " viscosity = "<< d->nu <<"\n"
	        << " lenght = "<< d->L <<"\n"
	        << " radius = "<< d->R <<"\n"
		<< " resistance = "<< d->resistance <<"\n"
		<<"######################\n \n";
    }

}

void
ResistanceProblem::run()

{
    Chrono chronoSet, chrono;

    chronoSet.start();
    typedef Oseen< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    // typedef boost::shared_ptr<BCVectorInterface>   bc_vector_interface;
    // Reading from data file

    GetPot dataFile( d->data_file_name.c_str() );

    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (d->comm->MyPID() == 0);


    // Boundary conditions
    BCHandler bcH( 5, BCHandler::HINT_BC_NONE );

    // fluid solver

    const RefFE*    refFE_vel;
    const QuadRule* qR_vel;
    const QuadRule* bdQr_vel;

    const RefFE*    refFE_press;
    const QuadRule* qR_press;
    const QuadRule* bdQr_press;

    DataNavierStokes<RegionMesh3D<LinearTetra> > dataNavierStokes;
    dataNavierStokes.setup( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(*dataNavierStokes.dataMesh()->mesh(), *d->comm);

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

    Dof uDof(*dataNavierStokes.dataMesh()->mesh(), *refFE_vel);

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
    if (verbose) std::cout << "Time discretization order " << dataNavierStokes.dataTime()->getBDF_order() << std::endl;

    dataNavierStokes.dataMesh()->setMesh(meshPart.mesh());

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

    Ensight<RegionMesh3D<LinearTetra> > ensight( dataFile, meshPart.mesh(), "poiseuille", d->comm->MyPID());

    vector_ptrtype velAndPressure ( new vector_type(fluid.solution(),Repeated ) );

    ensight.addVariable( ExporterData::Vector, "velocity", velAndPressure,
                         UInt(0), uFESpace.dof().numTotalDof() );

    ensight.addVariable( ExporterData::Scalar, "pressure", velAndPressure,
                         UInt(3*uFESpace.dof().numTotalDof() ),
                         UInt(  pFESpace.dof().numTotalDof() ) );


    // Initialization

    Real dt     = dataNavierStokes.dataTime()->getTimeStep();
    Real t0     = dataNavierStokes.dataTime()->getInitialTime();
    Real tFinal = dataNavierStokes.dataTime()->getEndTime();

    // bdf object to store the previous solutions

    BdfTNS<vector_type> bdf(dataNavierStokes.dataTime()->getBDF_order());

    // initialization with stokes solution
    if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

    dataNavierStokes.dataTime()->setTime(t0);

    vector_type beta( fullMap );
    vector_type rhs ( fullMap );

    BCFunctionBase uZero(zero_scalar);

    // parabolic profile

    BCFunctionBase  uPois( velocity );


    // Resistance condition:

    vector_type bcvector(fluid.getMap(),Repeated);

    bcvector.getEpetraVector().PutScalar(0.0);

    BCVector bcResistance(bcvector,uFESpace.dof().numTotalDof(),1);

    bcResistance.setResistanceCoef(d->resistance);

    // Boundary conditions

    bcH.addBC( "Wall",    1,   Essential,  Full,  uZero,  3 );
    bcH.addBC( "InFlow",  2,   Essential,  Full,  uPois  , 3 );
    bcH.addBC( "Edges",   20,  Essential,  Full,  uZero, 3 );
    bcH.addBC( "OutFlow", 3,   Resistance, Full,  bcResistance, 3);

    MPI_Barrier(MPI_COMM_WORLD);

    beta *= 0.;
    rhs  *= 0.;

    //initiliation fluid
    fluid.initialize(velocity, zero_scalar);
    //  fluid.initialize(zero_scalar, zero_scalar);

    bdf.bdf_u().initialize_unk(fluid.solution());

    fluid.setUp(dataFile);

    fluid.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);

    fluid.resetPrec();

    // Temporal loop

    int iter = 1, ii = 0;

    chronoSet.stop();

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {
        dataNavierStokes.dataTime()->setTime(time);

        chrono.start();

        double alpha = bdf.bdf_u().coeff_der( 0 ) / dataNavierStokes.dataTime()->getTimeStep();

        beta = bdf.bdf_u().extrap();
        rhs  = fluid.matrMass()*bdf.bdf_u().time_der( dataNavierStokes.dataTime()->getTimeStep() );

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdf_u().shift_right( fluid.solution() );
        *velAndPressure = fluid.solution();

        ensight.postProcess( time );

        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}





