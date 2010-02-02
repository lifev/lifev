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
   \author M.Pozzoli & C. Vergara
   \date 02-02-2010
 
   The test allows to impose the aortic flow rate at the inlet 
of a domain, through the introduction of a Lagrange multiplier.
The scheme consists in solving 2 Navier-Stokes problems at each time step.
For further details see 

Veneziani A., Vergara,  Flow rate defective Boundary Conditions in Haemodynamics Simulations, 
Int. Journ. Num. Meth. Fluids, 47, pp. 803--816, 2005,

and 

Vergara C., Ph.D. thesis.
*/

// #include <lifeconfig.h>

// #include <life/lifecore/life.hpp>
// #include <life/lifecore/GetPot.hpp>
// #include <life/lifecore/debug.hpp>

// #incLUDE <Life/lifefilters/importer.hpp>

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

#include <cylinder_1flux.hpp>
#include <iostream>

using namespace LifeV;
using namespace std;



Real lambda0 = -1.0;
Real lambda1 = -1.0;
Real lambdan = -1.0;
Real lambdan_old = lambdan;
Real pi = 3.14159265358979;



Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}



// Aortic flux: one chooses the frequency and the peak value
//
  Real a0 = 23.75;
  Real a1= 34.52;
  Real a2 = -30.81;
  Real a3 =- 12.2;
  Real a4= 4.452;
  Real f1 =-1.8;
  Real f2 = -0.49;
  Real f3 = 4.081;
  Real f4 = 5.77;
  Real ff = 1.25; // frequency
  Real w1 = 2 * ff * pi;
  Real Qmax = 0.165; // Flux at the systole (peak value)



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
    
  
  
};

Cylinder::Cylinder( int argc,
                    char** argv,
                    LifeV::AboutData const& /*ad*/,
                    LifeV::po::options_description const& /*od*/ )
    :
    d( new Private )
{
    GetPot command_line(argc, argv);
    //const char* data_file_name = command_line.follow("data", 2, "-f", "--file");
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
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
//    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
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

//    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (d->comm->MyPID() == 0);

    
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
    if (verbose) std::cout << "Time discretization order " << dataNavierStokes.getBDF_order() << std::endl;

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
                                              //bcH,
                                              *d->comm);
    EpetraMap fullMap(fluid.getMap());

    

    // Boundary conditions
    BCHandler bcH( 5, BCHandler::HINT_BC_NONE );
    BCFunctionBase uZero( zero_scalar );
    vector_type press_lagr(fluid.getMap(),Repeated); // Lagrange multiplier
    press_lagr.getEpetraVector().PutScalar(lambdan);
    BCVector uIn(press_lagr, fluid.velFESpace().dof().numTotalDof(), 1 );

    //cylinder
    //
    bcH.addBC( "Inlet",   1,   Natural,   Full,   uIn,   3 );
    bcH.addBC( "Outlet",   3,   Natural,   Full,    uZero, 3 );
    bcH.addBC( "Wall",     2,     Essential, Full,    uZero, 3 );
    bcH.addBC( "Slipwall", 4, Essential, Full, uZero, 3 );
    bcH.addBC( "Cylinder", 5, Essential, Full,      uZero, 3 );
    
    if (verbose) std::cout << "ok." << std::endl;

    fluid.setUp(dataFile);
    fluid.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);

    // Initialization

    Real dt     = dataNavierStokes.getTimeStep();
    Real t0     = dataNavierStokes.getInitialTime();
    Real tFinal = dataNavierStokes.getEndTime();


    // bdf object to store the previous solutions

    BdfTNS<vector_type> bdf(dataNavierStokes.getBDF_order());


    // initialization with stokes solution

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

    dataNavierStokes.setTime(t0);

    vector_type beta( fullMap );
    vector_type rhs ( fullMap );
    vector_type rhs0 ( fullMap );

    MPI_Barrier(MPI_COMM_WORLD);

    beta *= 0.;
    rhs  *= 0.;
    rhs0  *= 0.;

    fluid.initialize( zero_scalar, zero_scalar );

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

    Real QQn, QQ0;
    Real QQ;

    
    int jj = 1;
    int num_save = 1; // frequency in saving the results

    
   for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

      
      // Aortic flux  
      //
      QQ = -Qmax/105.35*(a0+a1*cos(w1*time+f1)+a2*cos(2*w1*time+f2)+a3*cos(3*w1*time+f3)+a4*cos(4*w1*time+f4)); 

      
      std::cout << "IMPOSED FLUX = " << -QQ << std::endl; 

        dataNavierStokes.setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << "We are now at time "<< dataNavierStokes.getTime() << " s. " << std::endl;
            std::cout << std::endl;
        }

        chrono.start();

	double alpha = bdf.bdf_u().coeff_der( 0 ) / dataNavierStokes.getTimeStep();

        beta = bdf.bdf_u().extrap();

	press_lagr.getEpetraVector().PutScalar(lambdan);
	rhs  = fluid.matrMass()*bdf.bdf_u().time_der( dataNavierStokes.getTimeStep() );

	fluid.setPrec(false);  // Compute new Preconditioner
        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

   	QQn = fluid.flux(1);
	
	vector_type X1(fluid.solution()); // Memorize solution of pb 1
	
        lambdan = lambda0;
	press_lagr.getEpetraVector().PutScalar(lambdan);
      
        fluid.setPrec(true); // Use the previous Preconditioner
        fluid.updateRHS(rhs0);
        fluid.iterate( bcH );
	
	QQ0 = fluid.flux(1);
	
	vector_type X2(fluid.solution()); // Memorize solution of pb 2

        // compute the variables to update the Lagrange multiplier
	//
        Real r0 = QQn - QQ;
	Real v;
	Real absr0;
	if(r0>0)
	  {
	    absr0=r0;
	    v=1;
	  }
	else
	  {
	    absr0=-r0;
	    v=-1;
	  }
	Real y = absr0 / QQ0;
	Real z = v * y;
	
	lambdan = lambdan_old + z; // Compute the Lagrange multplier
	std::cout << "lambda = " << lambdan << std::endl;

        lambdan_old = lambdan;
	
        vector_type X3(fluid.getMap(), Repeated); // Compute the true solution 
        X3 = X2;
        X3 *= (-z);
        X3 += X1;
               
        fluid.initialize(X3);

	QQn = fluid.flux(1); // Compute the numerical flux
        
	
	std::cout << " IMPOSED FLUX = " << QQ << std::endl;
	std::cout << " NUMERICAL FLUX = " << QQn << std::endl;

	bdf.bdf_u().shift_right( fluid.solution() );

//         if (((iter % save == 0) || (iter == 1 )))
//         {
        *velAndPressure = fluid.solution();

        if ( jj == num_save ){
        ensight.postProcess( time );
	//	        fluid.postProcess();
//         }
//	       postProcessFluxesPressures(fluid, time, verbose);
	jj = 0;
	}

        jj = jj + 1;

        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;

      
    }

}


//////////////////////


