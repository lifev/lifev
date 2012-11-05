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

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifesolver/OseenSolver.hpp>

#ifdef HAVE_HDF5
#include <life/lifefilters/ExporterHDF5.hpp>
#endif
#include <life/lifefilters/ExporterEnsight.hpp>

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

Real zero_scalar(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}


Real p0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 10;
}



Real velocity(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
    case 0:
        return 0;
        break;
    case 2:
        return   20*( 1-(x*x+y*y)/(0.5*0.5) );
        break;
    case 1:
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

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_Type;

    double Re;

    std::string data_file_name;

    double nu;  /**< viscosity (in cm^2/s) */
    double L;   /**< height and width of the domain (in cm) */
    double R;   /**< radius of the cylinder (in cm) */
    double resistance;
    boost::shared_ptr<Epetra_Comm>   comm;
};



ResistanceProblem::ResistanceProblem( int argc,
                                      char** argv )
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
    d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );

    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm = new Epetra_SerialComm();
#endif

    if (!d->comm->MyPID())
    {
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
    LifeChrono chronoSet, chrono;

    chronoSet.start();
    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef OseenSolver< mesh_Type >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
    // typedef boost::shared_ptr<BCVectorInterface>   bc_vector_interface;
    // Reading from data file

    GetPot dataFile( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);


    // Boundary conditions
    BCHandler bcH;

    // fluid solver

    const ReferenceFE*    refFE_vel;
    const QuadratureRule* qR_vel;
    const QuadratureRule* bdQr_vel;

    const ReferenceFE*    refFE_press;
    const QuadratureRule* qR_press;
    const QuadratureRule* bdQr_press;


    boost::shared_ptr<OseenData> oseenData(new OseenData());
    oseenData->setup( dataFile );

    MeshData meshData;
    meshData.setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr (new mesh_Type);
    readMesh(*fullMeshPtr, meshData);

    boost::shared_ptr<mesh_Type > localMeshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart(fullMeshPtr, d->comm);
        localMeshPtr = meshPart.meshPartition();
    }

    std::string uOrder =  dataFile( "fluid/discretization/vel_order", "P1");

    if ( uOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 velocity " << std::flush;
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( uOrder.compare("P1") == 0 )
    {
        if (verbose) std::cout << "P1 velocity ";
        refFE_vel = &feTetraP1;
        qR_vel    = &quadRuleTetra4pt;  // DoE 2
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( uOrder.compare("P1Bubble") == 0 )
    {
        if (verbose) std::cout << "P1-bubble velocity " << std::flush;
        refFE_vel = &feTetraP1bubble;
        qR_vel    = &quadRuleTetra64pt;  // DoE 2
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }

    std::string pOrder =  dataFile( "fluid/discretization/press_order", "P1");
    if ( pOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 pressure " << std::flush;
        refFE_press = &feTetraP2;
        qR_press    = &quadRuleTetra15pt; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else if ( pOrder.compare("P1") == 0 )
    {
        if (verbose) std::cout << "P1 pressure";
        refFE_press = &feTetraP1;
        qR_press    = &quadRuleTetra4pt;  // DoE 2
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }

    if (verbose) std::cout << std::endl;
    if (verbose) std::cout << "Time discretization order " << oseenData->dataTime()->orderBDF() << std::endl;



    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;
    feSpacePtr_Type uFESpacePtr( new feSpace_Type( localMeshPtr,
                                                   *refFE_vel,
                                                   *qR_vel,
                                                   *bdQr_vel,
                                                   3,
                                                   d->comm) );

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

    feSpacePtr_Type pFESpacePtr( new feSpace_Type( localMeshPtr,
                                                   *refFE_press,
                                                   *qR_press,
                                                   *bdQr_press,
                                                   1,
                                                   d->comm) );

    if (verbose)
        std::cout << "ok." << std::endl;

    UInt totalVelDof   = uFESpacePtr->map().map(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpacePtr->map().map(Unique)->NumGlobalElements();

    if (verbose) std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpacePtr,
                                    *pFESpacePtr,
                                    d->comm);
    MapEpetra fullMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;


#ifdef HAVE_HDF5
    ExporterHDF5<mesh_Type > ensight( dataFile, localMeshPtr, "resistance", d->comm->MyPID());
#else
    Ensight<mesh_Type > ensight( dataFile, localMeshPtr, "resistance", d->comm->MyPID());
#endif

    vectorPtr_Type velAndPressure ( new vector_Type(*fluid.solution(), ensight.mapType() ) );

    ensight.addVariable( ExporterData<mesh_Type>::VectorField, "velocity",
                         uFESpacePtr, velAndPressure, UInt(0) );

    ensight.addVariable( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                         velAndPressure, UInt(3*uFESpacePtr->dof().numTotalDof() ) );

    // Initialization

    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();

    // bdf object to store the previous solutions

    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup(oseenData->dataTime()->orderBDF());

    // initialization with stokes solution
    if (verbose) std::cout << "Computing the stokes solution ... " << std::endl << std::endl;

    oseenData->dataTime()->setTime(t0);

    vector_Type beta( fullMap );
    vector_Type rhs ( fullMap );

    BCFunctionBase uZero(zero_scalar);

    // parabolic profile

    BCFunctionBase  uPois( velocity );


    // Resistance condition:

    vector_Type bcvector(fluid.getMap(),Repeated);

    bcvector.epetraVector().PutScalar(0.0);

    BCVector bcResistance(bcvector,uFESpacePtr->dof().numTotalDof(),1);

    bcResistance.setResistanceCoeff(d->resistance);

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

    bdf.bdfVelocity().setInitialCondition(*fluid.solution());

    fluid.setUp(dataFile);

    fluid.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);

    fluid.resetPreconditioner();

    // Temporal loop

    int iter = 1;

    chronoSet.stop();

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {
        oseenData->dataTime()->setTime(time);

        chrono.start();

        double alpha = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / oseenData->dataTime()->timeStep();

        beta = bdf.bdfVelocity().extrapolation();

        bdf.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
        rhs  = fluid.matrixMass()*bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.updateSystem( alpha, beta, rhs );
        fluid.iterate( bcH );

        bdf.bdfVelocity().shiftRight( *fluid.solution() );
        *velAndPressure = *fluid.solution();

        ensight.postProcess( time );

        MPI_Barrier(MPI_COMM_WORLD);

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}





