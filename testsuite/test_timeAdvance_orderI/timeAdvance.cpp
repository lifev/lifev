/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): F. Nobile   <fabio.nobile@polimi.it>
             M. Pozzoli  <matteo1.pozzoli@mail.polimi.it>
             C. Vergara  <christian.vergara@polimi.it >
       Date: 2009-03-24

  Copyright (C) 2009 EPFL

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
   \file timeAdvance.cpp

   \date 2010-02-10
 */

/* ========================================================

Simple Problem test with Dirichlet Boundary condition

Solve the problem

              \frac{\partial u}{\partial t} - \Delta u = f

               u = u0 on the boundary

linear_function.hpp:

  uexact = exp(-sin(Pi/2*t))*(x+y+z);

  f = (3 \pi^2 + 1 )  exp(-t)  sin( \pi x) sin(\pi y) sin ( \pi z) on a cube

nonlinear_function.hpp:

   uexact = exp(-sin(Pi/2*t))*cos(x *Pi)*cos(y*Pi)*cos(z*Pi);

   f = Pi2/4*( sin(Pi/2*t)+cos(Pi/2*t)*cos(Pi/2*t) )*exp(-sin(Pi/2*t))*cos(x *Pi)*cos(y*Pi)*cos(z*Pi);
*/


// ===================================================
//! Includes
// ===================================================
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>

#include <life/lifesolver/SecondOrderSolver.hpp>
#include <life/lifesolver/dataSecondOrder.hpp>

#include <life/lifefem/dataTime.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifefem/timeAdvance_template.hpp>
#include <life/lifefem/newmark_template.hpp>
#include <life/lifefem/bdf_template.hpp>

#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifefilters/noexport.hpp>

#include <iostream>

#include "timeAdvance.hpp"

//choose of function.hpp
#include "linear_function.hpp"
//#include "nonlinear_function.hpp"

// ===================================================
//! Namespaces & Define
// ===================================================
using namespace LifeV;

const int TOP    = 6;
const int BOTTOM = 5;
const int LEFT   = 3;
const int RIGHT  = 4;
const int FRONT  = 2;
const int BACK   = 1;

typedef RegionMesh3D<LinearTetra> RegionMesh;


// ===================================================
//! Private members
// ===================================================
struct problem::Private
{
    Private():
            rho (1)
    {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;
    double rho;
    std::string    data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

    fct_type getUZero()
    {
        fct_type f;
        f = boost::bind(&UZero, _1, _2, _3, _4, _5);
        return f;
    }
};





// ===================================================
//! Constructors
// ===================================================
problem::problem( int          argc,
                  char**                argv,
                  boost::shared_ptr<Epetra_Comm>        structComm):
        members( new Private() )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );
    members->data_file_name = data_file_name;

    members->rho = dataFile( "problem/physics/density", 1. );

    std::cout << "density = " << members->rho << std::endl;

    members->comm = structComm;
    int ntasks = members->comm->NumProc();

    if (!members->comm->MyPID()) std::cout << "My PID = " << members->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
}

// ===================================================
//! Methods
// ===================================================
void
problem::run()
{
    typedef SecondOrderSolver< RegionMesh3D<LinearTetra> >::vector_type    vector_type;
    typedef boost::shared_ptr<vector_type>                                         vector_ptrtype;

    typedef boost::shared_ptr< TimeAdvance< vector_type > >                 TimeAdvance_type;

    bool verbose = (members->comm->MyPID() == 0);

    //
    // dataSecondOrder
    //

    GetPot dataFile( members->data_file_name.c_str() );

    DataSecondOrder dataProblem;
    dataProblem.setup(dataFile, "problem");

    DataMesh             dataMesh;
    dataMesh.setup(dataFile, "problem/space_discretization");

    boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);
    readMesh(*fullMeshPtr, dataMesh);

    partitionMesh< RegionMesh3D<LinearTetra> > meshPart( fullMeshPtr, members->comm );

    //
    // The Problem Solver
    //

    if (verbose) std::cout << "The Problem Solver" << std::flush;

    const RefFE*    refFE_beta  (0);
    const QuadRule* qR_beta     (0);
    const QuadRule* bdQr_beta   (0);

    refFE_beta = &feTetraP1;

    // Scalar Solution Space:

    std::string Order =  dataFile( "problem/space_discretization/order", "P1");
    const RefFE*    refFE(0);
    const QuadRule* qR(0);
    const QuadRule* bdQr(0);

    if ( Order.compare("P1") == 0 )
    {
        if (verbose) std::cout << "  Space order : P1" << std::flush;

        refFE = &feTetraP1;
        qR    = &quadRuleTetra15pt; // DoE 5
        bdQr  = &quadRuleTria4pt;   // DoE 2

    }
    else if ( Order.compare("P2") == 0 )
    {
        if (verbose) std::cout << " Space order : P2";

        refFE = &feTetraP2;
        qR    = &quadRuleTetra64pt; // DoE 6
        bdQr = &quadRuleTria4pt;   // DoE 2

    }
    if (verbose) std::cout << std::endl;


    // finite element space of the solution

    FESpace< RegionMesh, EpetraMap > FESpace(meshPart,
                                             *refFE,
                                             *qR,
                                             *bdQr,
                                             1,
                                             members->comm);

    // instantiation of the SecondOrderSolver class

    SecondOrderSolver< RegionMesh3D<LinearTetra> > problem                 (dataProblem,
                                                                            FESpace,
                                                                            members->comm);

// the boundary conditions

    BCFunctionBase uZero ( members->getUZero()  );
    BCFunctionBase uEx(uexact);

    BCHandler bcH( 6, BCHandler::HINT_BC_ONLY_ESSENTIAL );

    bcH.addBC( "Top",     TOP,    Essential, Full,      uEx, 1 );
    bcH.addBC( "Bottom",  BOTTOM, Essential, Full,      uEx, 1 );
    bcH.addBC( "Left",    LEFT,   Essential, Full,      uEx, 1 );
    bcH.addBC( "Right",   RIGHT,  Essential, Full,      uEx, 1 );
    bcH.addBC( "Front",   FRONT,  Essential, Full,      uEx, 1 );
    bcH.addBC( "Back",    BACK,   Essential, Full,      uEx, 1 );

    std::ofstream out_norm;
    if (verbose)
    {
        out_norm.open("norm.txt");
        out_norm << "  time   "
        <<"  L2_Error    "
        <<"  H1_Error    "
        <<"  L2_RelError "
        <<"  H1_RelError \n";
        out_norm.close();
    }

    Chrono chrono;

    std::string TimeAdvanceMethod =  dataFile( "problem/time_discretization/method", "Newmark");

    TimeAdvance_type  timeAdvance( TimeAdvanceFactory::instance().createObject( TimeAdvanceMethod ) );

    UInt OrderDev = 1;

    //! initialization of parameters of time Advance method:

    if (TimeAdvanceMethod =="Newmark")
        timeAdvance->setup( dataProblem.dataTime()->getNewmark_parameters() , OrderDev);

    if (TimeAdvanceMethod =="BDF")
        timeAdvance->setup(dataProblem.dataTime()->getBDF_order() , OrderDev);

    Real dt = dataProblem.dataTime()->getTimeStep();
    Real T  = dataProblem.dataTime()->getEndTime();

    chrono.start();
    problem.setUp(dataFile);

    problem.buildSystem();

    MPI_Barrier(MPI_COMM_WORLD);

    if (verbose ) std::cout << "ok." << std::endl;

    // building some vectors:
    EpetraMap uMap = problem.u().getMap();

    // computing the rhs
    vector_type rhs ( uMap, Unique );
    vector_type rhsV( uMap, Unique );

    // postProcess
    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
        exporter.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, "problem" ) );
    else
#endif
    {
        if (exporterType.compare("none") == 0)
            exporter.reset( new NoExport<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "problem", members ->comm->MyPID()) );
        else
            exporter.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "problem",   members->comm->MyPID()) );
    }

    exporter->setDirectory( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( meshPart.mesh(),  members->comm->MyPID() );

    vector_ptrtype U ( new vector_type(problem.u(), exporter->mapType() ) );
    vector_ptrtype V  ( new vector_type(problem.u(),  exporter->mapType() ) );
    vector_ptrtype Exact ( new vector_type(problem.u(), exporter->mapType() ) );
    vector_ptrtype vExact  ( new vector_type(problem.u(),  exporter->mapType() ) );
    vector_ptrtype RHS ( new vector_type(problem.u(), exporter->mapType() ) );
    exporter->addVariable( ExporterData::Scalar, "displacement", U,
                           UInt(0), FESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::Scalar, "velocity", V,
                           UInt(0), FESpace.dof().numTotalDof() );
    exporter->addVariable( ExporterData::Scalar, "uexact", Exact,
                           UInt(0), FESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::Scalar, "vexact", vExact,
                           UInt(0), FESpace.dof().numTotalDof() );

    exporter->postProcess( 0 );

    //initialization of unk

    //evaluate disp and vel as interpolate the bcFunction d0 and v0
    FESpace.interpolate(d0, *U, 0.0);
    FESpace.interpolate(v0, *V , 0.0);

//evaluate disp and vel as interpolate the bcFunction d0 and v0

    std::vector<vector_type> uv0;

    if (TimeAdvanceMethod =="Newmark")
    {
        uv0.push_back(*U);
        uv0.push_back(*V);
    }
    if (TimeAdvanceMethod =="BDF")
    {
        for ( int previousPass=0; previousPass < dataProblem.getBDF_order() ; previousPass++)
        {
            Real previousTimeStep = -previousPass*dt;
            FESpace.interpolate(uexact, *U, previousTimeStep );
            uv0.push_back(*U);
        }
    }

    timeAdvance->initialize_unk(uv0);

    timeAdvance-> setDeltaT(dataProblem.dataTime()->getTimeStep());

    timeAdvance->showMe();

    FESpace.interpolate(uexact, *Exact , 0);
    FESpace.interpolate(v0, *vExact , 0);

    *U = timeAdvance->unk(0);
    *V = timeAdvance->unk(1);

    for (Real time = dt; time <= T; time += dt)
    {
        dataProblem.setTime(time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << " P - Now we are at time " << dataProblem.getTime() << " s." << std::endl;
        }

        double alpha = timeAdvance->coeff_der( 0 ) / ( dataProblem.dataTime()->getTimeStep());

        //evaluate rhs
        rhs *=0;

        rhsV = timeAdvance->time_der(dt);
        FESpace.L2ScalarProduct(source_in, rhs, time);
        rhs += problem.matrMass() *rhsV;

        //update system
        problem.updateSystem(alpha, rhs );

        //solver system
        problem.iterate( bcH );    // Computes the matrices and solves the system

        //update unknowns of timeAdvance
        timeAdvance->shift_right(problem.u());

        //evaluate uexact solution
        FESpace.interpolate(uexact, *Exact , time);
        FESpace.interpolate(v0, *vExact , time);
        *U =  timeAdvance->unk(0);
        *V  =timeAdvance->vnk();

        //postProcess
        exporter->postProcess(time);

        // Error L2 and H1 Norms
        AnalyticalSol uExact;
        vector_type u (problem.u(), Repeated);

        Real H1_Error,H1_Error1, H1_RelError, L2_Error1, L2_Error, L2_RelError;

        L2_Error = FESpace.L2Error(uexact, u, time ,&L2_RelError);
        H1_Error = FESpace.H1Error(uExact, u, time ,&H1_RelError);


        //save the norm
        if (verbose)
        {
            out_norm.open("norm.txt", std::ios::app);
            out_norm << time             << "   "
            << L2_Error       << "   "
            << H1_Error      << "   "
            << L2_RelError << "   "
            << H1_RelError << "\n";
            out_norm.close();
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    chrono.stop();
    if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
}
