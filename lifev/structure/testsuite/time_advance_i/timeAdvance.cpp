//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER
/**
   \file timeAdvance.cpp

  Author(s): F. Nobile   <fabio.nobile@polimi.it>
             M. Pozzoli  <matteo1.pozzoli@mail.polimi.it>
             C. Vergara  <christian.vergara@polimi.it >

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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/structure/solver/VenantKirchhoffViscoelasticSolver.hpp>
#include <lifev/structure/solver/VenantKirchhoffViscoelasticData.hpp>

#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>

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


// ===================================================
//! Private members
// ===================================================
struct problem::Private
{
    Private() :
        rho (1)
    {}

    typedef std::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    double rho;
    std::string    data_file_name;


    std::shared_ptr<Epetra_Comm>      comm;

    fct_type getUZero()
    {
        fct_type f;
        f = std::bind (&UZero, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
        return f;
    }
};





// ===================================================
//! Constructors
// ===================================================
problem::problem ( int          argc,
                   char**                argv,
                   std::shared_ptr<Epetra_Comm>        structComm) :
    members ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    members->data_file_name = data_file_name;

    members->rho = dataFile ( "problem/physics/density", 1. );

    std::cout << "density = " << members->rho << std::endl;

    members->comm = structComm;
    int ntasks = members->comm->NumProc();

    if (!members->comm->MyPID() )
    {
        std::cout << "My PID = " << members->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}

// ===================================================
//! Methods
// ===================================================
void
problem::run()
{
    typedef RegionMesh<LinearTetra>                                   mesh_Type;
    typedef VenantKirchhoffViscoelasticSolver< mesh_Type >::vector_type vector_Type;
    typedef std::shared_ptr<vector_Type>                              vectorPtr_Type;

    typedef std::shared_ptr< TimeAdvance< vector_Type > >             TimeAdvance_type;

    typedef FESpace< mesh_Type, MapEpetra >                             FESpace_type;
    typedef  std::shared_ptr<FESpace_type>                            FESpace_ptrtype;

    bool verbose = (members->comm->MyPID() == 0);

    //
    // VenantKirchhoffViscoelasticData
    //

    GetPot dataFile ( members->data_file_name.c_str() );
    std::shared_ptr<VenantKirchhoffViscoelasticData> dataProblem (new VenantKirchhoffViscoelasticData( ) );
    dataProblem->setup (dataFile, "problem");


    MeshData             meshData;
    meshData.setup (dataFile, "problem/space_discretization");

    std::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( members->comm ) );
    readMesh (*fullMeshPtr, meshData);

    std::shared_ptr<mesh_Type > localMeshPtr;
    {
        MeshPartitioner< mesh_Type > meshPart ( fullMeshPtr, members->comm );
        localMeshPtr = meshPart.meshPartition();
    }

    //
    // The Problem Solver
    //

    if (verbose)
    {
        std::cout << "The Problem Solver" << std::flush;
    }

    // Scalar Solution Space:

    std::string Order =  dataFile ( "problem/space_discretization/order", "P1");

    if ( Order.compare ("P1") == 0 )
    {
        if (verbose)
        {
            std::cout << "  Space order : P1" << std::flush;
        }

    }
    else if ( Order.compare ("P2") == 0 )
    {
        if (verbose)
        {
            std::cout << " Space order : P2";
        }

    }
    if (verbose)
    {
        std::cout << std::endl;
    }


    // finite element space of the solution
    // finite element space of the solution
    std::string dOrder =  dataFile ( "problem/space_discretization/order", "P1");

    FESpace_ptrtype feSpace ( new FESpace_type (localMeshPtr, dOrder, 1, members->comm) );

    // instantiation of the VenantKirchhoffViscoelasticSolver class

    VenantKirchhoffViscoelasticSolver< mesh_Type > problem;

    problem.setup (dataProblem,  feSpace,
                   members->comm);

    problem.setDataFromGetPot (dataFile);
    // the boundary conditions

    BCFunctionBase uZero ( members->getUZero()  );
    BCFunctionBase uEx (uexact);

    BCHandler bcH;

    bcH.addBC ( "Top",     TOP,    Essential, Full,      uEx, 1 );
    bcH.addBC ( "Bottom",  BOTTOM, Essential, Full,      uEx, 1 );
    bcH.addBC ( "Left",    LEFT,   Essential, Full,      uEx, 1 );
    bcH.addBC ( "Right",   RIGHT,  Essential, Full,      uEx, 1 );
    bcH.addBC ( "Front",   FRONT,  Essential, Full,      uEx, 1 );
    bcH.addBC ( "Back",    BACK,   Essential, Full,      uEx, 1 );

    std::ofstream out_norm;
    if (verbose)
    {
        out_norm.open ("norm.txt");
        out_norm << "  time   "
                 << "  L2_Error    "
                 << "  H1_Error    "
                 << "  L2_RelError "
                 << "  H1_RelError \n";
        out_norm.close();
    }

    LifeChrono chrono;

    std::string TimeAdvanceMethod =  dataFile ( "problem/time_discretization/method", "Newmark");

    TimeAdvance_type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( TimeAdvanceMethod ) );

    UInt OrderDev = 1;

    //! initialization of parameters of time Advance method:

    if (TimeAdvanceMethod == "Newmark")
    {
        timeAdvance->setup ( dataProblem->dataTimeAdvance()->coefficientsNewmark() , OrderDev);
    }

    if (TimeAdvanceMethod == "BDF")
    {
        timeAdvance->setup (dataProblem->dataTimeAdvance()->orderBDF() , OrderDev);
    }

    Real dt = dataProblem->dataTime()->timeStep();
    Real T  = dataProblem->dataTime()->endTime();

    chrono.start();

    dataProblem->setup (dataFile, "problem");

    double alpha = timeAdvance->coefficientFirstDerivative ( 0 ) / ( dataProblem->dataTime()->timeStep() );

    problem.buildSystem (alpha);

    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose )
    {
        std::cout << "ok." << std::endl;
    }

    // building some vectors:
    MapEpetra uMap = problem.solution()->map();

    // computing the rhs
    vector_Type rhs ( uMap, Unique );
    vector_Type rhsV ( uMap, Unique );

    // postProcess
    std::shared_ptr< Exporter<mesh_Type > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<mesh_Type > ( dataFile, "problem" ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<mesh_Type > ( dataFile, localMeshPtr, "problem", members ->comm->MyPID() ) );
        }
        else
        {
            exporter.reset ( new ExporterEnsight<mesh_Type > ( dataFile, localMeshPtr, "problem",   members->comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId ( localMeshPtr,  members->comm->MyPID() );

    vectorPtr_Type U ( new vector_Type (*problem.solution(), exporter->mapType() ) );
    vectorPtr_Type V  ( new vector_Type (*problem.solution(),  exporter->mapType() ) );
    vectorPtr_Type Exact ( new vector_Type (*problem.solution(), exporter->mapType() ) );
    vectorPtr_Type vExact  ( new vector_Type (*problem.solution(),  exporter->mapType() ) );
    vectorPtr_Type RHS ( new vector_Type (*problem.solution(), exporter->mapType() ) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "displacement",
                            feSpace, U, UInt (0) );

    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "velocity",
                            feSpace, V, UInt (0) );

    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "uexact",
                            feSpace, Exact, UInt (0) );

    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "vexact",
                            feSpace, vExact, UInt (0) );

    exporter->postProcess ( 0 );

    //initialization of unk

    //evaluate disp and vel as interpolate the bcFunction d0 and v0
    feSpace->interpolate ( static_cast<FESpace_type::function_Type> ( d0 ), *U, 0.0 );
    feSpace->interpolate ( static_cast<FESpace_type::function_Type> ( v0 ), *V, 0.0 );

    //evaluate disp and vel as interpolate the bcFunction d0 and v0

    std::vector<vectorPtr_Type> uv0;

    if (TimeAdvanceMethod == "Newmark")
    {
        uv0.push_back (U);
        uv0.push_back (V);
    }
    if (TimeAdvanceMethod == "BDF")
    {
        for ( UInt previousPass = 0; previousPass < dataProblem->dataTimeAdvance()->orderBDF() ; previousPass++)
        {
            Real previousTimeStep = -previousPass * dt;
            feSpace->interpolate (static_cast<FESpace_type::function_Type> (uexact ), *U, previousTimeStep );
            uv0.push_back (U);
        }
    }

    //the uv0[0] should be the displacement
    //the uv0[1] should be the velocity

    //timeAdvance->setInitialCondition(uv0[0],uv0[1]);
    timeAdvance->setInitialCondition (uv0);

    timeAdvance-> setTimeStep (dataProblem->dataTime()->timeStep() );

    timeAdvance->showMe();

    feSpace->interpolate ( static_cast<FESpace_type::function_Type> ( uexact ), *Exact , 0);
    feSpace->interpolate ( static_cast<FESpace_type::function_Type> ( v0 ), *vExact , 0);

    *U = timeAdvance->solution();
    *V = timeAdvance->firstDerivative();

    for (Real time = dt; time <= T; time += dt)
    {
        dataProblem->dataTime()->setTime (time);

        if (verbose)
        {
            std::cout << std::endl;
            std::cout << " P - Now we are at time " << dataProblem->dataTime()->time() << " s." << std::endl;
        }

        //evaluate rhs
        rhs *= 0;
        timeAdvance->updateRHSContribution (dt);
        rhsV = timeAdvance->rhsContributionFirstDerivative();
        feSpace->l2ScalarProduct (source_in, rhs, time);
        rhs += problem.matrMass() * rhsV;

        //updateRHS
        problem.updateRHS (rhs);

        //solver system
        problem.iterate ( bcH );   // Computes the matrices and solves the system

        //update unknowns of timeAdvance
        timeAdvance->shiftRight (*problem.solution() );

        //evaluate uexact solution
        feSpace->interpolate ( static_cast<FESpace_type::function_Type> ( uexact ), *Exact , time);
        feSpace->interpolate ( static_cast<FESpace_type::function_Type> ( v0 ), *vExact , time);
        *U =  timeAdvance->solution();
        *V  = timeAdvance->firstDerivative();

        //postProcess
        exporter->postProcess (time);

        // Error L2 and H1 Norms
        AnalyticalSol uExact;
        vector_Type u (*problem.solution(), Repeated);

        Real H1_Error,  H1_RelError,  L2_Error, L2_RelError;

        L2_Error = feSpace->l2Error (uexact, u, time , &L2_RelError);
        H1_Error = feSpace->h1Error (uExact, u, time , &H1_RelError);

        //save the norm
        if (verbose)
        {
            out_norm.open ("norm.txt", std::ios::app);
            out_norm << time             << "   "
                     << L2_Error       << "   "
                     << H1_Error      << "   "
                     << L2_RelError << "   "
                     << H1_RelError << "\n";
            out_norm.close();
        }

        MPI_Barrier (MPI_COMM_WORLD);
    }
    chrono.stop();
    if (verbose)
    {
        std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }
}
