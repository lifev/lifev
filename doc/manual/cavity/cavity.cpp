/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/**
   \file cavity.cpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-25
 */

// LifeV definition files
#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifesolver/OseenSolver.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>
//#include <fstream> // To create an output for the flux

// Trilinos-MPI communication definitions
#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// Object type definitions
typedef LifeV::RegionMesh<LifeV::LinearTetra>       mesh_Type;
typedef LifeV::OseenSolver< mesh_Type >               fluid_Type;
typedef fluid_Type::vector_Type                       vector_Type;
typedef std::shared_ptr<vector_Type>                vectorPtr_Type;   //Pointer
typedef LifeV::FESpace< mesh_Type, LifeV::MapEpetra > feSpace_Type;
typedef std::shared_ptr<feSpace_Type>               feSpacePtr_Type;   //Pointer

// +-----------------------------------------------+
// | Data and functions for the boundary conditions|
// +-----------------------------------------------+
const int BACK   = 1;
const int FRONT  = 2;
const int LEFT   = 3;
const int RIGHT  = 4;
const int BOTTOM = 5;
const int TOP    = 6;

LifeV::Real lidBC (const LifeV::Real& /*t*/, const LifeV::Real& /*x*/, const LifeV::Real& /*y*/, const LifeV::Real& /*z*/, const LifeV::ID& i)
{
    switch (i)
    {
        case 1:
            return 1.0;
        default:
            return 0.0;
    }
}

LifeV::Real fZero ( const LifeV::Real& /* t */,
                    const LifeV::Real& /* x */,
                    const LifeV::Real& /* y */,
                    const LifeV::Real& /* z */,
                    const LifeV::ID& /* i */ )
{
    return 0.0;
}

int main (int argc, char** argv)
{

    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif

    std::shared_ptr<Epetra_Comm>   comm;
#ifdef EPETRA_MPI
    comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    int nproc;
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
#else
    comm.reset ( new Epetra_SerialComm() );
#endif

    bool verbose (false);
    if (comm->MyPID() == 0)
    {
        verbose = true;
        std::cout
                << " +-----------------------------------------------+" << std::endl
                << " |           Cavity example for LifeV            |" << std::endl
                << " +-----------------------------------------------+" << std::endl
                << std::endl
                << " +-----------------------------------------------+" << std::endl
                << " |           Author: Gwenol Grandperrin          |" << std::endl
                << " |             Date: October 25, 2010            |" << std::endl
                << " +-----------------------------------------------+" << std::endl
                << std::endl;

        std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
        std::cout << "Using MPI (" << nproc << " proc.)" << std::endl;
#else
        std::cout << "Using serial version" << std::endl;
#endif
    }

    // Chronometer
    LifeV::LifeChrono globalChrono;
    LifeV::LifeChrono initChrono;
    LifeV::LifeChrono iterChrono;
    globalChrono.start();
    initChrono.start();

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Loading the data]" << std::endl;
    }
    GetPot command_line (argc, argv);
    const std::string dataFileName = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

    // +-----------------------------------------------+
    // |               Loading the mesh                |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Loading the mesh]" << std::endl;
    }
    LifeV::MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");
    if (verbose)
    {
        std::cout << "Mesh file: " << meshData.meshDir() << meshData.meshFile() << std::endl;
    }
    std::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);
    LifeV::readMesh (*fullMeshPtr, meshData);
    // Split the mesh between processors
    std::shared_ptr< mesh_Type > localMeshPtr;
    {
        LifeV::MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, comm);
        localMeshPtr = meshPart.meshPartition();
    }

    // +-----------------------------------------------+
    // |            Creating the FE spaces             |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
    }
    std::string uOrder =  dataFile ( "fluid/space_discretization/vel_order",   "P2");
    std::string pOrder =  dataFile ( "fluid/space_discretization/press_order", "P1");
    if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                               << "FE for the pressure: " << pOrder << std::endl;

    if (verbose)
    {
        std::cout << "Building the velocity FE space... " << std::flush;
    }
    feSpacePtr_Type uFESpacePtr ( new feSpace_Type (localMeshPtr, uOrder, 3, comm) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Building the pressure FE space... " << std::flush;
    }
    feSpacePtr_Type pFESpacePtr ( new feSpace_Type (localMeshPtr, pOrder, 1, comm) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    // Total degrees of freedom (elements of matrix)
    LifeV::UInt totalVelDof   = uFESpacePtr->map().map (LifeV::Unique)->NumGlobalElements();
    LifeV::UInt totalPressDof = pFESpacePtr->map().map (LifeV::Unique)->NumGlobalElements();

    if (verbose)
    {
        std::cout << "Total Velocity Dof: " << totalVelDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total Pressure Dof: " << totalPressDof << std::endl;
    }

    // +-----------------------------------------------+
    // |             Boundary conditions               |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Boundary conditions]" << std::endl;
    }

    LifeV::BCFunctionBase uZero (fZero);
    LifeV::BCFunctionBase uLid (lidBC);

    std::vector<LifeV::ID> xComp (1);
    xComp[0] = 1;

    LifeV::BCHandler bcH;
    // A boundary condition in every face
    bcH.addBC ( "Top"   , TOP   , LifeV::Essential, LifeV::Full     , uLid , 3     );
    bcH.addBC ( "Left"  , LEFT  , LifeV::Essential, LifeV::Full     , uZero, 3     );
    bcH.addBC ( "Front" , FRONT , LifeV::Essential, LifeV::Component, uZero, xComp );
    bcH.addBC ( "Right" , RIGHT , LifeV::Essential, LifeV::Full     , uZero, 3     );
    bcH.addBC ( "Back"  , BACK  , LifeV::Essential, LifeV::Component, uZero, xComp );
    bcH.addBC ( "Bottom", BOTTOM, LifeV::Essential, LifeV::Full     , uZero, 3     );

    // Get the number of Lagrange Multiplyers (LM) and set the offsets
    std::vector<LifeV::bcName_Type> fluxVector = bcH.findAllBCWithType ( LifeV::Flux );
    LifeV::UInt numLM = static_cast<LifeV::UInt> ( fluxVector.size() );

    LifeV::UInt offset = uFESpacePtr->map().map (LifeV::Unique)->NumGlobalElements()
                         + pFESpacePtr->map().map (LifeV::Unique)->NumGlobalElements();

    for ( LifeV::UInt i = 0; i < numLM; ++i )
    {
        bcH.setOffset ( fluxVector[i], offset + i );
    }

    // +-----------------------------------------------+
    // |             Creating the problem              |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Creating the problem]" << std::endl;
    }
    std::shared_ptr<LifeV::OseenData> oseenData (new LifeV::OseenData() );
    oseenData->setup ( dataFile );

    if (verbose)
    {
        std::cout << "Time discretization order " << oseenData->dataTime()->orderBDF() << std::endl;
    }

    // The problem (matrix and rhs) is packed in an object called fluid
    LifeV::OseenSolver< mesh_Type > fluid (oseenData,
                                           *uFESpacePtr,
                                           *pFESpacePtr,
                                           comm,
                                           numLM);
    // Gets inputs from the data file
    fluid.setUp (dataFile);

    // Assemble the matrices
    fluid.buildSystem();

    // Communication map
    LifeV::MapEpetra fullMap (fluid.getMap() );

    // Synchronization
    MPI_Barrier (MPI_COMM_WORLD);

    // +-----------------------------------------------+
    // |       Initialization of the simulation        |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Initialization of the simulation]" << std::endl;
    }
    LifeV::Real dt     = oseenData->dataTime()->timeStep();
    LifeV::Real t0     = oseenData->dataTime()->initialTime();
    LifeV::Real tFinal = oseenData->dataTime()->endTime ();

    // bdf object to store the previous solutions
    LifeV::TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup (oseenData->dataTime()->orderBDF() );

    // Initialization with exact solution: either interpolation or "L2-NS"-projection
    t0 -= dt * bdf.bdfVelocity().order();

    vector_Type beta ( fullMap );
    vector_Type rhs ( fullMap );

    MPI_Barrier (MPI_COMM_WORLD);

    // We get the initial solution using a steady Stokes problem
    oseenData->dataTime()->setTime (t0);

    beta *= 0.;
    rhs  *= 0.;
    fluid.updateSystem (0.0, beta, rhs);
    fluid.iterate (bcH);
    bdf.bdfVelocity().setInitialCondition ( *fluid.solution() );

    LifeV::Real time = t0 + dt;
    for (  ; time <=  oseenData->dataTime()->initialTime() + dt / 2.; time += dt)
    {
        oseenData->dataTime()->setTime (time);

        fluid.updateSystem (0.0, beta, rhs);
        fluid.iterate (bcH);
        bdf.bdfVelocity().shiftRight ( *fluid.solution() );
    }

    // We erase the preconditioner build for Stokes
    // (A new one should be built for Navier-Stokes)
    fluid.resetPreconditioner();

    std::shared_ptr< LifeV::ExporterHDF5<mesh_Type> > exporter;

    vectorPtr_Type velAndPressure;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");

    exporter.reset ( new LifeV::ExporterHDF5<mesh_Type> ( dataFile, "cavity_example" ) );
    exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId ( localMeshPtr, comm->MyPID() );

    velAndPressure.reset ( new vector_Type (*fluid.solution(), exporter->mapType() ) );

    exporter->addVariable ( LifeV::ExporterData<mesh_Type>::VectorField, "velocity", uFESpacePtr,
                            velAndPressure, LifeV::UInt (0) );

    exporter->addVariable ( LifeV::ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                            velAndPressure, LifeV::UInt (3 * uFESpacePtr->dof().numTotalDof() ) );
    exporter->postProcess ( 0 );

    initChrono.stop();
    if (verbose)
    {
        std::cout << "Initialization time:  " << initChrono.diff() << " s." << std::endl;
    }

    // +-----------------------------------------------+
    // |             Solving the problem               |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Solving the problem]" << std::endl;
    }
    int iter = 1;

    for ( ; time <= tFinal + dt / 2.; time += dt, iter++)
    {

        oseenData->dataTime()->setTime (time);

        if (verbose)
        {
            std::cout << "[t = " << oseenData->dataTime()->time() << " s.]" << std::endl;
        }

        iterChrono.start();

        double alpha = bdf.bdfVelocity().coefficientFirstDerivative ( 0 ) / oseenData->dataTime()->timeStep();
        //Matteo
        //  beta =bdf.bdfVelocity().extrapolation();
        bdf.bdfVelocity().extrapolation ( beta ); // Extrapolation for the convective term
        bdf.bdfVelocity().updateRHSContribution (oseenData->dataTime()->timeStep() );
        rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.getDisplayer().leaderPrint ("alpha ", alpha);
        fluid.getDisplayer().leaderPrint ("\n");
        fluid.getDisplayer().leaderPrint ("norm beta ", beta.norm2() );
        fluid.getDisplayer().leaderPrint ("\n");
        fluid.getDisplayer().leaderPrint ("norm rhs  ", rhs.norm2() );
        fluid.getDisplayer().leaderPrint ("\n");

        fluid.updateSystem ( alpha, beta, rhs );
        fluid.iterate ( bcH );

        bdf.bdfVelocity().shiftRight ( *fluid.solution() );

        // Computation of the error
        vector_Type vel  (uFESpacePtr->map(), LifeV::Repeated);
        vector_Type press (pFESpacePtr->map(), LifeV::Repeated);
        vector_Type velpressure ( *fluid.solution(), LifeV::Repeated );

        velpressure = *fluid.solution();
        vel.subset (velpressure);
        press.subset (velpressure, uFESpacePtr->dim() *uFESpacePtr->fieldDim() );


        bool verbose = (comm->MyPID() == 0);


        // Exporting the solution
        *velAndPressure = *fluid.solution();
        exporter->postProcess ( time );


        MPI_Barrier (MPI_COMM_WORLD);

        iterChrono.stop();
        if (verbose)
        {
            std::cout << "Iteration time: " << iterChrono.diff() << " s." << std::endl << std::endl;
        }
    }

    globalChrono.stop();
    if (verbose)
    {
        std::cout << "Total simulation time:  " << globalChrono.diff() << " s." << std::endl;
    }

    exporter->closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}

