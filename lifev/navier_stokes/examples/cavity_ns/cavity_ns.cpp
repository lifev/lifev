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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/navier_stokes/solver/OseenSolver.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>



//#include <fstream> // To create an output for the flux

// Trilinos-MPI communication definitions
#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


using namespace LifeV;


// Object type definitions
typedef RegionMesh<LinearTetra>         mesh_Type;
typedef OseenSolver< mesh_Type >        fluid_Type;
typedef fluid_Type::vector_Type         vector_Type;
typedef std::shared_ptr<vector_Type>  vectorPtr_Type;   //Pointer
typedef FESpace< mesh_Type, MapEpetra > feSpace_Type;
typedef std::shared_ptr<feSpace_Type> feSpacePtr_Type;   //Pointer


// +-----------------------------------------------+
// | Data and functions for the boundary conditions|
// +-----------------------------------------------+
const int BACK   = 1;
const int FRONT  = 2;
const int LEFT   = 3;
const int RIGHT  = 4;
const int BOTTOM = 5;
const int TOP    = 6;

Real lidBC (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
        case 1:
            return 1.0;
        default:
            return 0.0;
    }
}

Real fZero ( const Real& /* t */,
             const Real& /* x */,
             const Real& /* y */,
             const Real& /* z */,
             const ID& /* i */ )
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
    LifeChrono globalChrono;
    LifeChrono initChrono;
    LifeChrono iterChrono;
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
    MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");
    if (verbose)
    {
        std::cout << "Mesh file: " << meshData.meshDir() << meshData.meshFile() << std::endl;
    }
    std::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);
    readMesh (*fullMeshPtr, meshData);
    // Split the mesh between processors
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, comm);

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
    feSpacePtr_Type uFESpacePtr ( new feSpace_Type (meshPart.meshPartition(), uOrder, 3, comm) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Building the pressure FE space... " << std::flush;
    }
    feSpacePtr_Type pFESpacePtr ( new feSpace_Type (meshPart.meshPartition(), pOrder, 1, comm) );
    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    // Total degrees of freedom (elements of matrix)
    UInt totalVelDof   = uFESpacePtr->map().map (Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpacePtr->map().map (Unique)->NumGlobalElements();

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

    BCFunctionBase uZero (fZero);
    BCFunctionBase uLid (lidBC);

    std::vector<ID> xComp (1);
    xComp[0] = 1;

    BCHandler bcH;
    // A boundary condition in every face
    bcH.addBC ( "Top"   , TOP   , Essential, Full     , uLid , 3     );
    bcH.addBC ( "Left"  , LEFT  , Essential, Full     , uZero, 3     );
    bcH.addBC ( "Front" , FRONT , Essential, Component, uZero, xComp );
    bcH.addBC ( "Right" , RIGHT , Essential, Full     , uZero, 3     );
    bcH.addBC ( "Back"  , BACK  , Essential, Component, uZero, xComp );
    bcH.addBC ( "Bottom", BOTTOM, Essential, Full     , uZero, 3     );

    // Get the number of Lagrange Multiplyers (LM) and set the offsets
    std::vector<bcName_Type> fluxVector = bcH.findAllBCWithType ( Flux );
    UInt numLM = static_cast<UInt> ( fluxVector.size() );

    UInt offset = uFESpacePtr->map().map (Unique)->NumGlobalElements()
                  + pFESpacePtr->map().map (Unique)->NumGlobalElements();

    for ( UInt i = 0; i < numLM; ++i )
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
    std::shared_ptr< OseenData > oseenData (new OseenData);
    oseenData->setup ( dataFile );


    if (verbose)
    {
        std::cout << "Time discretization order " << oseenData->dataTimeAdvance()->orderBDF() << std::endl;
    }

    // The problem (matrix and rhs) is packed in an object called fluid
    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpacePtr,
                                    *pFESpacePtr,
                                    comm,
                                    numLM);

    // Gets inputs from the data file
    fluid.setUp (dataFile);

    // Assemble the matrices
    fluid.buildSystem();

    // Communication map
    MapEpetra fullMap (fluid.getMap() );

    // +-----------------------------------------------+
    // |       Initialization of the simulation        |
    // +-----------------------------------------------+
    if (verbose)
    {
        std::cout << std::endl << "[Initialization of the simulation]" << std::endl;
    }
    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime ();

    // bdf object to store the previous solutions
    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup (oseenData->dataTimeAdvance()->orderBDF() );

    t0 -= dt * bdf.bdfVelocity().order();

    vector_Type beta ( fullMap );
    vector_Type rhs ( fullMap );
    double alpha = 0.;
    //MPI_Barrier(MPI_COMM_WORLD);

    // We get the initial solution using a steady Stokes problem
    oseenData->dataTime()->setTime (t0);

    beta *= 0.;
    rhs  *= 0.;
    fluid.updateSystem (alpha, beta, rhs);
    fluid.iterate (bcH);
    bdf.bdfVelocity().setInitialCondition ( *fluid.solution() );

    Real time = t0 + dt;
    for (  ; time <=  oseenData->dataTime()->initialTime() + dt / 2.; time += dt)
    {
        oseenData->dataTime()->setTime (time);

        fluid.updateSystem (alpha, beta, rhs);
        fluid.iterate (bcH);
        bdf.bdfVelocity().shiftRight ( *fluid.solution() );
    }

    // We erase the preconditioner build for Stokes
    // (A new one should be built for Navier-Stokes)
    fluid.resetPreconditioner();

    std::shared_ptr< ExporterHDF5<mesh_Type> > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");

    exporter.reset ( new ExporterHDF5<mesh_Type> ( dataFile, "cavity_example" ) );
    exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId ( meshPart.meshPartition(), comm->MyPID() );

    vectorPtr_Type velAndPressure;
    velAndPressure.reset ( new vector_Type (*fluid.solution(), exporter->mapType() ) );

    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpacePtr,
                            velAndPressure, UInt (0) );

    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                            velAndPressure, UInt (3 * uFESpacePtr->dof().numTotalDof() ) );
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

        bdf.bdfVelocity().extrapolation (beta); // Extrapolation for the convective term
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
        vector_Type vel  (uFESpacePtr->map(), Repeated);
        vector_Type press (pFESpacePtr->map(), Repeated);
        vector_Type velpressure ( *fluid.solution(), Repeated );

        velpressure = *fluid.solution();
        vel.subset (velpressure);
        press.subset (velpressure, uFESpacePtr->dim() *uFESpacePtr->fieldDim() );


        bool verbose = (comm->MyPID() == 0);


        // Exporting the solution
        *velAndPressure = *fluid.solution();
        exporter->postProcess ( time );

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

