/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Samuel Quinodoz <samuel.quinodoz@epfl.ch>
       Date: 2011-01-27

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
   \file main.cpp
   \author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   \date 2011-01-127


   Optimizations: keep the matrices, *= 0.0 only

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

#include <lifev/core/LifeV.hpp>


#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorBlockMonolithicEpetra.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/fem/GradientRecovery.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
//#include <lifev/core/mesh/RegionMesh3D.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/level_set/fem/LevelSetQRAdapter.hpp>

#include <iostream>

//#define BDF2_TIME 1

//#define FULL_CORRECTION 1
#define NORMAL_CORRECTION 1

using namespace LifeV;

// Physics
Real viscosityPlus (1e-3);
Real viscosityMinus (2e-5);
Real densityPlus (1000.0);
Real densityMinus (1.0);

Real liquidHeight (0.1);
Real cylinderRadius (0.144);
Real cylinderHeight (0.4);
Real alphaBL (3.0);

Real thetaInit (0.0);
Real omegaInit (0.0);
Real tRise (6.0);
Real omegaMax (2 * M_PI);
Real shakeRadius (0.1);

Real slipLength (0.05);
Real noSlipCoef (1e6);
Real viscosityRatio (50.0);
Real slipFriction (0.0);

// NS
Real NSSupg (1.0e-3); // 1e-3
Real NSPspg (1.0e-3); // 1e-3
Real NSdivdiv (5.0e2); // 5e2

bool NormalizeCorrection (true);

// LS advection
Real LSSupg (5e-2);

// HJ
UInt HJeach (1);
Real HJdt (2e-2);
Real HJtime (0.0);
Real HJfinalTime (HJdt);
Real HJsupg (5e-1);
Real HJpenal (1e5);

// Volume correction
Real volumeRelaxation (0.5);
Real epsilonArea (0.05);

// Export
UInt exportEach (1);



/* Level set initial position */
Real initLSFct ( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& z , const ID& /*i*/)
{
    Real distToSurface ( liquidHeight - z );

    return distToSurface;
};

/* Velocity initial */
Real initVelocity ( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& i)
{
    if (i == 2)
    {
        return 0.00;
    }
    return 0.0;
};


/* Velocity initial */
Real zeroFct ( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/)
{
    return 0.0;
};

/* Velocity initial */
Real oneFct ( const Real& /*t*/, const Real& /*x*/ , const Real& /*y*/, const Real& /*z*/ , const ID& /*i*/)
{
    return 1.0;
};


/* Velocity initial */
Real volumeForce ( const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    Real alpha;
    Real omega;
    Real theta;

    if (t < tRise)
    {
        theta = thetaInit + omegaInit * t + (omegaMax - omegaInit) * (t - tRise * std::cos (-M_PI / 2.0 + M_PI * t / tRise) / M_PI) / 2.0;
        omega = omegaInit + (omegaMax - omegaInit) * (std::sin (-M_PI / 2.0 + M_PI * t / tRise) + 1.0) / 2.0;
        alpha = (omegaMax - omegaInit) * M_PI * std::cos (-M_PI / 2.0 + M_PI * t / tRise) / (2 * tRise);
    }
    else
    {
        alpha = 0.0;
        omega = omegaMax;
        theta = thetaInit + (omegaMax - omegaInit) * (tRise / 2.0) + omegaMax * (t - tRise);
    }

    if (i == 0)
    {
        return shakeRadius * omega * omega * cos (theta) - shakeRadius * alpha * sin (theta);
    }
    if (i == 1)
    {
        return -shakeRadius * omega * omega * sin (theta) - shakeRadius * alpha * cos (theta);
    }
    if (i == 2)
    {
        return -9.806;
    }

    return 0.0;
};

void displayMotionParameters ( const Real& t)
{
    Real alpha;
    Real omega;
    Real theta;

    if (t < tRise)
    {
        theta = thetaInit + omegaInit * t + (omegaMax - omegaInit) * (t - tRise * std::cos (-M_PI / 2.0 + M_PI * t / tRise) / M_PI) / 2.0;
        omega = omegaInit + (omegaMax - omegaInit) * (std::sin (-M_PI / 2.0 + M_PI * t / tRise) + 1.0) / 2.0;
        alpha = (omegaMax - omegaInit) * M_PI * std::cos (-M_PI / 2.0 + M_PI * t / tRise) / (2 * tRise);
    }
    else
    {
        alpha = 0.0;
        omega = omegaMax;
        theta = thetaInit + (omegaMax - omegaInit) * (tRise / 2.0) + omegaMax * (t - tRise);
    }

    std::cout << " ###  Motion  ### " << std::endl;
    std::cout << " theta : " << theta << " (rad) " << std::endl;
    std::cout << " omega : " << 60.0 * omega / (2.0 * M_PI) << " (rpm) " << std::endl;
    std::cout << " alpha : " << alpha << " (rad/s2) " << std::endl;
}

Real volumeForce0 ( const Real& t, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return volumeForce (t, x, y, z, 0);
};
Real volumeForce1 ( const Real& t, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return volumeForce (t, x, y, z, 1);
};
Real volumeForce2 ( const Real& t, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return volumeForce (t, x, y, z, 2);
};

Real normalDirection ( const Real& /*t*/, const Real& x , const Real& y, const Real& /*z*/ , const ID& i)
{
    Real nnorm (x * x + y * y);
    if (i == 0)
    {
        return x / nnorm;
    }
    else if (i == 1)
    {
        return y / nnorm;
    }

    return 0.0;
};

Real slipFct ( Real ls )
{
    //Real frac(std::abs(ls)/(slipLength));

    if ( ls < -slipLength)
    {
        return noSlipCoef / viscosityRatio;
    }
    if ( ls < 0 )
    {
        return slipFriction / viscosityRatio;
    }
    if ( ls < slipLength )
    {
        return slipFriction;
    }
    return noSlipCoef;
};


void meshMap ( Real& x, Real& y, Real& z )
{

    if (alphaBL == 0.0)
    {
        x = x * cylinderRadius / 0.144;
        y = y * cylinderRadius / 0.144;

    }
    else
    {
        // Map to the unit disk
        x /= 0.144;
        y /= 0.144;

        // New radius there
        Real r (std::sqrt (x * x + y * y) );

        Real factor = (1 - std::exp (-alphaBL * r) ) / (1 - std::exp (-alphaBL) );
        x *= factor / r;
        y *= factor / r;

        // Map to the right radius
        x *= cylinderRadius;
        y *= cylinderRadius;
    }

    z = z * cylinderHeight / 0.4;
};


/* Heaviside function class */
class HeavisideFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& value)
    {
        if (value >= 0)
        {
            return 1.0;
        }
        return 0.0;
    }

    HeavisideFct() {}
    HeavisideFct (const HeavisideFct&) {}
    ~HeavisideFct() {}
};

/* Viscosity */
class ViscosityFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& value)
    {
        if (value >= 0)
        {
            return viscosityPlus;
        }
        return viscosityMinus;
    }

    ViscosityFct() {}
    ViscosityFct (const ViscosityFct&) {}
    ~ViscosityFct() {}
};

#define VISCOSITY eval(viscosity, value(ETlsFESpace,LSSolutionOld))

/* Density */
class DensityFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& value)
    {
        if (value >= 0)
        {
            return densityPlus;
        }
        return densityMinus;
    }

    DensityFct() {}
    DensityFct (const DensityFct&) {}
    ~DensityFct() {}
};

#define DENSITY eval(density, value(ETlsFESpace,LSSolutionOld))

/* NormalizeFct */
class NormalizeFct
{
public:
    typedef VectorSmall<3> return_Type;

    return_Type operator() (const VectorSmall<3>& value)
    {
        Real norm (sqrt ( value[0]*value[0] + value[1]*value[1] + value[2]*value[2]) );

        if (norm > 0)
        {
            return value * (1.0 / norm);
        }
        return value;
    }

    NormalizeFct() {}
    NormalizeFct (const NormalizeFct&) {}
    ~NormalizeFct() {}
};


/* NormalizeFct */
class NormFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<3> value)
    {
        Real norm (sqrt ( value[0]*value[0] + value[1]*value[1] + value[2]*value[2]) );

        return norm;
    }

    NormFct() {}
    NormFct (const NormFct&) {}
    ~NormFct() {}
};


/* SignFct */
class SignFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& value)
    {
        if (value > 0)
        {
            return 1.0;
        }
        else if (value < 0)
        {
            return -1.0;
        }
        return 0.0;
    }

    SignFct() {}
    SignFct (const SignFct&) {}
    ~SignFct() {}
};

/* SmoothSignFct */
class SmoothDeltaFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& value)
    {
        Real frac (value / epsilonArea);

        if (frac < -1)
        {
            return 0.0;
        }
        else if (frac > 1)
        {
            return 0.0;
        }
        else
        {
            return 0.5 * (1 + std::cos (M_PI * frac) ) / epsilonArea;
        }
    }

    SmoothDeltaFct() {}
    SmoothDeltaFct (const SmoothDeltaFct&) {}
    ~SmoothDeltaFct() {}
};




/* Register the preconditioner */
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}


/* Some typedef */
typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetraStructured<Real> matrix_block_type;
typedef VectorBlockMonolithicEpetra vector_block_type;
typedef MatrixEpetra<Real> matrix_type;
typedef VectorEpetra vector_type;


#define BOTTOM_LINE 1
#define TOP_LINE 2
#define BOTTOM 3
#define WALL 4
#define TOP 5

int main ( int argc, char** argv )
{


#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    // a flag to see who's the leader for output purposes
    bool verbose = (Comm->MyPID() == 0);

    // Open and read the data file
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    //GetPot dataFile( data_file_name );
    GetPot dataFile ( "data" );

    // Reading the data

    viscosityPlus = dataFile ("physics/viscosityPlus", 1e-3);
    viscosityMinus = dataFile ("physics/viscosityMinus", 2e-5);
    densityPlus = dataFile ("physics/densityPlus", 1000.0);
    densityMinus = dataFile ("physics/densityMinus", 1.0);

    NSSupg = dataFile ("oseen/stab/supg", 1e-3);
    NSPspg = dataFile ("oseen/stab/pspg", 1e-3);
    NSdivdiv = dataFile ("oseen/stab/divdiv", 5e2);

    NormalizeCorrection = dataFile ("oseen/pressure/normalize_correction", true);

    LSSupg = dataFile ("level_set/advection/supg", 5e-2);

    HJdt =  dataFile ("level_set/hj/dt", 2e-2);
    HJfinalTime = dataFile ("level_set/hj/final", HJdt);
    HJsupg = dataFile ("level_set/hj/supg", 5e-1);
    HJeach = dataFile ("level_set/hj/each", 1);
    HJpenal = dataFile ("level_set/hj/penalisation", 1e5);

    liquidHeight = dataFile ("physics/liquid_height", 0.1);
    cylinderRadius = dataFile ("physics/cylinder_radius", 0.144);
    cylinderHeight = dataFile ("physics/cylinder_height", 0.4);

    thetaInit = dataFile ("physics/theta_init", 0.0);
    omegaInit = 2 * M_PI * dataFile ("physics/omega_init", 0.0) / 60.0;
    tRise = dataFile ("physics/time_rise", 6.0);
    omegaMax = 2 * M_PI * dataFile ("physics/omega_final", 60.0) / 60.0;

    shakeRadius = dataFile ("physics/shake_radius", 0.1);

    slipLength = dataFile ("oseen/bc/slip_length", 0.05);
    noSlipCoef = dataFile ("oseen/bc/no_slip_coef", 1e6);
    viscosityRatio = dataFile ("oseen/bc/ratio", 50.0);
    slipFriction = dataFile ("oseen/bc/slip_friction", 0.0);

    exportEach = dataFile ("exporter/each", 1);

    alphaBL = dataFile ("mesh/alphaBL", 0.0);

    if (verbose)
    {
        std::cout << "alpha : " << alphaBL << std::endl;
    }


    Real exactVolume (M_PI * cylinderRadius * cylinderRadius * liquidHeight);
    if (verbose)
    {
        std::cout << "Exact : " << exactVolume << std::endl;
    }

    if (verbose)
    {
        std::cout << " Reading and partitioning the mesh " << std::endl;
    }


    // Load the mesh
    MeshData dataMesh;
    dataMesh.setup (dataFile, "mesh");

    std::shared_ptr < mesh_Type > fullMeshPtr (new mesh_Type);
    std::shared_ptr < mesh_Type > localMeshPtr (new mesh_Type);

    if (verbose)
    {
        std::cout << "Hello " << exactVolume << std::endl;
    }
    readMesh (*fullMeshPtr, dataMesh);

    // Scale the mesh
    //std::vector<Real> ScaleFactor(3,cylinderRadius/0.144);
    //ScaleFactor[2]=cylinderHeight/0.4;
    //std::vector<Real> DoNothing(3,0.0);
    //fullMeshPtr->transformMesh(ScaleFactor,DoNothing,DoNothing);

    if (verbose)
    {
        std::cout << " BL factor : " << alphaBL << std::endl;
    }
    //std::shared_ptr< MeshTransformer > meshTransformerObject ( new MeshTransfomer(fullMeshPtr) );
    //fullMeshPtr->transformMesh(meshMap);
    fullMeshPtr->meshTransformer().transformMesh (meshMap);

    // Partition the mesh
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
    localMeshPtr = meshPart.meshPartition();

    // Free the global mesh
    fullMeshPtr.reset();


    if (verbose)
    {
        std::cout << " Building FESpaces  " << std::endl;
    }


    //std::string uOrder("P1Bubble");
    std::string uOrder ("P1");
    std::string pOrder ("P1");
    std::string lsOrder ("P1");

    std::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace ( new FESpace< mesh_Type, MapEpetra > (localMeshPtr, uOrder, 3, Comm) );
    std::shared_ptr<FESpace< mesh_Type, MapEpetra > > pFESpace ( new FESpace< mesh_Type, MapEpetra > (localMeshPtr, pOrder, 1, Comm) );
    std::shared_ptr<FESpace< mesh_Type, MapEpetra > > lsFESpace ( new FESpace< mesh_Type, MapEpetra > (localMeshPtr, lsOrder, 1, Comm) );

    if (verbose)
    {
        std::cout << std::endl << " ### Dof Summary ###: " <<  std::endl;
    }
    if (verbose)
    {
        std::cout << " Velocity  : " << uFESpace->map().map (Unique)->NumGlobalElements() << std::endl;
    }
    if (verbose)
    {
        std::cout << " Pressure  : " << pFESpace->map().map (Unique)->NumGlobalElements() << std::endl;
    }
    if (verbose)
    {
        std::cout << " Level set : " << lsFESpace->map().map (Unique)->NumGlobalElements() << std::endl << std::endl;
    }


    if (verbose)
    {
        std::cout << " Building EA FESpaces  " << std::endl;
    }


    std::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuFESpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPart, & (uFESpace->refFE() ), Comm) );
    std::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETpFESpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, & (pFESpace->refFE() ), Comm) );
    std::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETlsFESpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, & (lsFESpace->refFE() ), Comm) );


    if (verbose)
    {
        std::cout << " Initial conditions " << std::endl;
    }

    vector_type LSSolution (ETlsFESpace->map(), Unique);

    lsFESpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (initLSFct), LSSolution, 0.0);

    vector_type LSSolutionOld (LSSolution, Repeated);
    vector_type HJSolutionOld (LSSolution, Repeated);
    vector_type HJProjSolution (LSSolution, Repeated);

    vector_type velocitySolution (ETuFESpace->map(), Unique);

    uFESpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (initVelocity), velocitySolution, 0.0);

    vector_type velocitySolutionOld (velocitySolution, Repeated);
#ifdef BDF2_TIME
    vector_type velocitySolutionOldOld (velocitySolution, Repeated);
#endif

    vector_block_type NSSolution (ETuFESpace->map() | ETpFESpace->map(), Unique);
    NSSolution *= 0.0;
    //NSSolution += velocitySolution; THIS DOES NOT WORK!

    vector_block_type NSSolutionOld (NSSolution , Repeated);
    NSSolutionOld *= 0.0; // JUST TO BE SURE


    if (verbose)
    {
        std::cout << " Building the solvers " << std::endl;
    }


    SolverAztecOO LSAdvectionSolver;
    LSAdvectionSolver.setCommunicator (Comm);
    LSAdvectionSolver.setDataFromGetPot (dataFile, "solver");
    LSAdvectionSolver.setupPreconditioner (dataFile, "prec");

    SolverAztecOO NSSolver;
    NSSolver.setCommunicator (Comm);
    NSSolver.setDataFromGetPot (dataFile, "solver");
    NSSolver.setupPreconditioner (dataFile, "prec");

    SolverAztecOO HJSolver;
    HJSolver.setCommunicator (Comm);
    HJSolver.setDataFromGetPot (dataFile, "solver");
    HJSolver.setupPreconditioner (dataFile, "prec");

    if (verbose)
    {
        std::cout << " Building the exporter " << std::endl;
    }


    ExporterHDF5<mesh_Type> exporter ( dataFile, meshPart.meshPartition(), "solution", Comm->MyPID() );
    exporter.setMultimesh (false);

    std::shared_ptr<vector_type> LSExported ( new vector_type (LSSolutionOld, Repeated) );
    std::shared_ptr<vector_type> NSExported ( new vector_type (NSSolutionOld, Repeated) );

    const UInt PressureOffset ( 3 * uFESpace->dof().numTotalDof() );

    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "level-set", lsFESpace, LSExported, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace, NSExported, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace, NSExported, PressureOffset);

    if (verbose)
    {
        std::cout << " Loading time stuff " << std::endl;
    }


    Real currentTime ( dataFile ("time/initial", 0.0) );
    Real finalTime ( dataFile ("time/final", 1.0) );
    Real dt ( dataFile ("time/dt", 0.1) );
    UInt niter (0);


    if (verbose)
    {
        std::cout << " Exporting the initial condition " << std::endl;
    }


    exporter.postProcess (currentTime);


    if (verbose)
    {
        std::cout << std::endl;
    }
    if (verbose)
    {
        std::cout << " ### Simulation times ### " << std::endl;
    }
    if (verbose)
    {
        std::cout << " From " << currentTime << " to " << finalTime << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time step: " << dt << std::endl;
    }
    if (verbose)
    {
        std::cout << std::endl;
    }




    while (currentTime < finalTime)
    {
        LifeChrono ChronoIteration;
        ChronoIteration.start();

        currentTime += dt;
        niter += 1;

        if (verbose)
        {
            std::cout << std::endl;
        }
        if (verbose)
        {
            std::cout << "----------------------------" << std::endl;
        }
        if (verbose)
        {
            std::cout << " Time : " << currentTime << std::endl;
        }
        if (verbose)
        {
            std::cout << " Iter : " << niter << std::endl;
        }
        if (verbose)
        {
            std::cout << "----------------------------" << std::endl;
        }
        if (verbose)
        {
            std::cout << std::endl;
        }

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Assembling the matrix ... " << std::flush;
        }

        LifeChrono ChronoItem;
        ChronoItem.start();

#define SUPG_TEST value(NSSupg)*h_K * (grad(phi_i)*eval(normalize,value(ETuFESpace,velocitySolutionOld)))

#define PSPG_TEST value(NSPspg)*h_K * grad(phi_i)

#define DIVDIV_TEST value(NSdivdiv)*h_K*div(phi_i)

#ifdef BDF2_TIME
        vector_type velocityExtrapolated (velocitySolutionOld * 2 - velocitySolutionOldOld, Repeated);
        Real alpha (1.5);
        vector_type velocityBdfRHS (velocitySolutionOld * 2 - velocitySolutionOldOld * 0.5, Repeated);
#else
        vector_type velocityExtrapolated (velocitySolutionOld, Repeated);
        Real alpha (1.0);
        vector_type velocityBdfRHS (velocitySolutionOld, Repeated);
#endif


        std::shared_ptr<matrix_block_type> NSMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
        *NSMatrix *= 0.0;

        {
            std::shared_ptr<DensityFct> density (new DensityFct);
            std::shared_ptr<ViscosityFct> viscosity (new ViscosityFct);
            std::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

            using namespace ExpressionAssembly;

            integrate (
                elements (ETuFESpace->mesh() ), // Mesh

                adapt (ETlsFESpace, LSSolutionOld, uFESpace->qr() ), // QR

                ETuFESpace,
                ETuFESpace,

                // NS
                DENSITY
                * ( value (alpha / dt) * dot (phi_i, phi_j)
                    + dot (grad (phi_j) *value (ETuFESpace, velocityExtrapolated), phi_i)
                  )

                + VISCOSITY * dot (grad (phi_i), grad (phi_j) )

                // Div div
                + DIVDIV_TEST
                *div (phi_j)

                // SUPG
                + dot (
                    grad (phi_j) *value (ETuFESpace, velocityExtrapolated) *DENSITY
                    + value (1.0 / dt) *DENSITY * phi_j
                    , SUPG_TEST)

            )
                    >> NSMatrix->block (0, 0);


            integrate (
                elements (ETuFESpace->mesh() ), // Mesh

                //uFESpace->qr(), // QR
                adapt (ETlsFESpace, LSSolutionOld, uFESpace->qr() ), // QR

                ETuFESpace,
                ETpFESpace,

                value (-1.0) *phi_j * div (phi_i)

                // SUPG
                + dot ( grad (phi_j) , SUPG_TEST)


            )
                    >> NSMatrix->block (0, 1);


            integrate (
                elements (ETuFESpace->mesh() ), // Mesh

                //uFESpace->qr(), // QR
                adapt (ETlsFESpace, LSSolutionOld, uFESpace->qr() ), // QR

                ETpFESpace,
                ETuFESpace,

                phi_i * div (phi_j)

                // PSPG
                + dot (
                    grad (phi_j) *value (ETuFESpace, velocityExtrapolated) *DENSITY
                    + value (1.0 / dt) *DENSITY * phi_j
                    , PSPG_TEST)


            )
                    >> NSMatrix->block (1, 0);


            integrate (
                elements (ETuFESpace->mesh() ), // Mesh

                //uFESpace->qr(), // QR
                adapt (ETlsFESpace, LSSolutionOld, uFESpace->qr() ), // QR

                ETpFESpace,
                ETpFESpace,

                // PSPG
                dot (
                    grad (phi_j)
                    , PSPG_TEST)

            )
                    >> NSMatrix->block (1, 1);


        }

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Boundary integral ... " << std::flush;
        }

        ChronoItem.start();

        //QuadratureBoundary myBDQR(buildTetraBDQR(uFESpace->bdQr()));
        QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );

        {
            using namespace ::LifeV::ExpressionAssembly;

            std::shared_ptr<DensityFct> density (new DensityFct);
            std::shared_ptr<ViscosityFct> viscosity (new ViscosityFct);
            std::shared_ptr<HeavisideFct> heaviside (new HeavisideFct);
            std::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

#ifdef NORMAL_CORRECTION
            integrate ( boundary (ETuFESpace->mesh(), WALL),
                        adapt (ETlsFESpace, LSSolutionOld, myBDQR),
                        //myBDQR,

                        ETuFESpace,
                        ETuFESpace,

                        value (-1.0) *VISCOSITY * dot ( (grad (phi_j) *Nface), Nface)
                        * dot ( phi_i, Nface )
                      ) >> NSMatrix->block (0, 0);

            integrate ( boundary (ETuFESpace->mesh(), WALL),
                        myBDQR,

                        ETuFESpace,
                        ETpFESpace,

                        value (1.0) *phi_j * dot ( phi_i, Nface )
                      ) >> NSMatrix->block (0, 1);
#endif

#ifdef FULL_CORRECTION
            integrate ( boundary (ETuFESpace->mesh(), WALL),
                        adapt (ETlsFESpace, LSSolutionOld, myBDQR),
                        //myBDQR,

                        ETuFESpace,
                        ETuFESpace,

                        value (-1.0) *VISCOSITY * dot ( (grad (phi_j) *Nface), phi_i)

                      ) >> NSMatrix->block (0, 0);

            integrate ( boundary (ETuFESpace->mesh(), WALL),
                        myBDQR,

                        ETuFESpace,
                        ETpFESpace,

                        value (1.0) *phi_j * dot ( phi_i, Nface )
                      ) >> NSMatrix->block (0, 1);
#endif

        }


        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Closing the matrix ... " << std::flush;
        }

        ChronoItem.start();

        NSMatrix->globalAssemble();

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Lifting the pressure ... " << std::flush;
        }

        ChronoItem.start();

        vector_type PressurePlus (ETpFESpace->map(), Repeated);
        vector_type PressureMinus (ETpFESpace->map(), Repeated);

        {
            vector_type N0 (GradientRecovery::ZZGradient (ETlsFESpace, LSSolutionOld, 0), Repeated);
            vector_type N1 (GradientRecovery::ZZGradient (ETlsFESpace, LSSolutionOld, 1), Repeated);
            vector_type N2 (GradientRecovery::ZZGradient (ETlsFESpace, LSSolutionOld, 2), Repeated);

            vector_type f0 (ETlsFESpace->map(), Unique);

            lsFESpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (volumeForce0), f0, currentTime);
            vector_type f1 (ETlsFESpace->map(), Unique);
            lsFESpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (volumeForce1), f1, currentTime);
            vector_type f2 (ETlsFESpace->map(), Unique);
            lsFESpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (volumeForce2), f2, currentTime);

            vector_type gn ( (densityPlus - densityMinus) * ( N0 * f0 + N1 * f1 + N2 * f2 ), Repeated);
            vector_type nNorm ( ( N0 * N0 + N1 * N1 + N2 * N2 ), Repeated);

            PressurePlus *= 0.0;
            PressureMinus *= 0.0;

            const UInt nbElements (ETpFESpace->mesh()->numElements() );
            const UInt nbLocalPDof (ETpFESpace->refFE().nbDof() );


            for (UInt iElement (0); iElement < nbElements; ++iElement)
            {
                for (UInt iDof (0); iDof < nbLocalPDof; ++iDof)
                {
                    UInt globalID ( ETlsFESpace->dof().localToGlobalMap (iElement, iDof) );

                    Real lsValue ( LSSolutionOld[globalID] );

                    if (lsValue >= 0.0)
                    {
                        if (NormalizeCorrection)
                        {
                            PressureMinus[globalID] = -lsValue * gn[globalID] / sqrt (nNorm[globalID]) ;
                        }
                        else
                        {
                            PressureMinus[globalID] = -lsValue * gn[globalID];
                        }
                    }
                    else
                    {
                        if (NormalizeCorrection)
                        {
                            PressurePlus[globalID] = lsValue * gn[globalID] / sqrt (nNorm[globalID]) ;
                        }
                        else
                        {
                            PressurePlus[globalID] = lsValue * gn[globalID];
                        }
                    }
                }
            }
        }

        ChronoItem.stop();

        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Assembling the rhs ... " << std::flush;
        }

        ChronoItem.start();

        vector_type forceRhs (ETuFESpace->map(), Unique);

        uFESpace->interpolate (static_cast<FESpace< mesh_Type, MapEpetra >::function_Type> (volumeForce), forceRhs, currentTime);

        vector_block_type NSRhs ( ETuFESpace->map() | ETpFESpace->map() , Repeated );
        NSRhs *= 0.0;

        {
            using namespace ExpressionAssembly;

            std::shared_ptr<DensityFct> density (new DensityFct);
            std::shared_ptr<ViscosityFct> viscosity (new ViscosityFct);
            std::shared_ptr<HeavisideFct> heaviside (new HeavisideFct);
            std::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

            integrate (
                elements (ETuFESpace->mesh() ), // Mesh

                adapt (ETlsFESpace, LSSolutionOld, uFESpace->qr() ), // QR

                ETuFESpace,

                DENSITY
                * dot ( value (1 / dt) *value (ETuFESpace, velocityBdfRHS) + value (ETuFESpace, forceRhs)
                        , phi_i)

                // Pressure correction
                + value (1.0)
                * ( eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) *value (ETpFESpace, PressurePlus)
                    + (value (1.0) - eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) ) *value (ETpFESpace, PressureMinus)  )
                *div (phi_i)

                // SUPG
                + dot (
                    DENSITY * value (1 / dt) *value (ETuFESpace, velocityBdfRHS)
                    + DENSITY * value (ETuFESpace, forceRhs)
                    - ( eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) *grad (ETpFESpace, PressurePlus)
                        + (value (1.0) - eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) ) *grad (ETpFESpace, PressureMinus) )
                    , SUPG_TEST)


            )
                    >> NSRhs.block (0);

            integrate (
                elements (ETuFESpace->mesh() ), // Mesh

                adapt (ETlsFESpace, LSSolutionOld, pFESpace->qr() ), // QR

                ETpFESpace,


                // PSPG
                dot (
                    DENSITY * value (1 / dt) *value (ETuFESpace, velocityBdfRHS)
                    + DENSITY * value (ETuFESpace, forceRhs)
                    - ( eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) *grad (ETpFESpace, PressurePlus)
                        + (value (1.0) - eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) ) *grad (ETpFESpace, PressureMinus) )
                    , PSPG_TEST)


            )
                    >> NSRhs.block (1);

#if defined(FULL_CORRECTION) || defined(NORMAL_CORRECTION)

            // BOUNDARY INTEGRAL FOR THE PRESSURE CORRECTION
            integrate ( boundary (ETuFESpace->mesh(), WALL),
                        //myBDQR,
                        adapt (ETlsFESpace, LSSolutionOld, myBDQR),
                        ETuFESpace,

                        value (-1.0) * ( eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) *value (ETpFESpace, PressurePlus)
                                         + (value (1.0) - eval (heaviside, value (ETlsFESpace, LSSolutionOld) ) ) *value (ETpFESpace, PressureMinus) )
                        * dot (phi_i, Nface)

                      )
                    >> NSRhs.block (0);
#endif

        }

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Closing the rhs ... " << std::flush;
        }

        ChronoItem.start();

        vector_block_type NSRhsUnique ( NSRhs, Unique );

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Applying boundary conditions ... " << std::flush;
        }

        ChronoItem.start();

        BCHandler NSBCHandler;
        BCFunctionBase ZeroBCFct (zeroFct);
        BCFunctionBase OneBCFct (oneFct);

        std::vector<UInt> compXY (2);
        compXY[0] = 0;
        compXY[1] = 1;

        // Make the Robin coefficient
        vector_type LSReplicated (ETlsFESpace->map() + ETlsFESpace->map() + ETlsFESpace->map() );
        LSReplicated *= 0.0;
        LSReplicated.add (LSSolution, 0);
        LSReplicated.add (LSSolution, ETlsFESpace->dof().numTotalDof() );
        LSReplicated.add (LSSolution, 2 * ETlsFESpace->dof().numTotalDof() );

        // Assume same space!
        LSReplicated.apply (slipFct);

        vector_type RobinVec (0.0 * velocitySolution, Repeated);
        vector_type RobinCoeff (LSReplicated, Repeated);

        BCVector uRobin (RobinVec, ETuFESpace->dof().numTotalDof(), 0);
        uRobin.setRobinCoeffVector (RobinCoeff);
        uRobin.setBetaCoeff (0.0);

        BCFunctionDirectional directionFct (zeroFct, normalDirection);

        NSBCHandler.addBC ("BottomL", BOTTOM_LINE, Essential, Full, ZeroBCFct , 3);
        NSBCHandler.addBC ("TopL", TOP_LINE, Essential, Full, ZeroBCFct , 3);
        NSBCHandler.addBC ("Bottom", BOTTOM, Essential, Full, ZeroBCFct , 3);

        NSBCHandler.addBC ("Wall", WALL, Essential, Normal, ZeroBCFct);
        //NSBCHandler.addBC("Wall", WALL, Essential, Component, ZeroBCFct ,compXY);
        //NSBCHandler.addBC("Wall", WALL, Essential, Directional, directionFct);

        NSBCHandler.addBC ("Wall", WALL, Robin, Full, uRobin , 3);
        //NSBCHandler.addBC("Wall", WALL, Natural, Full, OneBCFct ,3);

        NSBCHandler.addBC ("Top", TOP, Essential, Full, ZeroBCFct , 3);

        NSBCHandler.bcUpdate ( *meshPart.meshPartition(), uFESpace->feBd(), uFESpace->dof() );

        bcManage (*NSMatrix, NSRhsUnique,
                  *uFESpace->mesh(), uFESpace->dof(),
                  NSBCHandler, uFESpace->feBd(), 1.0, currentTime);

        //NSMatrix->diagonalize(3*ETuFESpace->dof().numTotalDof(),1.0,NSRhsUnique,0.0);

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        if (verbose)
        {
            std::cout << "[Navier-Stokes] Solving the system " << std::endl;
        }


        NSSolver.setMatrix (*NSMatrix);


        std::shared_ptr<matrix_type> NSMatrixNoBlock (new matrix_type ( NSMatrix->map() ) );
        *NSMatrixNoBlock += *NSMatrix;

        NSSolver.solveSystem (NSRhsUnique, NSSolution, NSMatrixNoBlock);


        if (verbose)
        {
            std::cout << "[Navier-Stokes] Time advancing ... " << std::flush;
        }

        ChronoItem.start();

        NSSolutionOld = NSSolution;
#ifdef BDF2_TIME
        velocitySolutionOldOld = velocitySolutionOld;
#endif
        velocitySolutionOld.subset (NSSolutionOld);

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << std::endl;
        }


        if (verbose)
        {
            std::cout << "[LS advection] Assembling the matrix ... " << std::flush;
        }

        ChronoItem.start();

        std::shared_ptr<matrix_type> LSAdvectionMatrix (new matrix_type ( ETlsFESpace->map() ) );
        *LSAdvectionMatrix *= 0.0;

        {
            using namespace ExpressionAssembly;

            std::shared_ptr<SignFct> sign (new SignFct);
            std::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

            integrate (
                elements (ETlsFESpace->mesh() ), // Mesh

                //adapt(ETlsFESpace,LSSolutionOld, lsFESpace->qr()), // QR
                lsFESpace->qr(),

                ETlsFESpace,
                ETlsFESpace,

                // Time derivative
                value (1 / dt) *phi_i * phi_j

                // Advection
                + dot ( value (ETuFESpace, velocitySolutionOld) , grad (phi_j) ) *phi_i

                // SUPG stabilization
                + ( value (LSSupg) *h_K )
                * ( value (1 / dt) *phi_j + dot ( value (ETuFESpace, velocitySolutionOld) , grad (phi_j) ) )
                //* dot( value(ETuFESpace,velocitySolutionOld) , grad(phi_i) )
                * dot ( eval (normalize, value (ETuFESpace, velocitySolutionOld) ) , grad (phi_i) )

            )
                    >> LSAdvectionMatrix;

        }

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[LS advection] Closing the matrix ... " << std::flush;
        }

        ChronoItem.start();

        LSAdvectionMatrix->globalAssemble();

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[LS advection] Assembling the rhs ... " << std::flush;
        }

        ChronoItem.start();

        vector_type LSAdvectionRhs ( ETlsFESpace->map(), Repeated ); // Repeated right?
        LSAdvectionRhs *= 0.0;

        {
            using namespace ExpressionAssembly;

            std::shared_ptr<SignFct> sign (new SignFct);
            std::shared_ptr<NormalizeFct> normalize (new NormalizeFct);

            integrate (
                elements (ETlsFESpace->mesh() ), // Mesh

                //adapt(ETlsFESpace,LSSolutionOld, lsFESpace->qr()), // QR
                lsFESpace->qr(), // QR

                ETlsFESpace,

                // Time derivative
                value (1 / dt) *value (ETlsFESpace, LSSolutionOld) *phi_i

                // SUPG stabilization
                + ( value (LSSupg) *h_K )
                * (value (1 / dt) *value (ETlsFESpace, LSSolutionOld) )
                //* dot( value(ETuFESpace,velocitySolutionOld) , grad(phi_i) )
                * dot ( eval (normalize, value (ETuFESpace, velocitySolutionOld) ) , grad (phi_i) )
            )
                    >> LSAdvectionRhs;
        }

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[LS advection] Closing the rhs ... " << std::flush;
        }

        ChronoItem.start();

        vector_type LSAdvectionRhsUnique ( LSAdvectionRhs, Unique );

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << "[LS advection] Solving the system " << std::endl;
        }


        LSAdvectionSolver.setMatrix (*LSAdvectionMatrix);
        LSAdvectionSolver.solveSystem (LSAdvectionRhsUnique, LSSolution, LSAdvectionMatrix);


        if (verbose)
        {
            std::cout << "[Hamilton-J.] Initializing ... " << std::flush;
        }

        ChronoItem.start();

        // To keep the original sign
        LSSolutionOld = LSSolution;
        HJSolutionOld = LSSolution;

        HJSolver.resetPreconditioner();

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }

        HJtime = 0.0;

        if (niter % HJeach == 0)
        {

            while ( HJtime < HJfinalTime)
            {
                HJtime += HJdt;

                if (verbose)
                {
                    std::cout << std::endl;
                }
                if (verbose)
                {
                    std::cout << "[Hamilton-J.] Time :  " << HJtime <<  std::endl;
                }
                if (verbose)
                {
                    std::cout << "[Hamilton-J.] Building the matrix ... " << std::flush;
                }


                ChronoItem.start();

                std::shared_ptr<matrix_type> HJMatrix (new matrix_type ( ETlsFESpace->map() ) );
                *HJMatrix *= 0.0;

#define BETA ( eval(sign, value(ETlsFESpace,LSSolutionOld) ) * eval(normalization, grad(ETlsFESpace,HJSolutionOld) ) )

                {
                    using namespace ExpressionAssembly;

                    std::shared_ptr<SignFct> sign (new SignFct);
                    std::shared_ptr<NormalizeFct> normalization (new NormalizeFct);

                    integrate (
                        elements (ETlsFESpace->mesh() ), // Mesh

                        adapt (ETlsFESpace, LSSolutionOld, lsFESpace->qr() ), // QR

                        ETlsFESpace,
                        ETlsFESpace,

                        // Time derivative
                        (
                            value (1 / HJdt) *phi_i * phi_j
                            + dot ( BETA , grad (phi_j) ) *phi_i

                            // SUPG stabilization
                            + ( value (HJsupg) *h_K )
                            * ( value (1 / HJdt) *phi_j + dot ( BETA , grad (phi_j) ) )
                            * dot ( BETA , grad (phi_i) )
                        )
                        * (value (1.0) - ifCrossed (ETlsFESpace, LSSolutionOld) )

                        +
                        value (HJpenal) * ( phi_i * phi_j  )
                        *ifCrossed (ETlsFESpace, LSSolutionOld)

                    )
                            >> HJMatrix;

                }

                ChronoItem.stop();
                if (verbose)
                {
                    std::cout << ChronoItem.diff() << " s" << std::endl;
                }


                if (verbose)
                {
                    std::cout << "[Hamilton-J.] Closing the matrix ... " << std::flush;
                }

                ChronoItem.start();

                HJMatrix->globalAssemble();

                ChronoItem.stop();
                if (verbose)
                {
                    std::cout << ChronoItem.diff() << " s" << std::endl;
                }


                if (verbose)
                {
                    std::cout << "[Hamilton-J.] Assembling the rhs " << std::flush;
                }

                ChronoItem.start();

                vector_type HJRhs ( ETlsFESpace->map(), Repeated );
                HJRhs *= 0.0;

                {
                    using namespace ExpressionAssembly;

                    std::shared_ptr<SignFct> sign (new SignFct);
                    std::shared_ptr<NormalizeFct> normalization (new NormalizeFct);
                    std::shared_ptr<NormFct> norm (new NormFct);

                    integrate (
                        elements (ETlsFESpace->mesh() ), // Mesh

                        adapt (ETlsFESpace, LSSolutionOld, lsFESpace->qr() ), // QR

                        ETlsFESpace,

                        (
                            value (1 / HJdt) *value (ETlsFESpace, HJSolutionOld) *phi_i
                            + eval (sign, value (ETlsFESpace, LSSolutionOld) ) *phi_i

                            + ( value (HJsupg) *h_K )
                            * ( value (1 / HJdt) *value (ETlsFESpace, HJSolutionOld) + eval (sign, value (ETlsFESpace, LSSolutionOld) ) )
                            * dot ( BETA, grad (phi_i) )
                        ) * (value (1.0) - ifCrossed (ETlsFESpace, LSSolutionOld) )

                        + value (HJpenal) * ( phi_i * value (ETlsFESpace, LSSolutionOld) / eval (norm, grad (ETlsFESpace, LSSolutionOld) ) ) * ifCrossed (ETlsFESpace, LSSolutionOld)
                    )

                            >> HJRhs;
                }

                ChronoItem.stop();
                if (verbose)
                {
                    std::cout << ChronoItem.diff() << " s" << std::endl;
                }


                if (verbose)
                {
                    std::cout << "[Hamilton-J.] Closing the rhs ... " << std::flush;
                }

                ChronoItem.start();

                vector_type HJRhsUnique ( HJRhs, Unique );

                ChronoItem.stop();
                if (verbose)
                {
                    std::cout << ChronoItem.diff() << " s" << std::endl;
                }


                if (verbose)
                {
                    std::cout << "[Hamilton-J.] Solving the system ... " << std::flush;
                }

                ChronoItem.start();

                HJSolver.setMatrix (*HJMatrix);
                HJSolver.solveSystem (HJRhsUnique, LSSolution, HJMatrix);

                ChronoItem.stop();
                if (verbose)
                {
                    std::cout << ChronoItem.diff() << " s" << std::endl;
                }


                if (verbose)
                {
                    std::cout << "[Hamilton-J.] Update the solution ... " << std::flush;
                }

                ChronoItem.start();

                HJSolutionOld = LSSolution;

                ChronoItem.stop();
                if (verbose)
                {
                    std::cout << ChronoItem.diff() << " s" << std::endl;
                }

                if (verbose)
                {
                    std::cout << std::endl;
                }
            } // end of the pseudo-time loop
        } // end if

        if (verbose)
        {
            std::cout << "[LS correction] Correction of the volume ... " << std::flush;
        }

        ChronoItem.start();

        Real localVolume = 0.0;
        Real globalVolume = 0.0;

        Real localArea = 0.0;
        Real globalArea = 0.0;

        {
            using namespace ExpressionAssembly;

            std::shared_ptr<HeavisideFct> heaviside (new HeavisideFct);
            std::shared_ptr<SmoothDeltaFct> delta (new SmoothDeltaFct);

            integrate (
                elements (ETlsFESpace->mesh() ), // Mesh
                adapt (ETlsFESpace, LSSolution, lsFESpace->qr() ), // QR
                eval (heaviside, value (ETlsFESpace, LSSolution) ) // Expression
            )
                    >> localVolume;

            integrate (
                elements (ETlsFESpace->mesh() ), // Mesh
                adapt (ETlsFESpace, LSSolution, lsFESpace->qr() ), // QR
                eval (delta, value (ETlsFESpace, LSSolution) ) // Expression
            )
                    >> localArea;

        }

        Comm->Barrier();
        Comm->SumAll (&localVolume, &globalVolume, 1);

        Comm->Barrier();
        Comm->SumAll (&localArea, &globalArea, 1);


        LSSolution += volumeRelaxation * (exactVolume - globalVolume) / globalArea;


        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }
        if (verbose)
        {
            std::cout << "[LS correction]     Rel. vol. diff:  " << (exactVolume - globalVolume) / exactVolume << std::endl;
        }
        if (verbose)
        {
            std::cout << "[LS correction]     Approx. area:  " << globalArea << std::endl;
        }


        if (verbose)
        {
            std::cout << "[LS advection] Time advancing ... " << std::flush;
        }

        ChronoItem.start();

        LSSolutionOld = LSSolution;

        ChronoItem.stop();
        if (verbose)
        {
            std::cout << ChronoItem.diff() << " s" << std::endl;
        }


        if (verbose)
        {
            std::cout << " Exporting " << std::endl;
        }


        *LSExported = LSSolutionOld;
        *NSExported = NSSolutionOld;

        if (niter % exportEach == 0)
        {
            exporter.postProcess (currentTime);
        }

        ChronoIteration.stop();
        if (verbose)
        {
            std::cout << std::endl << " Total iteration time : " << ChronoIteration.diff() << " s" << std::endl;
        }

        if (verbose)
        {
            displayMotionParameters (currentTime);
        }

    } // end time loop


    exporter.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif


    return ( EXIT_SUCCESS );
}
