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

//#include "life/lifesolver/NavierStokesSolver.hpp"
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/fem/FESpace.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEnsight.hpp>

#include <lifev/navier_stokes/solver/OseenSolver.hpp>

#include "flowConditions.hpp"
#include "resistance.hpp"
#include <iostream>


using namespace LifeV;

const int INLET       = 3;
const int WALL        = 200;
const int OUTLET      = 2;
const int RINGIN      = 30;
const int RINGOUT     = 20;


void
postProcessFluxesPressures ( OseenSolver< RegionMesh<LinearTetra> >& nssolver,
                             BCHandler& bcHandler,
                             const LifeV::Real& t, bool _verbose )
{
    LifeV::Real Q, P;
    UInt flag;

    for ( BCHandler::bcBaseIterator_Type it = bcHandler.begin();
            it != bcHandler.end(); ++it )
    {
        flag = it->flag();

        Q = nssolver.flux (flag);
        P = nssolver.pressure (flag);

        if ( _verbose )
        {
            std::ofstream outfile;
            std::stringstream filenamess;
            std::string filename;

            // file name contains the label
            filenamess << flag;
            // writing down fluxes
            filename = "flux_label" + filenamess.str() + ".m";
            outfile.open (filename.c_str(), std::ios::app);
            outfile << Q << " " << t << "\n";
            outfile.close();
            // writing down pressures
            filename = "pressure_label" + filenamess.str() + ".m";
            outfile.open (filename.c_str(), std::ios::app);
            outfile << P << " " << t << "\n";
            outfile.close();
            // reset ostringstream
            filenamess.str ("");
        }
    }

}


struct ResistanceTest::Private
{
    Private() :
        nu (1),
        H (1), D (1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_Type;

    double Re;

    std::string data_file_name;

    double      nu;  /**< viscosity (in m^2/s) */
    //const double rho; /**< density is constant (in kg/m^3) */
    double      H;   /**< height and width of the domain (in m) */
    double      D;   /**< diameter of the cylinder (in m) */
    bool        centered; /**< true if the cylinder is at the origin */

    std::string initial_sol;

    boost::shared_ptr<Epetra_Comm>   comm;

    // Static boost functions to impose boundary conditions
    // Inlet BCs (for this test Poiseuille)
    static Real fluxFunctionAneurysm (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
    {

        Real fluxFinal(0.0);
        Real rampAmpl (0.2);
        Real dt (0.001);

        if ( t <= rampAmpl )
        {
            fluxFinal = ( 0.09403 / rampAmpl) * t;
        }
        else
        {


            // We change the flux for our geometry
            const Real pi   = 3.141592653589793;
            const Real area = 0.1950; // BigMesh

            const Real areaFactor = area / ( 0.195 * 0.195 * pi);
            //const Real Average = (48.21 * pow (area, 1.84) ) * 60; //Mean Cebral's Flux per minut

            // Unit conversion from ml/min to cm^3/s
            const Real unitFactor = 1. / 60.;

            // T is the period of the cardiac cycle
            const Real T          = 0.8;

            // a0 is the average VFR (the value is taken from Karniadakis p970)
            const Real a0         = 255;
            //const Real volumetric = Average / a0; //VolumetricFactor per minut

            // Fourrier
            const Int M (7);
            const Real a[M] = { -0.152001, -0.111619, 0.043304, 0.028871, 0.002098, -0.027237, -0.000557};
            const Real b[M] = { 0.129013, -0.031435, -0.086106, 0.028263, 0.010177, 0.012160, -0.026303};

            Real flux (0);
            //      const Real xi(2*pi*t/T);
            const Real xi (2 * pi * (t - rampAmpl + dt) / T);

            flux = a0;
            Int k (1);
            for (; k <= M ; ++k)
            {
                flux += a0 * (a[k - 1] * cos (k * xi) + b[k - 1] * sin (k * xi) );
            }

            //return - (flux * areaFactor * unitFactor);
            fluxFinal =  (flux * areaFactor * unitFactor);
        }

        return fluxFinal;

    }

    static Real aneurismFluxInVectorial (const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i)
    {
        Real n1 (-0.67463);
        Real n2 (-0.19861);
        Real n3 (-0.71094);

        Real x0 (6.752113);
        Real y0 (9.398349);
        Real z0 (18.336601);

        Real flux (fluxFunctionAneurysm (t, x, y, z, i) );

        Real area (0.195);

        //Parabolic profile
        Real radius ( std::sqrt ( area / 3.14159265359 ) );

        Real radiusSquared = radius * radius;
        Real peak (0);
        peak = ( 2 * flux ) / ( area );

        switch (i)
        {
            case 0:
                // Flat profile: flux / area;
                // return n1 * flux / area;
                return n1 * std::max ( 0.0, peak * ( (radiusSquared - ( (x - x0) * (x - x0) + (y - y0) * (y - y0) ) ) / radiusSquared) );
            case 1:
                // Flat profile: flux / area;
                //return n2 * flux / area;
                return n2 * std::max ( 0.0 , peak * ( (radiusSquared - ( (x - x0) * (x - x0) + (y - y0) * (y - y0) ) ) / radiusSquared) );
            case 2:
                // Flat profile: flux / area;
                // return n3 * flux / area;
                return n3 * std::max ( 0.0, peak * ( (radiusSquared - ( (x - x0) * (x - x0) + (y - y0) * (y - y0) ) ) / radiusSquared) );
            default:
                return 0.0;
        }
    }

    // External BCs ( u = 0 )
    static Real zeroBCF (const Real&  /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
    {
        switch (i)
        {
            case 0:
                //Flat profile: flux / area;
                //return n1 * flux / area;
                return 0;
            case 1:
                //Flat profile: flux / area;
                //return n2 * flux / area;
                return 0;
            case 2:
                // Flat profile: flux / area;
                // return n3 * flux / area;
                return 0;
            default:
                return 0.0;
        }
    }

    // Outflow BC ( resistance BC - applied explicitly )
    // Defined by the class FlowConditions


};

ResistanceTest::ResistanceTest ( int argc,
                         char** argv )
    :
    parameters ( new Private )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );
    parameters->data_file_name = data_file_name;

    parameters->Re          = dataFile ( "fluid/problem/Re", 1. );
    parameters->nu          = dataFile ( "fluid/physics/viscosity", 1. ) /
                              dataFile ( "fluid/physics/density", 1. );
    parameters->H           = 20.;//dataFile( "fluid/problem/H", 20. );
    parameters->D           =               dataFile ( "fluid/problem/D", 1. );
    parameters->centered    = (bool)        dataFile ( "fluid/problem/centered", 0 );
    parameters->initial_sol = (std::string) dataFile ( "fluid/problem/initial_sol", "stokes");
    std::cout << parameters->initial_sol << std::endl;


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    int ntasks = 0;
    parameters->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );

#else
    parameters->comm.reset ( new Epetra_SerialComm() );
#endif

}

void
ResistanceTest::run()

{
    typedef RegionMesh<LinearTetra>               mesh_Type;
    typedef FESpace< mesh_Type, MapEpetra >       feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type>       feSpacePtr_Type;
    typedef OseenSolver< mesh_Type >::vector_Type vector_Type;
    typedef boost::shared_ptr<vector_Type>        vectorPtr_Type;
    // Reading from data file
    //
    GetPot dataFile ( parameters->data_file_name );

    //    int save = dataFile("fluid/miscellaneous/save", 1);

    bool verbose = (parameters->comm->MyPID() == 0);

    // Lagrange multiplier for flux BCs
    int numLM = 0;

    boost::shared_ptr<OseenData> oseenData (new OseenData() );
    oseenData->setup ( dataFile );

    MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");

    boost::shared_ptr<mesh_Type> fullMeshPtr ( new mesh_Type ( parameters->comm ) );
    readMesh (*fullMeshPtr, meshData);

    boost::shared_ptr<mesh_Type> meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, parameters->comm);
        meshPtr = meshPart.meshPartition();
    }

    //oseenData.meshData()->setMesh(meshPtr);

    std::string uOrder =  dataFile ( "fluid/space_discretization/vel_order", "P1");
    if (verbose)
    {
        std::cout << "Building the velocity FE space ... " << std::flush;
    }

    feSpacePtr_Type uFESpacePtr ( new feSpace_Type (meshPtr, uOrder, 3, parameters->comm) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }


    std::string pOrder =  dataFile ( "fluid/space_discretization/press_order", "P1");

    if (verbose)
    {
        std::cout << "Building the pressure FE space ... " << std::flush;
    }

    feSpacePtr_Type pFESpacePtr ( new feSpace_Type (meshPtr, pOrder, 1, parameters->comm) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    UInt totalVelDof   = uFESpacePtr->map().map (Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpacePtr->map().map (Unique)->NumGlobalElements();

    // Boundary conditions
    BCHandler bcH;
    BCFunctionBase uIn   (  Private::aneurismFluxInVectorial );
    BCFunctionBase uZero (  Private::zeroBCF );

    FlowConditions outFlowBC;

    // Read the resistance and hydrostatic pressure from data file
    Real resistance = dataFile ( "fluid/physics/resistance", 0.0 );
    Real hydrostatic = dataFile ( "fluid/physics/hydrostatic", 0.0 );

    // outFlowBC.initParameters( OUTLET, resistance, hydrostatic, "outlet-3" );

    vectorPtr_Type resistanceImplicit;
    resistanceImplicit.reset( new vector_Type( uFESpacePtr->map(), Repeated ) );
    resistanceImplicit->epetraVector().PutScalar(0.0);

    BCVector resistanceBCdefinition;
    resistanceBCdefinition.setRhsVector( *resistanceImplicit, uFESpacePtr->dof().numTotalDof(), 1 );
    resistanceBCdefinition.setResistanceCoeff( 600000.0 );


    // Explicit Resistance BC
    //BCFunctionBase resistanceBC( FlowConditions::outPressure0 );
    //cylinder

    bcH.addBC ( "Inlet",    INLET,    Essential,   Full,  uIn  , 3 );
    bcH.addBC ( "Wall",     WALL,     Essential,   Full,  uZero, 3 );

    // Explicit Resistance BC
    // bcH.addBC ( "Outlet",   OUTLET,   Natural,     Normal, resistanceBC );
    bcH.addBC ( "Outlet",   OUTLET,   Resistance, Full, resistanceBCdefinition, 3 );


    if (verbose)
    {
        std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;
    }

    if (verbose)
    {
        std::cout << "Calling the fluid constructor ... ";
    }

    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpacePtr,
                                    *pFESpacePtr,
                                    parameters->comm, numLM);
    MapEpetra fullMap (fluid.getMap() );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    fluid.setUp (dataFile);

    // Setting up the utility for post-processing
    fluid.setupPostProc( );

    std::cout << "Inlet area: " << fluid.area( INLET ) << std::endl;
    fluid.buildSystem();

    MPI_Barrier (MPI_COMM_WORLD);

    // Initialization

    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();


    // bdf object to store the previous solutions

    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup (oseenData->dataTimeAdvance()->orderBDF() );

    vector_Type beta ( fullMap );
    vector_Type rhs ( fullMap );

    std::string const exportFileName = dataFile ( "exporter/nameFile", "resistance");

#ifdef HAVE_HDF5
    ExporterHDF5<mesh_Type > exporter ( dataFile, meshPtr, exportFileName, parameters->comm->MyPID() );
#else
    ExporterEnsight<mesh_Type > exporter ( dataFile, meshPtr, exportFileName, parameters->comm->MyPID() );
#endif

    vectorPtr_Type velAndPressure ( new vector_Type (*fluid.solution(), exporter.mapType() ) );

    exporter.addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpacePtr,
                           velAndPressure, UInt (0) );

    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                           velAndPressure, UInt (3 * uFESpacePtr->dof().numTotalDof() ) );

    // Setting up the initial condition
    bdf.bdfVelocity().setInitialCondition ( *fluid.solution() );

    // Temporal loop

    LifeChrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt / 2.; time += dt, iter++)
    {

        std::cout << "Inlet area: " << fluid.area ( INLET ) << std::endl;

        // Updating the Neumann BC for resistance
        outFlowBC.renewParameters ( fluid, *velAndPressure );
        // if ( verbose )
        // {
        //     std::cout << std::endl;
        //     std::cout << "Name: "       << outFlowBC.name() << std::endl;
        //     std::cout << "Resistance: " << outFlowBC.resistance() << std::endl;
        //     std::cout << "Flow: "       << outFlowBC.flow() << std::endl;
        //     std::cout << "Hydrostatic: " << outFlowBC.hydrostatic() << std::endl;
        //     std::cout << "Total Pressure: " << outFlowBC.outP() << std::endl;
        //     std::cout << std::endl;
        // }

        oseenData->dataTime()->setTime (time);

        chrono.start();

        double alpha = bdf.bdfVelocity().coefficientFirstDerivative ( 0 ) / oseenData->dataTime()->timeStep();
        bdf.bdfVelocity().extrapolation (beta);
        bdf.bdfVelocity().updateRHSContribution ( oseenData->dataTime()->timeStep() );
        rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();

        fluid.updateSystem ( alpha, beta, rhs );
        fluid.iterate ( bcH );

        bdf.bdfVelocity().shiftRight ( *fluid.solution() );

        *velAndPressure = *fluid.solution();

        // if ( verbose )
        // {
        //     std::cout << "Post-processing!" << std::endl;
        // }
        exporter.postProcess ( time );

        MPI_Barrier (MPI_COMM_WORLD);

        chrono.stop();
        if (verbose)
        {
            std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
        }
    }

}


//////////////////////


