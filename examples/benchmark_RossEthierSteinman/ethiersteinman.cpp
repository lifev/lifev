/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Gilles Fourestey <gilles.fourestey@epfl.ch>
             Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2011-03-09

  Copyright (C) 2010 EPFL

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
   \file ethiersteiman.cpp
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2011-03-08
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
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifemesh/MeshPartitioner.hpp>
#include <life/lifemesh/MeshData.hpp>
#include <life/lifesolver/OseenData.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/TimeAdvanceBDFNavierStokes.hpp>
#include <life/lifefilters/ExporterEnsight.hpp>
#include <life/lifefilters/ExporterHDF5.hpp>
#include <life/lifefilters/ExporterEmpty.hpp>
#include <life/lifemesh/RegionMesh3DStructured.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "ethiersteinman.hpp"

using namespace LifeV;


Real zero_scalar( const Real& /* t */,
                  const Real& /* x */,
                  const Real& /* y */,
                  const Real& /* z */,
                  const ID& /* i */ )
{
    return 0.;
}



std::set<UInt> parseList( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    std::string::size_type commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );
        setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
        stringList = stringList.substr( commaPos+1 );
    }
    setList.insert( atoi( stringList.c_str() ) );
    return setList;
}



struct Ethiersteinman::Private
{
    Private() :
            nu    (1),
            steady(0)
    {}

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_Type;

    double         Re;

    std::string    data_file_name;

    double         nu;  /**< viscosity (in m^2/s) */
    //const double rho; /**< density is constant (in kg/m^3) */

    bool                             steady;
    boost::shared_ptr<Epetra_Comm>   comm;
};


Ethiersteinman::Ethiersteinman( int argc,
                                char** argv )
        :
        d( new Private )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    d->data_file_name = data_file_name;

    d->Re = dataFile( "fluid/problem/Re", 1. );
    d->nu = dataFile( "fluid/physics/viscosity", 1. ) /
            dataFile( "fluid/physics/density", 1. );

    // Test type
    string testType = dataFile("RossEthierSteinman/test", "none");
    if(testType == "none")
    {
        M_test = None;
    }
    else if(testType == "accuracy")
    {
        M_test = Accuracy;
    }
    else if(testType == "space_convergence")
    {
        M_test = SpaceConvergence;
    }
    else
    {
        std::cout << "[Error] Unknown test method" << std::endl;
        exit(1);
    }

    M_convTol     = dataFile("RossEthierSteinman/space_convergence_tolerance", 1.0);
    M_accuracyTol = dataFile("RossEthierSteinman/accuracy_tolerance", 1.0);

    // Method of initialization
    string initType = dataFile("RossEthierSteinman/initialization", "projection");
    if(initType == "projection")
    {
        M_initMethod = Projection;
    }
    else if(initType == "interpolation")
    {
        M_initMethod = Interpolation;
    }
    else
    {
        std::cout << "[Error] Unknown initialization method" << std::endl;
        exit(1);
    }

    M_exportNorms = dataFile("RossEthierSteinman/export_norms", false);
    M_exportExactSolutions = dataFile("RossEthierSteinman/export_exact_solutions", false);


    std::string meshSource =  dataFile( "RossEthierSteinman/mesh_source", "regular_mesh");
    if(meshSource == "regular_mesh")
    {
        M_meshSource = RegularMesh;
    }
    else if(meshSource == "file")
    {
        M_meshSource = File;
    }
    else
    {
        std::cout << "[Error] Unknown mesh source" << std::endl;
        exit(1);
    }

    //Checking the consistency of the data
    if(M_meshSource == File && M_test == SpaceConvergence)
    {
        std::cout << "[Error] You cannot use mesh files to test the space convergence." << std::endl;
        exit(1);
    }


#ifdef EPETRA_MPI

    //    MPI_Init(&argc,&argv);

    d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm.reset( new Epetra_SerialComm() );
#endif

}

void
Ethiersteinman::computeErrors(const vector_Type& velocityAndPressureSolution,
                              LifeV::Real& uL2Error, LifeV::Real& uRelError, fespace_Type& uFESpace,
                              LifeV::Real& pL2Error, LifeV::Real& pRelError, fespace_Type& pFESpace,
                              LifeV::Real time)
{
    // Computation of the error
    vector_Type vel  (uFESpace.map(), Repeated);
    vector_Type press(pFESpace.map(), Repeated);
    vector_Type velpressure ( velocityAndPressureSolution, Repeated );

    velpressure = velocityAndPressureSolution;
    vel.subset(velpressure);
    press.subset(velpressure, uFESpace.dim()*uFESpace.fieldDim());

    uL2Error = uFESpace.l2Error (Problem::uexact, vel  , time, &uRelError );
    pL2Error = pFESpace.l20Error(Problem::pexact, press, time, &pRelError );
}

bool
Ethiersteinman::checkConvergenceRate(const std::vector<std::string>& uFELabel,
                                     const std::vector<std::vector<LifeV::Real> >& uL2Error,
                                     const std::vector<UInt>& uConvergenceOrder,
                                     const std::vector<std::string>& pFELabel,
                                     const std::vector<std::vector<LifeV::Real> > pL2Error,
                                     const std::vector<UInt>& pConvergenceOrder,
                                     const std::vector<UInt>& meshDiscretizations,
                                     LifeV::Real convTolerance)
{
    // We want to check the convergence of the error and
    // see if it matches the theory.
    std::cout << "Checking the convergence:" << std::endl;

    // Test variable
    bool success(true); // Variable to keep trace of a previous error
    Real h1(0.0), h2(0.0); // Space discretization step
    Real uBound(0.0), pBound(0.0); // Velocity and pressure bounds
    Real uErrRatio(0.0), pErrRatio(0.0); // Ratio of the error E1/E2
    std::string status(""); // Information string

    UInt FEnumber(uFELabel.size());
    UInt discretizationNumber(meshDiscretizations.size());

    for (UInt iElem(0); iElem<FEnumber; ++iElem)
    {
        std::cout << "    - " << uFELabel[iElem] << "-" << pFELabel[iElem] << " ... ";

        // Everything is OK a priori
        status = "OK";

        for (UInt jDiscretization(0); jDiscretization<discretizationNumber-1; ++jDiscretization)
        {
            h1 = 1.0/meshDiscretizations[jDiscretization];
            h2 = 1.0/meshDiscretizations[jDiscretization+1];

            uBound = convTolerance*pow(h1/h2,int(uConvergenceOrder[iElem]));
            pBound = convTolerance*pow(h1/h2,int(pConvergenceOrder[iElem]));

            uErrRatio = uL2Error[iElem][jDiscretization]/uL2Error[iElem][jDiscretization+1]; // E1/E2
            pErrRatio = pL2Error[iElem][jDiscretization]/pL2Error[iElem][jDiscretization+1];

            if (uErrRatio < uBound)
            {
                status = "FAILED";
                success = false;
            }
            if (pErrRatio < pBound)
            {
                status = "FAILED";
                success = false;
            }
        }
        std::cout << status << std::endl;

    }

    return success;
}

void
Ethiersteinman::run()
{
    bool verbose = (d->comm->MyPID() == 0);
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    if(verbose){
        std::cout << " +-----------------------------------------------+" << std::endl
                  << " |    RossEthierSteinman benchmark for LifeV     |" << std::endl
                  << " +-----------------------------------------------+" << std::endl
                  << std::endl
                  << " +-----------------------------------------------+" << std::endl
                  << " |          Authors: Gwenol Grandperrin          |" << std::endl
                  << " |                   Gilles Fourestey            |" << std::endl
                  << " |                   Christophe Prud'homme       |" << std::endl
                  << " |             Date: 2010-03-09                  |" << std::endl
                  << " +-----------------------------------------------+" << std::endl
                  << std::endl;
        if (verbose) std::cout << "[[BEGIN_SIMULATION]]" << std::endl << std::endl;
           std::cout << "[Initilization of MPI]" << std::endl;
#ifdef HAVE_MPI
           std::cout << "Using MPI (" << nproc << " proc.)" << std::endl;
#else
           std::cout << "Using serial version" << std::endl;
#endif
    }

    // +-----------------------------------------------+
    // |             Begining of the test              |
    // +-----------------------------------------------+
    LifeChrono globalChrono;
    LifeChrono runChrono;
    LifeChrono initChrono;
    LifeChrono iterChrono;

    globalChrono.start();
    initChrono.start();

    // +-----------------------------------------------+
    // |               Loading the data                |
    // +-----------------------------------------------+
    if (verbose) std::cout << std::endl << "[Loading the data]" << std::endl;
    GetPot dataFile( d->data_file_name.c_str() );
    if (verbose)
    {
        switch(M_test)
        {
            case Accuracy:
                std::cout << "Test : checks the accuracy of the solution" << std::endl;
                break;
            case SpaceConvergence:
                std::cout << "Test : checks the convergence in space of the solution" << std::endl;
                break;
            case None:
                break;
        }
    }
    Problem::setParamsFromGetPot( dataFile );

    UInt discretizationNumber;
    if((M_test == SpaceConvergence) || (M_test == None))
    {
        // Loading the discretization to be tested
        discretizationNumber = dataFile( "fluid/space_discretization/mesh_number", 1 );
        for ( UInt i( 0 ); i < discretizationNumber; ++i )
        {
            meshDiscretization.push_back(dataFile( "fluid/space_discretization/mesh_size", 8, i ));
        }
    }
    else
    {
        meshDiscretization.push_back(0); // Just to be sure to have 1 element
        discretizationNumber = 1;
    }

    if(M_test == SpaceConvergence)
    {
        // Loading the convergence rate for the finite elements tested
        UInt FEnumber = dataFile( "fluid/space_discretization/FE_number", 1 );
        for ( UInt i( 0 ); i < FEnumber; ++i )
        {
            uConvergenceOrder.push_back(dataFile( "fluid/space_discretization/vel_conv_order_order", 2, i ));
        }
        for ( UInt i( 0 ); i < FEnumber; ++i )
        {
            pConvergenceOrder.push_back(dataFile( "fluid/space_discretization/press_conv_order", 2, i ));
        }
    }

    // Loading the Finite element to be tested
    UInt FEnumber = dataFile( "fluid/space_discretization/FE_number", 1 );
    for ( UInt i( 0 ); i < FEnumber; ++i )
    {
        uFE.push_back(dataFile( "fluid/space_discretization/vel_order", "P1", i ));
    }
    for ( UInt i( 0 ); i < FEnumber; ++i )
    {
        pFE.push_back(dataFile( "fluid/space_discretization/press_order", "P1", i ));
    }

    // Initialization of the errors array
    std::vector<std::vector<LifeV::Real> > uL2Error;
    std::vector<std::vector<LifeV::Real> > pL2Error;
    uL2Error.clear();
    pL2Error.clear();
    std::vector<LifeV::Real> tmpVec(discretizationNumber,0.0);
    for (UInt iElem(0); iElem<FEnumber; ++iElem)
    {
        uL2Error.push_back(tmpVec);
        pL2Error.push_back(tmpVec);
    }

    initChrono.stop();
    if (verbose) std::cout << "Initialization time (pre-run): " << initChrono.diff() << " s." << std::endl;

    // Loop on the mesh refinement
    for (UInt jDiscretization(0); jDiscretization<discretizationNumber; ++jDiscretization)
    {
        UInt mElem = meshDiscretization[jDiscretization];

        // Loop on the finite element
        for (UInt iElem(0); iElem<FEnumber; ++iElem)
        {
            if (verbose) std::cout << std::endl << "[[BEGIN_RUN_" << jDiscretization*FEnumber+iElem << "]]" << std::endl;
            runChrono.reset();
            runChrono.start();
            initChrono.reset();
            initChrono.start();

            if (verbose && M_exportNorms)
            {
                std::string fileName("norm_");
                std::ostringstream oss;
                oss << mElem;
                fileName.append(oss.str());
                fileName.append("_");
                fileName.append(uFE[iElem]);
                fileName.append(pFE[iElem]);
                fileName.append(".txt");
                out_norm.open(fileName.c_str());
                out_norm << "% time / u L2 error / L2 rel error   p L2 error / L2 rel error \n" << std::flush;
            }

            // +-----------------------------------------------+
            // |               Loading the mesh                |
            // +-----------------------------------------------+
            if (verbose) std::cout << "[Loading the mesh]" << std::endl;

            boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);

            // Building the mesh from the source
            if(M_meshSource == RegularMesh)
            {
                regularMesh3D( *fullMeshPtr,
                               1,
                               mElem, mElem, mElem,
                               false,
                               2.0,   2.0,   2.0,
                               -1.0,  -1.0,  -1.0);

                if (verbose) std::cout << "Mesh source: regular mesh("
                                       << mElem << "x" << mElem << "x" << mElem << ")" << std::endl;
            }
            else if(M_meshSource == File)
            {
                MeshData meshData;
                meshData.setup(dataFile, "fluid/space_discretization");
                readMesh(*fullMeshPtr, meshData);

                if (verbose) std::cout << "Mesh source: file("
                                       << meshData.meshDir() << meshData.meshFile() << ")" << std::endl;
            }
            else
            {
                if (verbose) std::cout << std::endl << "Error: Unknown source type for the mesh" << std::endl;
                exit(1);
            }

            if (verbose) std::cout << "Partitioning the mesh ... " << std::flush;
            MeshPartitioner< RegionMesh3D<LinearTetra> >   meshPart(fullMeshPtr, d->comm);
            fullMeshPtr.reset(); //Freeing the global mesh to save memory

            // +-----------------------------------------------+
            // |            Creating the FE spaces             |
            // +-----------------------------------------------+
            if (verbose) std::cout << std::endl << "[Creating the FE spaces]" << std::endl;
            std::string uOrder =  uFE[iElem];
            std::string pOrder =  pFE[iElem];

            if (verbose) std::cout << "FE for the velocity: " << uOrder << std::endl
                                   << "FE for the pressure: " << pOrder << std::endl;

            if (verbose) std::cout << "Building the velocity FE space ... " << std::flush;
            FESpace< RegionMesh3D<LinearTetra>, MapEpetra > uFESpace(meshPart, uOrder, 3, d->comm);
            if (verbose) std::cout << "ok." << std::endl;

            if (verbose) std::cout << "Building the pressure FE space ... " << std::flush;
            FESpace< RegionMesh3D<LinearTetra>, MapEpetra > pFESpace(meshPart, pOrder, 1, d->comm);
            if (verbose) std::cout << "ok." << std::endl;

            UInt totalVelDof   = uFESpace.dof().numTotalDof();
            UInt totalPressDof = pFESpace.dof().numTotalDof();

            if (verbose) std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
            if (verbose) std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;

            // +-----------------------------------------------+
            // |             Boundary conditions               |
            // +-----------------------------------------------+
            if (verbose) std::cout << std::endl << "[Boundary conditions]" << std::endl;
            std::string dirichletList = dataFile( "fluid/problem/dirichletList", "" );
            std::set<UInt> dirichletMarkers = parseList( dirichletList );
            std::string neumannList = dataFile( "fluid/problem/neumannList", "" );
            std::set<UInt> neumannMarkers = parseList( neumannList );

            BCHandler bcH;
            BCFunctionBase uWall( Problem::uexact );
            BCFunctionBase uNeumann( Problem::fNeumann );

            for (std::set<UInt>::const_iterator it = dirichletMarkers.begin();
                    it != dirichletMarkers.end(); ++it)
            {
                bcH.addBC( "Wall", *it, Essential, Full, uWall, 3 );
            }
            for (std::set<UInt>::const_iterator it = neumannMarkers.begin();
                    it != neumannMarkers.end(); ++it)
            {
                bcH.addBC( "Flux", *it, Natural, Full, uNeumann, 3 );
            }

            // If we change the FE we have to update the BCHandler (internal data)
            bcH.bcUpdate( *meshPart.meshPartition(), uFESpace.feBd(), uFESpace.dof());

            // +-----------------------------------------------+
            // |             Creating the problem              |
            // +-----------------------------------------------+
            if (verbose) std::cout<< std::endl << "[Creating the problem]" << std::endl;
            boost::shared_ptr<OseenData> oseenData(new OseenData());
            oseenData->setup( dataFile );

            if (verbose) std::cout << "Time discretization order " << oseenData->dataTime()->orderBDF() << std::endl;

            OseenSolver< RegionMesh3D<LinearTetra> > fluid (oseenData,
                                                      uFESpace,
                                                      pFESpace,
                                                      d->comm);
            MapEpetra fullMap(fluid.getMap());

            fluid.setUp(dataFile);
            fluid.buildSystem();

            MPI_Barrier(MPI_COMM_WORLD);

            // +-----------------------------------------------+
            // |       Initialization of the simulation        |
            // +-----------------------------------------------+
            if (verbose) std::cout<< std::endl << "[Initialization of the simulation]" << std::endl;
            Real dt     = oseenData->dataTime()->timeStep();
            Real t0     = oseenData->dataTime()->initialTime();
            Real tFinal = oseenData->dataTime()->endTime();


            // bdf object to store the previous solutions
            TimeAdvanceBDFNavierStokes<vector_Type> bdf;
            bdf.setup(oseenData->dataTime()->orderBDF());

            // initialization with exact solution: either interpolation or "L2-NS"-projection
            t0 -= dt*bdf.bdfVelocity().order();

            if (verbose) std::cout << "Computing the initial solution ... " << std::endl;

            vector_Type beta( fullMap );
            vector_Type rhs ( fullMap );

            MPI_Barrier(MPI_COMM_WORLD);

            oseenData->dataTime()->setTime(t0);
            fluid.initialize( Problem::uexact, Problem::pexact );

            bdf.bdfVelocity().setInitialCondition( *fluid.solution() );

            //
            // Initial solution loading (interpolation or projection)
            //
            Real time = t0 + dt;
            for (  ; time <=  oseenData->dataTime()->initialTime() + dt/2.; time += dt)
            {

                oseenData->dataTime()->setTime(time);

                beta *= 0.;
                rhs  *= 0.;

                fluid.initialize( Problem::uexact, Problem::pexact );

                beta = *fluid.solution();

                if (M_initMethod == Projection)
                {
                    uFESpace.interpolate(Problem::uderexact, rhs, time);
                    //uFESpace.l2ScalarProduct(Problem::uderexact, rhs, time);
                    rhs *= -1.;
                    rhs = fluid.matrixMass()*rhs;
                    fluid.updateSystem( 0., beta, rhs );
                    fluid.iterate(bcH);
                }

                // Computation of the error
                LifeV::Real urelerr, prelerr, ul2error, pl2error;

                computeErrors(*fluid.solution(),
                              ul2error, urelerr, uFESpace,
                              pl2error, prelerr, pFESpace,
                              time);

                if (verbose && M_exportNorms)
                {
                    out_norm << time  << " "
                    << ul2error << " "
                    << urelerr << " "
                    << pl2error << " "
                    << prelerr << "\n" << std::flush;
                }


                // Updating bdf
                bdf.bdfVelocity().shiftRight( *fluid.solution() );

            }

            fluid.resetPreconditioner();

            boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;

            vectorPtr_Type velAndPressure;
            // only for export -->
            vectorPtr_Type exactPressPtr;
            vector_Type exactPress(pFESpace.map(), Repeated);
            vectorPtr_Type exactVelPtr;
            vector_Type exactVel(uFESpace.map(), Repeated);
            // <--

            std::string const exporterType =  dataFile( "exporter/type", "ensight");

#ifdef HAVE_HDF5
            if (exporterType.compare("hdf5") == 0)
            {
                exporter.reset( new ExporterHDF5<RegionMesh3D<LinearTetra> > ( dataFile, "ethiersteinman" ) );
                exporter->setPostDir( "./" ); // This is a test to see if M_post_dir is working
                exporter->setMeshProcId( meshPart.meshPartition(), d->comm->MyPID() );
            }
            else
#endif
            {
                if (exporterType.compare("none") == 0)
                {
                    exporter.reset( new ExporterEmpty<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.meshPartition(), "ethiersteinman", d->comm->MyPID()) );
                }
                else
                {
                    exporter.reset( new ExporterEnsight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.meshPartition(), "ethiersteinman", d->comm->MyPID()) );
                }
            }

            velAndPressure.reset( new vector_Type(*fluid.solution(), exporter->mapType() ) );
            if(M_exportExactSolutions)
            {
                exactPressPtr.reset( new vector_Type(exactPress, exporter->mapType() ) );
                pFESpace.interpolate(Problem::pexact, *exactPressPtr, 0);
                exactVelPtr.reset( new vector_Type(exactVel, exporter->mapType() ) );
                uFESpace.interpolate(Problem::uexact, *exactVelPtr, 0);
            }

            exporter->addVariable( ExporterData::Vector, "velocity", velAndPressure,
                                   UInt(0), uFESpace.dof().numTotalDof() );

            exporter->addVariable( ExporterData::Scalar, "pressure", velAndPressure,
                                   UInt(3*uFESpace.dof().numTotalDof()),
                                   UInt(pFESpace.dof().numTotalDof()) );
            if(M_exportExactSolutions)
            {
                exporter->addVariable( ExporterData::Scalar, "exactPressure", exactPressPtr,
                                       UInt(0), UInt(pFESpace.dof().numTotalDof()) );
                exporter->addVariable( ExporterData::Vector, "exactVelocity", exactVelPtr,
                                       UInt(0), UInt(uFESpace.dof().numTotalDof()) );
            }
            exporter->postProcess( 0 );

            initChrono.stop();
            if (verbose) std::cout << "Initialization time: " << initChrono.diff() << " s." << std::endl;

            // +-----------------------------------------------+
            // |             Solving the problem               |
            // +-----------------------------------------------+
            if (verbose) std::cout<< std::endl << "[Solving the problem]" << std::endl;
            int iter = 1;

            for ( ; time <= tFinal + dt/2.; time += dt, iter++)
            {
                iterChrono.reset();
                iterChrono.start();

                oseenData->dataTime()->setTime(time);

                if (verbose) std::cout << "[t = "<< oseenData->dataTime()->time() << " s.]" << std::endl;

                double alpha = bdf.bdfVelocity().coefficientFirstDerivative( 0 ) / oseenData->dataTime()->timeStep();

                beta = bdf.bdfVelocity().extrapolation(); // Extrapolation for the convective term
                //beta *= 0;
                //uFESpace.interpolate(Problem::uexact, beta, time);
                bdf.bdfVelocity().updateRHSContribution( oseenData->dataTime()->timeStep());
                rhs  = fluid.matrixMass()*bdf.bdfVelocity().rhsContributionFirstDerivative();

                fluid.getDisplayer().leaderPrint("alpha ", alpha);
                fluid.getDisplayer().leaderPrint("\n");
                fluid.getDisplayer().leaderPrint("norm beta ", beta.norm2());
                fluid.getDisplayer().leaderPrint("\n");
                fluid.getDisplayer().leaderPrint("norm rhs  ", rhs.norm2());
                fluid.getDisplayer().leaderPrint("\n");

                fluid.updateSystem( alpha, beta, rhs );
                fluid.iterate( bcH );

                bdf.bdfVelocity().shiftRight( *fluid.solution() );

                // Computation of the error
                LifeV::Real urelerr, prelerr, ul2error, pl2error;

                computeErrors(*fluid.solution(),
                              ul2error, urelerr, uFESpace,
                              pl2error, prelerr, pFESpace,
                              time);

                if (verbose && M_exportNorms)
                {
                    out_norm << time  << " "
                    << ul2error << " "
                    << urelerr << " "
                    << pl2error << " "
                    << prelerr << "\n" << std::flush;
                }

                // Saving the errors for the final test
                uL2Error[iElem][jDiscretization] = ul2error;
                pL2Error[iElem][jDiscretization] = pl2error;

                // Exporting the solution
                *velAndPressure = *fluid.solution();
                if(M_exportExactSolutions)
                {
                    pFESpace.interpolate(Problem::pexact, *exactPressPtr, time);
                    uFESpace.interpolate(Problem::uexact, *exactVelPtr, time);
                }
                exporter->postProcess( time );


                MPI_Barrier(MPI_COMM_WORLD);

                iterChrono.stop();
                if (verbose) std::cout << "Iteration time: " << initChrono.diff() << " s." << std::endl << std::endl;
            }

            if (verbose && M_exportNorms)
            {
                out_norm.close();
            }

            // ** BEGIN Accuracy test **
            if(M_test == Accuracy)
            {
                // Computation of the error
                LifeV::Real urelerr, prelerr, ul2error, pl2error;

                computeErrors(*fluid.solution(),
                              ul2error, urelerr, uFESpace,
                              pl2error, prelerr, pFESpace,
                              time);

                if (verbose) std::cout << "Relative error: E(u)=" << urelerr << ", E(p)=" << prelerr << std::endl
                                       << "Tolerance=" << M_accuracyTol << std::endl;

                if (urelerr>M_accuracyTol || prelerr>M_accuracyTol)
                {
                    if (verbose) std::cout << "TEST_ROSSETHIERSTEINMAN STATUS: ECHEC" << std::endl;
                    throw Ethiersteinman::RESULT_CHANGED_EXCEPTION();
                }
            }
            // ** END Accuracy test **

            runChrono.stop();
            if (verbose) std::cout << "Total run time: " << runChrono.diff() << " s." << std::endl;
            if (verbose) std::cout << "[[END_RUN_" << jDiscretization*FEnumber+iElem << "]]" << std::endl;

        } // End of loop on the finite elements
    } // End of loop on the mesh refinement

    // ** BEGIN Space convergence test **
    if (verbose && (M_test == SpaceConvergence))
    {
        bool success;
        success = checkConvergenceRate(uFE, uL2Error, uConvergenceOrder,
                                       pFE, pL2Error, pConvergenceOrder,
                                       meshDiscretization,
                                       M_convTol);

        if (!success)
        {
            if (verbose) std::cout << "TEST_ROSSETHIERSTEINMAN STATUS: ECHEC" << std::endl;
            throw Ethiersteinman::RESULT_CHANGED_EXCEPTION();
        }
    }
    // ** END Space convergence test **
    globalChrono.stop();
    if (verbose) std::cout << std::endl << "Total simulation time:" << globalChrono.diff() << " s." << std::endl;
    if (verbose && (M_test != None)) std::cout << "TEST_ROSSETHIERSTEINMAN STATUS: SUCCESS" << std::endl;
    if (verbose) std::cout << std::endl << "[[END_SIMULATION]]" << std::endl;
}


//////////////////////


