/* -*- mode: c++ -*-

  This file is part of the LifeV library.
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
   \file main.cpp
   \author Davide Forti <davide.forti@epfl.ch>
   \date 2014-02-06
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/navier_stokes_blocks/solver/NavierStokesSolver.hpp>
#include <lifev/core/fem/TimeAndExtrapolationHandler.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "boundaryConditions.hpp"

using namespace LifeV;

int
main ( int argc, char** argv )
{
    bool verbose (false);
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( Comm->MyPID() == 0 )
    {
        verbose = true;
    }
#else
    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm () );
    verbose = true;
#endif

    {

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    // Reading the dataFile
    const std::string defaultDataName = "data";
    GetPot command_line (argc, argv);
    std::string data_file_name = command_line.follow (defaultDataName.c_str(), 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    // reading the mesh
    boost::shared_ptr<mesh_Type > fullMeshPtr ( new mesh_Type ( Comm ) );
    MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");
    readMesh (*fullMeshPtr, meshData);
    int numElementsTotal = fullMeshPtr->numElements();

    // mesh partitioner
    MeshPartitioner< mesh_Type >  meshPart (fullMeshPtr, Comm);
    boost::shared_ptr<mesh_Type > localMeshPtr ( new mesh_Type ( Comm ) );
    localMeshPtr = meshPart.meshPartition();
    fullMeshPtr.reset();

    // create the solver
    NavierStokesSolver ns( dataFile, Comm);
    ns.setup(localMeshPtr);
    ns.setParameters();
    ns.buildSystem();

    Real saveAfter = dataFile("fluid/save_after", 0.0);

    bool useStabilization = dataFile("fluid/stabilization/use", false);
    std::string stabilizationType = dataFile("fluid/stabilization/type", "none");
        
    int saveEvery = dataFile ( "fluid/save_every", 1 );
    int counterSaveEvery = saveEvery;

    // Time handler objects to deal with time advancing and extrapolation
    TimeAndExtrapolationHandler timeVelocity;
    Real dt       = dataFile("fluid/time_discretization/timestep",0.0);
    Real t0       = dataFile("fluid/time_discretization/initialtime",0.0);
    Real tFinal   = dataFile("fluid/time_discretization/endtime",0.0);
    UInt orderBDF = dataFile("fluid/time_discretization/BDF_order",2);

    // Order of BDF and extrapolation for the velocity
    timeVelocity.setBDForder(orderBDF);
    timeVelocity.setMaximumExtrapolationOrder(orderBDF);
    timeVelocity.setTimeStep(dt);

    // Initialize time advance
    vector_Type velocityInitial ( ns.uFESpace()->map() );
    std::vector<vector_Type> initialStateVelocity;
    velocityInitial *= 0 ;
    for ( UInt i = 0; i < orderBDF; ++i )
    	initialStateVelocity.push_back(velocityInitial);

    timeVelocity.initialize(initialStateVelocity);
        
    TimeAndExtrapolationHandler timePressure;
    if ( useStabilization && stabilizationType.compare("VMSLES_SEMI_IMPLICIT")==0 )
    {
        timePressure.setBDForder(orderBDF);
        timePressure.setMaximumExtrapolationOrder(orderBDF);
        timePressure.setTimeStep(dt);
        
        vector_Type pressureInitial ( ns.pFESpace()->map() );
        std::vector<vector_Type> initialStatePressure;
        pressureInitial.zero();
        for ( UInt i = 0; i < orderBDF; ++i )
            initialStatePressure.push_back(pressureInitial);
        
        timePressure.initialize(initialStatePressure);
    }
        
    // Exporter
    std::string outputName = dataFile ( "exporter/filename", "result");
    boost::shared_ptr< ExporterHDF5<mesh_Type > > exporter;

    exporter.reset ( new ExporterHDF5<mesh_Type > ( dataFile, outputName ) );
    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( localMeshPtr, Comm->MyPID() );
    vectorPtr_Type velocity( new vector_Type(ns.uFESpace()->map(), exporter->mapType() ) );
    vectorPtr_Type pressure( new vector_Type(ns.pFESpace()->map(), exporter->mapType() ) );
    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", ns.uFESpace(), velocity, UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", ns.pFESpace(), pressure, UInt (0) );

    if ( dataFile ( "fluid/stabilization/ode_fine_scale", false) )
        	ns.setExportFineScaleVelocity(*exporter, numElementsTotal);

    // exporter->postProcess ( t0 );

    // Boundary conditions
    boost::shared_ptr<BCHandler> bc ( new BCHandler (*BCh_fluid ()) );
    boost::shared_ptr<BCHandler> bc_drag ( new BCHandler (*BCh_drag ()) );
    boost::shared_ptr<BCHandler> bc_lift ( new BCHandler (*BCh_lift ()) );

    std::string preconditioner = dataFile("fluid/preconditionerType","none");

    // Time loop
    LifeChrono iterChrono;
    Real time = t0 + dt;

    vectorPtr_Type u_star( new vector_Type(ns.uFESpace()->map(), Unique ) );
    vectorPtr_Type p_star( new vector_Type(ns.pFESpace()->map(), Unique ) );
    vectorPtr_Type rhs_velocity( new vector_Type(ns.uFESpace()->map(), Unique ) );
    
    if ( useStabilization && stabilizationType.compare("VMSLES_SEMI_IMPLICIT")==0 )
        vectorPtr_Type p_star( new vector_Type(ns.pFESpace()->map(), Unique ) );
        
    ns.setAlpha(timeVelocity.alpha());
    ns.setTimeStep(dt);

    int time_step_count = 0;

    VectorSmall<2> AerodynamicCoefficients;
    Real S = 1.0*4.0;
    Real factor = 2.0/(1000.0*22.0*22.0*S); // 2/(rho*V^2*S)
    std::ofstream outputFile_coefficients;

    if ( verbose )
    {
    	outputFile_coefficients.open ("Aerodynamic_Coefficients.txt");
    	outputFile_coefficients << "% time / X force / Y force \n" << std::flush;
    }

    for ( ; time <= tFinal + dt / 2.; time += dt)
    {
    	time_step_count += 1;

    	if (verbose)
    		std::cout << "\nWe are at time " << time << " s\n\n";

    	iterChrono.reset();
    	iterChrono.start();

    	u_star->zero();
    	rhs_velocity->zero();
    	timeVelocity.extrapolate (orderBDF, *u_star);
    	timeVelocity.rhsContribution (*rhs_velocity);
        
        if ( useStabilization && stabilizationType.compare("VMSLES_SEMI_IMPLICIT")==0 )
        {
            timePressure.extrapolate (orderBDF,*p_star);
            ns.setExtrapolatedPressure(p_star);
        }
        
        ns.updateSystem ( u_star, rhs_velocity );
    	ns.iterate( bc, time );

    	AerodynamicCoefficients =  ns.computeForces(*bc_drag, *bc_lift);

    	if (verbose)
    	{
    		std::cout << "\nDrag coefficient " << AerodynamicCoefficients[0]*factor << "\n";
    		std::cout << "\nLift coefficient " << AerodynamicCoefficients[1]*factor << "\n";
    		outputFile_coefficients << time  << "  " << AerodynamicCoefficients[0]*factor    << "  " << AerodynamicCoefficients[1]*factor << "\n" << std::flush;
    	}

        ns.updateVelocity(velocity);
        ns.updatePressure(pressure);
        
    	iterChrono.stop();

    	if (verbose)
    		std::cout << "\nTimestep solved in " << iterChrono.diff() << " s\n";
        
    	// This part below handles the exporter of the solution.
    	// In particular, given a number of timesteps at which
    	// we ask to export the solution (from datafile), here
    	// the code takes care of exporting the solution also at
    	// the previous timesteps such that, if later a restart
    	// of the simulation is performed, it works correctly.
    	if ( orderBDF == 1 )
    	{
    		if ( time_step_count == (counterSaveEvery-1) )
    		{
    			exporter->postProcess ( time );
    		}
    		else if ( time_step_count == counterSaveEvery )
    		{
    			exporter->postProcess ( time );
    			counterSaveEvery += saveEvery;
    		}
    	}
    	else if ( orderBDF == 2 )
    	{
//    		if ( time_step_count == (counterSaveEvery-2) )
//    		{
//    			exporter->postProcess ( time );
//    		}
//    		else if ( time_step_count == (counterSaveEvery-1) )
//    		{
//    			exporter->postProcess ( time );
//    		}
    		if ( time_step_count == counterSaveEvery )
    		{
    			if ( time >= saveAfter )
    			{
    				exporter->postProcess ( time );
    			}
    			counterSaveEvery += saveEvery;
    		}
    	}
        
        timeVelocity.shift(*velocity);
        
        if ( useStabilization && stabilizationType.compare("VMSLES_SEMI_IMPLICIT")==0 )
            timePressure.shift(*pressure);
        
    }

    if (verbose )
    {
    	outputFile_coefficients.close();
    }

    exporter->closeFile();

	}

#ifdef HAVE_MPI
    if (verbose)
    {
        std::cout << "\nMPI Finalization" << std::endl;
    }
    MPI_Finalize();
#endif
    return ( EXIT_SUCCESS );
}


