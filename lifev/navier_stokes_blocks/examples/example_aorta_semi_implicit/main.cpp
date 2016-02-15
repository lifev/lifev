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

Real oneFunction (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
	return 1.0;
}

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
    boost::shared_ptr<BCHandler> bc_preprocessing ( new BCHandler (*BCh_preprocessing ()) );

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

    /* Preprocessing to define function phi used to define inflow of the aorta */

    UInt flag = 200;
    vectorPtr_Type laplacianSolution( new vector_Type( ns.uFESpace_scalar()->map() ) );
    laplacianSolution->zero();
    ns.solveLaplacian(flag, bc_preprocessing, laplacianSolution);

    BCFunctionBase one (oneFunction);

    // Inflow

    Real nx_inflow = 0.07780;
    Real ny_inflow = 0.0;
    Real nz_inflow = 0.99696;
    Real Q_hat_inflow = 0.0;
    vectorPtr_Type Phi_h_inflow;
    vectorPtr_Type V_hat_x_inflow;
    vectorPtr_Type V_hat_y_inflow;
    vectorPtr_Type V_hat_z_inflow;
    BCHandler bcH_laplacian_inflow;
    bcH_laplacian_inflow.addBC( "Inflow", 3, Essential, Full, one, 1 );
    bcH_laplacian_inflow.bcUpdate ( *ns.uFESpace_scalar()->mesh(), ns.uFESpace_scalar()->feBd(), ns.uFESpace_scalar()->dof() );

    ns.preprocessBoundary ( nx_inflow, ny_inflow, nz_inflow, bcH_laplacian_inflow, Q_hat_inflow, laplacianSolution, 3,
    		Phi_h_inflow, V_hat_x_inflow, V_hat_y_inflow, V_hat_z_inflow );

    if (verbose)
    	std::cout << "\tValue of inflow, Q_hat = " << Q_hat_inflow << std::endl;

    // Outflow 4
    Real nx_flag4 = 0.0;
    Real ny_flag4 = 0.206;
    Real nz_flag4 = 0.978;
    Real Q_hat_flag4 = 0.0;
    vectorPtr_Type Phi_h_outflow4;
    vectorPtr_Type V_hat_x_flag4;
    vectorPtr_Type V_hat_y_flag4;
    vectorPtr_Type V_hat_z_flag4;
    BCHandler bcH_laplacian_flag4;
    bcH_laplacian_flag4.addBC( "Outflow4", 4, Essential, Full, one, 1 );
    bcH_laplacian_flag4.bcUpdate ( *ns.uFESpace_scalar()->mesh(), ns.uFESpace_scalar()->feBd(), ns.uFESpace_scalar()->dof() );

    ns.preprocessBoundary ( nx_flag4, ny_flag4, nz_flag4, bcH_laplacian_flag4, Q_hat_flag4, laplacianSolution, 4,
    		Phi_h_outflow4, V_hat_x_flag4, V_hat_y_flag4, V_hat_z_flag4 );

    if (verbose)
    	std::cout << "\tValue of outflow 4, Q_hat = " << Q_hat_flag4 << std::endl;

    // Outflow 5
    Real nx_flag5 = 0.0;
    Real ny_flag5 = 0.0;
    Real nz_flag5 = 1.0;
    Real Q_hat_flag5 = 0.0;
    vectorPtr_Type Phi_h_outflow5;
    vectorPtr_Type V_hat_x_flag5;
    vectorPtr_Type V_hat_y_flag5;
    vectorPtr_Type V_hat_z_flag5;
    BCHandler bcH_laplacian_flag5;
    bcH_laplacian_flag5.addBC( "Outflow5", 5, Essential, Full, one, 1 );
    bcH_laplacian_flag5.bcUpdate ( *ns.uFESpace_scalar()->mesh(), ns.uFESpace_scalar()->feBd(), ns.uFESpace_scalar()->dof() );

    ns.preprocessBoundary ( nx_flag5, ny_flag5, nz_flag5, bcH_laplacian_flag5, Q_hat_flag5, laplacianSolution, 5,
    		Phi_h_outflow5, V_hat_x_flag5, V_hat_y_flag5, V_hat_z_flag5 );

    if (verbose)
    	std::cout << "\tValue of outflow 5, Q_hat = " << Q_hat_flag5 << std::endl;

    // Outflow 6
    Real nx_flag6 = 0.0;
    Real ny_flag6 = 0.51;
    Real nz_flag6 = 0.86;
    Real Q_hat_flag6 = 0.0;
    vectorPtr_Type Phi_h_outflow6;
    vectorPtr_Type V_hat_x_flag6;
    vectorPtr_Type V_hat_y_flag6;
    vectorPtr_Type V_hat_z_flag6;
    BCHandler bcH_laplacian_flag6;
    bcH_laplacian_flag6.addBC( "Outflow6", 6, Essential, Full, one, 1 );
    bcH_laplacian_flag6.bcUpdate ( *ns.uFESpace_scalar()->mesh(), ns.uFESpace_scalar()->feBd(), ns.uFESpace_scalar()->dof() );

    ns.preprocessBoundary ( nx_flag6, ny_flag6, nz_flag6, bcH_laplacian_flag6, Q_hat_flag6, laplacianSolution, 6,
    		Phi_h_outflow6, V_hat_x_flag6, V_hat_y_flag6, V_hat_z_flag6 );

    if (verbose)
    	std::cout << "\tValue of outflow 6, Q_hat = " << Q_hat_flag6 << std::endl;

    // Outflow 7
    Real nx_flag7 = 0.0;
    Real ny_flag7 = -0.807;
    Real nz_flag7 = -0.589;
    Real Q_hat_flag7 = 0.0;
    vectorPtr_Type Phi_h_outflow7;
    vectorPtr_Type V_hat_x_flag7;
    vectorPtr_Type V_hat_y_flag7;
    vectorPtr_Type V_hat_z_flag7;
    BCHandler bcH_laplacian_flag7;
    bcH_laplacian_flag7.addBC( "Outflow7", 7, Essential, Full, one, 1 );
    bcH_laplacian_flag7.bcUpdate ( *ns.uFESpace_scalar()->mesh(), ns.uFESpace_scalar()->feBd(), ns.uFESpace_scalar()->dof() );

    ns.preprocessBoundary ( nx_flag7, ny_flag7, nz_flag7, bcH_laplacian_flag7, Q_hat_flag7, laplacianSolution, 7,
    		Phi_h_outflow7, V_hat_x_flag7, V_hat_y_flag7, V_hat_z_flag7 );

    if (verbose)
    	std::cout << "\tValue of outflow 7, Q_hat = " << Q_hat_flag7 << std::endl;

    // Outflow 8
    Real nx_flag8 = 0.0;
    Real ny_flag8 = -0.185;
    Real nz_flag8 = 0.983;
    Real Q_hat_flag8 = 0.0;
    vectorPtr_Type Phi_h_outflow8;
    vectorPtr_Type V_hat_x_flag8;
    vectorPtr_Type V_hat_y_flag8;
    vectorPtr_Type V_hat_z_flag8;
    BCHandler bcH_laplacian_flag8;
    bcH_laplacian_flag8.addBC( "Outflow8", 8, Essential, Full, one, 1 );
    bcH_laplacian_flag8.bcUpdate ( *ns.uFESpace_scalar()->mesh(), ns.uFESpace_scalar()->feBd(), ns.uFESpace_scalar()->dof() );

    ns.preprocessBoundary ( nx_flag8, ny_flag8, nz_flag8, bcH_laplacian_flag8, Q_hat_flag8, laplacianSolution, 8,
    		Phi_h_outflow8, V_hat_x_flag8, V_hat_y_flag8, V_hat_z_flag8);

    if (verbose)
    	std::cout << "\tValue of outflow 8, Q_hat = " << Q_hat_flag8 << std::endl;

    // Outflow 9
    Real nx_flag9 = 0.0;
    Real ny_flag9 = 0.851;
    Real nz_flag9 = -0.525;
    Real Q_hat_flag9 = 0.0;
    vectorPtr_Type Phi_h_outflow9;
    vectorPtr_Type V_hat_x_flag9;
    vectorPtr_Type V_hat_y_flag9;
    vectorPtr_Type V_hat_z_flag9;
    BCHandler bcH_laplacian_flag9;
    bcH_laplacian_flag9.addBC( "Outflow9", 9, Essential, Full, one, 1 );
    bcH_laplacian_flag9.bcUpdate ( *ns.uFESpace_scalar()->mesh(), ns.uFESpace_scalar()->feBd(), ns.uFESpace_scalar()->dof() );

    ns.preprocessBoundary ( nx_flag9, ny_flag9, nz_flag9, bcH_laplacian_flag9, Q_hat_flag9, laplacianSolution, 9,
    		Phi_h_outflow9, V_hat_x_flag9, V_hat_y_flag9, V_hat_z_flag9);

    /* end preprocessing */

    /* Vectors needed to impose flowrates */

    vectorPtr_Type velAndPressure_flowrates;
    velAndPressure_flowrates.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    // Inflow
    vectorPtr_Type velAndPressure_inflow_reference;
    velAndPressure_inflow_reference.reset ( new vector_Type (ns.uFESpace()->map(), Unique ) );
    velAndPressure_inflow_reference->zero();
    velAndPressure_inflow_reference->subset(*V_hat_x_inflow,ns.uFESpace_scalar()->map(), 0, 0);
    velAndPressure_inflow_reference->subset(*V_hat_y_inflow,ns.uFESpace_scalar()->map(), 0, ns.uFESpace()->map().mapSize()/3 );
    velAndPressure_inflow_reference->subset(*V_hat_z_inflow,ns.uFESpace_scalar()->map(), 0, 2*ns.uFESpace()->map().mapSize()/3 );

    // Flag 4
    vectorPtr_Type velAndPressure_outflow4_reference;
    velAndPressure_outflow4_reference.reset ( new vector_Type (ns.uFESpace()->map(), Unique ) );
    velAndPressure_outflow4_reference->zero();
    velAndPressure_outflow4_reference->subset(*V_hat_x_flag4,ns.uFESpace_scalar()->map(), 0, 0);
    velAndPressure_outflow4_reference->subset(*V_hat_y_flag4,ns.uFESpace_scalar()->map(), 0, ns.uFESpace()->map().mapSize()/3 );
    velAndPressure_outflow4_reference->subset(*V_hat_z_flag4,ns.uFESpace_scalar()->map(), 0, 2*ns.uFESpace()->map().mapSize()/3 );

    // Flag 5
    vectorPtr_Type velAndPressure_outflow5_reference;
    velAndPressure_outflow5_reference.reset ( new vector_Type (ns.uFESpace()->map(), Unique ) );
    velAndPressure_outflow5_reference->zero();
    velAndPressure_outflow5_reference->subset(*V_hat_x_flag5,ns.uFESpace_scalar()->map(), 0, 0);
    velAndPressure_outflow5_reference->subset(*V_hat_y_flag5,ns.uFESpace_scalar()->map(), 0, ns.uFESpace()->map().mapSize()/3 );
    velAndPressure_outflow5_reference->subset(*V_hat_z_flag5,ns.uFESpace_scalar()->map(), 0, 2*ns.uFESpace()->map().mapSize()/3 );

    // Flag 6
    vectorPtr_Type velAndPressure_outflow6_reference;
    velAndPressure_outflow6_reference.reset ( new vector_Type (ns.uFESpace()->map(), Unique ) );
    velAndPressure_outflow6_reference->zero();
    velAndPressure_outflow6_reference->subset(*V_hat_x_flag6,ns.uFESpace_scalar()->map(), 0, 0);
    velAndPressure_outflow6_reference->subset(*V_hat_y_flag6,ns.uFESpace_scalar()->map(), 0, ns.uFESpace()->map().mapSize()/3 );
    velAndPressure_outflow6_reference->subset(*V_hat_z_flag6,ns.uFESpace_scalar()->map(), 0, 2*ns.uFESpace()->map().mapSize()/3 );

    // Flag 7
    vectorPtr_Type velAndPressure_outflow7_reference;
    velAndPressure_outflow7_reference.reset ( new vector_Type (ns.uFESpace()->map(), Unique ) );
    velAndPressure_outflow7_reference->zero();
    velAndPressure_outflow7_reference->subset(*V_hat_x_flag7,ns.uFESpace_scalar()->map(), 0, 0);
    velAndPressure_outflow7_reference->subset(*V_hat_y_flag7,ns.uFESpace_scalar()->map(), 0, ns.uFESpace()->map().mapSize()/3 );
    velAndPressure_outflow7_reference->subset(*V_hat_z_flag7,ns.uFESpace_scalar()->map(), 0, 2*ns.uFESpace()->map().mapSize()/3 );

    // Flag 8
    vectorPtr_Type velAndPressure_outflow8_reference;
    velAndPressure_outflow8_reference.reset ( new vector_Type (ns.uFESpace()->map(), Unique ) );
    velAndPressure_outflow8_reference->zero();
    velAndPressure_outflow8_reference->subset(*V_hat_x_flag8,ns.uFESpace_scalar()->map(), 0, 0);
    velAndPressure_outflow8_reference->subset(*V_hat_y_flag8,ns.uFESpace_scalar()->map(), 0, ns.uFESpace()->map().mapSize()/3 );
    velAndPressure_outflow8_reference->subset(*V_hat_z_flag8,ns.uFESpace_scalar()->map(), 0, 2*ns.uFESpace()->map().mapSize()/3 );

    // Flag 9
    vectorPtr_Type velAndPressure_outflow9_reference;
    velAndPressure_outflow9_reference.reset ( new vector_Type (ns.uFESpace()->map(), Unique ) );
    velAndPressure_outflow9_reference->zero();
    velAndPressure_outflow9_reference->subset(*V_hat_x_flag9,ns.uFESpace_scalar()->map(), 0, 0);
    velAndPressure_outflow9_reference->subset(*V_hat_y_flag9,ns.uFESpace_scalar()->map(), 0, ns.uFESpace()->map().mapSize()/3 );
    velAndPressure_outflow9_reference->subset(*V_hat_z_flag9,ns.uFESpace_scalar()->map(), 0, 2*ns.uFESpace()->map().mapSize()/3 );

    vectorPtr_Type velAndPressure_inflow;
    velAndPressure_inflow.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    vectorPtr_Type velAndPressure_outflow4;
    velAndPressure_outflow4.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    vectorPtr_Type velAndPressure_outflow5;
    velAndPressure_outflow5.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    vectorPtr_Type velAndPressure_outflow6;
    velAndPressure_outflow6.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    vectorPtr_Type velAndPressure_outflow7;
    velAndPressure_outflow7.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    vectorPtr_Type velAndPressure_outflow8;
    velAndPressure_outflow8.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    vectorPtr_Type velAndPressure_outflow9;
    velAndPressure_outflow9.reset ( new vector_Type ( ns.uFESpace()->map(), Unique ) );

    /* end preprocessing */

    Real i_HeartBeat = 0.0;

    int time_step_count = 0;

    for ( ; time <= tFinal + dt / 2.; time += dt)
    {
    	time_step_count += 1;

    	// ---------------------------------
    	// Evaluation of the inflow velocity
    	// ---------------------------------

    	Real Q = 0;
    	Real Q_inflow = 0;
    	Real Q_flag4  = 0;
    	Real Q_flag5  = 0;
    	Real Q_flag6  = 0;
    	Real Q_flag7  = 0;
    	Real Q_flag8  = 0;
    	Real Q_flag9  = 0;

    	Real T_heartbeat = 0.8;

    	if ( time < T_heartbeat )
    	{
    		i_HeartBeat = 0.0;
    	}
    	else if ( time >= T_heartbeat && time < 2*T_heartbeat )
    	{
    		i_HeartBeat = 1.0;
    	}
    	else if ( time >= 2*T_heartbeat && time < 3*T_heartbeat )
    	{
    		i_HeartBeat = 2.0;
    	}
    	else if ( time >= 3*T_heartbeat && time < 4*T_heartbeat )
    	{
    		i_HeartBeat = 3.0;
    	}

    	if ( (time >= 0.05 && time <= 0.42) || (time >= (0.05+T_heartbeat) && time <= (0.42+T_heartbeat) ) || (time >= (0.05+2*T_heartbeat) && time <= (0.42+2*T_heartbeat) ) || (time >= (0.05+3*T_heartbeat) && time <= (0.42+3*T_heartbeat) ) )
    	{
    		// nuova
    		Q = -2.314569820334801e+09*std::pow(time-i_HeartBeat*T_heartbeat,9) +
    				4.952537061974133e+09*std::pow(time-i_HeartBeat*T_heartbeat,8) -
    				4.532060231242586e+09*std::pow(time-i_HeartBeat*T_heartbeat,7) +
    				2.325743716202249e+09*std::pow(time-i_HeartBeat*T_heartbeat,6) -
    				7.387577876374097e+08*std::pow(time-i_HeartBeat*T_heartbeat,5) +
    				1.514516710083440e+08*std::pow(time-i_HeartBeat*T_heartbeat,4) -
    				2.018053394181958e+07*std::pow(time-i_HeartBeat*T_heartbeat,3) +
    				1.667954643625200e+06*std::pow(time-i_HeartBeat*T_heartbeat,2) -
    				7.160662399848596e+04*(time-i_HeartBeat*T_heartbeat) +
    				1.184312187078482e+03;
    		Q = Q/394;

    		// usata prima
    		//Q = 2.117637666632775e+04*std::pow(time-i_HeartBeat*T_heartbeat,6)-3.370930726888496e+04*std::pow(time-i_HeartBeat*T_heartbeat,5)+2.133377678002176e+04*std::pow(time-i_HeartBeat*T_heartbeat,4)-6.666366536069445e+03*std::pow(time-i_HeartBeat*T_heartbeat,3)+1.011772959679957e+03*std::pow(time-i_HeartBeat*T_heartbeat,2)-6.023975547926423e+01*(time-i_HeartBeat*T_heartbeat)+1.192718364532979e+00;
    	}
    	else
    	{
    		Q = 0.0;
    	}

    	Q_inflow = 394*Q;
    	Q_flag4  = 20.25*Q; // left_common_carotid
    	Q_flag5  = 21.82*Q; // right_common_carotid
    	Q_flag6  = 1.43*Q; // right_vertebral
    	Q_flag7  = 24.77*Q; // right_subclavian
    	Q_flag8  = 4.69*Q; // left_vertebral
    	Q_flag9  = 21.54*Q; // left_subclavian

    	Real pressureValue = 1500.0/2.51*(Q_inflow - Q_flag4 - Q_flag5 - Q_flag6 - Q_flag7 - Q_flag8 - Q_flag9);

    	if (verbose)
    	{
    		std::cout << "\nImposing outflow pressure of " << pressureValue << " dyne/cm^2\n\n";
    	}

    	Real alpha_flowrate_inflow = Q_inflow/Q_hat_inflow;
    	Real alpha_flowrate_flag4  = Q_flag4/Q_hat_flag4;
    	Real alpha_flowrate_flag5  = Q_flag5/Q_hat_flag5;
    	Real alpha_flowrate_flag6  = Q_flag6/Q_hat_flag6;
    	Real alpha_flowrate_flag7  = Q_flag7/Q_hat_flag7;
    	Real alpha_flowrate_flag8  = Q_flag8/Q_hat_flag8;
    	Real alpha_flowrate_flag9  = Q_flag9/Q_hat_flag9;

    	if (verbose)
    	{
    	    		std::cout << "Q_inflow: " << Q_inflow << std::endl << std::endl;
    	    		std::cout << "Q_flag4: "  << Q_flag4  << std::endl << std::endl;
    	    		std::cout << "Q_flag5: "  << Q_flag5  << std::endl << std::endl;
    	    		std::cout << "Q_flag6: "  << Q_flag6  << std::endl << std::endl;
    	    		std::cout << "Q_flag7: "  << Q_flag7  << std::endl << std::endl;
    	    		std::cout << "Q_flag8: "  << Q_flag8  << std::endl << std::endl;
    	    		std::cout << "Q_flag9: "  << Q_flag9  << std::endl << std::endl;
    	    		std::cout << "Q_outflow: "  << Q_inflow - Q_flag4 - Q_flag5 - Q_flag6 - Q_flag7 - Q_flag8 - Q_flag9  << std::endl << std::endl;
    	}

    	velAndPressure_inflow->zero();
    	*velAndPressure_inflow += *velAndPressure_inflow_reference;
    	*velAndPressure_inflow *= alpha_flowrate_inflow;

    	velAndPressure_outflow4->zero();
    	*velAndPressure_outflow4 += *velAndPressure_outflow4_reference;
    	*velAndPressure_outflow4 *= alpha_flowrate_flag4;

    	velAndPressure_outflow5->zero();
    	*velAndPressure_outflow5 += *velAndPressure_outflow5_reference;
    	*velAndPressure_outflow5 *= alpha_flowrate_flag5;

    	velAndPressure_outflow6->zero();
    	*velAndPressure_outflow6 += *velAndPressure_outflow6_reference;
    	*velAndPressure_outflow6 *= alpha_flowrate_flag6;

    	velAndPressure_outflow7->zero();
    	*velAndPressure_outflow7 += *velAndPressure_outflow7_reference;
    	*velAndPressure_outflow7 *= alpha_flowrate_flag7;

    	velAndPressure_outflow8->zero();
    	*velAndPressure_outflow8 += *velAndPressure_outflow8_reference;
    	*velAndPressure_outflow8 *= alpha_flowrate_flag8;

    	velAndPressure_outflow9->zero();
    	*velAndPressure_outflow9 += *velAndPressure_outflow9_reference;
    	*velAndPressure_outflow9 *= alpha_flowrate_flag9;

    	velAndPressure_flowrates->zero();
    	*velAndPressure_flowrates += *velAndPressure_inflow;
    	*velAndPressure_flowrates += *velAndPressure_outflow4;
    	*velAndPressure_flowrates += *velAndPressure_outflow5;
    	*velAndPressure_flowrates += *velAndPressure_outflow6;
    	*velAndPressure_flowrates += *velAndPressure_outflow7;
    	*velAndPressure_flowrates += *velAndPressure_outflow8;
    	*velAndPressure_flowrates += *velAndPressure_outflow9;

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
    	ns.iterate( bc, time, velAndPressure_flowrates );

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


