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
    @file multirate.cpp
   @author L. Oldani <luca.oldani@mail.polimi.it>
   @date 2013-12-27
*/

/*!
    Simple 2D multirate test.
*/

// ===================================================
//! Includes
// ===================================================

#include "multirate.hpp"
#include "user_fun.hpp"

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

// ===================================================
//!              Standard functions
// ===================================================

// Standard functions
Real UOne( const Real& /* t */,
           const Real& /* x */,
           const Real& /* y */,
           const Real& /* z */,
           const ID&   /* icomp */)
{
    return 1.;
}

// ===================================================
//! Private Members
// ===================================================

struct multirate_2d::Private
{
    Private() {}

    // Policy for scalar functions
    typedef boost::function<Real ( const Real&, const Real&,
                                   const Real&, const Real&, const ID& )>
    fct_type;

    // Policy for the flux function
    typedef boost::function<Vector ( const Real&, const Real&,
                                     const Real&, const Real&,
                                     const std::vector<Real>& )>
    vectorFct_type;

    std::string    data_file_name;
    std::string    discretization_section;

    boost::shared_ptr<Epetra_Comm>   comm;

    // Function Types

    fct_type getUOne()
    {
        fct_type f;
        f = boost::bind( &UOne, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getAnalyticalSolution()
    {
        fct_type f;
        f = boost::bind( &DataProblem::analyticalSolution, _1, _2, _3, _4, _5 );
        return f;
    }

};

// ===================================================
//! Constructors
// ===================================================

multirate_2d::multirate_2d ( int argc, char* argv[] )
        : Members ( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "multirate";

#ifdef EPETRA_MPI
    Members->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    Int ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    Members->comm.reset( new Epetra_SerialComm() );
#endif

}

// ===================================================
//! Methods
// ===================================================

Real
multirate_2d::run()
{
    using namespace DataProblem;

    // Life chonos
    LifeChrono chronoTotal;
    LifeChrono chronoReadAndPartitionMesh;
    LifeChrono chronoBoundaryCondition;
    LifeChrono chronoFiniteElementSpace;
    LifeChrono chronoProblem;
    LifeChrono chronoProcess;
    LifeChrono chronoTimeStep;
    LifeChrono chronoError;

    // Start chronoTotal for measure the total time for the computation.
    chronoTotal.start();

    // Displayer to print on screen.
    boost::shared_ptr < Displayer > displayer ( new Displayer ( Members->comm ) );

    // Reading from data file.
    GetPot dataFile( Members->data_file_name.c_str() );

    // The multirate Solver.
    displayer->leaderPrint ( "The multirate solver" );

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes.
    chronoReadAndPartitionMesh.start();

    // Create the data file. 
    multirateDataPtr_Type multirateData ( new multirateData_Type );

    // Set up the data.
    multirateData->setup ( dataFile );

    // Create the mesh file handler.
    MeshData meshData;

    // Set up the mesh file.
    meshData.setup ( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the mesh.
    regionMeshPtr_Type fullMeshPtr ( new regionMesh_Type ( Members->comm ) );

    // Select if the mesh is structured or not.
    if ( meshData.meshType() != "structured" )
    {
        // Set up the mesh.
        // This method is not implemented for quadrilateral meshes.
        // readMesh ( *fullMeshPtr, meshData );
        ERROR_MSG ( "Filter for quadrangular mesh not implemented" );
    }
    else
    {
        // Section of the structured mesh.
        const std::string structuredSection = Members->discretization_section + "/space_discretization/";

        // Set up the structured mesh.
        regularMesh2D ( *fullMeshPtr, 0,
                        dataFile( ( structuredSection + "nx" ).data(), 4 ),
                        dataFile( ( structuredSection + "ny" ).data(), 4 ),
                        dataFile( ( structuredSection + "verbose" ).data(), false ),
                        dataFile( ( structuredSection + "lx" ).data(), 1. ),
                        dataFile( ( structuredSection + "ly" ).data(), 1. ) );
    }

    // Local mesh part.
    regionMeshPtr_Type meshPtr ( new regionMesh_Type ( Members->comm ) );

    // Create the mesh partitioned.
    meshPartitionerPtr_Type meshPart ( new meshPartitioner_Type ); ////////////////////

    // Do the partition using ParMetis.
    meshPart->doPartition ( fullMeshPtr, Members->comm );////////////////////

    // Assign the local part of the mesh.
    meshPtr = meshPart->meshPartition();////////////////////

    // Clear the partitioner, now useless.
    meshPart.reset();////////////////////

    // Stop chronoReadAndPartitionMesh.
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh.
    displayer->leaderPrint ( "Time for read and partition the mesh ",
                             chronoReadAndPartitionMesh.diff(), "\n" );

    // Create the boundary conditions.

    // Start chronoBoundaryCondition for measure the total time for create the boundary conditions
    chronoBoundaryCondition.start();

    bcHandlerPtr_Type bcMultirate ( new bcHandler_Type );

    setBoundaryConditions ( bcMultirate );

    // Stop chronoBoundaryCondition
    chronoBoundaryCondition.stop();

    // The leader process print chronoBoundaryCondition.
    displayer->leaderPrint ( "Time for create the boundary conditions handler ",
                             chronoBoundaryCondition.diff(), "\n" );

    // Create the solution spaces.

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces.
    chronoFiniteElementSpace.start();

    // Primal solution parameters
    const ReferenceFE* refFE ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* qR ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* bdQr ( static_cast<QuadratureRule*>(NULL) );

    refFE = &feQuadQ0;
    qR = &quadRuleQuad4pt;
    bdQr = &quadRuleSeg1pt;

    // Interpolate of velocity field parameters.
    const ReferenceFE* velocity_refFE ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* velocity_qR ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* velocity_bdQr ( static_cast<QuadratureRule*>(NULL) );

    velocity_refFE = &feQuadQ0;
    velocity_qR = &quadRuleQuad4pt;
    velocity_bdQr = &quadRuleSeg1pt;

    // Finite element space of the interpolation of dual variable.
    fESpacePtr_Type velocity_FESpacePtr ( new fESpace_Type ( meshPtr, *velocity_refFE, *velocity_qR,
                                                             *velocity_bdQr, 2, Members->comm ) ); ////////////////////

    // Vector for the interpolated dual solution.
    vectorPtr_Type velocityVector ( new vector_Type ( velocity_FESpacePtr->map(), Repeated ) ); ////////////////////

    // Finite element space
    fESpacePtr_Type feSpacePtr( new fESpace_Type ( meshPtr, *refFE, *qR, *bdQr, 1,
                                                  Members->comm ) ); ////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    // Parameters Of TimeStepField
    const ReferenceFE* P0_FE ( static_cast<ReferenceFE*>(NULL) );
    const QuadratureRule* P0_qR ( static_cast<QuadratureRule*>(NULL) );
    const QuadratureRule* P0_bdQr ( static_cast<QuadratureRule*>(NULL) );

    P0_FE = &feQuadQ0;
    P0_qR = &quadRuleQuad4pt;
    P0_bdQr = &quadRuleSeg1pt;

// Finite element space of the interpolation of dual variable.
    fESpacePtr_Type P0_fESpacePtr ( new fESpace_Type ( meshPtr, *P0_FE, *P0_qR,
                                                             *P0_bdQr, 1, Members->comm ) ); ////////////////////

// Vector for the interpolated dual solution.
    vectorPtr_Type P0Vector ( new vector_Type ( P0_fESpacePtr->map(), Repeated ) ); ////////////////////

// Create the ghost map for the subdomain interface elements.
    GhostHandler<regionMesh_Type> P0_ghost( fullMeshPtr, meshPtr, P0_fESpacePtr->mapPtr(), Members->comm ); ////////////////////

*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create the ghost map for the subdomain interface elements.
    GhostHandler<regionMesh_Type> ghost( fullMeshPtr, meshPtr, feSpacePtr->mapPtr(), Members->comm ); ////////////////////

    // Clear the global mesh.
    fullMeshPtr.reset(); ////////////////////

    // Solution field.
    scalarFieldPtr_Type solutionField ( new scalarField_Type ( feSpacePtr ) ); ////////////////////

	// CFL field
	scalarFieldPtr_Type CFLField ( new scalarField_Type ( feSpacePtr ) ); ////////////////////

	// region field
	scalarFieldPtr_Type regionField ( new scalarField_Type ( feSpacePtr ) ); ////////////////////

    // Stop chronoFiniteElementSpace
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace
    displayer->leaderPrint ( "Time for create the finite element spaces ",
                             chronoFiniteElementSpace.diff(), "\n" );

    // Start chronoProblem for measure the total time for create the problem
    chronoProblem.start();

    // Instantiation of the multirateSolver class
	multirateSolver_Type multirateSolver;

    // Stop chronoProblem
    chronoProblem.stop();

    // The leader process print chronoProblem
    displayer->leaderPrint ( "Time for create the problem ", chronoProblem.diff(), "\n" );

    // Process the problem

    // Start chronoProcess for measure the total time for the simulation
    chronoProcess.start();

    // Set the data for the solver.
    multirateSolver.setData ( multirateData );

    // Set the displayer.
    multirateSolver.setDisplayer ( displayer );

    // Set the solution field.
    multirateSolver.setSolutionField ( solutionField, ghost.ghostMapOnElementsP0() );

	//Set the TimeStep field
	multirateSolver.setTimeStepField ( CFLField, ghost.ghostMapOnElementsP0() );

	//Set the region field
	multirateSolver.setRegionField ( regionField, ghost.ghostMapOnElementsP0() );

    // Clear the ghost structure.
    ghost.clean();

    // Set the initial solution.
    scalarFctPtr_Type initialConditionFct ( new initialCondition );
    multirateSolver.setInitialSolution ( initialConditionFct );

    // Set the mass function.
    scalarFctPtr_Type massFct ( new massTerm );
    multirateSolver.setMassTerm ( massFct );

    /*    
    //// Create the physical flux
    physicalFluxPtr_Type physicalFlux ( new HyperbolicUserFlux );

    // Fill the physical flux with the user defined flux and flux prime
    physicalFlux->setup( dataFile, Members->discretization_section + "/flux" );

    // Create the numerical flux
    numericalFluxPtr_Type numericalFlux ( new HyperbolicGodunovFlux );

    // Save the physical flux in the numerical flux
    numericalFlux->setPhysicalFlux ( physicalFlux );

    //// Set the numerical flux usign the physical flux.
    multirateSolver.setNumericalFlux( numericalFlux );
    */
    
    // Set the boudary conditions
    multirateSolver.setBoundaryCondition( bcMultirate );

    // Set the exporter for the solution
    boost::shared_ptr< Exporter< regionMesh_Type > > exporter;

    // Shared pointer used in the exporter for the solution
    vectorPtr_Type exporterSolution;

	// Set the shared pointer for the CFL
	vectorPtr_Type exporterCFL;

	// Set the shared pointer for the region
	vectorPtr_Type exporterRegion;

    // Shared pointer for the exact solution, error
    vectorPtr_Type exact;
    vectorPtr_Type exporterError;

    // Type of the exporter
    const std::string exporterType =  dataFile( "exporter/type", "none");

    // The name of the file
    const std::string exporterFileName = dataFile( "exporter/file_name", "Concentration" );

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare("hdf5") == 0 )
    {
        exporter.reset( new ExporterHDF5< regionMesh_Type > ( dataFile, exporterFileName ) );
    }
    else
#endif
    {
        exporter.reset( new ExporterEmpty< regionMesh_Type > ( dataFile, exporterFileName ) );
    }

    // Set directory where to save the solution
    exporter->setPostDir( dataFile( "exporter/folder", "./" ) );

    // Export the partitioning
    exporter->exportPID( meshPtr, Members->comm );

    // Set the mesh.
    exporter->setMeshProcId( meshPtr, Members->comm->MyPID() );

    // Set the exporter solution
    exporterSolution.reset( new vector_Type ( solutionField->getVector(),
                                              exporter->mapType() ) );

    // Add the solution to the exporter
    exporter->addVariable ( ExporterData<regionMesh_Type>::ScalarField,
                            dataFile ( "exporter/name_field", "Concentration" ),
                            feSpacePtr,
                            exporterSolution,
                            static_cast<UInt>( 0 ),
                            ExporterData<regionMesh_Type>::UnsteadyRegime,
                            ExporterData<regionMesh_Type>::Cell );

		

	// Set the exporter CFL
    	exporterCFL.reset( new vector_Type ( CFLField->getVector(),
                                              exporter->mapType() ) );

    	// Add the CFL to the exporter
    	exporter->addVariable ( ExporterData<regionMesh_Type>::ScalarField,
                            dataFile ( "exporter/name_field", "CFL" ),
                            feSpacePtr,
                            exporterCFL,
                            static_cast<UInt>( 0 ),
                            ExporterData<regionMesh_Type>::UnsteadyRegime,
                            ExporterData<regionMesh_Type>::Cell );

	// Set the exporter Region
    	exporterRegion.reset( new vector_Type ( regionField->getVector(),
                                              exporter->mapType() ) );

    	// Add the Region to the exporter
    	exporter->addVariable ( ExporterData<regionMesh_Type>::ScalarField,
                            dataFile ( "exporter/name_field", "region" ),
                            feSpacePtr,
                            exporterRegion,
                            static_cast<UInt>( 0 ),
                            ExporterData<regionMesh_Type>::UnsteadyRegime,
                            ExporterData<regionMesh_Type>::Cell );


    // set the exact pointer
    exact.reset (new vector_Type (solutionField->getVector(), exporter->mapType() ) );

    // set the error pointer
    exporterError.reset (new vector_Type (solutionField->getVector(), exporter->mapType() ) );

    // Add the exact solution to the exporter
    exporter->addVariable ( ExporterData < regionMesh_Type >::ScalarField,
                            dataFile( "exporter/name_field", "Concentration" ) + std::string("Error"),
                            feSpacePtr,
                            exporterError,
                            static_cast<UInt>( 0 ),
                            ExporterData < regionMesh_Type >::UnsteadyRegime,
                            ExporterData < regionMesh_Type >::Cell );


    // Display the total number of unknowns
    displayer->leaderPrint( "Number of unknowns : ",
                            feSpacePtr->map().map(Unique)->NumGlobalElements(), "\n" );

    // Solve the problem

    // Copy the initial solution to the exporter
    *exporterSolution = solutionField->getVector();

    // Save the initial solution into the exporter
    exporter->postProcess( multirateData->dataTimePtr()->initialTime() );

    // Define the end time
    const Real endTime = multirateData->dataTimePtr()->endTime();

    // Flag for the last time step that does not coincide with the last advance
    bool isLastTimeStep ( false );
	
    // Current tiem step
    Real timeStep;
	Real timeStepNew;

    // A loop for the simulation, it starts from \Delta t and end in N \Delta t = T
    while ( multirateData->dataTimePtr()->time() < endTime && !isLastTimeStep )
    {
	// Start chronoTimeStep for measure the time for the current time step
        chronoTimeStep.start();

	////////////////////////////////////////////
     
	timeStep = multirateSolver.computeTimeStep ();
 	    
	// Set the last time step for the simulation.
        multirateData->dataTimePtr()->setTimeStep( timeStep );

	multirateSolver.setupTimeSteps();
	
	//compute two regions setting m
	multirateSolver.computeRegionAccordingTimeStep(4);
 	        
	if (multirateSolver.isSlow())
	{	

		UInt i=0;
	      	while (i<multirateSolver.getLevel() && !isLastTimeStep)
		{
				 
			// Check if the time step is consistent.
        		if ( multirateData->dataTimePtr()->isLastTimeStep() )
        		{
	    			// Compute the last time step.
            			timeStepNew = multirateData->dataTimePtr()->leftTime();

				// This is the last time step in the simulation
            			isLastTimeStep = true;

				// Set the last time step for the simulation.
        			multirateData->dataTimePtr()->setTimeStep( timeStepNew );
	
				// Set timestep for slow region and fast region				
				multirateSolver.setFastTimeStep(timeStepNew);
				multirateSolver.setSlowTimeStep(timeStepNew+(i)*timeStep);
        		}
			//////////////////////////ACTUNG//////////////////////////
			
			++i;

			// Update time
		        multirateData->dataTimePtr()->updateTime();

			// The leader process prints the temporal data.
		        if ( displayer->isLeader() )
		        {
		            multirateData->dataTimePtr()->showMe();
		        }

		        // solve one step of the multirate problem.
		        multirateSolver.solveOneFastTimeStep();

			if (i==multirateSolver.getLevel() || isLastTimeStep==true)
			{				
				multirateSolver.solveOneSlowTimeStep();
			}		

			// Save the solution
	 
        	        // Copy the solution to the exporter
		        *exporterSolution = solutionField->getVector ();
		
		        // Interpolate the exact solution
		        feSpacePtr->interpolate ( Members->getAnalyticalSolution(), *exact, multirateData->dataTimePtr()->time() );

			// Copy the CFL to the exporter
			*exporterCFL = CFLField->getVector ();

			// Copy the region to the exporter
			*exporterRegion = regionField->getVector ();
	 
 		        // Compute error in exporter
		        *exporterError = *exact - solutionField->getVector();
		        exporterError->abs();

			// Save the solution into the exporter
			exporter->postProcess( multirateData->dataTimePtr()->time() );

			// Stop chronoTimeStep
			chronoTimeStep.stop();

        		// The leader process print chronoTimeStep
			displayer->leaderPrint( "Time for current time step ", chronoTimeStep.diff(), "\n" );
		
		}
	      
	}
	else
	{	
		
		// Check if the time step is consistent.
	        if ( multirateData->dataTimePtr()->isLastTimeStep() )
	        {
			// Compute the last time step.
	            	timeStep = multirateData->dataTimePtr()->leftTime();
	
	            	// This is the last time step in the simulation
	            	isLastTimeStep = true;
			
			// Set the last time step for the simulation.
	        	multirateData->dataTimePtr()->setTimeStep( timeStep );
	        }	
		
		// Update time
	        multirateData->dataTimePtr()->updateTime();

	        // The leader process prints the temporal data.
	        if ( displayer->isLeader() )
	        {
	            multirateData->dataTimePtr()->showMe();
	        }

	        // solve one step of the multirate problem.
	        multirateSolver.solveOneTimeStep();

	        // Save the solution

	        // Copy the solution to the exporter
	        *exporterSolution = solutionField->getVector ();
	
	        // Interpolate the exact solution
	        feSpacePtr->interpolate ( Members->getAnalyticalSolution(), *exact, multirateData->dataTimePtr()->time() );

		// Copy the CFL to the exporter
	        *exporterCFL = CFLField->getVector ();

		// Copy the regionfield to the exporter
	        *exporterRegion = regionField->getVector ();
	
	        // Compute error in exporter
	        *exporterError = *exact - solutionField->getVector();
	        exporterError->abs();

	        // Save the solution into the exporter
	        exporter->postProcess( multirateData->dataTimePtr()->time() );

	        // Stop chronoTimeStep
	        chronoTimeStep.stop();

        	// The leader process print chronoTimeStep
	        displayer->leaderPrint( "Time for current time step ", chronoTimeStep.diff(), "\n" );
	}

    }

    // Stop chronoProcess
    chronoProcess.stop();

    // The leader process print chronoProcess
    displayer->leaderPrint( "Time for process ", chronoProcess.diff(), "\n" );

    // Compute the errors

    // Start chronoError for measure the total time for computing the errors.
    chronoError.start();

    // Compute the error L2 norms
    Real L2Norm(0), exactL2Norm(0), L2Error(0), L2RelativeError(0);

    // Norms and errors for the pressure
    displayer->leaderPrint( "\nERROR\n" );

    // Compute the L2 norm for the solution
    L2Norm = feSpacePtr->l2Norm( solutionField->getVector() );

    // Display the L2 norm for the solution
    displayer->leaderPrint( " L2 norm of solution:            ", L2Norm, "\n" );

    // Compute the L2 norm for the analytical solution
    exactL2Norm = feSpacePtr->l2NormFunction( Members->getAnalyticalSolution(),
                                              multirateData->dataTimePtr()->endTime() );

    // Display the L2 norm for the analytical solution
    displayer->leaderPrint( " L2 norm of exact solution:      ", exactL2Norm, "\n" );

    // Compute the L2 error for the solution
    L2Error = feSpacePtr->l2ErrorWeighted( Members->getAnalyticalSolution(),
                                           solutionField->getVector(),
                                           Members->getUOne(),
                                           multirateData->dataTimePtr()->endTime() );

    // Display the L2 error for the solution
    displayer->leaderPrint( " L2 error:                       ", L2Error, "\n" );

    // Compute the L2 realative error for the solution
    L2RelativeError = L2Error / L2Norm;

    // Display the L2 relative error for the solution
    displayer->leaderPrint( " L2 relative error:              ", L2RelativeError, "\n" );

    // Stop chronoError
    chronoError.stop();

    // The leader process print chronoError
    displayer->leaderPrint( "Time for compute errors ", chronoError.diff(), "\n" );

    // Stop chronoTotal
    chronoTotal.stop();

    // The leader process print chronoTotal
    displayer->leaderPrint( "Total time for the computation ", chronoTotal.diff(), "\n" );

    // Return the error, needed for the succes/failure of the test
    return L2Error;


}
