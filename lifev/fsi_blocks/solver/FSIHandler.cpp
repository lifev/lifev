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

/*!
 @brief Settings - File handling the solution of the FSI problem

 @author Davide Forti <davide.forti@epfl.ch>
 @date 28-10-2014

@maintainer Simone Deparis <simone.deparis@epfl.ch>
*/

#include <lifev/core/LifeV.hpp>
#include <lifev/fsi_blocks/solver/FSIHandler.hpp>

namespace LifeV
{

FSIHandler::FSIHandler( const commPtr_Type& communicator ) :
M_comm ( communicator ),
M_displayer ( communicator ),
M_dt (0.0),
M_t_zero (0.0),
M_t_end (0.0),
M_exporterFluid (),
M_exporterStructure (),
M_applyOperatorResidual(new Operators::FSIApplyOperator),
M_applyOperatorJacobian(new Operators::FSIApplyOperator),
M_applyOperatorJacobianNonConforming(new Operators::FSIApplyOperatorNonConforming),
M_prec(new Operators::DirichletNeumannPreconditioner),
M_invOper(),
M_useShapeDerivatives( false ),
M_printResiduals ( false ),
M_printSteps ( false ),
M_NewtonIter ( 0 ),
M_extrapolateInitialGuess ( false ),
M_orderExtrapolationInitialGuess ( 3 ),
M_usePartitionedMeshes ( false ),
M_subiterateFluidDirichlet ( false ),
M_gravity ( 0.0 ),
M_considerGravity ( false ),
M_moveMesh ( true ),
M_nonconforming ( false ),
M_lambda_num_structure ( false ),
M_precPtrBuilt ( false ),
M_linearElasticity ( true ),
M_disregardRestart ( false ),
M_prescribeInflowFlowrate ( false ),
M_useBDF ( false )
{
}

FSIHandler::~FSIHandler( )
{
}

void
FSIHandler::setGravity ( const Real& gravity, const Real& gravity_direction)
{
	M_gravity = gravity;
	M_gravityDirection = gravity_direction;
	M_considerGravity = true;
}

void
FSIHandler::setDatafile( const GetPot& dataFile)
{
    M_datafile = dataFile;
    setParameterLists( );
}

void
FSIHandler::setParameterLists( )
{
	Teuchos::RCP<Teuchos::ParameterList> solversOptions = Teuchos::getParametersFromXmlFile ("solversOptionsFast.xml");
	std::string precType = M_datafile ("fluid/preconditionerType", "none");
	M_prec->setFluidPreconditioner(precType);
	M_prec->setOptions(*solversOptions);
	setSolversOptions(*solversOptions);
}

void
FSIHandler::setSolversOptions(const Teuchos::ParameterList& solversOptions)
{
    std::shared_ptr<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("MonolithicOperator")) );
    M_pListLinSolver = monolithicOptions;
}

void
FSIHandler::readMeshes( )
{
    M_fluidMesh.reset ( new mesh_Type ( M_comm ) );
    M_meshDataFluid.reset( new MeshData ( ) );
    M_meshDataFluid->setup (M_datafile, "fluid/space_discretization");
    readMesh (*M_fluidMesh, *M_meshDataFluid);
    int severityLevelFluid = M_fluidMesh->check(1,true,true);

    M_displayer.leaderPrintMax ( "\tSeverity Level fluid mesh is: ", severityLevelFluid ) ;

    M_structureMesh.reset ( new mesh_Type ( M_comm ) );
    M_meshDataStructure.reset( new MeshData ( ) );
    M_meshDataStructure->setup (M_datafile, "solid/space_discretization");
    readMesh (*M_structureMesh, *M_meshDataStructure);
    int severityLevelSolid = M_structureMesh->check(1,true,true);

    M_displayer.leaderPrintMax ( "\tSeverity Level solid mesh is: ", severityLevelSolid ) ;

}

void
FSIHandler::partitionMeshes( )
{
    M_fluidPartitioner.reset ( new MeshPartitioner< mesh_Type > (M_fluidMesh, M_comm) );
    M_fluidLocalMesh.reset ( new mesh_Type ( M_comm ) );
    M_fluidLocalMesh = M_fluidPartitioner->meshPartition();

    M_structurePartitioner.reset ( new MeshPartitioner< mesh_Type > (M_structureMesh, M_comm) );
    M_structureLocalMesh.reset ( new mesh_Type ( M_comm ) );
    M_structureLocalMesh = M_structurePartitioner->meshPartition();
}

void
FSIHandler::readPartitionedMeshes( )
{
	const std::string fluidHdf5File (M_datafile ("offlinePartioner/fluidPartitionedMesh", "fluid.h5") );
	const std::string solidHdf5File (M_datafile ("offlinePartioner/solidPartitionedMesh", "solid.h5") );

	std::shared_ptr<Epetra_MpiComm> comm = std::dynamic_pointer_cast<Epetra_MpiComm>(M_comm);

	// Load fluid mesh part from HDF5
	M_displayer.leaderPrint ( "\tReading the fluid mesh parts\n" ) ;
	PartitionIO<mesh_Type > partitionIO (fluidHdf5File, comm);
	partitionIO.read (M_fluidLocalMesh);

	// Load fluid mesh part from HDF5
	M_displayer.leaderPrint ( "\tReading the solid mesh parts\n" ) ;
	PartitionIO<mesh_Type > partitionIOstructure (solidHdf5File, comm);
	partitionIOstructure.read (M_structureLocalMesh);
}

void FSIHandler::setup ( )
{
	M_useBDF = M_datafile ( "solid/time_discretization/useBDF", false );

	M_linearElasticity = M_datafile ( "solid/linear_elasticity", true );

	M_saveEvery = M_datafile ( "exporter/save_every", 1 );

	M_disregardRestart = M_datafile ( "exporter/disregardRestart", false );

	M_counterSaveEvery = M_saveEvery;

	M_usePartitionedMeshes = M_datafile ( "offlinePartioner/readPartitionedMeshes", false );

    M_restart = M_datafile ( "importer/restart", false );
    
	// Fluid
	M_fluid.reset ( new NavierStokesSolverBlocks ( M_datafile, M_comm ) );
	M_fluid->setup ( M_fluidLocalMesh );

    // Temporary fix for Shape derivatives
    M_fluid->pFESpace()->setQuadRule(M_fluid->uFESpace()->qr());

	// Structure
    if ( M_linearElasticity )
    {
    	M_structure.reset ( new LinearElasticity ( M_comm ) );
    }
    else
    {
    	M_structureNeoHookean.reset ( new NeoHookean ( M_comm ) );
    }

	// setup of the structural class
	setupStructure();

	// This beacuse the ale solver requires that the FESpace is given from outside
	createAleFESpace();

	updateBoundaryConditions();

	M_extrapolateInitialGuess = M_datafile ( "newton/extrapolateInitialGuess", false );
	M_orderExtrapolationInitialGuess = M_datafile ( "newton/orderExtrapolation", 3 );

    // Exporters
    setupExporters( );
    
	initializeTimeAdvance ( );

	// Setting auxiliary variables for the fluid solver
	M_fluid->setAlpha(M_fluidTimeAdvance->alpha());
	M_fluid->setTimeStep(M_dt);
	M_fluid->buildSystem();
	M_fluid->setupPostProc();

	// Ale
	M_ale.reset( new ALESolver ( *M_aleFESpace, M_comm ) );
	M_ale->setUp( M_datafile );

	// Structure
	if ( M_linearElasticity )
	{
		if ( M_useBDF )
		{
			M_structure->assemble_matrices(M_dt, M_structureTimeAdvanceBDF->massCoefficient(), M_structureBC, M_useBDF);
		}
		else
		{
			M_structure->assemble_matrices(M_dt, M_structureTimeAdvance->get_beta(), M_structureBC, M_useBDF);
		}
	}

	// Data needed by the Newton algorithm
	M_relativeTolerance = M_datafile ( "newton/reltol", 1.e-4);
	M_absoluteTolerance = M_datafile ( "newton/abstol", 1.e-4);
	M_etaMax = M_datafile ( "newton/etamax", 1e-4);
	M_maxiterNonlinear = M_datafile ( "newton/maxiter", 10);

	M_displayer.leaderPrintMax ( " Maximum Newton iterations = ", M_maxiterNonlinear ) ;

	M_nonLinearLineSearch = M_datafile ( "newton/NonLinearLineSearch", 0);
	if (M_comm->MyPID() == 0)
    {
		M_out_res.open ("residualsNewton");
        M_outputLinearIterations.open("linearIterations.dat");
        M_outputPreconditionerComputation.open("preconditionerComputation.dat");
        M_outputTimeLinearSolver.open ("solutionLinearSystem.dat");
    }
    
    M_useShapeDerivatives = M_datafile ( "newton/useShapeDerivatives", false);
    M_subiterateFluidDirichlet = M_datafile ( "fluid/subiterateFluidDirichlet", false);

	M_printResiduals = M_datafile ( "newton/output_Residuals", false);
	M_printSteps = M_datafile ( "newton/output_Steps", false);

	// setting also in the preconditioner if using shape derivatives
	M_prec->setUseShapeDerivatives(M_useShapeDerivatives);

	// setting also in the preconditioner if using shape derivatives
	M_prec->setSubiterateFluidDirichlet(M_subiterateFluidDirichlet);

	// Solver
	std::string solverType(M_pListLinSolver->get<std::string>("Linear Solver Type"));
	M_invOper.reset(Operators::InvertibleOperatorFactory::instance().createObject(solverType));
	M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));

	// Preconditioner
	M_prec->setComm ( M_comm );

	// Output Files
	if(M_comm->MyPID()==0)
	{
		M_outputTimeStep << std::scientific;
		M_outputTimeStep.open ("TimeStep.dat" );

		if (M_printResiduals)
			M_outputResiduals.open ("Residuals.txt" );

		if (M_printSteps)
			M_outputSteps.open ("Steps.txt" );
	}

	M_moveMesh = M_datafile ( "fluid/mesh/move_mesh", true);
}

void FSIHandler::setupExporters( )
{
	// Exporters
	std::string outputNameFluid = M_datafile ( "exporter/fluid_filename", "fluid");
	std::string outputNameStructure = M_datafile ( "exporter/structure_filename", "fluid");
	instantiateExporter(M_exporterFluid, M_fluidLocalMesh, outputNameFluid);
	instantiateExporter(M_exporterStructure, M_structureLocalMesh, outputNameStructure);

	M_fluidVelocity.reset ( new VectorEpetra ( M_fluid->uFESpace()->map(), M_exporterFluid->mapType() ) );
	M_fluidPressure.reset ( new VectorEpetra ( M_fluid->pFESpace()->map(), M_exporterFluid->mapType() ) );
	M_fluidDisplacement.reset ( new VectorEpetra ( M_aleFESpace->map(), M_exporterFluid->mapType() ) );
    M_Lagrange.reset ( new VectorEpetra ( M_fluid->uFESpace()->map(), M_exporterFluid->mapType() ) );
	M_structureDisplacement.reset ( new VectorEpetra ( M_displacementFESpace->map(), M_exporterStructure->mapType() ) );
	M_structureVelocity.reset ( new VectorEpetra ( M_displacementFESpace->map(), M_exporterStructure->mapType() ) );
	M_structureAcceleration.reset ( new VectorEpetra ( M_displacementFESpace->map(), M_exporterStructure->mapType() ) );

    M_fluidVelocity->zero();
    M_fluidPressure->zero();
    M_fluidDisplacement->zero();
    M_Lagrange->zero(); // vector needed for a correct restart of the FSI simulation
    M_structureDisplacement->zero();
    M_structureVelocity->zero();
    M_structureAcceleration->zero();
    
	M_exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f - velocity", M_fluid->uFESpace(), M_fluidVelocity, UInt (0) );
	M_exporterFluid->addVariable ( ExporterData<mesh_Type>::ScalarField, "f - pressure", M_fluid->pFESpace(), M_fluidPressure, UInt (0) );
	M_exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f - displacement", M_aleFESpace, M_fluidDisplacement, UInt (0) );
    M_exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f - lagrange", M_fluid->uFESpace(), M_Lagrange, UInt (0) );
    
	M_exporterStructure->addVariable ( ExporterData<mesh_Type>::VectorField, "s - displacement", M_displacementFESpace, M_structureDisplacement, UInt (0) );
	if ( !M_useBDF )
	{
		M_exporterStructure->addVariable ( ExporterData<mesh_Type>::VectorField, "s - velocity", M_displacementFESpace, M_structureVelocity, UInt (0) );
		M_exporterStructure->addVariable ( ExporterData<mesh_Type>::VectorField, "s - acceleration", M_displacementFESpace, M_structureAcceleration, UInt (0) );
	}
}

void
FSIHandler::instantiateExporter( std::shared_ptr< Exporter<mesh_Type > >& exporter,
							     const meshPtr_Type& localMesh,
							     const std::string& outputName)
{
	std::string const exporterType =  M_datafile ( "exporter/type", "ensight");
#ifdef HAVE_HDF5
	if (exporterType.compare ("hdf5") == 0)
	{
		exporter.reset ( new ExporterHDF5<mesh_Type > ( M_datafile, outputName ) );
		exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
		exporter->setMeshProcId ( localMesh, M_comm->MyPID() );
	}
	else if(exporterType.compare ("vtk") == 0)
#endif
	{
		exporter.reset ( new ExporterVTK<mesh_Type > ( M_datafile, outputName ) );
		exporter->setPostDir ( "./" ); // This is a test to see if M_post_dir is working
		exporter->setMeshProcId ( localMesh, M_comm->MyPID() );
	}
}

void FSIHandler::setupStructure ( )
{
	const std::string dOrder = M_datafile ( "solid/space_discretization/order", "P2");

	Real density = M_datafile("solid/physics/density",0.0);
	Real young = M_datafile("solid/model/young",0.0);
	Real poisson = M_datafile("solid/model/poisson",0.0);

	if ( M_linearElasticity )
	{
		M_structure->setCoefficients(density, young, poisson);

		if ( M_datafile("solid/thin_layer/use_thin", false ) && M_datafile("solid/linear_elasticity", true ) )
		{
			M_structure->setCoefficientsThinLayer ( M_datafile("solid/thin_layer/density",0.0),
										  	  	  	M_datafile("solid/thin_layer/young_thin",0.0),
										  	  	  	M_datafile("solid/thin_layer/poisson_thin",0.0),
										  	  	  	M_datafile("solid/thin_layer/h_thin",0.0),
										  	  	  	M_datafile("interface/flag",0.0) );
		}

		M_structure->setup(M_structureLocalMesh, dOrder);
		M_displacementFESpace = M_structure->fespace();
		M_displacementETFESpace = M_structure->et_fespace();
		M_displacementFESpaceScalar.reset ( new FESpace_Type (M_structureLocalMesh, dOrder, 1, M_comm) );
	}
	else
	{
		M_structureNeoHookean->setCoefficients(density, young, poisson);
		M_structureNeoHookean->setup(M_structureLocalMesh, dOrder);
		M_displacementFESpace = M_structureNeoHookean->fespace();
		M_displacementETFESpace = M_structureNeoHookean->et_fespace();
		M_displacementFESpaceScalar.reset ( new FESpace_Type (M_structureLocalMesh, dOrder, 1, M_comm) );
	}

	if ( !M_usePartitionedMeshes )
		M_displacementFESpaceSerial.reset ( new FESpace_Type (M_structureMesh, dOrder, 3, M_comm) );

	M_displayer.leaderPrintMax ( " Number of DOFs for the structure = ", M_displacementFESpace->dof().numTotalDof()*3 ) ;

}

void FSIHandler::createAleFESpace()
{
	const std::string aleOrder = M_datafile ( "ale/space_discretization/order", "P2");
	M_aleFESpace.reset ( new FESpace_Type (M_fluidLocalMesh, aleOrder, 3, M_comm) );
	M_displayer.leaderPrintMax ( " Number of DOFs for the ale = ", M_aleFESpace->dof().numTotalDof()*3 ) ;
}
    
void FSIHandler::setBoundaryConditions ( const bcPtr_Type& fluidBC,
                                         const bcPtr_Type& structureBC,
                                         const bcPtr_Type& aleBC )
{
    M_fluidBC.reset ( new BCHandler ( ) );
    M_fluidBC_residual_natural.reset ( new BCHandler ( ) );
    M_fluidBC_residual_essential.reset ( new BCHandler ( ) );
    
    for ( std::vector<BCBase>::iterator it = fluidBC->begin(); it != fluidBC->end(); it++ )
    {
        if ( it->type() > 3 ) // meaning esssential
        {
            M_fluidBC->addBC( *it );
            M_fluidBC_residual_essential->addBC( *it );
        }
        else if ( it->type() == 0 ) // meaning natural
        {
            M_fluidBC_residual_natural->addBC( *it );
        }
    }
    
    M_structureBC.reset ( new BCHandler ( ) );
    M_structureBC_residual_natural.reset ( new BCHandler ( ) );
    M_structureBC_residual_essential.reset ( new BCHandler ( ) );
    
    for ( std::vector<BCBase>::iterator it = structureBC->begin(); it != structureBC->end(); it++ )
    {
        if ( it->type() > 3 ) // meaning esssential
        {
            M_structureBC->addBC( *it );
            M_structureBC_residual_essential->addBC( *it );
        }
        else if ( it->type() == 0 ) // meaning natural
        {
            M_structureBC_residual_natural->addBC( *it );
        }
    }
    
    M_aleBC.reset ( new BCHandler ( ) );
    M_aleBC_residual_natural.reset ( new BCHandler ( ) );
    M_aleBC_residual_essential.reset ( new BCHandler ( ) );
    
    for ( std::vector<BCBase>::iterator it = aleBC->begin(); it != aleBC->end(); it++ )
    {
        if ( it->type() > 3 ) // meaning esssential
        {
            M_aleBC->addBC( *it );
            M_aleBC_residual_essential->addBC( *it );
        }
        else if ( it->type() == 0 ) // meaning natural
        {
            M_aleBC_residual_natural->addBC( *it );
        }
    }
    
    for ( std::vector<BCBase>::iterator it = M_interfaceFluidBC->begin(); it != M_interfaceFluidBC->end(); it++ )
    {
        M_aleBC->addBC( *it );
    }
}

void FSIHandler::updateBoundaryConditions ( )
{
	M_fluidBC->bcUpdate ( *M_fluid->uFESpace()->mesh(), M_fluid->uFESpace()->feBd(), M_fluid->uFESpace()->dof() );
	M_structureBC->bcUpdate ( *M_displacementFESpace->mesh(), M_displacementFESpace->feBd(), M_displacementFESpace->dof() );
	M_aleBC->bcUpdate ( *M_aleFESpace->mesh(), M_aleFESpace->feBd(), M_aleFESpace->dof() );
	M_aleBC_residual_essential->bcUpdate ( *M_aleFESpace->mesh(), M_aleFESpace->feBd(), M_aleFESpace->dof() );
    M_aleBC_residual_natural->bcUpdate ( *M_aleFESpace->mesh(), M_aleFESpace->feBd(), M_aleFESpace->dof() );
}

void FSIHandler::initializeTimeAdvance ( )
{
	// Fluid
	M_fluidTimeAdvance.reset( new TimeAndExtrapolationHandler ( ) );
	M_dt       = M_datafile("fluid/time_discretization/timestep",0.0);
	M_t_zero   = M_datafile("fluid/time_discretization/initialtime",0.0);
	M_t_end    = M_datafile("fluid/time_discretization/endtime",0.0);
	M_orderBDF = M_datafile("fluid/time_discretization/BDF_order",2);
    
	// Initialize and create objects for the fluid time advance
	M_fluidTimeAdvance->setBDForder(M_orderBDF);
	M_fluidTimeAdvance->setMaximumExtrapolationOrder(M_orderBDF);
	M_fluidTimeAdvance->setTimeStep(M_dt);
    vector_Type velocityInitial ( M_fluid->uFESpace()->map() );
    std::vector<vector_Type> initialStateVelocity;
    
    // Initialize and create objects for the ALE time advance
    M_aleTimeAdvance.reset( new TimeAndExtrapolationHandler ( ) );
    M_aleTimeAdvance->setMaximumExtrapolationOrder(M_orderBDF);
    M_aleTimeAdvance->setBDForder(M_orderBDF);
    M_aleTimeAdvance->setTimeStep(M_dt);
    vector_Type disp_mesh_initial ( M_aleFESpace->map() );
    std::vector<vector_Type> initialStateALE;
    std::vector<vector_Type> initialStateSolid;
    
    if ( M_useBDF )
    {
    	M_structureTimeAdvanceBDF.reset ( new BDFSecondOrderDerivative );
    	M_structureTimeAdvanceBDF->setTimeStep(M_dt);
    	M_orderBDFSolid = M_datafile("solid/time_discretization/BDF_order",2);
    	M_structureTimeAdvanceBDF->setBDForder( M_orderBDFSolid );
    }
    else
    {
		// Structure time advance
		M_structureTimeAdvance.reset ( new Newmark );
		M_structureTimeAdvance->set_beta(0.49); // 0.25
		M_structureTimeAdvance->set_gamma(0.9); // 0.5
		M_structureTimeAdvance->set_timestep(M_dt);
    }

    if ( !M_restart )
    {
        // Fluid
        velocityInitial.zero() ;
        for ( UInt i = 0; i < M_orderBDF; ++i )
            initialStateVelocity.push_back(velocityInitial);
    
        // ALE
        disp_mesh_initial.zero() ;
        for ( UInt i = 0; i < M_orderBDF; ++i )
            initialStateALE.push_back(disp_mesh_initial);

        // Solid
        if ( M_useBDF )
        {
        	vector_Type ds_initial ( M_displacementFESpace->map() );
        	ds_initial.zero();
        	for ( UInt i = 0; i < (M_orderBDFSolid+1); ++i )
        		initialStateSolid.push_back(ds_initial);
        }
        else
        {
        	vectorPtr_Type ds_initial ( new vector_Type ( M_displacementFESpace->map() ) );
        	vectorPtr_Type vs_initial ( new vector_Type ( M_displacementFESpace->map() ) );
        	vectorPtr_Type as_initial ( new vector_Type ( M_displacementFESpace->map() ) );
        	ds_initial->zero();
        	vs_initial->zero();
        	as_initial->zero();
        	M_structureTimeAdvance->initialize(ds_initial, vs_initial, as_initial);
        }
    }
    else
    {
    	std::string const importerType       =  M_datafile ( "importer/type", "hdf5");
    	std::string const fileNameFluid      =  M_datafile ( "importer/fluid_filename", "SolutionRestarted");
    	std::string const fileNameStructure  =  M_datafile ( "importer/structure_filename", "SolutionRestarted");
    	std::string const initialLoaded      =  M_datafile ( "importer/initSol", "NO_DEFAULT_VALUE");

    	M_importerFluid.reset ( new  ExporterHDF5<mesh_Type > ( M_datafile, fileNameFluid) );
    	M_importerFluid->setMeshProcId (M_fluid->uFESpace()->mesh(), M_comm->MyPID() );

    	M_importerStructure.reset ( new  ExporterHDF5<mesh_Type > ( M_datafile, fileNameStructure) );
    	M_importerStructure->setMeshProcId (M_displacementFESpace->mesh(), M_comm->MyPID() );

    	std::string iterationString;
    	iterationString = initialLoaded;

    	if ( M_useBDF )
    	{
    		for (UInt iterInit = 0; iterInit < (M_orderBDF+1); iterInit++ )
    		{
    			vectorPtr_Type velocityRestart;
    			velocityRestart.reset ( new vector_Type (M_fluid->uFESpace()->map(),  Unique ) );
    			velocityRestart->zero();

    			vectorPtr_Type pressureRestart;
    			pressureRestart.reset ( new vector_Type (M_fluid->pFESpace()->map(),  Unique ) );
    			pressureRestart->zero();

    			vectorPtr_Type aleRestart;
    			aleRestart.reset ( new vector_Type (M_aleFESpace->map(),  Unique ) );
    			aleRestart->zero();

    			vectorPtr_Type lagrangeRestart;
    			lagrangeRestart.reset ( new vector_Type (M_fluid->uFESpace()->map(),  Unique ) );
    			lagrangeRestart->zero();

    			vectorPtr_Type structureRestart;
    			structureRestart.reset ( new vector_Type (M_displacementFESpace->map(),  Unique ) );
    			structureRestart->zero();

    			LifeV::ExporterData<mesh_Type> velocityReader (LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("f - velocity." + iterationString),
    					M_fluid->uFESpace(),
    					velocityRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> pressureReader (LifeV::ExporterData<mesh_Type>::ScalarField,
    					std::string ("f - pressure." + iterationString),
    					M_fluid->pFESpace(),
    					pressureRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> aleReader      (LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("f - displacement." + iterationString),
    					M_aleFESpace,
    					aleRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> lagrangeReader (LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("f - lagrange." + iterationString),
    					M_fluid->uFESpace(),
    					lagrangeRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> structureReader(LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("s - displacement." + iterationString),
    					M_displacementFESpace,
    					structureRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			M_importerFluid->readVariable (velocityReader);
    			M_importerFluid->readVariable (pressureReader);
    			M_importerFluid->readVariable (aleReader);
    			M_importerFluid->readVariable (lagrangeReader);
    			M_importerStructure->readVariable (structureReader);

    			int iterations = std::atoi (iterationString.c_str() );

    			if ( iterInit == 0 )
    			{
    				*M_fluidVelocity = *velocityRestart;
    				*M_fluidPressure = *pressureRestart;
    				*M_fluidDisplacement = *aleRestart;
    				*M_Lagrange = *lagrangeRestart;
    				*M_structureDisplacement = *structureRestart;
    			}

    			if ( M_datafile ( "importer/change_dt", false ) )
    			{
    				Real previousDt = M_datafile ( "importer/previous_dt", 0.0001 );
    				iterations = iterations - M_dt/previousDt;
    			}
    			else
    			{
    				iterations--;
    			}

    			std::ostringstream iter;
    			iter.fill ( '0' );
    			iter << std::setw (5) << ( iterations );
    			iterationString = iter.str();

    			if ( iterInit < M_orderBDFSolid )
    			{
    				initialStateVelocity.push_back(*velocityRestart);
    				initialStateALE.push_back(*aleRestart);
    			}
    			initialStateSolid.push_back(*structureRestart);
    		}

    		std::reverse(initialStateVelocity.begin(),initialStateVelocity.end());
    		std::reverse(initialStateALE.begin(),initialStateALE.end());
    		std::reverse(initialStateSolid.begin(),initialStateSolid.end());
    	}
    	else
    	{
    		for (UInt iterInit = 0; iterInit < M_orderBDF; iterInit++ )
    		{
    			vectorPtr_Type velocityRestart;
    			velocityRestart.reset ( new vector_Type (M_fluid->uFESpace()->map(),  Unique ) );
    			velocityRestart->zero();

    			vectorPtr_Type pressureRestart;
    			pressureRestart.reset ( new vector_Type (M_fluid->pFESpace()->map(),  Unique ) );
    			pressureRestart->zero();

    			vectorPtr_Type aleRestart;
    			aleRestart.reset ( new vector_Type (M_aleFESpace->map(),  Unique ) );
    			aleRestart->zero();

    			vectorPtr_Type lagrangeRestart;
    			lagrangeRestart.reset ( new vector_Type (M_fluid->uFESpace()->map(),  Unique ) );
    			lagrangeRestart->zero();

    			vectorPtr_Type structureRestart;
    			structureRestart.reset ( new vector_Type (M_displacementFESpace->map(),  Unique ) );
    			structureRestart->zero();

    			vectorPtr_Type structureVelocityRestart;
    			structureVelocityRestart.reset ( new vector_Type (M_displacementFESpace->map(),  Unique ) );
    			structureVelocityRestart->zero();

    			vectorPtr_Type structureAccelerationRestart;
    			structureAccelerationRestart.reset ( new vector_Type (M_displacementFESpace->map(),  Unique ) );
    			structureAccelerationRestart->zero();

    			LifeV::ExporterData<mesh_Type> velocityReader (LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("f - velocity." + iterationString),
    					M_fluid->uFESpace(),
    					velocityRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> pressureReader (LifeV::ExporterData<mesh_Type>::ScalarField,
    					std::string ("f - pressure." + iterationString),
    					M_fluid->pFESpace(),
    					pressureRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> aleReader      (LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("f - displacement." + iterationString),
    					M_aleFESpace,
    					aleRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> lagrangeReader (LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("f - lagrange." + iterationString),
    					M_fluid->uFESpace(),
    					lagrangeRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> structureReader(LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("s - displacement." + iterationString),
    					M_displacementFESpace,
    					structureRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> structureVelocityReader(LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("s - velocity." + iterationString),
    					M_displacementFESpace,
    					structureVelocityRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );

    			LifeV::ExporterData<mesh_Type> structureAccelerationReader(LifeV::ExporterData<mesh_Type>::VectorField,
    					std::string ("s - acceleration." + iterationString),
    					M_displacementFESpace,
    					structureAccelerationRestart,
    					UInt (0),
    					LifeV::ExporterData<mesh_Type>::UnsteadyRegime );


    			M_importerFluid->readVariable (velocityReader);
    			M_importerFluid->readVariable (pressureReader);
    			M_importerFluid->readVariable (aleReader);
    			M_importerFluid->readVariable (lagrangeReader);
    			M_importerStructure->readVariable (structureReader);
    			M_importerStructure->readVariable (structureVelocityReader);
    			M_importerStructure->readVariable (structureAccelerationReader);

    			int iterations = std::atoi (iterationString.c_str() );

    			if ( iterInit == 0 )
    			{
    				*M_fluidVelocity = *velocityRestart;
    				*M_fluidPressure = *pressureRestart;
    				*M_fluidDisplacement = *aleRestart;
    				*M_Lagrange = *lagrangeRestart;
    				*M_structureDisplacement = *structureRestart;
    				*M_structureVelocity = *structureVelocityRestart;
    				*M_structureAcceleration = *structureAccelerationRestart;
    			}

    			iterations--;

    			std::ostringstream iter;
    			iter.fill ( '0' );
    			iter << std::setw (5) << ( iterations );
    			iterationString = iter.str();

    			initialStateVelocity.push_back(*velocityRestart);
    			initialStateALE.push_back(*aleRestart);

    			if ( iterInit == 0 )
    			{
    				M_structureTimeAdvance->restart(M_structureDisplacement, M_structureVelocity, M_structureAcceleration);
    			}
    		}

    		// For BDF 1 it does not change anything, for BDF2 it is necessary. We do the std::reverse
    		// just to the fluid and ALE stencil because the one of the structure is ordered in the
    		// other way round.
    		std::reverse(initialStateVelocity.begin(),initialStateVelocity.end());
    		std::reverse(initialStateALE.begin(),initialStateALE.end());
    	}
    }
    
    // Fluid
	M_fluidTimeAdvance->initialize(initialStateVelocity);

    // ALE
    M_aleTimeAdvance->initialize(initialStateALE);

    // Solid - BDF case
    if ( M_useBDF )
    {
    	M_structureTimeAdvanceBDF->initialize(initialStateSolid);
    }

    // Post-Processing of the initial solution
	M_exporterFluid->postProcess(M_t_zero);
	M_exporterStructure->postProcess(M_t_zero);
}

void FSIHandler::buildInterfaceMaps ()
{
	M_nonconforming = M_datafile("interface/nonconforming", false);
	M_lambda_num_structure = M_datafile("interface/lambda_num_structure", false);
	markerID_Type interface = M_datafile("interface/flag", 1);
	Real tolerance = M_datafile("interface/tolerance", 1.0);
	Int flag = M_datafile("interface/fluid_vertex_flag", 123);
	M_useMasses = M_datafile("interface/useMasses", true);

	M_displayer.leaderPrintMax ( " Flag of the interface = ", interface ) ;
	M_displayer.leaderPrintMax ( " Tolerance for dofs on the interface = ", tolerance ) ;

	if ( !M_nonconforming )
	{
		if ( !M_usePartitionedMeshes )
		{
			M_dofStructureToFluid.reset ( new DOFInterface3Dto3D );
			M_dofStructureToFluid->setup ( M_fluid->uFESpace()->refFE(), M_fluid->uFESpace()->dof(), M_displacementFESpaceSerial->refFE(), M_displacementFESpaceSerial->dof() );
			M_dofStructureToFluid->update ( *M_fluid->uFESpace()->mesh(), interface, *M_displacementFESpaceSerial->mesh(), interface, tolerance, &flag);
			M_localDofMap.reset(new std::map<UInt, UInt> ( M_dofStructureToFluid->localDofMap ( ) ) );
			M_displacementFESpaceSerial.reset();
		}
		else
		{
			const std::string interfaceHdf5File (M_datafile ("offlinePartioner/interfacePartitioned", "interface.h5") );
			std::shared_ptr<Epetra_MpiComm> comm = std::dynamic_pointer_cast<Epetra_MpiComm>(M_comm);
			DOFInterfaceIO interfaceIO (interfaceHdf5File, comm);
			interfaceIO.read (M_localDofMap);
		}

		createInterfaceMaps ( *M_localDofMap );

		if ( M_lambda_num_structure )
			constructInterfaceMap ( *M_localDofMap, M_displacementFESpace->map().map(Unique)->NumGlobalElements()/nDimensions );
		else
			constructInterfaceMap ( *M_localDofMap, M_fluid->uFESpace()->map().map(Unique)->NumGlobalElements()/nDimensions );

		M_displayer.leaderPrintMax ( " Number of DOFs on the interface = ", M_lagrangeMap->mapSize() ) ;
	}
	else
	{
		M_displayer.leaderPrint ( " Using nonconforming meshes " ) ;

		// Reading xml parameter file for solver of the interpolation
		Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
		belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList_rbf3d.xml" );

		const std::string uOrder = M_datafile ( "fluid/space_discretization/vel_order", "P1");
		FESpacePtr_Type fluid_FESpace_whole ( new FESpace_Type ( M_fluidMesh, uOrder, 3, M_comm) );

		const std::string dOrder = M_datafile ( "solid/space_discretization/order", "P1");
		FESpacePtr_Type solid_FESpace_whole ( new FESpace_Type ( M_structureMesh, dOrder, 3, M_comm) );

		vectorPtr_Type fluid_known ( new vector_Type ( M_fluidVelocity->map() ) );
		vectorPtr_Type structure_known ( new vector_Type ( M_structureDisplacement->map() ) );
		fluid_known->zero();
		structure_known->zero();

		M_FluidToStructureInterpolant.reset( new interpolation_Type);
		M_FluidToStructureInterpolant->setup(M_datafile, belosList);
		M_FluidToStructureInterpolant->setMeshSize(MeshUtility::MeshStatistics::computeSize(*M_fluidMesh).maxH);
		M_FluidToStructureInterpolant->setFlag( M_datafile("interface/flag", 1) );
		M_FluidToStructureInterpolant->setVectors(fluid_known, structure_known);
		M_FluidToStructureInterpolant->buildTableDofs_known ( fluid_FESpace_whole );
		M_FluidToStructureInterpolant->buildTableDofs_unknown ( solid_FESpace_whole );
		M_FluidToStructureInterpolant->identifyNodes_known ( );
		M_FluidToStructureInterpolant->identifyNodes_unknown ( );
		M_FluidToStructureInterpolant->buildKnownInterfaceMap();
		M_FluidToStructureInterpolant->buildUnknownInterfaceMap();
		M_FluidToStructureInterpolant->buildOperators();
		M_FluidToStructureInterpolant->getKnownInterfaceMap(M_fluidInterfaceMap);

		M_StructureToFluidInterpolant.reset( new interpolation_Type);
		M_StructureToFluidInterpolant->setup(M_datafile, belosList);
		M_StructureToFluidInterpolant->setMeshSize(MeshUtility::MeshStatistics::computeSize(*M_structureMesh).maxH);
		M_StructureToFluidInterpolant->setFlag( M_datafile("interface/flag", 1) );
		M_StructureToFluidInterpolant->setVectors(structure_known, fluid_known);
		M_StructureToFluidInterpolant->buildTableDofs_known ( solid_FESpace_whole );
		M_StructureToFluidInterpolant->buildTableDofs_unknown ( fluid_FESpace_whole );
		M_StructureToFluidInterpolant->identifyNodes_known ( );
		M_StructureToFluidInterpolant->identifyNodes_unknown ( );
		M_StructureToFluidInterpolant->buildKnownInterfaceMap();
		M_StructureToFluidInterpolant->buildUnknownInterfaceMap();
		M_StructureToFluidInterpolant->buildOperators();
		M_StructureToFluidInterpolant->getKnownInterfaceMap(M_structureInterfaceMap);

		// Getting a scalar map that will be used for the lagrange multipliers
		M_FluidToStructureInterpolant->getInterpolationOperatorMap( M_lagrangeMapScalar );

		// Building map for the lagrange multipliers, vectorial
		M_lagrangeMap.reset ( new MapEpetra ( *M_lagrangeMapScalar ) );
		*M_lagrangeMap += *M_lagrangeMapScalar;
		*M_lagrangeMap += *M_lagrangeMapScalar;

		// Setting the vectors for the numeration of the interface in the fluid and the structure
		M_FluidToStructureInterpolant->getNumerationInterfaceKnown(M_numerationInterfaceFluid);
		M_StructureToFluidInterpolant->getNumerationInterfaceKnown(M_numerationInterfaceStructure);
	}
}

void FSIHandler::createInterfaceMaps(std::map<ID, ID> const& locDofMap)
{
	std::vector<int> dofInterfaceFluid;
	dofInterfaceFluid.reserve ( locDofMap.size() );

	std::map<ID, ID>::const_iterator i;

	for (UInt dim = 0; dim < nDimensions; ++dim)
		for ( i = locDofMap.begin(); i != locDofMap.end(); ++i )
		{
			dofInterfaceFluid.push_back (i->first + dim * M_fluid->uFESpace()->dof().numTotalDof() );
		}

	int* pointerToDofs (0);
	if (dofInterfaceFluid.size() > 0)
	{
		pointerToDofs = &dofInterfaceFluid[0];
	}

	M_fluidInterfaceMap.reset ( new MapEpetra ( -1, static_cast<int> (dofInterfaceFluid.size() ), pointerToDofs, M_fluid->uFESpace()->map().commPtr() ) );

	M_fluid->uFESpace()->map().commPtr()->Barrier();

	std::vector<int> dofInterfaceSolid;
	dofInterfaceSolid.reserve ( locDofMap.size() );

	for (UInt dim = 0; dim < nDimensions; ++dim)
		for ( i = locDofMap.begin(); i != locDofMap.end(); ++i )
		{
			dofInterfaceSolid.push_back (i->second + dim * M_displacementFESpace->dof().numTotalDof() );
		}

	pointerToDofs = 0;
	if (dofInterfaceSolid.size() > 0)
	{
		pointerToDofs = &dofInterfaceSolid[0];
	}

	M_structureInterfaceMap.reset ( new MapEpetra ( -1, static_cast<int> (dofInterfaceSolid.size() ), pointerToDofs, M_displacementFESpace->map().commPtr() ) );
}

void FSIHandler::constructInterfaceMap ( const std::map<ID, ID>& locDofMap,
										 const UInt subdomainMaxId )
{
	std::map<ID, ID>::const_iterator ITrow;

	Int numtasks = M_comm->NumProc();
	int* numInterfaceDof (new int[numtasks]);
	int pid = M_comm->MyPID();
	int numMyElements;

	if ( M_lambda_num_structure )
	{
		numMyElements = M_structureInterfaceMap->map (Unique)->NumMyElements();
	}
	else
	{
		numMyElements = M_fluidInterfaceMap->map (Unique)->NumMyElements();
	}

	numInterfaceDof[pid] = numMyElements;

	mapPtr_Type subMap;

	if ( M_lambda_num_structure )
	{
		subMap.reset ( new map_Type ( *M_structureInterfaceMap->map (Unique), (UInt) 0, subdomainMaxId) );
	}
	else
	{
		subMap.reset ( new map_Type (*M_fluidInterfaceMap->map (Unique), (UInt) 0, subdomainMaxId) );
	}

	M_numerationInterface.reset (new VectorEpetra (*subMap, Unique) );

	for (int j = 0; j < numtasks; ++j)
	{
		M_comm->Broadcast ( &numInterfaceDof[j], 1, j);
	}

	for (int j = numtasks - 1; j > 0 ; --j)
	{
		numInterfaceDof[j] = numInterfaceDof[j - 1];
	}
	numInterfaceDof[0] = 0;
	for (int j = 1; j < numtasks ; ++j)
	{
		numInterfaceDof[j] += numInterfaceDof[j - 1];
	}

	UInt l = 0;

	if ( M_lambda_num_structure )
	{
		Real M_interface = (UInt) M_structureInterfaceMap->map (Unique)->NumGlobalElements() / nDimensions;
		for (l = 0, ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
		{
			if (M_structureInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second) ) >= 0)
			{
				(*M_numerationInterface) [ITrow->second ] = l + (int) (numInterfaceDof[pid] / nDimensions);
				if ( (int) (*M_numerationInterface) (ITrow->second ) != floor (l + numInterfaceDof[pid] / nDimensions + 0.2) )
				{
					std::cout << "ERROR! the numeration of the coupling map is not correct" << std::endl;
				}
				++l;
			}
		}

		std::vector<int> couplingVector;
		couplingVector.reserve ( (int) (M_structureInterfaceMap->map (Unique)->NumMyElements() ) );

		for (UInt dim = 0; dim < nDimensions; ++dim)
		{
			for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
			{
				if (M_structureInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second) ) >= 0)
				{
					couplingVector.push_back ( (*M_numerationInterface) (ITrow->second ) + dim * M_interface );
				}
			}
		}

		M_lagrangeMap.reset (new MapEpetra (-1, static_cast< Int> ( couplingVector.size() ), &couplingVector[0], M_comm) );
	}
	else
	{
		Real M_interface = (UInt) M_fluidInterfaceMap->map (Unique)->NumGlobalElements() / nDimensions;
		for (l = 0, ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
		{
			if (M_fluidInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first) ) >= 0)
			{
				(*M_numerationInterface) [ITrow->first ] = l + (int) (numInterfaceDof[pid] / nDimensions);
				if ( (int) (*M_numerationInterface) (ITrow->first ) != floor (l + numInterfaceDof[pid] / nDimensions + 0.2) )
				{
					std::cout << "ERROR! the numeration of the coupling map is not correct" << std::endl;
				}
				++l;
			}
		}

		std::vector<int> couplingVector;
		couplingVector.reserve ( (int) (M_fluidInterfaceMap->map (Unique)->NumMyElements() ) );

		for (UInt dim = 0; dim < nDimensions; ++dim)
		{
			for ( ITrow = locDofMap.begin(); ITrow != locDofMap.end() ; ++ITrow)
			{
				if (M_fluidInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first) ) >= 0)
				{
					couplingVector.push_back ( (*M_numerationInterface) (ITrow->first ) + dim * M_interface );
				}
			}
		}

		M_lagrangeMap.reset (new MapEpetra (-1, static_cast< Int> ( couplingVector.size() ), &couplingVector[0], M_comm) );
	}

	delete [] numInterfaceDof;
}

void
FSIHandler::assembleCoupling ( )
{
	M_coupling.reset ( new FSIcouplingCE ( M_comm ) );

	if ( M_lambda_num_structure )
	{
		if ( M_useBDF )
		{
			M_coupling->setUp ( M_dt, M_structureInterfaceMap->mapSize()/3.0 , M_structureTimeAdvanceBDF->coefficientFirstDerivative(),
							    M_lagrangeMap, M_fluid->uFESpace(), M_displacementFESpace, M_numerationInterface );
		}
		else
		{
			M_coupling->setUp ( M_dt, M_structureInterfaceMap->mapSize()/3.0 , M_structureTimeAdvance->get_beta(), M_structureTimeAdvance->get_gamma(),
								M_lagrangeMap, M_fluid->uFESpace(), M_displacementFESpace, M_numerationInterface );
		}
	}
	else
	{
		if ( M_useBDF )
		{
			M_coupling->setUp ( M_dt, M_fluidInterfaceMap->mapSize()/3.0 , M_structureTimeAdvanceBDF->coefficientFirstDerivative(),
								M_lagrangeMap, M_fluid->uFESpace(), M_displacementFESpace, M_numerationInterface );
		}
		else
		{
			M_coupling->setUp ( M_dt, M_fluidInterfaceMap->mapSize()/3.0 , M_structureTimeAdvance->get_beta(), M_structureTimeAdvance->get_gamma(),
								M_lagrangeMap, M_fluid->uFESpace(), M_displacementFESpace, M_numerationInterface );
		}
	}

	M_coupling->buildBlocks ( *M_localDofMap, M_lambda_num_structure, M_useBDF );
}

void
FSIHandler::buildMonolithicMap ( )
{
	// Fluid velocity
	M_monolithicMap.reset( new map_Type ( M_fluid->uFESpace()->map() ) );

	// Fluid pressure
	*M_monolithicMap += M_fluid->pFESpace()->map();

	// Structure displacement
	*M_monolithicMap += M_displacementFESpace->map();

	// Weak stresses
	M_displayer.leaderPrintMax ("\nNumber of DOFs of the Lagrange multipliers: ", M_lagrangeMap->map(Unique)->NumGlobalElements());
	*M_monolithicMap += *M_lagrangeMap;

	// ALE
	*M_monolithicMap += M_aleFESpace->map();
}

void
FSIHandler::getMatrixStructure ( )
{
	M_matrixStructure.reset (new matrix_Type ( M_displacementFESpace->map() ) );

	M_matrixStructure->zero();

	M_interface_mass_structure_robin.reset (new matrix_Type ( M_displacementFESpace->map() ) );
	M_interface_mass_structure_robin->zero();

	if ( M_datafile("solid/robin_external_wall", false ) )
	{
		M_displayer.leaderPrint ("\nUsing Robin BC at the external wall of the structure\n");
		Real alpha_robin = M_datafile("solid/robin_elastic", 100000.0 );

		// ASSEMBLE ROBIN BC MATRIX AT THE INTERFACE
		QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
		{
			using namespace ExpressionAssembly;
			integrate ( boundary (M_displacementETFESpace->mesh(), M_datafile("solid/externalWallFlag", 210 ) ),
					myBDQR,
					M_displacementETFESpace,
					M_displacementETFESpace,
					value ( alpha_robin ) * dot(phi_i, phi_j)
			)
			>>M_interface_mass_structure_robin;
		}

	}

	M_interface_mass_structure_robin->globalAssemble();

	*M_matrixStructure += *M_interface_mass_structure_robin;
	*M_matrixStructure += *(M_structure->jacobian());

	M_matrixStructure->globalAssemble();
}

void
FSIHandler::get_structure_coupling_velocities ( )
{
	if ( M_useBDF )
	{
		M_rhsStructureVelocity.reset ( new vector_Type ( M_displacementFESpace->map(), Unique ) );

		M_structureTimeAdvanceBDF->first_der_old_dts( *M_rhsStructureVelocity );

		*M_rhsStructureVelocity *= ( -1.0 / M_dt );

	}
	else
	{
		M_structureTimeAdvance->compute_csi( );

		M_rhsStructureVelocity.reset ( new vector_Type ( M_displacementFESpace->map(), Unique ) );

		M_rhsStructureVelocity->zero();

		*M_rhsStructureVelocity = ( ( M_dt * M_structureTimeAdvance->get_gamma() ) * (*M_structureTimeAdvance->get_csi() )   -
									( M_dt * ( 1.0 - M_structureTimeAdvance->get_gamma() ) ) * ( *M_structureTimeAdvance->old_second_derivative() ) -
									*M_structureTimeAdvance->old_first_derivative() );
	}
}

void
FSIHandler::updateRhsCouplingVelocities ( )
{
	if ( M_lambda_num_structure )
	{
		vector_Type lambda (*M_structureInterfaceMap, Unique);
		structureToInterface ( lambda, *M_rhsStructureVelocity);
		M_rhsCouplingVelocities.reset( new VectorEpetra ( *M_lagrangeMap ) );
		M_rhsCouplingVelocities->zero();

		std::map<ID, ID>::const_iterator ITrow;

		UInt interface (M_structureInterfaceMap->mapSize()/3.0 );
		UInt totalDofs (M_displacementFESpace->dof().numTotalDof() );

		for (UInt dim = 0; dim < 3; ++dim)
		{
			for ( ITrow = M_localDofMap->begin(); ITrow != M_localDofMap->end() ; ++ITrow)
			{
				if (M_structureInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->second ) ) >= 0 )
				{
						(*M_rhsCouplingVelocities) [  (int) (*M_numerationInterface) [ITrow->second ] + dim * interface ] = lambda ( ITrow->second + dim * totalDofs );
				}
			}
		}
	}
	else
	{
		vector_Type lambda (*M_structureInterfaceMap, Unique);
		structureToInterface ( lambda, *M_rhsStructureVelocity);
		M_rhsCouplingVelocities.reset( new VectorEpetra ( *M_lagrangeMap ) );
		M_rhsCouplingVelocities->zero();

		std::map<ID, ID>::const_iterator ITrow;

		UInt interface (M_fluidInterfaceMap->mapSize()/3.0 );
		UInt totalDofs (M_displacementFESpace->dof().numTotalDof() );

		for (UInt dim = 0; dim < 3; ++dim)
		{
			for ( ITrow = M_localDofMap->begin(); ITrow != M_localDofMap->end() ; ++ITrow)
			{
				if (M_fluidInterfaceMap->map (Unique)->LID ( static_cast<EpetraInt_Type> (ITrow->first ) ) >= 0 )
				{
					(*M_rhsCouplingVelocities) [  (int) (*M_numerationInterface) [ITrow->first ] + dim * interface ] = lambda ( ITrow->second + dim * totalDofs );
				}
			}
		}
	}
}

void
FSIHandler::updateRhsCouplingVelocities_nonconforming ( )
{
	// Update known field for the interpolation
	M_StructureToFluidInterpolant->updateRhs ( M_rhsStructureVelocity );
	M_StructureToFluidInterpolant->interpolate();

	// Get the vector on the fluid interface
	M_StructureToFluidInterpolant->getSolutionOnGamma (M_rhsCouplingVelocities);
}

void
FSIHandler::structureToInterface (vector_Type& VectorOnGamma, const vector_Type& VectorOnStructure)
{
    if (VectorOnStructure.mapType() == Repeated)
    {
        vector_Type const  VectorOnStructureUnique (VectorOnStructure, Unique);
        structureToInterface (VectorOnGamma, VectorOnStructureUnique);
        return;
    }
    if (VectorOnGamma.mapType() == Repeated)
    {
        vector_Type  VectorOnGammaUn (VectorOnGamma.map(), Unique);
        structureToInterface ( VectorOnGammaUn, VectorOnStructure);
        VectorOnGamma = VectorOnGammaUn;
        return;
    }

    MapEpetra subMap (*VectorOnStructure.map().map (Unique), 0, VectorOnStructure.map().map (Unique)->NumGlobalElements() );
    vector_Type subVectorOnStructure (subMap, Unique);
    subVectorOnStructure.subset (VectorOnStructure, 0);
    VectorOnGamma = subVectorOnStructure;
}

void
FSIHandler::solveFSIproblem ( )
{
	LifeChrono iterChrono;
	LifeChrono smallThingsChrono;
	M_time = M_t_zero + M_dt;

	// Clean unpartitioned mesh
	M_fluidMesh.reset();
	M_structureMesh.reset();

	// Build monolithic map
	buildMonolithicMap ( );

    M_solution.reset ( new VectorEpetra ( *M_monolithicMap ) );
    M_solution->zero();

    if ( M_restart )
    {
        M_solution.reset ( new VectorEpetra ( *M_monolithicMap ) );
        M_solution->zero();

        M_solution->subset(*M_fluidVelocity, M_fluid->uFESpace()->map(), 0, 0);
        M_solution->subset(*M_fluidPressure, M_fluid->pFESpace()->map(), 0, M_fluid->uFESpace()->map().mapSize());
        M_solution->subset(*M_structureDisplacement, M_displacementFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize());
        M_LagrangeRestart.reset( new vector_Type (*M_lagrangeMap, Unique ) );
        if ( !M_nonconforming )
        {
        	*M_LagrangeRestart = ( *(M_coupling->fluidVelocityToLambda()) * (*M_Lagrange) );
        	M_solution->subset(*M_LagrangeRestart, *M_lagrangeMap, 0, M_fluid->uFESpace()->map().mapSize() +
        			M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() );
        }
        else
        {
        	vectorPtr_Type Lagrange ( new vector_Type (*M_lagrangeMap, Unique ) );
        	M_FluidToStructureInterpolant->restrictOmegaToGamma_Known( M_Lagrange, Lagrange );
        	*M_LagrangeRestart = *Lagrange;
        	M_solution->subset(*M_LagrangeRestart, *M_lagrangeMap, 0, M_fluid->uFESpace()->map().mapSize() +
        	        			M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() );
        }
        M_solution->subset(*M_fluidDisplacement, M_aleFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize() );
    }

	M_ale->applyBoundaryConditions ( *M_aleBC );

	// Apply boundary conditions for the structure problem (the matrix will not change during the simulation, it is linear elasticity)
	if ( M_linearElasticity )
	{
		getMatrixStructure ( );
		M_displayer.leaderPrint ( "\t Set and approximate structure block in the preconditioner.. " ) ;
		smallThingsChrono.start();
		M_prec->setStructureBlock ( M_matrixStructure );
		M_prec->updateApproximatedStructureMomentumOperator ( );
		smallThingsChrono.stop();
		M_displayer.leaderPrintMax ( "done in ", smallThingsChrono.diff() ) ;
	}

	// Set blocks for the preconditioners: geometry
	smallThingsChrono.reset();
	M_displayer.leaderPrint ( "\t Set and approximate geometry block in the preconditioner... " ) ;
	smallThingsChrono.start();
	M_prec->setGeometryBlock ( M_ale->matrix() );
	M_prec->updateApproximatedGeometryOperator ( );
	smallThingsChrono.stop();
	M_displayer.leaderPrintMax ( "done in ", smallThingsChrono.diff() ) ;

	if ( M_nonconforming )
	{
		M_prec->setCouplingOperators_nonconforming(M_FluidToStructureInterpolant, M_StructureToFluidInterpolant, M_lagrangeMap);
		if ( M_useMasses)
		{
			M_displayer.leaderPrint ( "[FSI] - Assemble interface mass of the structure for coupling of stresses\n" ) ;
			assembleStructureInterfaceMass ( );
		}
	}
	else
	{
		// Set the coupling blocks in the preconditioner
		M_prec->setCouplingBlocks ( M_coupling->lambdaToFluidMomentum(),
									M_coupling->lambdaToStructureMomentum(),
									M_coupling->structureDisplacementToLambda(),
									M_coupling->fluidVelocityToLambda(),
									M_coupling->structureDisplacementToFluidDisplacement() );
	}

	M_prec->setMonolithicMap ( M_monolithicMap );

	int time_step_count = 0;

	for ( ; M_time <= M_t_end + M_dt / 2.; M_time += M_dt)
	{
        if ( M_comm->MyPID()==0 )
        {
            M_outputLinearIterations << std::endl;
            M_outputPreconditionerComputation << std::endl;
            M_outputTimeLinearSolver << std::endl;
        }

        time_step_count += 1;

		M_displayer.leaderPrint ( "\n-----------------------------------\n" ) ;
		M_displayer.leaderPrintMax ( "FSI - solving now for time ", M_time ) ;
		M_displayer.leaderPrint ( "\n" ) ;
		iterChrono.start();

		updateSystem ( );

		if ( M_extrapolateInitialGuess && M_time == (M_t_zero + M_dt) )
		{
			M_displayer.leaderPrint ( "FSI - initializing extrapolation of initial guess\n" ) ;
			initializeExtrapolation ( );
		}

		if ( M_extrapolateInitialGuess )
		{
			M_displayer.leaderPrint ( "FSI - Extrapolating initial guess for Newton\n" ) ;
			M_extrapolationSolution->extrapolate (M_orderExtrapolationInitialGuess, *M_solution);
		}

		// Apply current BC to the solution vector
		applyBCsolution ( M_solution );

		M_maxiterNonlinear = 10;
        
		NonLinearRichardson ( *M_solution, *this, M_absoluteTolerance, M_relativeTolerance, M_maxiterNonlinear, M_etaMax, M_nonLinearLineSearch, 0, 2, M_out_res, M_time);

		iterChrono.stop();

		M_displayer.leaderPrint ( "\n" ) ;
		M_displayer.leaderPrintMax ( "FSI - timestep solved in ", iterChrono.diff() ) ;
		M_displayer.leaderPrint ( "-----------------------------------\n\n" ) ;

		// Writing the norms into a file
		if ( M_comm->MyPID()==0 )
		{
			M_outputTimeStep << M_time << ", " << iterChrono.diff() << std::endl;
		}

		iterChrono.reset();

		if ( M_extrapolateInitialGuess )
			M_extrapolationSolution->shift(*M_solution);

		// Export the solution obtained at the current timestep
		M_fluidVelocity->subset(*M_solution, M_fluid->uFESpace()->map(), 0, 0);
		M_fluidPressure->subset(*M_solution, M_fluid->pFESpace()->map(), M_fluid->uFESpace()->map().mapSize(), 0);
        
        vectorPtr_Type lagrange ( new vector_Type (*M_lagrangeMap, Unique) );
        vectorPtr_Type Lagrange;
        lagrange->zero();
        lagrange->subset(*M_solution, *M_lagrangeMap, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() +
                        M_displacementFESpace->map().mapSize(), 0);

        if ( !M_nonconforming )
        {
        	*M_Lagrange = ( *(M_coupling->lambdaToFluidMomentum()) * (*lagrange) );
        }
        else
        {

        	M_FluidToStructureInterpolant->expandGammaToOmega_Known( lagrange, Lagrange );
        	*M_Lagrange = *Lagrange;
        }

        M_fluidDisplacement->subset(*M_solution, M_aleFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize(), 0);
		M_structureDisplacement->subset(*M_solution, M_displacementFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize(), 0);

		// Updating all the time-advance objects
		M_fluidTimeAdvance->shift(*M_fluidVelocity);

		if ( M_useBDF )
		{
			M_structureTimeAdvanceBDF->shift(*M_structureDisplacement);
		}
		else
		{
			M_structureTimeAdvance->shift(M_structureDisplacement);
		}

		M_aleTimeAdvance->shift(*M_fluidDisplacement);

		if ( !M_useBDF )
		{
			*M_structureVelocity  = *M_structureTimeAdvance->old_first_derivative();

			*M_structureAcceleration = *M_structureTimeAdvance->old_second_derivative();
		}

		// This part below handles the exporter of the solution.
		// In particular, given a number of timesteps at which
		// we ask to export the solution (from datafile), here
		// the code takes care of exporting the solution also at
		// the previous timesteps such that, if later a restart
		// of the simulation is performed, it works correctly.

		if ( M_orderBDF == 1 )
		{
			if ( time_step_count == M_counterSaveEvery )
			{
				M_exporterFluid->postProcess(M_time);
				M_exporterStructure->postProcess(M_time);
				M_counterSaveEvery += M_saveEvery;
			}
		}
		else if ( M_orderBDF == 2 )
		{
			if ( time_step_count == (M_counterSaveEvery-1) && !M_disregardRestart )
			{
				M_exporterFluid->postProcess(M_time);
				M_exporterStructure->postProcess(M_time);
			}
			else if ( time_step_count == M_counterSaveEvery )
			{
				M_exporterFluid->postProcess(M_time);
				M_exporterStructure->postProcess(M_time);
				M_counterSaveEvery += M_saveEvery;
			}
		}
	}

	M_exporterFluid->closeFile();
	M_exporterStructure->closeFile();
}

void
FSIHandler::intializeTimeLoop ( )
{
	LifeChrono smallThingsChrono;

	// Clean unpartitioned mesh
	M_fluidMesh.reset();
	M_structureMesh.reset();

	// Build monolithic map
	buildMonolithicMap ( );

	M_solution.reset ( new VectorEpetra ( *M_monolithicMap ) );
	M_solution->zero();

	if ( M_restart )
	{
		M_solution.reset ( new VectorEpetra ( *M_monolithicMap ) );
		M_solution->zero();

		M_solution->subset(*M_fluidVelocity, M_fluid->uFESpace()->map(), 0, 0);
		M_solution->subset(*M_fluidPressure, M_fluid->pFESpace()->map(), 0, M_fluid->uFESpace()->map().mapSize());
		M_solution->subset(*M_structureDisplacement, M_displacementFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize());
        *M_Lagrange = ( *(M_coupling->fluidVelocityToLambda()) * (*M_LagrangeRestart) );
		M_solution->subset(*M_fluidDisplacement, M_aleFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize()+M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize() );
	}

	// Apply boundary conditions for the ale problem (the matrix will not change during the simulation)
	M_ale->applyBoundaryConditions ( *M_aleBC );

	// Apply boundary conditions for the structure problem (the matrix will not change during the simulation, it is linear elasticity)
	if ( M_linearElasticity )
	{
		getMatrixStructure ( );
		M_displayer.leaderPrint ( "\t Set and approximate structure block in the preconditioner.. " ) ;
		smallThingsChrono.start();
		M_prec->setStructureBlock ( M_matrixStructure );
		M_prec->updateApproximatedStructureMomentumOperator ( );
		smallThingsChrono.stop();
		M_displayer.leaderPrintMax ( "done in ", smallThingsChrono.diff() ) ;
	}

	// Set blocks for the preconditioners: geometry
	smallThingsChrono.reset();
	M_displayer.leaderPrint ( "\t Set and approximate geometry block in the preconditioner... " ) ;
	smallThingsChrono.start();
	M_prec->setGeometryBlock ( M_ale->matrix() );
	M_prec->updateApproximatedGeometryOperator ( );
	smallThingsChrono.stop();
	M_displayer.leaderPrintMax ( "done in ", smallThingsChrono.diff() ) ;

	if ( M_nonconforming )
	{
		M_prec->setCouplingOperators_nonconforming(M_FluidToStructureInterpolant, M_StructureToFluidInterpolant, M_lagrangeMap);
		if ( M_useMasses)
		{
			M_displayer.leaderPrint ( "[FSI] - Assemble interface mass of the structure for coupling of stresses\n" ) ;
			assembleStructureInterfaceMass ( );
		}
	}
	else
	{
		// Set the coupling blocks in the preconditioner
		M_prec->setCouplingBlocks ( M_coupling->lambdaToFluidMomentum(),
									M_coupling->lambdaToStructureMomentum(),
									M_coupling->structureDisplacementToLambda(),
									M_coupling->fluidVelocityToLambda(),
									M_coupling->structureDisplacementToFluidDisplacement() );
	}

	M_prec->setMonolithicMap ( M_monolithicMap );
}

void
FSIHandler::solveTimeStep ( )
{
	if ( M_comm->MyPID()==0 )
	{
		M_outputLinearIterations << std::endl;
		M_outputPreconditionerComputation << std::endl;
		M_outputTimeLinearSolver << std::endl;
	}

	M_displayer.leaderPrint ( "\n-----------------------------------\n" ) ;
	M_displayer.leaderPrintMax ( "FSI - solving now for time ", M_time ) ;
	M_displayer.leaderPrint ( "\n" ) ;
	LifeChrono iterChrono;
	iterChrono.start();

	updateSystem ( );

	if ( M_extrapolateInitialGuess && M_time == (M_t_zero + M_dt) )
	{
		M_displayer.leaderPrint ( "FSI - initializing extrapolation of initial guess\n" ) ;
		initializeExtrapolation ( );
	}

	if ( M_extrapolateInitialGuess )
	{
		M_displayer.leaderPrint ( "FSI - Extrapolating initial guess for Newton\n" ) ;
		M_extrapolationSolution->extrapolate (M_orderExtrapolationInitialGuess, *M_solution);
	}

	// Apply current BC to the solution vector
	applyBCsolution ( M_solution );

	M_maxiterNonlinear = 10;

	// Using the solution at the previous timestep as initial guess -> TODO: extrapolation
	NonLinearRichardson ( *M_solution, *this, M_absoluteTolerance, M_relativeTolerance, M_maxiterNonlinear, M_etaMax,
			M_nonLinearLineSearch, 0, 2, M_out_res, M_time);

	iterChrono.stop();
	M_displayer.leaderPrint ( "\n" ) ;
	M_displayer.leaderPrintMax ( "FSI - timestep solved in ", iterChrono.diff() ) ;
	M_displayer.leaderPrint ( "-----------------------------------\n\n" ) ;

	// Writing the norms into a file
	if ( M_comm->MyPID()==0 )
	{
		// M_outputTimeStep << "Time = " << M_time << " solved in " << iterChrono.diff() << " seconds" << std::endl;
		M_outputTimeStep << M_time << ", " << iterChrono.diff() << std::endl;
	}

	iterChrono.reset();

	if ( M_extrapolateInitialGuess )
		M_extrapolationSolution->shift(*M_solution);

	// Export the solution obtained at the current timestep
	M_fluidVelocity->subset(*M_solution, M_fluid->uFESpace()->map(), 0, 0);
	M_fluidPressure->subset(*M_solution, M_fluid->pFESpace()->map(), M_fluid->uFESpace()->map().mapSize(), 0);
    M_Lagrange->subset(*M_solution, *M_lagrangeMap, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() +
                        M_displacementFESpace->map().mapSize(), 0);
    M_fluidDisplacement->subset(*M_solution, M_aleFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize(), 0);
	M_structureDisplacement->subset(*M_solution, M_displacementFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize(), 0);

	// Updating all the time-advance objects
	M_fluidTimeAdvance->shift(*M_fluidVelocity);

	if ( M_useBDF )
	{
		M_structureTimeAdvanceBDF->shift(*M_structureDisplacement);
	}
	else
	{
		M_structureTimeAdvance->shift(M_structureDisplacement);
	}

	M_aleTimeAdvance->shift(*M_fluidDisplacement);

	if ( !M_useBDF )
	{
		*M_structureVelocity  = *M_structureTimeAdvance->old_first_derivative();

		*M_structureAcceleration = *M_structureTimeAdvance->old_second_derivative();
	}
}

void
FSIHandler::postprocessResults( const int& time_step_count )
{
	if ( M_orderBDF == 1 )
	{
		if ( time_step_count == M_counterSaveEvery )
		{
			M_exporterFluid->postProcess(M_time);
			M_exporterStructure->postProcess(M_time);
			M_counterSaveEvery += M_saveEvery;
		}
	}
	else if ( M_orderBDF == 2 )
	{
		if ( time_step_count == (M_counterSaveEvery-1) && !M_disregardRestart )
		{
			M_exporterFluid->postProcess(M_time);
			M_exporterStructure->postProcess(M_time);
		}
		else if ( time_step_count == M_counterSaveEvery )
		{
			M_exporterFluid->postProcess(M_time);
			M_exporterStructure->postProcess(M_time);
			M_counterSaveEvery += M_saveEvery;
		}
	}
}

void
FSIHandler::finalizeExporters ( )
{
	M_exporterFluid->closeFile();
	M_exporterStructure->closeFile();
}

void
FSIHandler::initializeExtrapolation( )
{
	// Initialize vector in the extrapolation object for the initial guess of Newton
	vector_Type solutionInitial ( *M_monolithicMap );
	std::vector<vector_Type> initialStateSolution;
	solutionInitial.zero();
	for ( UInt i = 0; i < M_orderExtrapolationInitialGuess; ++i )
		initialStateSolution.push_back(solutionInitial);

	M_extrapolationSolution.reset( new TimeAndExtrapolationHandler ( ) );
	M_extrapolationSolution->setBDForder(M_orderExtrapolationInitialGuess);
	M_extrapolationSolution->setMaximumExtrapolationOrder(M_orderExtrapolationInitialGuess);
	M_extrapolationSolution->initialize(initialStateSolution);
}

void
FSIHandler::applyBCsolution(vectorPtr_Type& M_solution)
{
	//! Extract each component of the input vector
	VectorEpetra velocity(M_fluid->uFESpace()->map(), Unique);
	velocity.subset(*M_solution, M_fluid->uFESpace()->map(), 0, 0);

	VectorEpetra displacement(M_displacementFESpace->map(), Unique);
	displacement.subset(*M_solution, M_displacementFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize(), 0 );

	VectorEpetra geometry(M_aleFESpace->map(), Unique);
	geometry.subset(*M_solution, M_aleFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() +
			                                   	   	  M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize(), 0 );

	//! Apply BC on each component
	if ( !M_fluidBC->bcUpdateDone() )
		M_fluidBC->bcUpdate ( *M_fluid->uFESpace()->mesh(), M_fluid->uFESpace()->feBd(), M_fluid->uFESpace()->dof() );

	bcManageRhs ( velocity, *M_fluid->uFESpace()->mesh(), M_fluid->uFESpace()->dof(), *M_fluidBC, M_fluid->uFESpace()->feBd(), 1.0, M_time );

	if ( !M_structureBC->bcUpdateDone() )
		M_structureBC->bcUpdate ( *M_displacementFESpace->mesh(), M_displacementFESpace->feBd(), M_displacementFESpace->dof() );

	bcManageRhs ( displacement, *M_displacementFESpace->mesh(), M_displacementFESpace->dof(), *M_structureBC, M_displacementFESpace->feBd(), 1.0, M_time );

	if ( !M_aleBC_residual_essential->bcUpdateDone() )
		M_aleBC_residual_essential->bcUpdate ( *M_aleFESpace->mesh(), M_aleFESpace->feBd(), M_aleFESpace->dof() );

	bcManageRhs ( geometry, *M_aleFESpace->mesh(), M_aleFESpace->dof(), *M_aleBC_residual_essential, M_aleFESpace->feBd(), 1.0, M_time );

	//! Push local contributions into the global one
	M_solution->subset(velocity, M_fluid->uFESpace()->map(), 0, 0);
	M_solution->subset(displacement, M_displacementFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() );
	M_solution->subset(geometry, M_aleFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() +
     	   	  	  	  	  	  	  	  	  	  	  	  	 M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize() );
}

void
FSIHandler::applyBCresidual(VectorEpetra& r_u, VectorEpetra& r_ds, VectorEpetra& r_df)
{
	//! Extract each component of the input vector
	VectorEpetra ru_copy(r_u, Unique);
	ru_copy *= -1;

	//! Apply BC on each component
	if ( !M_fluidBC_residual_natural->bcUpdateDone() )
		M_fluidBC_residual_natural->bcUpdate ( *M_fluid->uFESpace()->mesh(), M_fluid->uFESpace()->feBd(), M_fluid->uFESpace()->dof() );

	bcManageRhs ( ru_copy,
                 *M_fluid->uFESpace()->mesh(),
                 M_fluid->uFESpace()->dof(),
                 *M_fluidBC_residual_natural,
                 M_fluid->uFESpace()->feBd(),
                 1.0,
                 M_time );

	ru_copy *= -1;

    if ( !M_fluidBC_residual_essential->bcUpdateDone() )
        M_fluidBC_residual_essential->bcUpdate ( *M_fluid->uFESpace()->mesh(), M_fluid->uFESpace()->feBd(), M_fluid->uFESpace()->dof() );

    bcManageRhs ( ru_copy,
                 *M_fluid->uFESpace()->mesh(),
                 M_fluid->uFESpace()->dof(),
                 *M_fluidBC_residual_essential,
                 M_fluid->uFESpace()->feBd(),
                 0.0, // Homogeneous condition for Dirichlet type on the residual
                 M_time );
    
	r_u.zero();
	r_u = ru_copy;

	if ( !M_structureBC_residual_natural->bcUpdateDone() )
		M_structureBC_residual_natural->bcUpdate ( *M_displacementFESpace->mesh(), M_displacementFESpace->feBd(), M_displacementFESpace->dof() );

	if ( M_datafile("solid/robin_external_wall", false ) )
	{
		VectorEpetra rds_robin( M_displacementFESpace->map(), Unique);
		rds_robin.zero();
		rds_robin = *M_interface_mass_structure_robin * (*M_dsk);
		r_ds += rds_robin;
	}

	bcManageRhs ( r_ds,
                 *M_displacementFESpace->mesh(),
                 M_displacementFESpace->dof(),
                 *M_structureBC_residual_natural,
                 M_displacementFESpace->feBd(),
                 1.0,
                 M_time );

    if ( !M_structureBC_residual_essential->bcUpdateDone() )
        M_structureBC_residual_essential->bcUpdate ( *M_displacementFESpace->mesh(), M_displacementFESpace->feBd(), M_displacementFESpace->dof() );

    bcManageRhs ( r_ds,
                 *M_displacementFESpace->mesh(),
                 M_displacementFESpace->dof(),
                 *M_structureBC_residual_essential,
                 M_displacementFESpace->feBd(),
                 0.0, // Homogeneous condition for Dirichlet type on the residual
                 M_time );
    
	if ( !M_aleBC_residual_natural->bcUpdateDone() )
		M_aleBC_residual_natural->bcUpdate ( *M_aleFESpace->mesh(), M_aleFESpace->feBd(), M_aleFESpace->dof() );
    
	bcManageRhs ( r_df, *M_aleFESpace->mesh(), M_aleFESpace->dof(), *M_aleBC_residual_natural, M_aleFESpace->feBd(), 1.0, M_time );
    
    if ( !M_aleBC_residual_essential->bcUpdateDone() )
        M_aleBC_residual_essential->bcUpdate ( *M_aleFESpace->mesh(), M_aleFESpace->feBd(), M_aleFESpace->dof() );
    
    bcManageRhs ( r_df, *M_aleFESpace->mesh(), M_aleFESpace->dof(), *M_aleBC_residual_essential, M_aleFESpace->feBd(), 0.0, M_time );
}

void
FSIHandler::assembleStructureInterfaceMass()
{
	// INITIALIZE MATRIX WITH THE MAP OF THE INTERFACE
	matrixPtr_Type structure_interfaceMass( new matrix_Type ( M_displacementFESpace->map(), 50 ) );
	structure_interfaceMass->zero();

	// ASSEMBLE MASS MATRIX AT THE INTERFACE
	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
	{
		using namespace ExpressionAssembly;
		integrate ( boundary (M_displacementETFESpace->mesh(), M_datafile("interface/flag", 1) ),
					myBDQR,
					M_displacementETFESpace,
					M_displacementETFESpace,
					//Boundary Mass
					dot(phi_i,phi_j)
		)
		>> structure_interfaceMass;
	}
	structure_interfaceMass->globalAssemble();

	mapPtr_Type M_mapStuctureGammaVectorial;
	M_StructureToFluidInterpolant->getVectorialInterpolationMap(M_mapStuctureGammaVectorial);

	// RESTRICT MATRIX TO INTERFACE DOFS ONLY
	M_interface_mass_structure.reset(new matrix_Type ( *M_mapStuctureGammaVectorial, 50 ) );
	structure_interfaceMass->restrict ( M_mapStuctureGammaVectorial, M_numerationInterfaceStructure,
										M_structureDisplacement->size()/3, M_interface_mass_structure );
}

void
FSIHandler::applyInverseFluidMassOnGamma ( const vectorPtr_Type& weakLambda, vectorPtr_Type& strongLambda )
{
	Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
	belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList_rbf3d.xml" );

	// Preconditioner
	if ( !M_precPtrBuilt )
	{
		prec_Type* precRawPtr;
		precRawPtr = new prec_Type;
		precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
		M_precPtr.reset ( precRawPtr );
		M_precPtrBuilt = true;
	}

	// Linear Solver
	LinearSolver solverOne;
	solverOne.setCommunicator ( M_comm );
	solverOne.setParameters ( *belosList );
	solverOne.setPreconditioner ( M_precPtr );

	// Solution
	strongLambda->zero();

	// Solve system
	solverOne.setOperator ( M_interface_mass_fluid );
	solverOne.setRightHandSide ( weakLambda );
	solverOne.solve ( strongLambda );
}

void
FSIHandler::evalResidual(vector_Type& residual, const vector_Type& solution, const UInt iter_newton)
{
	residual.zero();
	M_NewtonIter = iter_newton;

	// Create empty vectors to store each component of the FSI residual

	vectorPtr_Type res_u ( new vector_Type ( M_fluid->uFESpace()->map(), Unique ) );
	vectorPtr_Type res_p ( new vector_Type ( M_fluid->pFESpace()->map(), Unique ) );
	vectorPtr_Type res_ds ( new vector_Type ( M_displacementFESpace->map(), Unique ) );
	vectorPtr_Type res_lambda ( new vector_Type ( *M_lagrangeMap, Unique ) );
	vectorPtr_Type res_df ( new vector_Type (M_aleFESpace->map(), Unique ) );

	res_u->zero();
	res_p->zero();
	res_ds->zero();
	res_lambda->zero();
	res_df->zero();

	// Create empty vectors to store each component of the FSI solution from previous newton step

	vectorPtr_Type u_k ( new vector_Type ( M_fluid->uFESpace()->map() ) );
	vectorPtr_Type p_k ( new vector_Type ( M_fluid->pFESpace()->map() ) );
	vectorPtr_Type ds_k ( new vector_Type ( M_displacementFESpace->map() ) );
	vectorPtr_Type lambda_k ( new vector_Type ( *M_lagrangeMap ) );
	vectorPtr_Type df_k ( new vector_Type ( M_aleFESpace->map() ) );

	u_k->zero();
	p_k->zero();
	ds_k->zero();
	lambda_k->zero();
	df_k->zero();

	// Copy in the vectors u_k, p_k, ds_k, lambda_k and df_k each component of the solution at the previous newton step k

	u_k->subset (solution,  M_fluid->uFESpace()->map(), 0, 0);
	p_k->subset (solution,  M_fluid->pFESpace()->map(), M_fluid->uFESpace()->map().mapSize(), 0);
	ds_k->subset (solution,  M_displacementFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize(), 0);
	lambda_k->subset(solution,  *M_lagrangeMap, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize(), 0);
	df_k->subset(solution,  M_aleFESpace->map(), M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize(), 0);

	M_dsk.reset( new vector_Type ( *ds_k ) );

	//---------------//
	// Move the mesh //
	//---------------//

	vectorPtr_Type mmRep ( new vector_Type (*df_k, Repeated ) );
	moveMesh ( *mmRep );

	//---------------------------------//
	// Compute the fluid mesh velocity //
	//---------------------------------//

    vectorPtr_Type meshVelocity ( new vector_Type (M_aleFESpace->map() ) );
    meshVelocity->zero();
    vector_Type meshVelocity_bdf ( M_aleFESpace->map() ); // velocity from previous timesteps
    meshVelocity_bdf.zero();
    M_aleTimeAdvance->rhsContribution( meshVelocity_bdf );
    *meshVelocity = ( M_aleTimeAdvance->alpha()*(*df_k)/(M_dt) - meshVelocity_bdf );

    //---------------------------------//
    // Compute the convective velocity //
    //---------------------------------//

	M_beta.reset ( new VectorEpetra ( M_fluid->uFESpace()->map ( ) ) );
	M_beta->zero();
	*M_beta = ( *u_k - *meshVelocity );

	//--------------------------------------------------------------------//
	// Re-assemble the fluid constant blocks since we have moved the mesh //
	//--------------------------------------------------------------------//

	M_fluid->buildSystem();
	M_fluid->setupPostProc();

	if ( M_nonconforming )
	{
		//----------------------------------------------//
		// Residual, part 1: compute the fluid residual //
		//----------------------------------------------//

		M_displayer.leaderPrint ( "[FSI] - Computing residual of the fluid \n" ) ;

		// Evaluate the residual coming from the fluid
		M_fluid->evaluateResidual( M_beta, u_k, p_k, M_rhs_velocity, res_u, res_p);

		vectorPtr_Type lambda_km1_omegaF ( new vector_Type ( M_fluid->uFESpace()->map() ) );
		M_FluidToStructureInterpolant->expandGammaToOmega_Known(lambda_k, lambda_km1_omegaF);

		*res_u += *lambda_km1_omegaF;

		//--------------------------------------------------//
		// Residual, part 2: compute the structure residual //
		//--------------------------------------------------//

		M_displayer.leaderPrint ( "[FSI] - Computing residual structure \n" ) ;

		// Vector containing the residual of the structural part
		res_ds->zero();

        if ( M_linearElasticity )
        {
            if ( M_useBDF )
            {
                vectorPtr_Type old_second_der_terms ( new vector_Type ( M_displacementFESpace->map() ) );
                old_second_der_terms->zero();
                M_structureTimeAdvanceBDF->second_der_old_dts(*old_second_der_terms);

                *res_ds = ( ( ( 1.0 / ( M_dt * M_dt ) * M_structureTimeAdvanceBDF->massCoefficient() ) * ( (*M_structure->mass_matrix_no_bc() ) * (*ds_k) ) ) +
                          ( (*M_structure->stiffness_matrix_no_bc() ) * (*ds_k) ) +
                          ( ( 1.0 / ( M_dt * M_dt ) ) * ( ( *M_structure->mass_matrix_no_bc() ) * ( *old_second_der_terms ) ) ) );
            }
            else
            {
                *res_ds = ( ( ( 1.0/ ( M_dt * M_dt * M_structureTimeAdvance->get_beta() ) ) * ( (*M_structure->mass_matrix_no_bc() ) * (*ds_k) ) ) +
                          ( (*M_structure->stiffness_matrix_no_bc() ) * (*ds_k) ) -
                          ( (*M_structure->mass_matrix_no_bc() ) * (*M_structureTimeAdvance->get_csi()) ) );
            }
        }
        else
        {
            if ( M_useBDF )
            {
                Real coefficient = 1.0/ ( M_dt * M_dt ) * M_structureTimeAdvanceBDF->massCoefficient();
                vectorPtr_Type old_second_der_terms ( new vector_Type ( M_displacementFESpace->map() ) );
                old_second_der_terms->zero();
                M_structureTimeAdvanceBDF->second_der_old_dts(*old_second_der_terms);
                *old_second_der_terms *= ( 1.0 / ( M_dt * M_dt ) );
                M_structureNeoHookean->evaluate_residual(ds_k, coefficient, old_second_der_terms, res_ds );
            }
            else
            {
                Real coefficient = 1.0/ ( M_dt * M_dt * M_structureTimeAdvance->get_beta() );
                M_structureNeoHookean->evaluate_residual(ds_k, coefficient, M_structureTimeAdvance->get_csi(), res_ds );
            }

        }

        if ( M_useMasses )
		{
			M_displayer.leaderPrint ( "[FSI] - Assemble interface mass of the fluid \n" ) ;
			M_fluid->assembleInterfaceMass ( M_interface_mass_fluid, M_lagrangeMap, M_datafile("interface/flag", 1),
											 M_numerationInterfaceFluid, M_fluid->uFESpace()->map().mapSize()/3 );

			M_displayer.leaderPrint ( "[FSI] - From weak to strong residual - fluid \n" ) ;
			vectorPtr_Type tmp_gamma_f ( new vector_Type ( M_lagrangeMap ) );
			tmp_gamma_f->zero();
			applyInverseFluidMassOnGamma( lambda_k, tmp_gamma_f );

			M_displayer.leaderPrint ( "[FSI] - Interpolating strong residual \n" ) ;
			vectorPtr_Type tmp_omega_f ( new vector_Type ( M_fluid->uFESpace()->map() ) );
			tmp_omega_f->zero ( );

			M_FluidToStructureInterpolant->expandGammaToOmega_Known( tmp_gamma_f, tmp_omega_f );

			M_FluidToStructureInterpolant->updateRhs(tmp_omega_f);

			M_FluidToStructureInterpolant->interpolate();

			vectorPtr_Type tmp_omega_s ( new vector_Type ( M_displacementFESpace->map() ) );

			tmp_omega_s->zero();

			M_FluidToStructureInterpolant->solution(tmp_omega_s);

			vectorPtr_Type tmp_gamma_s;

			M_StructureToFluidInterpolant->restrictOmegaToGamma_Known(tmp_omega_s, tmp_gamma_s);

			vectorPtr_Type out_gamma_s ( new vector_Type ( tmp_gamma_s->map( ) ) );
			out_gamma_s->zero();

			M_displayer.leaderPrint ( "[FSI] - From strong to weak residual - structure \n" ) ;
			*out_gamma_s = (*M_interface_mass_structure) * ( *tmp_gamma_s );

			vectorPtr_Type structure_weak_residual ( new vector_Type ( M_displacementFESpace->map() ) );
			structure_weak_residual->zero();

			M_StructureToFluidInterpolant->expandGammaToOmega_Known( out_gamma_s, structure_weak_residual );

			*res_ds -= *structure_weak_residual;
		}
		else
		{
			M_FluidToStructureInterpolant->updateRhs(lambda_km1_omegaF);

			M_FluidToStructureInterpolant->interpolate();

			vectorPtr_Type structure_weak_residual ( new vector_Type ( M_displacementFESpace->map() ) );
			structure_weak_residual->zero();

			M_FluidToStructureInterpolant->solution(structure_weak_residual);

			*res_ds -= *structure_weak_residual;
		}

		//-------------------------------------------------//
		// Residual, part 3: compute the coupling residual //
		//-------------------------------------------------//

		M_displayer.leaderPrint ( "[FSI] - Computing residual coupling velocities \n" ) ;

		res_lambda->zero();

		vectorPtr_Type velocity_km1_gamma( new vector_Type ( *M_lagrangeMap ) );
		velocity_km1_gamma->zero();
		M_FluidToStructureInterpolant->restrictOmegaToGamma_Known(u_k, velocity_km1_gamma);

		vectorPtr_Type structure_vel( new vector_Type ( M_displacementFESpace->map() ) );
		structure_vel->zero();

        if ( M_useBDF )
        {
            *structure_vel = ( ( M_structureTimeAdvanceBDF->coefficientFirstDerivative()/ M_dt ) * (*ds_k) );
        }
        else
        {
            *structure_vel = ( ( M_structureTimeAdvance->get_gamma()/(M_dt * M_structureTimeAdvance->get_beta() ) ) * (*ds_k) );
        }
        
		vectorPtr_Type res_couplingVel_omega_f ( new vector_Type ( M_fluid->uFESpace()->map() ) );
		res_couplingVel_omega_f->zero();

		M_StructureToFluidInterpolant->updateRhs ( structure_vel );
		M_StructureToFluidInterpolant->interpolate();
		M_StructureToFluidInterpolant->solution(res_couplingVel_omega_f);

		vectorPtr_Type res_couplingVel_gamma_f ( new vector_Type ( *M_lagrangeMap ) );
		res_couplingVel_gamma_f->zero();

		M_FluidToStructureInterpolant->restrictOmegaToGamma_Known(res_couplingVel_omega_f, res_couplingVel_gamma_f);

		*res_lambda = (*velocity_km1_gamma - *res_couplingVel_gamma_f + *M_rhsCouplingVelocities );

		//--------------------------------------------//
		// Residual, part 4: compute the ALE residual //
		//--------------------------------------------//

		M_displayer.leaderPrint ( "[FSI] - Computing residual ALE \n" ) ;

		M_StructureToFluidInterpolant->updateRhs ( ds_k );

		M_StructureToFluidInterpolant->interpolate();

		vectorPtr_Type res_ds_ale ( new vector_Type ( M_fluid->uFESpace()->map() ) );

		res_ds_ale->zero();

		M_StructureToFluidInterpolant->solution(res_ds_ale);

		vectorPtr_Type res_ale_ale ( new vector_Type ( M_aleFESpace->map() ) );

		res_ale_ale->zero();

		*res_ale_ale = ( *M_ale->matrix ( ) ) * ( *df_k );

		res_df->zero();

		*res_df -= *res_ds_ale;

		*res_df += *res_ale_ale;

		//---------------------------------------//
		// Residual, part 5: monolithic residual //
		//---------------------------------------//

		M_displayer.leaderPrint ( "[FSI] - Gathering all components of the residual  \n" ) ;

	}
	else
	{
		//------------------------//
		// Residual for the fluid //
		//------------------------//

        // Evaluate the residual coming from the fluid ( without coupling for the moment )
		M_fluid->evaluateResidual( M_beta, u_k, p_k, M_rhs_velocity, res_u, res_p);

        // Adding the coupling contribution
		*res_u += ( *(M_coupling->lambdaToFluidMomentum()) * (*lambda_k) );

		//------------------------//
		// Residual for the solid //
		//------------------------//

		if ( M_linearElasticity )
		{
			if ( M_useBDF )
			{
				vectorPtr_Type old_second_der_terms ( new vector_Type ( M_displacementFESpace->map() ) );
				old_second_der_terms->zero();
				M_structureTimeAdvanceBDF->second_der_old_dts(*old_second_der_terms);

				*res_ds = ( ( ( 1.0 / ( M_dt * M_dt ) * M_structureTimeAdvanceBDF->massCoefficient() ) * ( (*M_structure->mass_matrix_no_bc() ) * (*ds_k) ) ) +
						    ( (*M_structure->stiffness_matrix_no_bc() ) * (*ds_k) ) +
						    ( (*M_coupling->lambdaToStructureMomentum() ) * (*lambda_k) ) +
						    ( ( 1.0 / ( M_dt * M_dt ) ) * ( ( *M_structure->mass_matrix_no_bc() ) * ( *old_second_der_terms ) ) ) );
			}
			else
			{
				*res_ds = ( ( ( 1.0/ ( M_dt * M_dt * M_structureTimeAdvance->get_beta() ) ) * ( (*M_structure->mass_matrix_no_bc() ) * (*ds_k) ) ) +
							( (*M_structure->stiffness_matrix_no_bc() ) * (*ds_k) ) +
							( (*M_coupling->lambdaToStructureMomentum() ) * (*lambda_k) ) -
							( (*M_structure->mass_matrix_no_bc() ) * (*M_structureTimeAdvance->get_csi()) ) );
			}
		}
		else
		{
			if ( M_useBDF )
			{
				Real coefficient = 1.0/ ( M_dt * M_dt ) * M_structureTimeAdvanceBDF->massCoefficient();
				vectorPtr_Type old_second_der_terms ( new vector_Type ( M_displacementFESpace->map() ) );
				old_second_der_terms->zero();
				M_structureTimeAdvanceBDF->second_der_old_dts(*old_second_der_terms);
				*old_second_der_terms *= ( 1.0 / ( M_dt * M_dt ) );
				M_structureNeoHookean->evaluate_residual(ds_k, coefficient, old_second_der_terms, res_ds );
			}
			else
			{
				Real coefficient = 1.0/ ( M_dt * M_dt * M_structureTimeAdvance->get_beta() );
				M_structureNeoHookean->evaluate_residual(ds_k, coefficient, M_structureTimeAdvance->get_csi(), res_ds );
			}
			*res_ds += ( (*M_coupling->lambdaToStructureMomentum() ) * (*lambda_k) );
		}
		//---------------------------------------//
		// Residual for the lagrange multipliers //
		//---------------------------------------//

		*res_lambda = ( ( (*M_coupling->fluidVelocityToLambda()) * (*u_k) ) +
						( (*M_coupling->structureDisplacementToLambda()) * (*ds_k) ) +
						(  *M_rhsCouplingVelocities )  );

		//----------------------//
		// Residual for the ALE //
		//----------------------//

		*res_df = ( ( (*M_coupling->structureDisplacementToFluidDisplacement()) * (*ds_k) ) +
					( (*M_ale->matrix()) * (*df_k) ) );

	}

	applyBCresidual ( *res_u, *res_ds, *res_df );

	residual.subset(*res_u, M_fluid->uFESpace()->map(), 0, 0 );
	residual.subset(*res_p, M_fluid->pFESpace()->map(), 0, M_fluid->uFESpace()->map().mapSize() );
	residual.subset(*res_ds, M_displacementFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() );
	residual.subset(*res_lambda, *M_lagrangeMap, 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() );
	residual.subset(*res_df, M_aleFESpace->map(), 0, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize() );

	if (M_printResiduals)
	{
		// Computing the norms
		Real norm_full_residual = residual.normInf();
		Real norm_u_residual = res_u->normInf();
		Real norm_p_residual = res_p->normInf();
		Real norm_ds_residual = res_ds->normInf();
		Real norm_lambda_residual = res_lambda->normInf();
		Real norm_df_residual = res_df->normInf();

		// Writing the norms into a file
		if ( M_comm->MyPID()==0 && iter_newton == 0)
		{
			M_outputResiduals << "------------------------------------------------------------------------------------" << std::endl;
			M_outputResiduals << "# time = " << M_time << std::endl;
			M_outputResiduals << "Initial norms: " << std::endl;
			M_outputResiduals << "||r|| =" << norm_full_residual << ", ||r_u|| =" << norm_u_residual << ", ||r_p|| =" << norm_p_residual;
			M_outputResiduals << ", ||r_ds|| =" << norm_ds_residual << ", ||r_lambda|| =" << norm_lambda_residual << ", ||r_df|| =" << norm_df_residual <<  std::endl;
			M_outputResiduals << "iter    ||r||    ||r_u||        ||r_p||    ||r_ds||      ||r_lambda||    ||r_df||" << std::endl;
		}
		else if ( M_comm->MyPID()==0 && iter_newton > 0 )
		{
			M_outputResiduals << iter_newton << "  " << norm_full_residual << "  " << norm_u_residual << "  " << norm_p_residual << "  " <<
					norm_ds_residual << "  " << norm_lambda_residual << "  " << norm_df_residual << std::endl;
		}
	}

	//---------------------------------------------//
	// Post-step: update the Jacobian of the fluid //
	//---------------------------------------------//

	M_displayer.leaderPrint ( "[FSI] - Update Jacobian terms: \n" ) ;
	M_fluid->updateConvectiveTerm(M_beta);
	M_fluid->updateJacobian(*u_k);
	if ( M_fluid->useStabilization() )
		M_fluid->updateStabilization(*M_beta, *u_k, *p_k, *M_rhs_velocity);

	M_fluid->applyBoundaryConditionsJacobian ( M_fluidBC );

	if ( !M_linearElasticity )
	{
		M_displayer.leaderPrint ( "\nAssembly Jacobian structure" ) ;
		M_matrixStructure.reset (new matrix_Type ( M_displacementFESpace->map() ) );
		M_matrixStructure->zero();

		Real coefficient;
		if ( M_useBDF )
		{
			coefficient = 1.0/ ( M_dt * M_dt ) * M_structureTimeAdvanceBDF->massCoefficient();
		}
		else
		{
			coefficient = 1.0/ ( M_dt * M_dt * M_structureTimeAdvance->get_beta() );
		}

		M_structureNeoHookean->update_jacobian( ds_k, coefficient, M_matrixStructure);
		bcManageMatrix( *M_matrixStructure, *M_displacementFESpace->mesh(), M_displacementFESpace->dof(), *M_structureBC_residual_essential,
				        M_displacementFESpace->feBd(), 1.0, M_time);
	}

	//----------------------------------------------------//
	// Post-step bis: Compute the shape derivatives block //
	//----------------------------------------------------//

	if (M_useShapeDerivatives)
	{

		Real alpha = M_fluidTimeAdvance->alpha() / M_dt;
		Real density = M_fluid->density();
		Real viscosity = M_fluid->viscosity();

		vector_Type un (M_fluid->uFESpace()->map() );
		vector_Type uk (M_fluid->uFESpace()->map() + M_fluid->pFESpace()->map() );
		//vector_Type pk (M_fluid->pFESpace()->map() );

		vector_Type meshVelRep (  *meshVelocity, Repeated ) ;

		//When this class is used, the convective term is used implictly
		un.subset ( solution, 0 );

		uk.subset ( solution, 0 );
		//vector_Type veloFluidMesh ( M_uFESpace->map(), Repeated );
		//this->transferMeshMotionOnFluid ( meshVelRep, veloFluidMesh );


		// Simone check with Davide
		M_ale->updateShapeDerivatives ( alpha,
										density,
										viscosity,
										un,
										uk,
										// pk
										meshVelRep, // or veloFluidMesh
										*M_fluid->uFESpace(),
										*M_fluid->pFESpace(),
										true, //This flag tells the method to consider the velocity of the domain implicitly,
										true, //This flag tells the method to consider the convective term implicitly ,
										*M_fluidBC);
	}

	//------------------------------------------------------------//
	// Post-step tris: update the apply operator for the Jacobian //
	//------------------------------------------------------------//

	if ( !M_nonconforming )
		initializeApplyOperatorJacobian();
}

void
FSIHandler::solveJac( vector_Type& increment, const vector_Type& residual, const Real /*linearRelTol*/ )
{
	increment.zero();

	if ( M_nonconforming )
	{
		M_applyOperatorJacobianNonConforming->setComm ( M_comm );
		M_applyOperatorJacobianNonConforming->setTimeStep ( M_dt );
		M_applyOperatorJacobianNonConforming->setDatafile ( M_datafile );
		M_applyOperatorJacobianNonConforming->setMonolithicMap ( M_monolithicMap );
		M_applyOperatorJacobianNonConforming->setMaps ( M_fluid->uFESpace()->mapPtr(),
														M_fluid->pFESpace()->mapPtr(),
														M_displacementFESpace->mapPtr(),
														M_lagrangeMap,
														M_aleFESpace->mapPtr());

		M_applyOperatorJacobianNonConforming->setInterfaceMassMatrices (  M_interface_mass_fluid, M_interface_mass_structure );

		if ( M_fluid->useStabilization() )
			M_applyOperatorJacobianNonConforming->setFluidBlocks ( M_fluid->block00(), M_fluid->block01(), M_fluid->block10(), M_fluid->block11());
		else
			M_applyOperatorJacobianNonConforming->setFluidBlocks ( M_fluid->block00(), M_fluid->block01(), M_fluid->block10() );

		M_applyOperatorJacobianNonConforming->setStructureBlock ( M_matrixStructure );

		M_applyOperatorJacobianNonConforming->setALEBlock ( M_ale->matrix() );

		M_applyOperatorJacobianNonConforming->setUseShapeDerivatives ( M_useShapeDerivatives );

		M_applyOperatorJacobianNonConforming->setTimeStep(M_dt);

		M_applyOperatorJacobianNonConforming->setInterpolants ( M_FluidToStructureInterpolant, M_StructureToFluidInterpolant, M_useMasses );

        if ( M_useBDF )
        {
            M_applyOperatorJacobianNonConforming->setBDFcoeff ( M_structureTimeAdvanceBDF->coefficientFirstDerivative() );
        }
        else
        {
            M_applyOperatorJacobianNonConforming->setGamma ( M_structureTimeAdvance->get_gamma() );
            
            M_applyOperatorJacobianNonConforming->setBeta ( M_structureTimeAdvance->get_beta() );
        }
        
		if ( M_useShapeDerivatives )
		{
			M_applyOperatorJacobianNonConforming->setShapeDerivativesBlocks ( M_ale->shapeDerivativesVelocity(),
																			  M_ale->shapeDerivativesPressure() );
		}

		M_invOper->setOperator(M_applyOperatorJacobianNonConforming);
	}
	else
	{
		M_invOper->setOperator(M_applyOperatorJacobian);
	}

	//---------------------------------------------------//
	// First: set the fluid blocks in the preconditioner //
	//---------------------------------------------------//

	if ( !M_fluid->useStabilization() )
		M_prec->setFluidBlocks( M_fluid->block00(), M_fluid->block01(), M_fluid->block10() );
	else
		M_prec->setFluidBlocks( M_fluid->block00(), M_fluid->block01(), M_fluid->block10(), M_fluid->block11() );

	if (M_useShapeDerivatives)
	{
		M_prec->setShapeDerivativesBlocks(M_ale->shapeDerivativesVelocity(), M_ale->shapeDerivativesPressure());
	}

	if ( !M_linearElasticity )
	{
		LifeChrono smallThingsChrono;
		M_displayer.leaderPrint ( "\t Set and approximate structure block in the preconditioner.. " ) ;
		smallThingsChrono.start();
		M_prec->setStructureBlock ( M_matrixStructure );
		M_prec->updateApproximatedStructureMomentumOperator ( );
		smallThingsChrono.stop();
		M_displayer.leaderPrintMax ( "done in ", smallThingsChrono.diff() ) ;
	}

    if ( M_nonconforming)
    {
        if ( M_useBDF )
        {
            M_prec->setBDFcoeff ( M_structureTimeAdvanceBDF->coefficientFirstDerivative() );
        }
        else
        {
            M_prec->setGamma ( M_structureTimeAdvance->get_gamma() );
            M_prec->setBeta ( M_structureTimeAdvance->get_beta() );
        }
        M_prec->setVelocityFESpace ( M_fluid->uFESpace() );
    	M_prec->setBCInterface ( M_interfaceFluidBC );
    	M_prec->setTimeStep ( M_dt );
    	M_prec->setMonolithicMap ( M_monolithicMap );
    }
    else
    {
    	M_prec->setVelocityFESpace ( M_fluid->uFESpace() );
    	M_prec->setBCInterface ( M_interfaceFluidBC );
        M_prec->setDomainMap(M_applyOperatorJacobian->OperatorDomainBlockMapPtr());
        M_prec->setRangeMap(M_applyOperatorJacobian->OperatorRangeBlockMapPtr());
    }

	//--------------------------------------------------------------------------------//
	// Second: update the operators associated to shur complements and fluid momentum //
	//--------------------------------------------------------------------------------//

	LifeChrono smallThingsChrono;
	M_displayer.leaderPrint ( "\n Set preconditioner for the fluid momentum and the shur complements\n" ) ;
	M_displayer.leaderPrint ( "\t Set and approximate fluid momentum in the preconditioner.. " ) ;
	smallThingsChrono.start();

    // To be changed
    M_prec->updateApproximatedFluidOperator();
	smallThingsChrono.stop();
	M_displayer.leaderPrintMax ( "done in ", smallThingsChrono.diff() ) ;

	if ( M_comm->MyPID()==0 )
	{
		M_outputPreconditionerComputation << " " << smallThingsChrono.diff();
	}

	//------------------------------------------------------//
	// Third: set the preconditioner of the jacobian system //
	//------------------------------------------------------//

    M_invOper->setPreconditioner(M_prec);

	//-------------------------//
	// Forth: solve the system //
	//-------------------------//

	M_invOper->ApplyInverse(residual.epetraVector() , increment.epetraVector());

    if ( M_comm->MyPID()==0 )
    {
        M_outputLinearIterations << " " << M_invOper->NumIter();
        M_outputTimeLinearSolver << " " << M_invOper->TimeSolver();
    }

	M_displayer.leaderPrint (" FSI-  End of solve Jac ...                      ");

	if (M_printSteps)
	{
		// Defining vectors for the single components of the residual
		vectorPtr_Type s_fluid_vel ( new vector_Type (M_fluid->uFESpace()->map() ) );
		vectorPtr_Type s_fluid_press ( new vector_Type (M_fluid->pFESpace()->map() ) );
		vectorPtr_Type s_structure_displ ( new vector_Type (M_displacementFESpace->map() ) );
		vectorPtr_Type s_coupling ( new vector_Type (*M_lagrangeMap ) );
		vectorPtr_Type s_fluid_displ ( new vector_Type (M_aleFESpace->map() ) );

		// Getting each single contribution
		s_fluid_vel->subset(increment);
		s_fluid_press->subset(increment, M_fluid->uFESpace()->dof().numTotalDof() * 3);
		s_structure_displ->subset(increment, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() );
		s_coupling->subset(increment, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() );
		s_fluid_displ->subset(increment, M_fluid->uFESpace()->map().mapSize() + M_fluid->pFESpace()->map().mapSize() + M_displacementFESpace->map().mapSize() + M_lagrangeMap->mapSize() );

		// Computing the norms
		Real norm_full_step = increment.normInf();
		Real norm_u_step = s_fluid_vel->normInf();
		Real norm_p_step = s_fluid_press->normInf();
		Real norm_ds_step = s_structure_displ->normInf();
		Real norm_lambda_step = s_coupling->normInf();
		Real norm_df_step = s_fluid_displ->normInf();

		// Writing the norms into a file
		if ( M_comm->MyPID()==0 && M_NewtonIter == 0)
		{
			M_outputSteps << "------------------------------------------------------------------------------------" << std::endl;
			M_outputSteps << "# time = " << M_time << std::endl;
			M_outputSteps << "iter    ||s||        ||s_u||        ||s_p||        ||s_ds||     ||s_lambda||      ||s_df||" << std::endl;
			M_outputSteps << M_NewtonIter+1 << std::setw (15) << norm_full_step << std::setw (15) << norm_u_step << std::setw (15) << norm_p_step << std::setw (15) <<
					norm_ds_step << std::setw (15) << norm_lambda_step << std::setw (15) << norm_df_step << std::endl;
		}
		else if ( M_comm->MyPID()==0 && M_NewtonIter > 0 )
		{
			M_outputSteps << M_NewtonIter+1 << std::setw (15) << norm_full_step << std::setw (15) << norm_u_step << std::setw (15) << norm_p_step << std::setw (15) <<
					norm_ds_step << std::setw (15) << norm_lambda_step << std::setw (15) << norm_df_step << std::endl;
		}
	}

}

void
FSIHandler::moveMesh ( const VectorEpetra& displacement )
{
    M_displayer.leaderPrint (" FSI-  Moving the mesh ...                      ");
    M_fluidLocalMesh->meshTransformer().moveMesh (displacement,  M_aleFESpace->dof().numTotalDof() );
    M_displayer.leaderPrint ( "done\n" );
}

void
FSIHandler::updateSystem ( )
{
	M_rhs_velocity.reset ( new VectorEpetra ( M_fluid->uFESpace()->map ( ) ) );
	M_rhs_velocity->zero();

	// Compute rhs contribution due to the time derivative
	M_fluidTimeAdvance->rhsContribution (*M_rhs_velocity);

	// Get the vector needed for the coupling of the velocities
	get_structure_coupling_velocities ( );

	if ( M_nonconforming )
		updateRhsCouplingVelocities_nonconforming ( );
	else
		updateRhsCouplingVelocities ( );
}

void
FSIHandler::initializeApplyOperatorJacobian ( )
{
	Operators::FSIApplyOperator::operatorPtrContainer_Type operDataJacobian(5,5);
	operDataJacobian(0,0) = M_fluid->block00()->matrixPtr();
	operDataJacobian(0,1) = M_fluid->block01()->matrixPtr();
    operDataJacobian(0,3) = M_coupling->lambdaToFluidMomentum()->matrixPtr();
	operDataJacobian(1,0) = M_fluid->block10()->matrixPtr();
	if ( M_fluid->useStabilization() )
		operDataJacobian(1,1) = M_fluid->block11()->matrixPtr();
	operDataJacobian(2,2) = M_matrixStructure->matrixPtr();
	operDataJacobian(2,3) = M_coupling->lambdaToStructureMomentum()->matrixPtr();
	operDataJacobian(3,0) = M_coupling->fluidVelocityToLambda()->matrixPtr();
	operDataJacobian(3,2) = M_coupling->structureDisplacementToLambda()->matrixPtr();
	operDataJacobian(4,2) = M_coupling->structureDisplacementToFluidDisplacement()->matrixPtr();
	operDataJacobian(4,4) = M_ale->matrix()->matrixPtr();

    if (M_useShapeDerivatives)
    {
        operDataJacobian(0,4) = M_ale->shapeDerivativesVelocity()->matrixPtr();  // shape derivatives
        operDataJacobian(1,4) = M_ale->shapeDerivativesPressure()->matrixPtr();  // shape derivatives
    }
    M_applyOperatorJacobian->setMonolithicMap(M_monolithicMap);
	M_applyOperatorJacobian->setUp(operDataJacobian, M_comm);
}

}
