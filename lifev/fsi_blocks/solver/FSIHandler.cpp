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
M_applyOperator(new Operators::FSIApplyOperator)
{
}

FSIHandler::~FSIHandler( )
{
}

void
FSIHandler::setDatafile( const GetPot& dataFile)
{
    M_datafile = dataFile;
}
    
void
FSIHandler::readMeshes( )
{
    M_fluidMesh.reset ( new mesh_Type ( M_comm ) );
    M_meshDataFluid.reset( new MeshData ( ) );
    M_meshDataFluid->setup (M_datafile, "fluid/space_discretization");
    readMesh (*M_fluidMesh, *M_meshDataFluid);
    
    M_structureMesh.reset ( new mesh_Type ( M_comm ) );
    M_meshDataStructure.reset( new MeshData ( ) );
    M_meshDataStructure->setup (M_datafile, "solid/space_discretization");
    readMesh (*M_structureMesh, *M_meshDataStructure);
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
    
void FSIHandler::setup ( )
{
	// Fluid
	M_fluid.reset ( new NavierStokesSolver ( M_datafile, M_comm ) );
	M_fluid->setup ( M_fluidLocalMesh );

	// Structure data
	M_dataStructure.reset ( new StructuralConstitutiveLawData ( ) );
	M_dataStructure->setup ( M_datafile );

	// This beacuse the structural solver requires that the FESpaces are given from outside
	createStructureFESpaces();

	// This beacuse the ale solver requires that the FESpace is given from outside
	createAleFESpace();

	updateBoundaryConditions();

	initializeTimeAdvance ( );

	// Setting auxiliary variables for the fluid solver
	M_fluid->setAlpha(M_fluidTimeAdvance->alpha());
	M_fluid->setTimeStep(M_dt);
	M_fluid->buildSystem();

	// Structure
	M_structure.reset (new StructuralOperator<mesh_Type> ( ) );
	M_structure->setup ( M_dataStructure, M_displacementFESpace, M_displacementETFESpace, M_structureBC, M_comm);
	double timeAdvanceCoefficient = M_structureTimeAdvance->coefficientSecondDerivative ( 0 ) / ( M_dt * M_dt );
	M_structure->buildSystem (timeAdvanceCoefficient);

	// Ale
	M_ale.reset( new ALESolver ( *M_aleFESpace, M_comm ) );
	M_ale->setUp( M_datafile );

	// Exporters
	setupExporters( );

	// Data needed by the Newton algorithm
	M_relativeTolerance = M_datafile ( "newton/abstol", 1.e-4);
	M_absoluteTolerance = M_datafile ( "newton/reltol", 1.e-4);
	M_etaMax = M_datafile ( "newton/etamax", 1e-4);
	M_maxiterNonlinear = M_datafile ( "newton/maxiter", 10);
	M_nonLinearLineSearch = M_datafile ( "newton/NonLinearLineSearch", 0);
	if (M_comm->MyPID() == 0)
		M_out_res.open ("residualsNewton");
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
	M_structureDisplacement.reset ( new VectorEpetra ( M_displacementFESpace->map(), M_exporterStructure->mapType() ) );

	*M_fluidVelocity *= 0;
	*M_fluidPressure *= 0;
	*M_fluidDisplacement *= 0;
	*M_structureDisplacement *= 0;

	M_exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f - velocity", M_fluid->uFESpace(), M_fluidVelocity, UInt (0) );
	M_exporterFluid->addVariable ( ExporterData<mesh_Type>::ScalarField, "f - pressure", M_fluid->pFESpace(), M_fluidPressure, UInt (0) );
	M_exporterFluid->addVariable ( ExporterData<mesh_Type>::VectorField, "f - displacement", M_aleFESpace, M_fluidDisplacement, UInt (0) );
	M_exporterStructure->addVariable ( ExporterData<mesh_Type>::VectorField, "s - displacement", M_displacementFESpace, M_structureDisplacement, UInt (0) );

	M_exporterFluid->postProcess(M_t_zero);
	M_exporterStructure->postProcess(M_t_zero);
}

void
FSIHandler::instantiateExporter( boost::shared_ptr< Exporter<mesh_Type > >& exporter,
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

void FSIHandler::createStructureFESpaces ( )
{
	const std::string dOrder = M_datafile ( "solid/space_discretization/order", "P2");
	M_displacementFESpace.reset ( new FESpace_Type (M_structureLocalMesh, dOrder, 3, M_comm) );
	M_displacementETFESpace.reset ( new solidETFESpace_Type (*M_structurePartitioner, & (M_displacementFESpace->refFE() ), & (M_displacementFESpace->fe().geoMap() ), M_comm) );
	M_displacementFESpaceSerial.reset ( new FESpace_Type (M_structureMesh, dOrder, 3, M_comm) );

	M_displayer.leaderPrintMax ( " Number of DOFs for the structure = ", M_displacementFESpace->dof().numTotalDof()*3 ) ;
}

void FSIHandler::createAleFESpace()
{
	const std::string aleOrder = M_datafile ( "ale/space_discretization/order", "P2");
	M_aleFESpace.reset ( new FESpace_Type (M_fluidLocalMesh, aleOrder, 3, M_comm) );
	M_displayer.leaderPrintMax ( " Number of DOFs for the ale = ", M_aleFESpace->dof().numTotalDof()*3 ) ;
}

void FSIHandler::setBoundaryConditions ( const bcPtr_Type& fluidBC, const bcPtr_Type& structureBC, const bcPtr_Type& aleBC)
{
	M_fluidBC 	  = fluidBC;
	M_structureBC = structureBC;
	M_aleBC       = aleBC;
}

void FSIHandler::updateBoundaryConditions ( )
{
	M_fluidBC->bcUpdate ( *M_fluid->uFESpace()->mesh(), M_fluid->uFESpace()->feBd(), M_fluid->uFESpace()->dof() );
	M_structureBC->bcUpdate ( *M_displacementFESpace->mesh(), M_displacementFESpace->feBd(), M_displacementFESpace->dof() );
	M_aleBC->bcUpdate ( *M_aleFESpace->mesh(), M_aleFESpace->feBd(), M_aleFESpace->dof() );
}

void FSIHandler::initializeTimeAdvance ( )
{
	// Fluid
	M_fluidTimeAdvance.reset( new TimeAndExtrapolationHandler ( ) );
	M_dt       = M_datafile("fluid/time_discretization/timestep",0.0);
	M_t_zero   = M_datafile("fluid/time_discretization/initialtime",0.0);
	M_t_end    = M_datafile("fluid/time_discretization/endtime",0.0);
	M_orderBDF = M_datafile("fluid/time_discretization/BDF_order",2);

	// Order of BDF and extrapolation for the velocity
	M_fluidTimeAdvance->setBDForder(M_orderBDF);
	M_fluidTimeAdvance->setMaximumExtrapolationOrder(M_orderBDF);
	M_fluidTimeAdvance->setTimeStep(M_dt);

	// Initialize time advance
	vector_Type velocityInitial ( M_fluid->uFESpace()->map() );
	std::vector<vector_Type> initialStateVelocity;
	velocityInitial *= 0 ;
	for ( UInt i = 0; i < M_orderBDF; ++i )
		initialStateVelocity.push_back(velocityInitial);

	M_fluidTimeAdvance->initialize(initialStateVelocity);

	// Structure
	const std::string timeAdvanceMethod =  M_datafile ( "solid/time_discretization/method", "Newmark");
	M_structureTimeAdvance.reset ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );
	UInt OrderDev = 2;
	M_structureTimeAdvance->setup (M_dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);
	M_structureTimeAdvance->setTimeStep ( M_dt );
	vectorPtr_Type disp (new VectorEpetra (M_displacementFESpace->map(), Unique) );
	*disp *= 0;
	std::vector<vectorPtr_Type> uv0;
	for ( UInt previousPass = 0; previousPass < M_structureTimeAdvance->size() ; previousPass++)
	{
		Real previousTimeStep = M_t_zero - previousPass * M_dt;
		uv0.push_back (disp);
	}
	M_structureTimeAdvance->setInitialCondition (uv0);
	M_structureTimeAdvance->updateRHSContribution ( M_dt );

	// Ale - the same method as for the structure
	M_aleTimeAdvance.reset ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );
	M_aleTimeAdvance->setup (M_dataStructure->dataTimeAdvance()->orderBDF() , 1);
	M_aleTimeAdvance->setTimeStep ( M_dt );
	vectorPtr_Type fluidDisp (new VectorEpetra (M_aleFESpace->map(), Unique) );
	*fluidDisp *= 0;
	std::vector<vectorPtr_Type> df0;
	for ( UInt previousPass = 0; previousPass < M_aleTimeAdvance->size() ; previousPass++)
	{
		Real previousTimeStep = M_t_zero - previousPass * M_dt;
		df0.push_back (fluidDisp);
	}
	M_aleTimeAdvance->setInitialCondition (df0);
	M_aleTimeAdvance->updateRHSContribution ( M_dt );
}

void FSIHandler::buildInterfaceMaps ()
{
	markerID_Type interface = M_datafile("interface/flag", 1);
	Real tolerance = M_datafile("interface/tolerance", 1);
	Int flag = M_datafile("interface/fluid_vertex_flag", 0);

	M_dofStructureToFluid.reset ( new DOFInterface3Dto3D );
	M_dofStructureToFluid->setup ( M_fluid->uFESpace()->refFE(), M_fluid->uFESpace()->dof(), M_displacementFESpaceSerial->refFE(), M_displacementFESpaceSerial->dof() );
	M_dofStructureToFluid->update ( *M_fluid->uFESpace()->mesh(), interface, *M_displacementFESpaceSerial->mesh(),  interface, tolerance, &flag);

	createInterfaceMaps ( M_dofStructureToFluid->localDofMap ( ) );

	constructInterfaceMap ( M_dofStructureToFluid->localDofMap ( ), M_displacementFESpace->map().map(Unique)->NumGlobalElements()/nDimensions );
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
	int numMyElements = M_structureInterfaceMap->map (Unique)->NumMyElements();
	numInterfaceDof[pid] = numMyElements;
	MapEpetra subMap (*M_structureInterfaceMap->map (Unique), (UInt) 0, subdomainMaxId);

	M_numerationInterface.reset (new VectorEpetra (subMap, Unique) );

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
	}// so the map for the coupling part of the matrix is just Unique

	M_lagrangeMap.reset (new MapEpetra (-1, static_cast< Int> ( couplingVector.size() ), &couplingVector[0], M_comm) );
}

void
FSIHandler::assembleCoupling ( )
{
	M_coupling.reset ( new FSIcouplingCE ( M_comm ) );

	M_coupling->setUp ( M_dt, M_structureInterfaceMap->mapSize()/3.0 , M_structureTimeAdvance->coefficientFirstDerivative ( 0 ),
						M_lagrangeMap, M_fluid->uFESpace(), M_displacementFESpace, M_numerationInterface );

	M_coupling->buildBlocks ( M_dofStructureToFluid->localDofMap ( ) );

}

void
FSIHandler::buildMonolithicMap ( )
{
	M_monolithicMap.reset( new map_Type ( M_fluid->uFESpace()->map() ) );
	*M_monolithicMap += M_fluid->pFESpace()->map();
	*M_monolithicMap += M_displacementFESpace->map();
	*M_monolithicMap += *M_lagrangeMap;
	*M_monolithicMap += M_aleFESpace->map();
}

void
FSIHandler::solveFSIproblem ( )
{
	LifeChrono iterChrono;
	M_time = M_t_zero + M_dt;

	buildMonolithicMap ( );
	M_solution.reset ( new VectorEpetra ( *M_monolithicMap ) );
	*M_solution *= 0;

	// Apply boundary conditions for the ale problem (the matrix will not change during the simulation)
	M_ale->applyBoundaryConditions ( *M_aleBC );

	// Apply boundary conditions for the structure problem (the matrix will not change during the simulation, it is linear elasticity)
	getMatrixStructure ( );
	getRhsStructure ( );
	applyBCstructure ( );

	for ( ; M_time <= M_t_end + M_dt / 2.; M_time += M_dt)
	{
		M_displayer.leaderPrint ( "\n-----------------------------------\n" ) ;
		M_displayer.leaderPrintMax ( "FSI - solving now for time ", M_time ) ;
		M_displayer.leaderPrint ( "\n" ) ;
		iterChrono.start();

		updateSystem ( );

		// Using the solution at the previous timestep as initial guess -> TODO: extrapolation
		UInt status = NonLinearRichardson ( *M_solution, *this, M_absoluteTolerance, M_relativeTolerance, M_maxiterNonlinear, M_etaMax,
											M_nonLinearLineSearch, 0, 2, M_out_res, M_time);

		iterChrono.stop();
		M_displayer.leaderPrint ( "\n" ) ;
		M_displayer.leaderPrintMax ( "FSI - timestep solved in ", iterChrono.diff() ) ;
		iterChrono.reset();
		M_displayer.leaderPrint ( "-----------------------------------\n\n" ) ;

		// Updating all the time-advance objects


		// Export the solution obtained at the current timestep
		M_exporterFluid->postProcess(M_time);
		M_exporterStructure->postProcess(M_time);
	}

	M_exporterFluid->closeFile();
	M_exporterStructure->closeFile();
}

void
FSIHandler::evalResidual(vector_Type& residual, const vector_Type& solution, const UInt iter_newton)
{
	residual = 0.;

	//---------------------------------------------------------//
	// First: extract the fluid displacement and move the mesh //
	//---------------------------------------------------------//

	UInt offset ( M_fluid->uFESpace()->map().mapSize() +
				  M_fluid->pFESpace()->map().mapSize() +
				  M_displacementFESpace->map().mapSize() +
				  M_lagrangeMap->mapSize() );

	vectorPtr_Type meshDisplacement ( new vector_Type (M_aleFESpace->map() ) );
	meshDisplacement->subset (solution, offset);
	vectorPtr_Type mmRep ( new vector_Type (*meshDisplacement, Repeated ) );
	moveMesh ( *mmRep );

	//-------------------------------------------------------------------//
	// Second: re-assemble the fluid blocks since we have moved the mesh //
	//-------------------------------------------------------------------//

	M_fluid->buildSystem();
	M_fluid->updateSystem ( M_beta_star, M_rhs_velocity );
	M_fluid->applyBoundaryConditions ( M_fluidBC, M_time );

	//--------------------------------------//
	// Third: initialize the apply operator //
	//--------------------------------------//

	initializeApplyOperator ( );

	//-----------------------------//
	// Forth: compute the residual //
	//-----------------------------//

	VectorEpetra tmp(*M_monolithicMap);
	tmp.zero();
	M_applyOperator->Apply(solution.epetraVector(),residual.epetraVector());

	M_displayer.leaderPrint (" END OF EVAL RESIDUAL ");
	int kkk;
	std::cin >> kkk;
}

void
FSIHandler::solveJac( vector_Type& increment, const vector_Type& residual, const Real linearRelTol )
{

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
	// Fluid update - initialize vectors
	M_u_star.reset ( new VectorEpetra ( M_fluid->uFESpace()->map ( ) ) );
	M_w_star.reset ( new VectorEpetra ( M_aleFESpace->map ( ) ) );
	M_beta_star.reset ( new VectorEpetra ( M_fluid->uFESpace()->map ( ) ) );
	M_rhs_velocity.reset ( new VectorEpetra ( M_fluid->uFESpace()->map ( ) ) );

	*M_u_star *= 0;
	*M_w_star *= 0;
	*M_beta_star *= 0;
	*M_rhs_velocity *= 0;

	// Compute extrapolation of the velocity for the fluid and rhs contribution due to the time derivative
	M_fluidTimeAdvance->extrapolate (M_orderBDF, *M_u_star);
	M_fluidTimeAdvance->rhsContribution (*M_rhs_velocity);

	// Extrapolate the mesh velocity
	M_aleTimeAdvance->extrapolation ( *M_w_star );

	// compute beta* = u*-w*
	*M_beta_star += *M_u_star;
	*M_beta_star -= *M_w_star;
}

void
FSIHandler::initializeApplyOperator ( )
{
	Operators::FSIApplyOperator::operatorPtrContainer_Type operData(5,5);
	operData(0,0) = M_fluid->getF()->matrixPtr();
	operData(0,1) = M_fluid->getBtranspose()->matrixPtr();
	operData(0,3) = M_coupling->lambdaToFluidMomentum()->matrixPtr();
	operData(1,0) = M_fluid->getB()->matrixPtr();
	operData(2,2) = M_matrixStructure->matrixPtr();
	operData(2,3) = M_coupling->lambdaToStructureMomentum()->matrixPtr();
	operData(3,0) = M_coupling->fluidVelocityToLambda()->matrixPtr();
	operData(3,2) = M_coupling->structureDisplacementToLambda()->matrixPtr();
	operData(4,2) = M_coupling->structureDisplacementToFluidDisplacement()->matrixPtr();
	operData(4,4) = M_ale->matrix()->matrixPtr();
	M_applyOperator->setUp(operData, M_comm);
}

void
FSIHandler::getMatrixStructure ( )
{
	M_matrixStructure.reset (new matrix_Type ( M_displacementFESpace->map(), 1 ) );
	*M_matrixStructure *= 0.0;

	vectorPtr_Type solidPortion ( new vector_Type ( M_displacementFESpace->map() ) );
	*solidPortion *= 0;

	M_structure->material()->updateJacobianMatrix ( *solidPortion, M_dataStructure, M_structure->mapMarkersVolumes(), M_structure->mapMarkersIndexes(), M_structure->displayerPtr() );
	*M_matrixStructure += *M_structure->massMatrix();
	*M_matrixStructure += * (M_structure->material()->jacobian() );
}

void
FSIHandler::getRhsStructure ( )
{
	Real timeAdvanceCoefficient = M_structureTimeAdvance->coefficientSecondDerivative ( 0 ) / ( M_dt * M_dt );
	M_structureTimeAdvance->updateRHSContribution ( M_dt );

	M_rhsStructure.reset ( new VectorEpetra ( M_displacementFESpace->map() ) );
	*M_rhsStructure *= 0;
	*M_rhsStructure += *M_structure->massMatrix() * M_structureTimeAdvance->rhsContributionSecondDerivative() / timeAdvanceCoefficient;
}

void
FSIHandler::applyBCstructure ( )
{
	if ( !M_structureBC->bcUpdateDone() )
		M_structureBC->bcUpdate ( *M_displacementFESpace->mesh(), M_displacementFESpace->feBd(), M_displacementFESpace->dof() );



	bcManage ( *M_matrixStructure, *M_rhsStructure, *M_displacementFESpace->mesh(), M_displacementFESpace->dof(), *M_structureBC, M_displacementFESpace->feBd(), 1.0, M_time );

	M_matrixStructure->globalAssemble();
}

}
