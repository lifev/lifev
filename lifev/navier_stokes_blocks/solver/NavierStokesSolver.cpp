#include <lifev/navier_stokes_blocks/solver/NavierStokesSolver.hpp>


namespace LifeV
{

NavierStokesSolver::NavierStokesSolver(const dataFile_Type dataFile, const commPtr_Type& communicator):
		M_comm(communicator),
		M_dataFile(dataFile),
		M_displayer(communicator),
		M_graphIsBuilt(false),
        M_oper(new Operators::NavierStokesOperator),
        M_invOper(),
        M_fullyImplicit(false),
        M_graphPCDisBuilt(false),
        M_steady ( dataFile("fluid/miscellaneous/steady", false) ),
        M_density ( dataFile("fluid/physics/density", 1.0 ) ),
        M_viscosity ( dataFile("fluid/physics/viscosity", 0.035 ) )
{
	M_prec.reset ( Operators::NSPreconditionerFactory::instance().createObject (dataFile("fluid/preconditionerType","none")));
}

NavierStokesSolver::~NavierStokesSolver()
{
}

void NavierStokesSolver::setParameters( )
{
    std::string optionsPrec = M_dataFile("fluid/options_preconditioner","solverOptionsFast");
    optionsPrec += ".xml";
	Teuchos::RCP<Teuchos::ParameterList> solversOptions = Teuchos::getParametersFromXmlFile (optionsPrec);
	M_prec->setOptions(*solversOptions);
	setSolversOptions(*solversOptions);
}


void NavierStokesSolver::setup(const meshPtr_Type& mesh)
{
	std::string uOrder = M_dataFile("fluid/space_discretization/vel_order","P1");
	std::string pOrder = M_dataFile("fluid/space_discretization/pres_order","P1");

	Int geoDimensions = mesh_Type::S_geoDimensions;

	M_velocityFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, uOrder, geoDimensions, M_comm) );
	M_pressureFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, pOrder, 1, M_comm) );
	M_velocityFESpaceScalar.reset (new FESpace<mesh_Type, map_Type> (mesh, uOrder, 1, M_comm) );

	M_fespaceUETA.reset( new ETFESpace_velocity(M_velocityFESpace->mesh(), &(M_velocityFESpace->refFE()), M_comm));
	M_fespacePETA.reset( new ETFESpace_pressure(M_pressureFESpace->mesh(), &(M_pressureFESpace->refFE()), M_comm));

	M_uExtrapolated.reset( new vector_Type ( M_velocityFESpace->map(), Repeated ) );
	M_uExtrapolated->zero();

	M_stiffStrain = M_dataFile("fluid/space_discretization/stiff_strain", true);

    M_velocity.reset( new vector_Type(M_velocityFESpace->map()) );
    M_pressure.reset( new vector_Type(M_pressureFESpace->map()) );

    M_block00.reset( new matrix_Type(M_velocityFESpace->map() ) );
    M_block01.reset( new matrix_Type(M_velocityFESpace->map() ) );
    M_block10.reset( new matrix_Type(M_pressureFESpace->map() ) );
    M_block11.reset( new matrix_Type(M_pressureFESpace->map() ) );

    M_block00->zero();
    M_block10->zero();
    M_block01->zero();
    M_block11->zero();

    M_fullyImplicit = M_dataFile ( "newton/convectiveImplicit", false);

    if ( M_fullyImplicit )
    {
    	M_displayer.leaderPrint ( " F - solving nonlinear Navier-Stokes\n");
    	M_relativeTolerance = M_dataFile ( "newton/abstol", 1.e-4);
    	M_absoluteTolerance = M_dataFile ( "newton/reltol", 1.e-4);
    	M_etaMax = M_dataFile ( "newton/etamax", 1e-4);
    	M_maxiterNonlinear = M_dataFile ( "newton/maxiter", 10);

    	M_nonLinearLineSearch = M_dataFile ( "newton/NonLinearLineSearch", 0);

    	if (M_comm->MyPID() == 0)
    		M_out_res.open ("residualsNewton");

    	M_monolithicMap.reset( new map_Type ( M_velocityFESpace->map() ) );
    	*M_monolithicMap += M_pressureFESpace->map();

    	M_solution.reset( new vector_Type ( *M_monolithicMap ) );
    }

    M_displayer.leaderPrintMax ( " Number of DOFs for the velocity = ", M_velocityFESpace->dof().numTotalDof()*3 ) ;
    M_displayer.leaderPrintMax ( " Number of DOFs for the pressure = ", M_pressureFESpace->dof().numTotalDof() ) ;
}

void NavierStokesSolver::setSolversOptions(const Teuchos::ParameterList& solversOptions)
{
    boost::shared_ptr<Teuchos::ParameterList> monolithicOptions;
    monolithicOptions.reset(new Teuchos::ParameterList(solversOptions.sublist("MonolithicOperator")) );
    M_pListLinSolver = monolithicOptions;
}

void NavierStokesSolver::buildGraphs()
{
	M_displayer.leaderPrint ( " F - Pre-building the graphs... ");
	LifeChrono chrono;
	chrono.start();

	{
		using namespace ExpressionAssembly;

		if ( !M_steady )
		{
			// Graph velocity mass -> block (0,0)
			M_Mu_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 M_density * dot ( phi_i, phi_j )
			) >> M_Mu_graph;
			M_Mu_graph->GlobalAssemble();
			M_Mu_graph->OptimizeStorage();
		}

		// Graph block (0,1) of NS
		M_Btranspose_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespacePETA,
					 value(-1.0) * phi_j * div(phi_i)
		) >> M_Btranspose_graph;
		M_Btranspose_graph->GlobalAssemble( *(M_pressureFESpace->map().map (Unique)), *(M_velocityFESpace->map().map (Unique)) );
        M_Btranspose_graph->OptimizeStorage();

		// Graph block (1,0) of NS
		M_B_graph.reset (new Epetra_FECrsGraph (Copy, *(M_pressureFESpace->map().map (Unique)), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespaceUETA,
					 phi_i * div(phi_j)
		) >> M_B_graph;
		M_B_graph->GlobalAssemble( *(M_velocityFESpace->map().map (Unique)), *(M_pressureFESpace->map().map (Unique)) );
        M_B_graph->OptimizeStorage();

		// Graph convective term, block (0,0)
		M_C_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespaceUETA,
					 dot( M_density*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i)
		) >> M_C_graph;
		M_C_graph->GlobalAssemble();
        M_C_graph->OptimizeStorage();

		// Graph stiffness, block (0,0)
		M_A_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		if (M_stiffStrain)
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 value( 0.5 * M_viscosity ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A_graph;
		}
		else
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 M_viscosity * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A_graph;
		}
		M_A_graph->GlobalAssemble();
        M_A_graph->OptimizeStorage();

		// Graph of entire block (0,0)
		M_F_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		if (M_stiffStrain)
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 dot ( phi_i, phi_j ) + // mass
						 dot( M_density*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 value( 0.5 * M_viscosity ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
			) >> M_F_graph;


			M_Jacobian_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 dot ( phi_i, phi_j ) + // mass
						 dot( M_density*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 dot( M_density * phi_j * grad(M_fespaceUETA, *M_uExtrapolated), phi_i ) + // part of the Jacobian
						 value( 0.5 * M_viscosity ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
						) >> M_Jacobian_graph;
			M_Jacobian_graph->GlobalAssemble();
			M_Jacobian_graph->OptimizeStorage();

		}
		else
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 dot ( phi_i, phi_j ) + // mass
						 dot( M_density*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 M_viscosity * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
			) >> M_F_graph;


			M_Jacobian_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 dot ( phi_i, phi_j ) + // mass
						 dot( M_density*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 dot( M_density * phi_j * grad(M_fespaceUETA, *M_uExtrapolated), phi_i ) + // part of the Jacobian
						 M_viscosity * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
					  ) >> M_Jacobian_graph;
			M_Jacobian_graph->GlobalAssemble();
			M_Jacobian_graph->OptimizeStorage();

		}
		M_F_graph->GlobalAssemble();
        M_F_graph->OptimizeStorage();
	}

	M_graphIsBuilt = true;

	chrono.stop();
	M_displayer.leaderPrintMax ( "   done in ", chrono.diff() ) ;
}

void NavierStokesSolver::updatePCD(const vectorPtr_Type& velocity)
{
	// Note: for semi-implicit treatment of the convective term velocity is an extrapolation. When the convective term is treated in a fully-implicit
	// manner, velocoty is given from the previous newton step. Pay attantion when calling updatePCD, give in input the correct vector. In this way the
	// method is flexible and can handle both scenarios.

	vectorPtr_Type velocityRepeated( new vector_Type ( *velocity, Repeated ) );

	if ( !M_graphPCDisBuilt )
		buildPCDGraphs();

	M_displayer.leaderPrint ( " F - Assembling PCD terms... ");
	LifeChrono chrono;
	chrono.start();

	{
		using namespace ExpressionAssembly;

		M_Mp.reset (new matrix_Type ( M_pressureFESpace->map(), *M_Mp_graph ) );
		M_Mp->zero();
		integrate ( elements (M_fespaceUETA->mesh() ),
					M_pressureFESpace->qr(),
					M_fespacePETA,
					M_fespacePETA,
					phi_i * phi_j
				  ) >> M_Mp;
		M_Mp->globalAssemble();

		// Graph pressure mass
		M_Fp.reset (new matrix_Type ( M_pressureFESpace->map(), *M_Fp_graph ) );
		integrate ( elements (M_fespacePETA->mesh() ),
					M_pressureFESpace->qr(),
					M_fespacePETA,
					M_fespacePETA,
					value( M_density*M_alpha/M_timeStep ) * phi_i * phi_j
					+ value( 0.5 * M_viscosity ) * dot( grad(phi_i) + grad(phi_i) , grad(phi_j) + grad(phi_j) )
					+ value( M_density ) * dot( value(M_fespaceUETA, *velocityRepeated),grad(phi_j)) * phi_i
				  ) >> M_Fp;
	}

	bcManageMatrix( *M_Fp, *M_pressureFESpace->mesh(), M_pressureFESpace->dof(), *M_bcPCD, M_pressureFESpace->feBd(), 1.0, 0.0);
	M_Fp->globalAssemble();

	chrono.stop();
	M_displayer.leaderPrintMax ( "   done in ", chrono.diff() ) ;
}

void NavierStokesSolver::buildPCDGraphs()
{
	M_displayer.leaderPrint ( " F - Pre-building the graphs for PCD... ");
	LifeChrono chrono;
	chrono.start();

	{
		using namespace ExpressionAssembly;

		// Graph pressure mass
		M_Mp_graph.reset (new Epetra_FECrsGraph (Copy, * (M_pressureFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespacePETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespacePETA,
					 phi_i * phi_j
				   ) >> M_Mp_graph;
		M_Mp_graph->GlobalAssemble();
		M_Mp_graph->OptimizeStorage();

		// Graph pressure mass
		M_Fp_graph.reset (new Epetra_FECrsGraph (Copy, * (M_pressureFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespacePETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespacePETA,
					 value( M_density*M_alpha/M_timeStep ) * phi_i * phi_j
					 + value( 0.5 * M_viscosity ) * dot( grad(phi_i) + grad(phi_i) , grad(phi_j) + grad(phi_j) )
					 + value( M_density ) * dot( value(M_fespaceUETA, *M_uExtrapolated),grad(phi_j)) * phi_i
				   ) >> M_Fp_graph;
		M_Fp_graph->GlobalAssemble();
		M_Fp_graph->OptimizeStorage();

	}

	M_graphPCDisBuilt = true;

	chrono.stop();
	M_displayer.leaderPrintMax ( "   done in ", chrono.diff() ) ;
}

void NavierStokesSolver::buildSystem()
{
	if ( !M_graphIsBuilt )
		buildGraphs();

	M_displayer.leaderPrint ( " F - Assembling constant terms... ");
	LifeChrono chrono;
	chrono.start();

	{
		using namespace ExpressionAssembly;

		if ( !M_steady )
		{
			// Graph velocity mass -> block (0,0)
			M_Mu.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Mu_graph ) );
			M_Mu->zero();
			integrate ( elements (M_fespaceUETA->mesh() ),
						M_velocityFESpace->qr(),
						M_fespaceUETA,
						M_fespaceUETA,
						M_density * dot ( phi_i, phi_j )
					  ) >> M_Mu;
			M_Mu->globalAssemble();
		}

		M_Btranspose.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Btranspose_graph ) );
		M_Btranspose->zero();
		integrate( elements(M_fespaceUETA->mesh()),
				   M_velocityFESpace->qr(),
				   M_fespaceUETA,
				   M_fespacePETA,
				   value(-1.0) * phi_j * div(phi_i)
		) >> M_Btranspose;
		M_Btranspose->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );

		M_B.reset (new matrix_Type ( M_pressureFESpace->map(), *M_B_graph ) );
		M_B->zero();
		integrate( elements(M_fespaceUETA->mesh()),
				   M_pressureFESpace->qr(),
				   M_fespacePETA,
				   M_fespaceUETA,
				   phi_i * div(phi_j)
		) >> M_B;
		M_B->globalAssemble( M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());

		M_A.reset (new matrix_Type ( M_velocityFESpace->map(), *M_A_graph ) );
		M_A->zero();
		if ( M_stiffStrain )
		{
			integrate( elements(M_fespaceUETA->mesh()),
					   M_velocityFESpace->qr(),
					   M_fespaceUETA,
					   M_fespaceUETA,
					   value( 0.5 * M_viscosity ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A;
		}
		else
		{
			integrate( elements(M_fespaceUETA->mesh()),
					   M_velocityFESpace->qr(),
					   M_fespaceUETA,
				 	   M_fespaceUETA,
					   M_viscosity * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A;
		}
		M_A->globalAssemble();

		M_C.reset (new matrix_Type ( M_velocityFESpace->map(), *M_C_graph ) );
		M_C->zero();
		M_C->globalAssemble();

		M_F.reset (new matrix_Type ( M_velocityFESpace->map(), *M_F_graph ) );
		M_F->zero();
		M_F->globalAssemble();

		M_Jacobian.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Jacobian_graph ) );
		M_Jacobian->zero();
		M_Jacobian->globalAssemble();

		M_block01->zero();
		M_block10->zero();
		*M_block01 += *M_Btranspose;
		*M_block10 += *M_B;
		M_block01->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );
		M_block10->globalAssemble(M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());
	}

	chrono.stop();
	M_displayer.leaderPrintMax ( " done in ", chrono.diff() ) ;
}

void NavierStokesSolver::updateSystem( const vectorPtr_Type& u_star, const vectorPtr_Type& rhs_velocity )
{
	// Note that u_star HAS to extrapolated from outside. Hence it works also for FSI in this manner.
	M_uExtrapolated.reset( new vector_Type ( *u_star, Repeated ) );

	// Update convective term
	M_C->zero();
	{
		using namespace ExpressionAssembly;
		integrate( elements(M_fespaceUETA->mesh()),
				M_velocityFESpace->qr(),
				M_fespaceUETA,
				M_fespaceUETA,
				dot( M_density * value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) // semi-implicit treatment of the convective term
		)
		>> M_C;
	}
	M_C->globalAssemble();

	// Get the matrix corresponding to the block (0,0)
	M_block00->zero();
	*M_block00 += *M_Mu;
	*M_block00 *= M_alpha/M_timeStep;
	*M_block00 += *M_A;
	*M_block00 += *M_C;
	M_block00->globalAssemble();

	// Get the right hand side with inertia contribution
	M_rhs.reset( new vector_Type ( M_velocityFESpace->map(), Unique ) );
	M_rhs->zero();

	if ( !M_steady )
		*M_rhs = *M_Mu* (*rhs_velocity);
}

void NavierStokesSolver::applyGravityForce ( const Real& gravity, const Real& gravityDirection)
{
	vectorPtr_Type gravity_vector ( new vector_Type ( M_velocityFESpace->map(), Unique ) );
	vectorPtr_Type gravity_component ( new vector_Type ( M_velocityFESpaceScalar->map(), Unique ) );

	gravity_component->zero();
	*gravity_component += gravity;
	gravity_vector->subset(*gravity_component, M_velocityFESpaceScalar->map(), 0, gravityDirection*M_velocityFESpaceScalar->dof().numTotalDof() );

	*M_rhs += *M_Mu* (*gravity_vector);
}

void NavierStokesSolver::applyBoundaryConditions ( bcPtr_Type & bc, const Real& time )
{
	updateBCHandler(bc);
	bcManage ( *M_block00, *M_rhs, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 1.0, time );

	for (BCHandler::bcBaseIterator_Type it = bc->begin(); it != bc->end(); ++it)
	{
		if ( it->type() != Essential)
			continue;

		bcManageMatrix( *M_block01, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 0.0, 0.0);
	}
}

void NavierStokesSolver::iterate( bcPtr_Type & bc, const Real& time )
{
	applyBoundaryConditions ( bc, time );

    //(1) Set up the OseenOperator
    M_displayer.leaderPrint( "\tNS operator - set up the block operator...");
    LifeChrono chrono;
    chrono.start();

    Operators::NavierStokesOperator::operatorPtrContainer_Type operData(2,2);
    operData(0,0) = M_block00->matrixPtr();
    operData(0,1) = M_block01->matrixPtr();
    operData(1,0) = M_block10->matrixPtr();
    M_oper->setUp(operData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(2) Set the data for the preconditioner

	M_displayer.leaderPrint( "\tPreconditioner operator - set up the block operator...");
    chrono.reset();
    chrono.start();

    if ( std::strcmp(M_prec->Label(),"aSIMPLEOperator")==0 )
    {
    	M_prec->setUp(M_block00, M_block10, M_block01);
    	M_prec->setDomainMap(M_oper->OperatorDomainBlockMapPtr());
    	M_prec->setRangeMap(M_oper->OperatorRangeBlockMapPtr());
    	M_prec->updateApproximatedMomentumOperator();
    	M_prec->updateApproximatedSchurComplementOperator();
    }
    else if ( std::strcmp(M_prec->Label(),"aPCDOperator")==0 )
    {
    	updatePCD(M_uExtrapolated);
    	M_prec->setUp(M_block00, M_block10, M_block01, M_Fp, M_Mp, M_Mu);
    	M_prec->setDomainMap(M_oper->OperatorDomainBlockMapPtr());
    	M_prec->setRangeMap(M_oper->OperatorRangeBlockMapPtr());
    	M_prec->updateApproximatedMomentumOperator();
    	M_prec->updateApproximatedSchurComplementOperator();
    	M_prec->updateApproximatedPressureMassOperator();
    }

    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(3) Set the solver for the linear system
    M_displayer.leaderPrint( "\tset up the Trilinos solver...");
    chrono.start();
    std::string solverType(M_pListLinSolver->get<std::string>("Linear Solver Type"));
    M_invOper.reset(Operators::InvertibleOperatorFactory::instance().createObject(solverType));

    M_invOper->setOperator(M_oper);
    M_invOper->setPreconditioner(M_prec);
    M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));

    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    // Solving the system
    BlockEpetra_Map upMap;
    upMap.setUp ( M_velocityFESpace->map().map(Unique), M_pressureFESpace->map().map(Unique));

    BlockEpetra_MultiVector up(upMap, 1), rhs(upMap, 1);
    rhs.block(0).Update(1.0, M_rhs->epetraVector(), 0.);

    // Solving the linear system
    M_invOper->ApplyInverse(rhs,up);

    M_velocity->epetraVector().Update(1.0,up.block(0),0.0);
    M_pressure->epetraVector().Update(1.0,up.block(1),0.0);

}

void NavierStokesSolver::iterate_nonlinear( bcPtr_Type & bc, const Real& time )
{
	// Initialize the solution
	M_solution->zero();

	applyBoundaryConditionsSolution ( bc, time ); // the second argument is zero since the problem is steady

	// Call Newton
	UInt status = NonLinearRichardson ( *M_solution, *this, M_absoluteTolerance, M_relativeTolerance, M_maxiterNonlinear, M_etaMax,
			M_nonLinearLineSearch, 0, 2, M_out_res, 0.0);

	M_velocity->subset ( *M_solution, M_velocityFESpace->map(), 0, 0 );
	M_pressure->subset ( *M_solution, M_pressureFESpace->map(), M_velocityFESpace->map().mapSize(), 0 );
}

void NavierStokesSolver::iterate_steady( bcPtr_Type & bc )
{
	// Initialize the solution
	M_solution->zero();

	applyBoundaryConditionsSolution ( bc, 0.0 ); // the second argument is zero since the problem is steady

	// Call Newton
	UInt status = NonLinearRichardson ( *M_solution, *this, M_absoluteTolerance, M_relativeTolerance, M_maxiterNonlinear, M_etaMax,
			M_nonLinearLineSearch, 0, 2, M_out_res, 0.0);

	M_velocity->subset ( *M_solution, M_velocityFESpace->map(), 0, 0 );
	M_pressure->subset ( *M_solution, M_pressureFESpace->map(), M_velocityFESpace->map().mapSize(), 0 );
}

void NavierStokesSolver::applyBoundaryConditionsSolution ( bcPtr_Type & bc, const Real& time )
{
	updateBCHandler(bc);
	bcManageRhs ( *M_solution, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 1.0, time );
}

void NavierStokesSolver::updateBCHandler( bcPtr_Type & bc )
{
	bc->bcUpdate ( *M_velocityFESpace->mesh(), M_velocityFESpace->feBd(), M_velocityFESpace->dof() );
}

void NavierStokesSolver::evaluateResidual( const vectorPtr_Type& convective_velocity,
					   	   	   	   	   	   const vectorPtr_Type& velocity_km1,
					   	   	   	   	   	   const vectorPtr_Type& pressure_km1,
					   	   	   	   	   	   const vectorPtr_Type& rhs_velocity,
					   	   	   	   	   	   vectorPtr_Type& residual)
{
	residual->zero();

	// Residual vector for the velocity and pressure components
	vectorPtr_Type res_velocity ( new vector_Type ( M_velocityFESpace->map(), Unique ) );
	vectorPtr_Type res_pressure ( new vector_Type ( M_pressureFESpace->map(), Unique ) );
	res_velocity->zero();
	res_pressure->zero();

	// Get repeated versions of input vectors for the assembly
	vectorPtr_Type convective_velocity_repeated ( new vector_Type (*convective_velocity, Repeated) );
	vectorPtr_Type u_km1_repeated ( new vector_Type (*velocity_km1, Repeated) );
	vectorPtr_Type p_km1_repeated ( new vector_Type (*pressure_km1, Repeated) );
	vectorPtr_Type rhs_velocity_repeated ( new vector_Type (*rhs_velocity, Repeated) );

	{
		using namespace ExpressionAssembly;

		if ( M_stiffStrain )
		{
			integrate ( elements ( M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_density * value ( M_alpha / M_timeStep ) * dot ( value ( M_fespaceUETA, *u_km1_repeated ), phi_i ) -
					M_density * dot ( value ( M_fespaceUETA, *rhs_velocity_repeated ), phi_i ) +
					value ( 0.5 ) * M_viscosity * dot ( grad ( phi_i )  + transpose ( grad ( phi_i ) ), grad ( M_fespaceUETA, *u_km1_repeated ) + transpose ( grad ( M_fespaceUETA, *u_km1_repeated ) ) ) +
					M_density * dot ( value ( M_fespaceUETA, *convective_velocity_repeated ) * grad ( M_fespaceUETA, *u_km1_repeated ), phi_i ) +
					value ( -1.0 ) * value ( M_fespacePETA, *p_km1_repeated ) * div ( phi_i )
			) >> res_velocity;
		}
		else
		{
			integrate ( elements ( M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_density * value ( M_alpha / M_timeStep ) * dot ( value ( M_fespaceUETA, *u_km1_repeated ), phi_i ) -
					M_density * dot ( value ( M_fespaceUETA, *rhs_velocity_repeated ), phi_i ) +
					M_viscosity * dot ( grad ( phi_i ), grad ( M_fespaceUETA, *u_km1_repeated ) + transpose ( grad ( M_fespaceUETA, *u_km1_repeated ) ) ) +
					M_density * dot ( value ( M_fespaceUETA, *convective_velocity_repeated ) * grad ( M_fespaceUETA, *u_km1_repeated ), phi_i ) +
					value ( -1.0 ) * value ( M_fespacePETA, *p_km1_repeated ) * div ( phi_i )
			) >> res_velocity;
		}

		integrate ( elements ( M_fespaceUETA->mesh() ),
				M_pressureFESpace->qr(),
				M_fespacePETA,
				trace ( grad ( M_fespaceUETA, *u_km1_repeated ) ) * phi_i
		) >> res_pressure;
	}

	res_velocity->globalAssemble();
	res_pressure->globalAssemble();

	residual->subset ( *res_velocity, M_velocityFESpace->map(), 0, 0 );
	residual->subset ( *res_pressure, M_pressureFESpace->map(), 0, M_velocityFESpace->map().mapSize() );

}

void NavierStokesSolver::evalResidual(vector_Type& residual, const vector_Type& solution, const UInt iter_newton)
{
	// Residual to zero
	residual.zero();

	// Extract the velocity and the pressure at the previous Newton step ( index k minus 1, i.e. km1)
	vectorPtr_Type u_km1 ( new vector_Type ( M_velocityFESpace->map(), Repeated ) );
	vectorPtr_Type p_km1 ( new vector_Type ( M_pressureFESpace->map(), Repeated ) );
	u_km1->zero();
	p_km1->zero();
	u_km1->subset ( solution, M_velocityFESpace->map(), 0, 0 );
	p_km1->subset ( solution, M_pressureFESpace->map(), M_velocityFESpace->map().mapSize(), 0 );

	// Residual vector for the velocity and pressure components
	vectorPtr_Type res_velocity ( new vector_Type ( M_velocityFESpace->map(), Unique ) );
	vectorPtr_Type res_pressure ( new vector_Type ( M_pressureFESpace->map(), Unique ) );
	res_velocity->zero();
	res_pressure->zero();

	vectorPtr_Type rhs_velocity_repeated;
	if ( !M_steady )
	{
		rhs_velocity_repeated.reset( new vector_Type ( *M_velocityRhs, Repeated ) );
	}

	if ( M_steady )
	{
		// Assemble the residual vector:
		using namespace ExpressionAssembly;

		if ( M_stiffStrain )
		{
			integrate ( elements ( M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					value ( 0.5 ) * M_viscosity * dot ( grad ( phi_i )  + transpose ( grad ( phi_i ) ), grad ( M_fespaceUETA, *u_km1 ) + transpose ( grad ( M_fespaceUETA, *u_km1 ) ) ) +
					M_density * dot ( value ( M_fespaceUETA, *u_km1 ) * grad ( M_fespaceUETA, *u_km1 ), phi_i ) +
					value ( -1.0 ) * value ( M_fespacePETA, *p_km1 ) * div ( phi_i )

			) >> res_velocity;
		}
		else
		{
			integrate ( elements ( M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_viscosity * dot ( grad ( phi_i ), grad ( M_fespaceUETA, *u_km1 ) + transpose ( grad ( M_fespaceUETA, *u_km1 ) ) ) +
					M_density * dot ( value ( M_fespaceUETA, *u_km1 ) * grad ( M_fespaceUETA, *u_km1 ), phi_i ) +
					value ( -1.0 ) * value ( M_fespacePETA, *p_km1 ) * div ( phi_i )
			) >> res_velocity;
		}

		integrate ( elements ( M_fespaceUETA->mesh() ),
				M_pressureFESpace->qr(),
				M_fespacePETA,
				trace ( grad ( M_fespaceUETA, *u_km1 ) ) * phi_i
		) >> res_pressure;
	}
	else
	{
		using namespace ExpressionAssembly;

		if ( M_stiffStrain )
		{
			integrate ( elements ( M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_density * value ( M_alpha / M_timeStep ) * dot ( value ( M_fespaceUETA, *u_km1 ), phi_i ) -
					M_density * dot ( value ( M_fespaceUETA, *rhs_velocity_repeated ), phi_i ) +
					value ( 0.5 ) * M_viscosity * dot ( grad ( phi_i )  + transpose ( grad ( phi_i ) ), grad ( M_fespaceUETA, *u_km1 ) + transpose ( grad ( M_fespaceUETA, *u_km1 ) ) ) +
					M_density * dot ( value ( M_fespaceUETA, *u_km1 ) * grad ( M_fespaceUETA, *u_km1 ), phi_i ) +
					value ( -1.0 ) * value ( M_fespacePETA, *p_km1 ) * div ( phi_i )
			) >> res_velocity;
		}
		else
		{
			integrate ( elements ( M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_density * value ( M_alpha / M_timeStep ) * dot ( value ( M_fespaceUETA, *u_km1 ), phi_i ) -
					M_density * dot ( value ( M_fespaceUETA, *rhs_velocity_repeated ), phi_i ) +
					M_viscosity * dot ( grad ( phi_i ), grad ( M_fespaceUETA, *u_km1 ) + transpose ( grad ( M_fespaceUETA, *u_km1 ) ) ) +
					M_density * dot ( value ( M_fespaceUETA, *u_km1 ) * grad ( M_fespaceUETA, *u_km1 ), phi_i ) +
					value ( -1.0 ) * value ( M_fespacePETA, *p_km1 ) * div ( phi_i )
			) >> res_velocity;
		}

		integrate ( elements ( M_fespaceUETA->mesh() ),
				M_pressureFESpace->qr(),
				M_fespacePETA,
				trace ( grad ( M_fespaceUETA, *u_km1 ) ) * phi_i
		) >> res_pressure;
	}

    res_velocity->globalAssemble();
    res_pressure->globalAssemble();

    bcManageRhs ( *res_velocity, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *M_bc, M_velocityFESpace->feBd(), 0.0, 0.0 );

    residual.subset ( *res_velocity, M_velocityFESpace->map(), 0, 0 );
    residual.subset ( *res_pressure, M_pressureFESpace->map(), 0, M_velocityFESpace->map().mapSize() );

    // We need to update the linearized terms in the jacobian coming from the convective part

    updateConvectiveTerm(u_km1);
    updateJacobian(*u_km1);

}

void NavierStokesSolver::updateConvectiveTerm ( const vectorPtr_Type& velocity)
{
    // Note that u_star HAS to extrapolated from outside. Hence it works also for FSI in this manner.
    vectorPtr_Type velocity_repeated;
    velocity_repeated.reset ( new vector_Type ( *velocity, Repeated ) );

    // Update convective term
    M_C->zero();
    {
        using namespace ExpressionAssembly;
        integrate ( elements (M_fespaceUETA->mesh() ),
                    M_velocityFESpace->qr(),
                    M_fespaceUETA,
                    M_fespaceUETA,
                    dot ( M_density * value (M_fespaceUETA, *velocity_repeated) *grad (phi_j), phi_i)
                  )
                >> M_C;
    }

    // Assuming steady NS
    M_block00->zero();
    if ( !M_steady )
    {
    	*M_block00 += *M_Mu;
    	*M_block00 *= M_alpha / M_timeStep;
    }
    *M_block00 += *M_A;
    *M_block00 += *M_C;
    if ( !M_fullyImplicit )
    	M_block00->globalAssemble( );
}

void NavierStokesSolver::updateJacobian( const vector_Type& u_k )
{
	vector_Type uk_rep ( u_k, Repeated );
	M_Jacobian->zero();

	M_displayer.leaderPrint ( "[F] - Update Jacobian convective term\n" ) ;

	if ( M_fullyImplicit )
	{
		using namespace ExpressionAssembly;
		integrate( elements(M_fespaceUETA->mesh()),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_fespaceUETA,
					dot( M_density * phi_j * grad(M_fespaceUETA, uk_rep), phi_i )
		)
		>> M_Jacobian;
	}

	*M_block00 += *M_Jacobian;
	M_block00->globalAssemble( );
}

void NavierStokesSolver::solveJac( vector_Type& increment, const vector_Type& residual, const Real linearRelTol )
{
	// Apply BCs on the jacobian matrix
	applyBoundaryConditionsJacobian ( M_bc );

	//(1) Set up the OseenOperator
	M_displayer.leaderPrint( "\tNS operator - set up the block operator...");
	LifeChrono chrono;
	chrono.start();

	Operators::NavierStokesOperator::operatorPtrContainer_Type operDataJacobian ( 2, 2 );
	operDataJacobian ( 0, 0 ) = M_block00->matrixPtr();
	operDataJacobian ( 0, 1 ) = M_block01->matrixPtr();
	operDataJacobian ( 1, 0 ) = M_block10->matrixPtr();

	M_oper->setUp(operDataJacobian, M_displayer.comm());
	chrono.stop();
	M_displayer.leaderPrintMax(" done in " , chrono.diff() );

	M_displayer.leaderPrint( "\tPreconditioner operator - set up the block operator...");
	chrono.reset();
	chrono.start();

	//(2) Set the data for the preconditioner
	if ( std::strcmp(M_prec->Label(),"aSIMPLEOperator")==0 )
	{
		M_prec->setUp(M_block00, M_block10, M_block01);
		M_prec->setDomainMap(M_oper->OperatorDomainBlockMapPtr());
		M_prec->setRangeMap(M_oper->OperatorRangeBlockMapPtr());
		M_prec->updateApproximatedMomentumOperator();
		M_prec->updateApproximatedSchurComplementOperator();
	}
	else if ( std::strcmp(M_prec->Label(),"aPCDOperator")==0 )
	{
		updatePCD(M_uExtrapolated);
		M_prec->setUp(M_block00, M_block10, M_block01, M_Fp, M_Mp, M_Mu);
		M_prec->setDomainMap(M_oper->OperatorDomainBlockMapPtr());
		M_prec->setRangeMap(M_oper->OperatorRangeBlockMapPtr());
		M_prec->updateApproximatedMomentumOperator();
		M_prec->updateApproximatedSchurComplementOperator();
		M_prec->updateApproximatedPressureMassOperator();
	}
	chrono.stop();
	M_displayer.leaderPrintMax(" done in " , chrono.diff() );

	//(3) Set the solver for the linear system
	//(3) Set the solver for the linear system
	M_displayer.leaderPrint( "\tset up the Trilinos solver...");
	chrono.start();
	std::string solverType(M_pListLinSolver->get<std::string>("Linear Solver Type"));
	M_invOper.reset(Operators::InvertibleOperatorFactory::instance().createObject(solverType));
	chrono.stop();
	M_displayer.leaderPrintMax(" done in " , chrono.diff() );

	M_invOper->setOperator(M_oper);
	M_invOper->setPreconditioner(M_prec);
	M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));

	increment.zero();

	// Solving the linear system
	M_invOper->ApplyInverse(residual.epetraVector(), increment.epetraVector());
}

void NavierStokesSolver::applyBoundaryConditionsJacobian ( bcPtr_Type & bc )
{
	bcManageMatrix( *M_block00, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 1.0, 0.0);
	bcManageMatrix( *M_block01, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 0.0, 0.0);
}

}
