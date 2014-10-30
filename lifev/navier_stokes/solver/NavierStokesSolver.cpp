#include <lifev/navier_stokes/solver/NavierStokesSolver.hpp>

namespace LifeV
{

NavierStokesSolver::NavierStokesSolver(const dataFile_Type dataFile, const commPtr_Type& communicator):
		M_comm(communicator),
		M_dataFile(dataFile),
		M_displayer(communicator),
		M_graphIsBuilt(false),
        M_oper(new Operators::NavierStokesOperator),
        M_prec(new Operators::aSIMPLEOperator),
        M_invOper()
{
}

NavierStokesSolver::~NavierStokesSolver()
{
}


void NavierStokesSolver::setup(const meshPtr_Type& mesh)
{
	M_fluidData.reset( new OseenData() );
	M_fluidData->setup( M_dataFile );

	std::string uOrder = M_dataFile("fluid/space_discretization/vel_order","P1");
	std::string pOrder = M_dataFile("fluid/space_discretization/pres_order","P1");

	Int geoDimensions = mesh_Type::S_geoDimensions;

	M_velocityFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, uOrder, geoDimensions, M_comm) );
	M_pressureFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, pOrder, 1, M_comm) );

	M_fespaceUETA.reset( new ETFESpace_velocity(M_velocityFESpace->mesh(), &(M_velocityFESpace->refFE()), M_comm));
	M_fespacePETA.reset( new ETFESpace_pressure(M_pressureFESpace->mesh(), &(M_pressureFESpace->refFE()), M_comm));

	M_uExtrapolated.reset( new vector_Type ( M_velocityFESpace->map(), Repeated ) );
	*M_uExtrapolated *= 0;

	M_stiffStrain = M_dataFile("fluid/space_discretization/stiff_strain", true);

    Teuchos::RCP<Teuchos::ParameterList> solversOptions = Teuchos::getParametersFromXmlFile ("solversOptionsFast.xml");
    M_prec->setOptions(*solversOptions);
    setSolversOptions(*solversOptions);

    M_velocity.reset( new vector_Type(M_velocityFESpace->map()) );
    M_pressure.reset( new vector_Type(M_pressureFESpace->map()) );

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

		// Graph velocity mass -> block (0,0)
		M_Mu_graph.reset (new Epetra_FECrsGraph (Copy, * (M_velocityFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_fespaceUETA->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespaceUETA,
					 M_fluidData->density() * dot ( phi_i, phi_j )
		) >> M_Mu_graph;
		M_Mu_graph->GlobalAssemble();
        M_Mu_graph->OptimizeStorage();

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
					 dot( M_fluidData->density()*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i)
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
						 value( 0.5 * M_fluidData->viscosity() ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A_graph;
		}
		else
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 M_fluidData->viscosity() * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) )
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
						 dot( M_fluidData->density()*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 value( 0.5 * M_fluidData->viscosity() ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
			) >> M_F_graph;
		}
		else
		{
			buildGraph ( elements (M_fespaceUETA->mesh() ),
						 quadRuleTetra4pt,
						 M_fespaceUETA,
						 M_fespaceUETA,
						 dot ( phi_i, phi_j ) + // mass
						 dot( M_fluidData->density()*value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) + // convective term
						 M_fluidData->viscosity() * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) ) // stiffness
			) >> M_F_graph;
		}
		M_F_graph->GlobalAssemble();
        M_F_graph->OptimizeStorage();
	}

	M_graphIsBuilt = true;

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

		// Graph velocity mass -> block (0,0)
		M_Mu.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Mu_graph ) );
		*M_Mu *= 0;
		integrate ( elements (M_fespaceUETA->mesh() ),
					M_velocityFESpace->qr(),
					M_fespaceUETA,
					M_fespaceUETA,
					M_fluidData->density() * dot ( phi_i, phi_j )
		) >> M_Mu;
		M_Mu->globalAssemble();

		M_Btranspose.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Btranspose_graph ) );
		*M_Btranspose *= 0;
		integrate( elements(M_fespaceUETA->mesh()),
				   M_velocityFESpace->qr(),
				   M_fespaceUETA,
				   M_fespacePETA,
				   value(-1.0) * phi_j * div(phi_i)
		) >> M_Btranspose;

		M_B.reset (new matrix_Type ( M_pressureFESpace->map(), *M_B_graph ) );
		*M_B *= 0;
		integrate( elements(M_fespaceUETA->mesh()),
				   M_pressureFESpace->qr(),
				   M_fespacePETA,
				   M_fespaceUETA,
				   phi_i * div(phi_j)
		) >> M_B;
		M_B->globalAssemble( M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());

		M_A.reset (new matrix_Type ( M_velocityFESpace->map(), *M_A_graph ) );
		*M_A *= 0;
		if ( M_stiffStrain )
		{
			integrate( elements(M_fespaceUETA->mesh()),
					   M_velocityFESpace->qr(),
					   M_fespaceUETA,
					   M_fespaceUETA,
					   value( 0.5 * M_fluidData->viscosity() ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A;
		}
		else
		{
			integrate( elements(M_fespaceUETA->mesh()),
					   M_velocityFESpace->qr(),
					   M_fespaceUETA,
				 	   M_fespaceUETA,
					   M_fluidData->viscosity() * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) )
			) >> M_A;
		}
		M_A->globalAssemble();

		M_C.reset (new matrix_Type ( M_velocityFESpace->map(), *M_C_graph ) );
		*M_C *= 0;
		M_C->globalAssemble();

		M_F.reset (new matrix_Type ( M_velocityFESpace->map(), *M_F_graph ) );
		*M_F *= 0;
		M_F->globalAssemble();


	}

	chrono.stop();
	M_displayer.leaderPrintMax ( " done in ", chrono.diff() ) ;
}

void NavierStokesSolver::updateSystem( const vectorPtr_Type& u_star, const vectorPtr_Type& rhs_velocity )
{
	// Note that u_star HAS to extrapolated from outside. Hence it works also for FSI in this manner.
	M_uExtrapolated.reset( new vector_Type ( *u_star, Repeated ) );

	// Update convective term
	*M_C *= 0;
	{
		using namespace ExpressionAssembly;
		integrate( elements(M_fespaceUETA->mesh()),
				   M_velocityFESpace->qr(),
				   M_fespaceUETA,
				   M_fespaceUETA,
				   dot( M_fluidData->density() * value(M_fespaceUETA, *M_uExtrapolated)*grad(phi_j), phi_i) // semi-implicit treatment of the convective term
		)
		>> M_C;
	}
	M_C->globalAssemble();

	// Get the matrix corresponding to the block (0,0)
	*M_F *= 0;
	*M_F += *M_Mu;
	*M_F *= M_alpha/M_timeStep;
	*M_F += *M_A;
	*M_F += *M_C;

	// Get the right hand side with inertia contribution
	M_rhs.reset( new vector_Type ( M_velocityFESpace->map(), Unique ) );
	*M_rhs = *M_Mu* (*rhs_velocity);
}

void NavierStokesSolver::applyBoundaryConditions ( bcPtr_Type & bc, const Real& time )
{
	updateBCHandler(bc);
	bcManage ( *M_F, *M_rhs, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 1.0, time );

	for (BCHandler::bcBaseIterator_Type it = bc->begin(); it != bc->end(); ++it)
	{
		if ( it->type() != Essential)
			continue;

		bcManageMatrix( *M_Btranspose, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 0.0, 0.0);
	}

	M_Btranspose->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );
}

void NavierStokesSolver::iterate( bcPtr_Type & bc, const Real& time )
{
	applyBoundaryConditions ( bc, time );

    //(1) Set up the OseenOperator
    M_displayer.leaderPrint( "\tNS operator - set up the block operator...");
    LifeChrono chrono;
    chrono.start();

    Operators::NavierStokesOperator::operatorPtrContainer_Type operData(2,2);
    operData(0,0) = M_F->matrixPtr();
    operData(0,1) = M_Btranspose->matrixPtr();
    operData(1,0) = M_B->matrixPtr();
    M_oper->setUp(operData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(2) Set the data for the preconditioner
    M_displayer.leaderPrint( "\tPreconditioner operator - set up the block operator...");
    chrono.reset();
    chrono.start();
    M_prec->setUp(M_F, M_B, M_Btranspose, M_comm);
    M_prec->setDomainMap(M_oper->OperatorDomainBlockMapPtr());
    M_prec->setRangeMap(M_oper->OperatorRangeBlockMapPtr());
    M_prec->updateApproximatedMomentumOperator();
    M_prec->updateApproximatedSchurComplementOperator();
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

void NavierStokesSolver::updateBCHandler( bcPtr_Type & bc )
{
	bc->bcUpdate ( *M_velocityFESpace->mesh(), M_velocityFESpace->feBd(), M_velocityFESpace->dof() );
}

}
