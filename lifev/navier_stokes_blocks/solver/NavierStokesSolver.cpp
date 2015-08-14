#include <lifev/navier_stokes_blocks/solver/NavierStokesSolver.hpp>


namespace LifeV
{

NavierStokesSolver::NavierStokesSolver(const dataFile_Type dataFile, const commPtr_Type& communicator):
		M_comm(communicator),
		M_dataFile(dataFile),
		M_displayer(communicator),
		M_graphIsBuilt(false),
        M_oper(new Operators::NavierStokesOperator),
        M_operLoads(new Operators::NavierStokesOperator),
        M_invOper(),
        M_fullyImplicit(false),
        M_graphPCDisBuilt(false),
        M_steady ( dataFile("fluid/miscellaneous/steady", false) ),
        M_density ( dataFile("fluid/physics/density", 1.0 ) ),
        M_viscosity ( dataFile("fluid/physics/viscosity", 0.035 ) ),
        M_useStabilization ( dataFile("fluid/stabilization/use", false) ),
        M_stabilizationType ( dataFile("fluid/stabilization/type", "none") ),
        M_computeAerodynamicLoads ( dataFile("fluid/compute_loads", false) ),
        M_nonconforming ( dataFile("fluid/nonconforming", false) )
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

void NavierStokesSolver::setup(const meshPtr_Type& mesh, const int& id_domain)
{
	std::string uOrder;
	std::string pOrder;

	if ( !M_nonconforming )
	{
		uOrder = M_dataFile("fluid/space_discretization/vel_order","P1");
		pOrder = M_dataFile("fluid/space_discretization/pres_order","P1");
		M_stiffStrain = M_dataFile("fluid/space_discretization/stiff_strain", true);

	}
	else
	{
		if ( id_domain == 0 )
		{
			uOrder = M_dataFile("fluid/master_space_discretization/vel_order","P1");
			pOrder = M_dataFile("fluid/master_space_discretization/pres_order","P1");
			M_stiffStrain = M_dataFile("fluid/master_space_discretization/stiff_strain", true);
		}
		else if ( id_domain == 1 )
		{
			uOrder = M_dataFile("fluid/master_space_discretization/vel_order","P1");
			pOrder = M_dataFile("fluid/slave_space_discretization/pres_order","P1");
			M_stiffStrain = M_dataFile("fluid/slave_space_discretization/stiff_strain", true);
		}
	}

	M_useGraph = M_dataFile("fluid/use_graph", true);

	Int geoDimensions = mesh_Type::S_geoDimensions;

	M_velocityFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, uOrder, geoDimensions, M_comm) );
	M_pressureFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, pOrder, 1, M_comm) );
	M_velocityFESpaceScalar.reset (new FESpace<mesh_Type, map_Type> (mesh, uOrder, 1, M_comm) );

	M_fespaceUETA.reset( new ETFESpace_velocity(M_velocityFESpace->mesh(), &(M_velocityFESpace->refFE()), M_comm));
	M_fespacePETA.reset( new ETFESpace_pressure(M_pressureFESpace->mesh(), &(M_pressureFESpace->refFE()), M_comm));

	M_uExtrapolated.reset( new vector_Type ( M_velocityFESpace->map(), Repeated ) );
	M_uExtrapolated->zero();
    
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

    M_monolithicMap.reset( new map_Type ( M_velocityFESpace->map() ) );
    *M_monolithicMap += M_pressureFESpace->map();
    
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

    	M_solution.reset( new vector_Type ( *M_monolithicMap ) );
    	M_solution->zero();
    }

    if ( M_useStabilization )
    {
    	M_stabilization.reset ( StabilizationFactory::instance().createObject ( M_dataFile("fluid/stabilization/type","none") ) );
    	M_stabilization->setVelocitySpace(M_velocityFESpace);
    	M_stabilization->setPressureSpace(M_pressureFESpace);
    	M_stabilization->setDensity(M_density);
    	M_stabilization->setViscosity(M_viscosity);
    	int vel_order = M_dataFile ( "fluid/stabilization/vel_order", 1 );
    	M_stabilization->setConstant( vel_order );
    	M_stabilization->setBDForder ( M_dataFile ( "fluid/time_discretization/BDF_order", 1 ) );
    	M_stabilization->setCommunicator ( M_velocityFESpace->map().commPtr() );
    	M_stabilization->setTimeStep ( M_dataFile ( "fluid/time_discretization/timestep", 0.001 ) );
    	M_stabilization->setETvelocitySpace ( M_fespaceUETA );
    	M_stabilization->setETpressureSpace ( M_fespacePETA );
        M_stabilization->setUseGraph( M_useGraph );
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

        if ( M_useStabilization )
        {
        	M_stabilization->buildGraphs();
        }

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
	if ( M_useGraph && !M_graphIsBuilt )
		buildGraphs();

	M_displayer.leaderPrint ( " F - Assembling constant terms... ");
	LifeChrono chrono;
	chrono.start();
    
    if ( M_useGraph )
    {
        M_Mu.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Mu_graph ) );
        M_Mu->zero();
        
        M_Btranspose.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Btranspose_graph ) );
        M_Btranspose->zero();
        
        M_B.reset (new matrix_Type ( M_pressureFESpace->map(), *M_B_graph ) );
        M_B->zero();
        
        M_A.reset (new matrix_Type ( M_velocityFESpace->map(), *M_A_graph ) );
        M_A->zero();
        
        M_C.reset (new matrix_Type ( M_velocityFESpace->map(), *M_C_graph ) );
        M_C->zero();
        
        M_F.reset (new matrix_Type ( M_velocityFESpace->map(), *M_F_graph ) );
        M_F->zero();
        
        M_Jacobian.reset (new matrix_Type ( M_velocityFESpace->map(), *M_Jacobian_graph ) );
        M_Jacobian->zero();
    }
    else
    {
        M_Mu.reset (new matrix_Type ( M_velocityFESpace->map() ) );
        M_Mu->zero();
        
        M_Btranspose.reset (new matrix_Type ( M_velocityFESpace->map() ) );
        M_Btranspose->zero();
        
        M_B.reset (new matrix_Type ( M_pressureFESpace->map() ) );
        M_B->zero();
        
        M_A.reset (new matrix_Type ( M_velocityFESpace->map() ) );
        M_A->zero();
        
        M_C.reset (new matrix_Type ( M_velocityFESpace->map() ) );
        M_C->zero();
        
        M_F.reset (new matrix_Type ( M_velocityFESpace->map() ) );
        M_F->zero();
        
        M_Jacobian.reset (new matrix_Type ( M_velocityFESpace->map() ) );
        M_Jacobian->zero();
    }

	{
		using namespace ExpressionAssembly;

		if ( !M_steady )
		{
			// Velocity mass -> block (0,0)
			integrate ( elements (M_fespaceUETA->mesh() ),
						M_velocityFESpace->qr(),
						M_fespaceUETA,
						M_fespaceUETA,
						M_density * dot ( phi_i, phi_j )
					  ) >> M_Mu;
			M_Mu->globalAssemble();
		}

		integrate( elements(M_fespaceUETA->mesh()),
				   M_velocityFESpace->qr(),
				   M_fespaceUETA,
				   M_fespacePETA,
				   value(-1.0) * phi_j * div(phi_i)
		         ) >> M_Btranspose;
        
		M_Btranspose->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );

		integrate( elements(M_fespaceUETA->mesh()),
				   M_pressureFESpace->qr(),
				   M_fespacePETA,
				   M_fespaceUETA,
				   phi_i * div(phi_j)
                  ) >> M_B;
		
        M_B->globalAssemble( M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());

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
		
		M_block01->zero();
		M_block10->zero();
		*M_block01 += *M_Btranspose;
		*M_block10 += *M_B;
        
        if ( !M_useStabilization )
        {
            M_block01->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );
            M_block10->globalAssemble(M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());
        }
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

		QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

		integrate( boundary(M_fespaceUETA->mesh(), 6),
				   myBDQR,
				   M_fespaceUETA,
				   M_fespaceUETA,
				    value(-1.0)*value(M_viscosity)*dot( phi_i, ( grad(phi_j) + transpose( grad(phi_j) ) ) * Nface )
				   -value(M_viscosity)*dot( ( grad(phi_i) + transpose( grad(phi_i) ) ) * Nface, phi_j )
				   -value(M_density)*dot( phi_i, dot( value(M_fespaceUETA, *M_uExtrapolated), Nface ) * phi_j )
				   +value(30.0)*value(M_viscosity)/value(0.005)*dot( phi_i - dot( phi_i, Nface)*Nface, phi_j - dot( phi_j, Nface)*Nface )
				   +value(30.0)*value(M_viscosity)/value(0.005)*( dot( phi_i, Nface)*dot( phi_j, Nface) )
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
    if ( !M_useStabilization )
        M_block00->globalAssemble();

	// Get the right hand side with inertia contribution
	M_rhs.reset( new vector_Type ( M_velocityFESpace->map(), Unique ) );
	M_rhs->zero();

	if ( !M_steady )
		*M_rhs = *M_Mu* (*rhs_velocity);
    
    if ( !M_fullyImplicit && M_useStabilization )
    {
        M_displayer.leaderPrint ( "\tF - Assembling semi-implicit stabilization terms... ");
        LifeChrono chrono;
        chrono.start();
        
        if ( M_stabilizationType.compare("SUPG_SEMI_IMPLICIT")==0  )
            M_stabilization->apply_matrix( *u_star );
        else if ( M_stabilizationType.compare("VMSLES_SEMI_IMPLICIT")==0 )
            M_stabilization->apply_matrix( *u_star, *M_pressure_extrapolated, *rhs_velocity);
        
        *M_block00 += *M_stabilization->block_00();
        M_block00->globalAssemble();
        
        M_block01->zero();
        *M_block01 += *M_Btranspose;
        *M_block01 += *M_stabilization->block_01();

        {
        	using namespace ExpressionAssembly;

        	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

        	integrate( boundary(M_fespaceUETA->mesh(), 6),
        			   myBDQR,
        			   M_fespaceUETA,
        			   M_fespacePETA,
        			   dot( phi_i, phi_j * Nface)
        	)
        	>> M_block01;
        }

        M_block01->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );
        
        M_block10->zero();
        *M_block10 += *M_B;
        *M_block10 += *M_stabilization->block_10();

        {
        	using namespace ExpressionAssembly;

        	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

        	integrate( boundary(M_fespaceUETA->mesh(), 6),
        			   myBDQR,
        			   M_fespacePETA,
        			   M_fespaceUETA,
        			   value(-1.0)*dot(phi_i*Nface, phi_j)
        	)
        	>> M_block10;
        }

        M_block10->globalAssemble(M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());
        
        M_block11->zero();
        *M_block11 += *M_stabilization->block_11();
        M_block11->globalAssemble();
        
        M_rhs_pressure.reset( new vector_Type ( M_pressureFESpace->map(), Unique ) );
        M_rhs_pressure->zero();

        M_stabilization->apply_vector( M_rhs, M_rhs_pressure, *u_star, *rhs_velocity);
        
        M_rhs->globalAssemble();
        M_rhs_pressure->globalAssemble();
        
        chrono.stop();
        M_displayer.leaderPrintMax ( " done in ", chrono.diff() ) ;

    }
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
	if ( M_computeAerodynamicLoads )
	{
		M_displayer.leaderPrint( "\tNS operator - Compute Loads: TRUE");

		M_block00_noBC.reset( new matrix_Type(M_velocityFESpace->map() ) );
		M_block01_noBC.reset( new matrix_Type(M_velocityFESpace->map() ) );

		M_block00_noBC->zero();
		M_block01_noBC->zero();

		*M_block00_noBC += *M_block00;
		*M_block01_noBC += *M_block01;

        M_block00_noBC->globalAssemble();
        M_block01_noBC->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );
        
		// Right hand side
		M_rhs_noBC.reset( new vector_Type ( M_velocityFESpace->map(), Unique ) );
		M_rhs_noBC->zero();
		*M_rhs_noBC += *M_rhs;
    }

	updateBCHandler(bc);
	bcManage ( *M_block00, *M_rhs, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 1.0, time );
    bcManageMatrix( *M_block01, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 0.0, 0.0);
    //M_block01->spy("blocco01");
    //M_block01->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );
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
    if ( M_useStabilization )
        operData(1,1) = M_block11->matrixPtr();
    
    M_oper->setUp(operData, M_displayer.comm());
    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    //(2) Set the data for the preconditioner

	M_displayer.leaderPrint( "\tPreconditioner operator - set up the block operator...");
    chrono.reset();
    chrono.start();

    if ( std::strcmp(M_prec->Label(),"aSIMPLEOperator")==0 )
    {
        if (M_useStabilization)
        {
            M_prec->setUp(M_block00, M_block10, M_block01, M_block11);
        }
        else
        {
            M_prec->setUp(M_block00, M_block10, M_block01);
        }

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

    M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));
    M_invOper->setOperator(M_oper);
    M_invOper->setPreconditioner(M_prec);

    chrono.stop();
    M_displayer.leaderPrintMax(" done in " , chrono.diff() );

    // Solving the system
    BlockEpetra_Map upMap;
    upMap.setUp ( M_velocityFESpace->map().map(Unique), M_pressureFESpace->map().map(Unique));

    if ( M_useStabilization )
    {
        vector_Type rhs ( *M_monolithicMap, Unique );
        vector_Type sol ( *M_monolithicMap, Unique );
     
        rhs.zero();
        sol.zero();
        
        rhs.subset ( *M_rhs, M_velocityFESpace->map(), 0, 0 );
        rhs.subset ( *M_rhs_pressure, M_pressureFESpace->map(), 0, M_velocityFESpace->map().mapSize() );
        
        M_invOper->ApplyInverse(rhs.epetraVector(), sol.epetraVector());
        
        M_velocity->subset ( sol, M_velocityFESpace->map(), 0, 0 );
        M_pressure->subset ( sol, M_pressureFESpace->map(), M_velocityFESpace->map().mapSize(), 0 );

        if ( M_computeAerodynamicLoads )
        {
        	// Computation of aerodynamic forces

        	/* Residual form, see [Forti, Dede 2015] */
        	/*

			M_forces.reset ( new vector_Type ( M_velocityFESpace->map(), Unique ) );
        	M_forces->zero();

            *M_forces += *M_rhs_noBC;
            
            *M_forces -= *M_block00_noBC * ( *M_velocity );
            
            *M_forces -= *M_block01_noBC * ( *M_pressure );
            */

        	/* Direct evaluation by integrals */

        	vector_Type forces_rep( M_velocityFESpace->map(), Repeated );
        	forces_rep.zero();

        	vector_Type velo_rep ( *M_velocity, Repeated);
        	vector_Type pres_rep ( *M_pressure, Repeated);

        	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
        	MatrixSmall<3, 3> Eye;
        	Eye *= 0.0;
        	Eye[0][0] = 1;
        	Eye[1][1] = 1;
        	Eye[2][2] = 1;

        	{
        		using namespace ::LifeV::ExpressionAssembly;

        		integrate ( boundary (M_fespaceUETA->mesh(), 6), // Note that 6 is for the toy problem, hardcoded! TODO: generalize to generic flag
        		            myBDQR,
        		            M_fespaceUETA,
        		             value(M_viscosity)* dot( phi_i, ( grad(M_fespaceUETA, velo_rep) + transpose( grad(M_fespaceUETA, velo_rep) ) ) * Nface )
        		            -dot( phi_i, value(M_fespacePETA, pres_rep) * Nface )
        		          ) >> forces_rep;
        	}

        	forces_rep.globalAssemble();

        	M_forces.reset ( new vector_Type ( forces_rep, Unique ) );
        	*M_forces *= -1;
        }
    }
    else
    {
        BlockEpetra_MultiVector up(upMap, 1), rhs(upMap, 1);
        rhs.block(0).Update(1.0, M_rhs->epetraVector(), 0.);
    
        // Solving the linear system
        M_invOper->ApplyInverse(rhs,up);

        M_velocity->epetraVector().Update(1.0,up.block(0),0.0);
        M_pressure->epetraVector().Update(1.0,up.block(1),0.0);
    }
}

VectorSmall<2> NavierStokesSolver::computeForces( BCHandler& bcHDrag,
												  BCHandler& bcHLift )
{
	bcHDrag.bcUpdate ( *M_velocityFESpace->mesh(), M_velocityFESpace->feBd(), M_velocityFESpace->dof() );
	bcHLift.bcUpdate ( *M_velocityFESpace->mesh(), M_velocityFESpace->feBd(), M_velocityFESpace->dof() );

	vector_Type onesOnBodyDrag(M_velocityFESpace->map(), Unique);
	onesOnBodyDrag.zero();

	vector_Type onesOnBodyLift(M_velocityFESpace->map(), Unique);
	onesOnBodyLift.zero();

	bcManageRhs ( onesOnBodyDrag, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),  bcHDrag, M_velocityFESpace->feBd(), 1., 0.);
	bcManageRhs ( onesOnBodyLift, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(),  bcHLift, M_velocityFESpace->feBd(), 1., 0.);

	Real drag (0.0);
	Real lift (0.0);

	drag = M_forces->dot(onesOnBodyDrag);
	lift = M_forces->dot(onesOnBodyLift);

	VectorSmall<2> Forces;
	Forces[0] = drag;
	Forces[1] = lift;

	return Forces;
}

void NavierStokesSolver::iterate_nonlinear( bcPtr_Type & bc, const Real& time )
{
	applyBoundaryConditionsSolution ( bc, time ); // the second argument is zero since the problem is steady

	// Call Newton
	UInt status = NonLinearRichardson ( *M_solution, *this, M_absoluteTolerance, M_relativeTolerance, M_maxiterNonlinear, M_etaMax,
			M_nonLinearLineSearch, 0, 2, M_out_res, 0.0);

    if ( M_computeAerodynamicLoads )
    {
      M_forces.reset ( new vector_Type ( *M_monolithicMap, Unique ) );
      M_forces->zero();
      computeForcesNonLinear(M_forces, M_solution);
    }

    if (status == EXIT_FAILURE)
    {
        M_displayer.leaderPrint(" WARNING: Newton  failed " );
    }
    
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
					   	   	   	   	   	   vectorPtr_Type& residualVelocity,
					   	   	   	   	   	   vectorPtr_Type& residualPressure)
{
	residualVelocity->zero();
	residualPressure->zero();

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
			) >> residualVelocity;
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
			) >> residualVelocity;
		}

		integrate ( elements ( M_fespaceUETA->mesh() ),
				M_pressureFESpace->qr(),
				M_fespacePETA,
				trace ( grad ( M_fespaceUETA, *u_km1_repeated ) ) * phi_i
		) >> residualPressure;
	}

	if ( M_useStabilization )
	{
		M_displayer.leaderPrint ( "[F] - Assembly residual of the stabilization \n" ) ;
		M_stabilization->apply_vector(residualVelocity, residualPressure, *convective_velocity, *velocity_km1, *pressure_km1, *rhs_velocity);
	}

	residualVelocity->globalAssemble();
	residualPressure->globalAssemble();
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

	if ( M_useStabilization )
	{
		M_displayer.leaderPrint ( "[F] - Assembly residual of the stabilization \n" ) ;
		M_stabilization->apply_vector(res_velocity, res_pressure, *convective_velocity, *velocity_km1, *pressure_km1, *rhs_velocity);
	}

	res_velocity->globalAssemble();
	res_pressure->globalAssemble();

	residual->subset ( *res_velocity, M_velocityFESpace->map(), 0, 0 );
	residual->subset ( *res_pressure, M_pressureFESpace->map(), 0, M_velocityFESpace->map().mapSize() );

}

void NavierStokesSolver::assembleInterfaceMass( matrixPtr_Type& mass_interface, const mapPtr_Type& interface_map,
												markerID_Type interfaceFlag, const vectorPtr_Type& numerationInterface,
												const UInt& offset)
{
	// INITIALIZE MATRIX WITH THE MAP OF THE INTERFACE
	matrixPtr_Type fluid_interfaceMass( new matrix_Type ( M_velocityFESpace->map(), 50 ) );
	fluid_interfaceMass->zero();

	// ASSEMBLE MASS MATRIX AT THE INTERFACE
	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
	{
		using namespace ExpressionAssembly;
		integrate ( boundary (M_fespaceUETA->mesh(), interfaceFlag ),
					myBDQR,
					M_fespaceUETA,
					M_fespaceUETA,
					//Boundary Mass
					dot(phi_i,phi_j)
				  )
				  >> fluid_interfaceMass;
	}
	fluid_interfaceMass->globalAssemble();

	// RESTRICT MATRIX TO INTERFACE DOFS ONLY
	mass_interface.reset(new matrix_Type ( *interface_map, 50 ) );
	fluid_interfaceMass->restrict ( interface_map, numerationInterface, offset, mass_interface );
}

void NavierStokesSolver::computeForcesNonLinear(vectorPtr_Type& force, const vectorPtr_Type& solution)
{
  // Forces to zero                                                                                                                                   
  force->zero();

  // Extract the velocity and the pressure from the solution vector                                                     
  vectorPtr_Type u_km1 ( new vector_Type ( M_velocityFESpace->map(), Repeated ) );
  vectorPtr_Type p_km1 ( new vector_Type ( M_pressureFESpace->map(), Repeated ) );
  u_km1->zero();
  p_km1->zero();
  u_km1->subset ( *solution, M_velocityFESpace->map(), 0, 0 );
  p_km1->subset ( *solution, M_pressureFESpace->map(), M_velocityFESpace->map().mapSize(), 0 );

  // Force vector for the velocity and pressure components                                                                                           
  vectorPtr_Type res_velocity ( new vector_Type ( M_velocityFESpace->map(), Unique ) );
  vectorPtr_Type res_pressure ( new vector_Type ( M_pressureFESpace->map(), Unique ) );
  res_velocity->zero();
  res_pressure->zero();

  vectorPtr_Type rhs_velocity_repeated;
  rhs_velocity_repeated.reset( new vector_Type ( *M_velocityRhs, Repeated ) );

  using namespace ExpressionAssembly;

  if ( M_stiffStrain )
  {
      integrate ( elements ( M_fespaceUETA->mesh() ),
		  M_velocityFESpace->qr(),
		  M_fespaceUETA,
		  M_density * value ( M_alpha / M_timeStep ) * dot ( value ( M_fespaceUETA, *u_km1 ), phi_i ) -
		  M_density * dot ( value ( M_fespaceUETA, *rhs_velocity_repeated ), phi_i ) +
		  value ( 0.5 ) * M_viscosity * dot ( grad ( phi_i )  + transpose ( grad ( phi_i ) ), grad ( M_fespaceUETA, *u_km1 ) + \
						      transpose ( grad ( M_fespaceUETA, *u_km1 ) ) ) +
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

  if ( M_useStabilization )
  {
      M_stabilization->apply_vector(res_velocity, res_pressure, *u_km1, *p_km1, *rhs_velocity_repeated);
  }

  res_velocity->globalAssemble();
  res_pressure->globalAssemble();

  force->subset ( *res_velocity, M_velocityFESpace->map(), 0, 0 );
  force->subset ( *res_pressure, M_pressureFESpace->map(), 0, M_velocityFESpace->map().mapSize() );
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
    
	if ( M_useStabilization )
	{
		M_displayer.leaderPrint ( "[F] - Assembly residual of the stabilization \n" ) ;
		M_stabilization->apply_vector(res_velocity, res_pressure, *u_km1, *p_km1, *rhs_velocity_repeated);
	}

    res_velocity->globalAssemble();
    res_pressure->globalAssemble();

    bcManageRhs ( *res_velocity, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *M_bc, M_velocityFESpace->feBd(), 0.0, 0.0 );

    residual.subset ( *res_velocity, M_velocityFESpace->map(), 0, 0 );
    residual.subset ( *res_pressure, M_pressureFESpace->map(), 0, M_velocityFESpace->map().mapSize() );

    // We need to update the linearized terms in the jacobian coming from the convective part

    updateConvectiveTerm(u_km1);
    updateJacobian(*u_km1);

    if ( M_useStabilization )
    {
    	M_stabilization->apply_matrix(*u_km1, *p_km1, *M_velocityRhs);

    	*M_block00 += *M_stabilization->block_00();
    	M_block00->globalAssemble();

    	M_block01->zero();
    	*M_block01 += *M_Btranspose;
    	*M_block01 += *M_stabilization->block_01();
    	M_block01->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );

    	M_block10->zero();
    	*M_block10 += *M_B;
    	*M_block10 += *M_stabilization->block_10();
    	M_block10->globalAssemble(M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());

    	M_block11->zero();
    	*M_block11 += *M_stabilization->block_11();
    	M_block11->globalAssemble();

//    	M_block00->spy("block00");
//    	M_block01->spy("block01");
//    	M_block10->spy("block10");
//    	M_block11->spy("block11");
    }

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
	if ( !M_useStabilization )
		M_block00->globalAssemble( );
}

void NavierStokesSolver::updateStabilization( const vector_Type& convective_velocity_previous_newton_step,
	 	 	 	 	 	 	 	 	 	 	  const vector_Type& velocity_previous_newton_step,
	 	 	 	 	 	 	 	 	 	 	  const vector_Type& pressure_previous_newton_step,
	 	 	 	 	 	 	 	 	 	 	  const vector_Type& velocity_rhs )
{
	M_displayer.leaderPrint ( "[F] - Update Jacobian stabilization terms\n" ) ;
	M_stabilization->apply_matrix(convective_velocity_previous_newton_step,
								  velocity_previous_newton_step,
								  pressure_previous_newton_step,
								  velocity_rhs);

	*M_block00 += *M_stabilization->block_00();
	M_block00->globalAssemble();

	M_block01->zero();
	*M_block01 += *M_Btranspose;
	*M_block01 += *M_stabilization->block_01();
	M_block01->globalAssemble( M_pressureFESpace->mapPtr(), M_velocityFESpace->mapPtr() );

	M_block10->zero();
	*M_block10 += *M_B;
	*M_block10 += *M_stabilization->block_10();
	M_block10->globalAssemble(M_velocityFESpace->mapPtr(), M_pressureFESpace->mapPtr());

	M_block11->zero();
	*M_block11 += *M_stabilization->block_11();
	M_block11->globalAssemble();
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
	if ( M_useStabilization )
		operDataJacobian ( 1, 1 ) = M_block11->matrixPtr();

	M_oper->setUp(operDataJacobian, M_displayer.comm());
	chrono.stop();
	M_displayer.leaderPrintMax(" done in " , chrono.diff() );

	M_displayer.leaderPrint( "\tPreconditioner operator - set up the block operator...");
	chrono.reset();
	chrono.start();

	//(2) Set the data for the preconditioner
	if ( std::strcmp(M_prec->Label(),"aSIMPLEOperator")==0 )
	{
		if (M_useStabilization)
		{
			M_prec->setUp(M_block00, M_block10, M_block01, M_block11);
		}
		else
		{
			M_prec->setUp(M_block00, M_block10, M_block01);
		}

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
    M_invOper->setParameterList(M_pListLinSolver->sublist(solverType));

    chrono.stop();
	M_displayer.leaderPrintMax(" done in " , chrono.diff() );

	M_invOper->setOperator(M_oper);
	M_invOper->setPreconditioner(M_prec);
	
	increment.zero();

	// Solving the linear system
	M_invOper->ApplyInverse(residual.epetraVector(), increment.epetraVector());
}

void NavierStokesSolver::applyBoundaryConditionsJacobian ( bcPtr_Type & bc )
{
	bcManageMatrix( *M_block00, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 1.0, 0.0);
	bcManageMatrix( *M_block01, *M_velocityFESpace->mesh(), M_velocityFESpace->dof(), *bc, M_velocityFESpace->feBd(), 0.0, 0.0);
}

void NavierStokesSolver::integrateForces ( const vectorPtr_Type & velocity, const vectorPtr_Type & pressure)
{
	vector_Type forces_rep( M_velocityFESpace->map(), Repeated );
	forces_rep.zero();

	vector_Type velo_rep ( *velocity, Repeated);
	vector_Type pres_rep ( *pressure, Repeated);

	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
	MatrixSmall<3, 3> Eye;
	Eye *= 0.0;
	Eye[0][0] = 1;
	Eye[1][1] = 1;
	Eye[2][2] = 1;

	{
		using namespace ::LifeV::ExpressionAssembly;

		integrate ( boundary (M_fespaceUETA->mesh(), 6), // Note that 6 is for the toy problem, hardcoded! TODO: generalize to generic flag
				myBDQR,
				M_fespaceUETA,
				value(M_viscosity)* dot( phi_i, ( grad(M_fespaceUETA, velo_rep) + transpose( grad(M_fespaceUETA, velo_rep) ) ) * Nface )
				-dot( phi_i, value(M_fespacePETA, pres_rep) * Nface )
		) >> forces_rep;
	}

	forces_rep.globalAssemble();

	M_forces.reset ( new vector_Type ( forces_rep, Unique ) );
	*M_forces *= -1;
}


}
