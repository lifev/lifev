#include <lifev/navier_stokes_blocks/solver/StabilizationSUPG_semi_implicit_ale.hpp>

// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1.0)/( eval(squareroot,TAU_M_DEN) )
#define TAU_M_DEN      TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC
#define TAU_M_DEN_DT   value(M_density*M_density)*value(M_bdfOrder*M_bdfOrder)/value(M_timestep * M_timestep)
#define TAU_M_DEN_VEL  value(M_density*M_density)*dot(value(M_fespaceUETA, beta_rep), value(M_fespaceUETA, beta_rep))/(h_K*h_K)
#define TAU_M_DEN_VISC value(M_C_I)*value(M_viscosity*M_viscosity)/(h_K*h_K*h_K*h_K)

#define TAU_M_TILDE 	       value(1.0)/( TAU_M_DEN_DT_NOSQUARED + eval(squareroot,TAU_M_DEN_VEL+TAU_M_DEN_VISC) )
#define TAU_M_DEN_DT_NOSQUARED value(M_density)*value(M_alpha)/value(M_timestep)

// MACRO TO DEFINE TAU_C
#define TAU_C (h_K*h_K)/(TAU_M)

// MACRO BELOW FOR FINE SCALE
#define TAU_M_NO_DT value(1.0)/( eval( squareroot, TAU_M_DEN_VEL + TAU_M_DEN_VISC ) )
#define TAU_C_NO_DT (h_K*h_K)/(TAU_M_NO_DT)

// MACRO FOR EVALUATION OF THE FINE SCALE VELOCITY
#define TAU_M_UPRIME          value(1.0)/( TAU_M_DEN_DT_NOSQUARED + eval(squareroot,TAU_M_DEN_UPRIME) )
#define TAU_M_DEN_UPRIME      TAU_M_DEN_VEL_UPRIME + TAU_M_DEN_VISC_UPRIME
#define TAU_M_DEN_VEL_UPRIME  value(M_density*M_density)*dot(value(M_fespaceUETA, *velocity_repeated), value(M_fespaceUETA, *velocity_repeated))/(h_K*h_K)
#define TAU_M_DEN_VISC_UPRIME value(M_C_I)*value(M_viscosity*M_viscosity)/(h_K*h_K*h_K*h_K)

// MACRO FOR EVALUATION OF THE FINE SCALE PRESSURE
#define TAU_M_NO_DT_UPRIME value(1.0)/( eval( squareroot, TAU_M_DEN_VEL_UPRIME + TAU_M_DEN_VISC_UPRIME ) )
#define TAU_C_NO_DT_PPRIME (h_K*h_K)/(TAU_M_NO_DT_UPRIME)

namespace LifeV
{

//=============================================================================
// Constructor
//=============================================================================

StabilizationSUPG_semi_implicit_ale::StabilizationSUPG_semi_implicit_ale():
		M_label("SUPG_semi_implicit_ale"),
		M_useODEfineScale ( false )
{
}

//=============================================================================
// Methods
//=============================================================================

void StabilizationSUPG_semi_implicit_ale::setConstant(const int & value)
{
	if ( value == 1 )
		M_C_I = 30;
	else if ( value == 2 )
		M_C_I = 60;
	else
		ASSERT(0!=0, "Please implement a suitable value for M_C_I for your velocity FE order");
}

void StabilizationSUPG_semi_implicit_ale::setUseODEfineScale ( const bool& useODEfineScale )
{
	M_useODEfineScale = useODEfineScale;
	setupODEfineScale();
}

void StabilizationSUPG_semi_implicit_ale::setupODEfineScale ( )
{
	int numVolumes = M_uFESpace->mesh()->numVolumes();
	UInt numQuadraturePointsVelocity = M_uFESpace->qr().nbQuadPt();

	std::vector<std::vector<VectorSmall<3>>> velocity_fine;
	velocity_fine.resize(numVolumes);

	M_fineScaleVelocity.resize(numVolumes);
	M_fineScalePressure.resize(numVolumes);
	M_fineScaleVelocityRhs.resize(numVolumes);

	for ( int i = 0; i < numVolumes; ++i )
	{
		velocity_fine[i].resize( numQuadraturePointsVelocity );
		M_fineScaleVelocityRhs[i].resize( numQuadraturePointsVelocity );
		M_fineScaleVelocity[i].resize( numQuadraturePointsVelocity );
		M_fineScalePressure[i].resize( numQuadraturePointsVelocity );

		// Here we assume that same quadrature rule
		// is used for pressure and velocity.
		for ( int j = 0; j < numQuadraturePointsVelocity; ++j )
		{
			velocity_fine[i][j](0) = 0.0;
			velocity_fine[i][j](1) = 0.0;
			velocity_fine[i][j](2) = 0.0;

			M_fineScaleVelocity[i][j](0) = 0.0;
			M_fineScaleVelocity[i][j](1) = 0.0;
			M_fineScaleVelocity[i][j](2) = 0.0;

			M_fineScalePressure[i][j](0) = 0.0;

			M_fineScaleVelocityRhs[i][j](0) = 0.0;
			M_fineScaleVelocityRhs[i][j](1) = 0.0;
			M_fineScaleVelocityRhs[i][j](2) = 0.0;
		}
	}

	M_handlerFineScaleVelocity.reset( new TimeAndExtrapolationHandlerQuadPts<3> ( ) );
	M_handlerFineScaleVelocity->setBDForder(M_bdfOrder);
	M_handlerFineScaleVelocity->setTimeStep(M_timestep);
	std::vector<std::vector<std::vector<VectorSmall<3>>>> initial_state_velocity;
	for ( int i = 0; i < M_bdfOrder; ++i )
	{
		initial_state_velocity.push_back ( velocity_fine );
	}
	M_handlerFineScaleVelocity->initialize(initial_state_velocity);

}

void StabilizationSUPG_semi_implicit_ale::updateODEfineScale ( const vectorPtr_Type& velocity, const vectorPtr_Type& pressure )
{
	computeFineScales(velocity, pressure );
	computeFineScalesForVisualization ( velocity, pressure );
	M_handlerFineScaleVelocity->shift(M_fineScaleVelocity);
	M_handlerFineScaleVelocity->rhsContribution(M_fineScaleVelocityRhs);
}

void StabilizationSUPG_semi_implicit_ale::computeFineScales ( const vectorPtr_Type& velocity, const vectorPtr_Type& pressure )
{
	vectorPtr_Type velocity_repeated( new vector_Type( *velocity, Repeated ) );
	vectorPtr_Type pressure_repeated( new vector_Type( *pressure, Repeated ) );
	vectorPtr_Type velocity_rhs_repeated( new vector_Type( *M_rhsVelocity, Repeated) );

	using namespace ExpressionAssembly;

	boost::shared_ptr<SquareRoot_supg_semi_implicit_ale> squareroot(new SquareRoot_supg_semi_implicit_ale());

	EvaluateAtQuadrature ( elements (  M_uFESpace->mesh() ),
						   M_uFESpace->qr(),
						   M_fespaceUETA,
						   TAU_M_UPRIME * (
	                              value(-M_density*M_alpha/M_timestep)*value(M_fespaceUETA, *velocity_repeated)
	                              +value(M_density)*value(M_fespaceUETA, *velocity_rhs_repeated)
	                              -value(M_density)*value(M_fespaceUETA, *velocity_repeated)*grad(M_fespaceUETA, *velocity_repeated)
	                              -grad(M_fespacePETA, *pressure_repeated)
	                              +value(M_viscosity)*laplacian(M_fespaceUETA, *velocity_repeated)
	                              +value(M_density)*quadpts(M_fespaceUETA, M_fineScaleVelocityRhs )
	                              )
	) >> M_fineScaleVelocity;

	EvaluateAtQuadrature ( elements (  M_uFESpace->mesh() ),
						   M_uFESpace->qr(),
						   M_fespacePETA,
						   value(-1.0) * TAU_C_NO_DT_PPRIME * trace( grad ( M_fespaceUETA, *velocity_repeated ) )
	) >> M_fineScalePressure;
}

//=============================================================================================
void StabilizationSUPG_semi_implicit_ale::computeFineScalesForVisualization ( const vectorPtr_Type& velocity, const vectorPtr_Type& pressure )
{
	 QuadratureRule qr ( quadRuleTetra1pt );

	 using namespace ExpressionAssembly;

	 vectorPtr_Type velocity_repeated( new vector_Type( *velocity, Repeated ) );
	 vectorPtr_Type pressure_repeated( new vector_Type( *pressure, Repeated ) );
	 vectorPtr_Type velocity_rhs_repeated( new vector_Type( *M_rhsVelocity, Repeated) );

	 boost::shared_ptr<SquareRoot_supg_semi_implicit_ale> squareroot(new SquareRoot_supg_semi_implicit_ale());

	 ComputeFineScaleVel ( elements (  M_fespaceUETA->mesh() ),
			 	 	 	 	 	qr,
			 	 	 	 	    M_fespaceUETA,
			 	 	 	 	    TAU_M_UPRIME * (
			 	 	 	 	    		value(-M_density*M_alpha/M_timestep)*value(M_fespaceUETA, *velocity_repeated)
			 	 	 	 	    		+value(M_density)*value(M_fespaceUETA, *velocity_rhs_repeated)
			 	 	 	 	    		-value(M_density)*value(M_fespaceUETA, *velocity_repeated)*grad(M_fespaceUETA, *velocity_repeated)
			 	 	 	 	    		-grad(M_fespacePETA, *pressure_repeated)
			 	 	 	 	    		+value(M_viscosity)*laplacian(M_fespaceUETA, *velocity_repeated)
			 	 	 	 	    		+value(M_density)*quadpts(M_fespaceUETA, M_fineScaleVelocityRhs )
			 	 	 	 	    )
	 ) >> *M_fineVelocity;

	 ComputeFineScalePres ( elements (  M_pFESpace->mesh() ),
			 	 	 	    qr,
			 	 	 	    M_fespacePETA,
			 	 	 	    value(-1.0) * TAU_C_NO_DT_PPRIME * trace( grad ( M_fespaceUETA, *velocity_repeated ) )
	 ) >> *M_finePressure;

}

void StabilizationSUPG_semi_implicit_ale::setExportFineScaleVelocity ( ExporterHDF5<mesh_Type> & exporter, const int& numElementsTotal )
{
	int numVolumes = M_uFESpace->mesh()->numVolumes();
	std::vector<int> id_elem;
	std::vector<int> id_elem_scalar;
	for ( int i = 0; i < numVolumes; ++i )
	{
		id_elem_scalar.push_back ( M_uFESpace->mesh()->element(i).id() );

		for ( int j = 0; j < 3; ++j )
		{
			id_elem.push_back ( M_uFESpace->mesh()->element(i).id() + numElementsTotal * j );
		}
	}

	int* pointerToDofs_scalar (0);
	pointerToDofs_scalar = &id_elem_scalar[0];
	boost::shared_ptr<MapEpetra> map_scalar ( new MapEpetra ( -1, static_cast<int> (id_elem_scalar.size() ), pointerToDofs_scalar, M_uFESpace->map().commPtr() ) );

	int* pointerToDofs (0);
	pointerToDofs = &id_elem[0];
	boost::shared_ptr<MapEpetra> map ( new MapEpetra ( -1, static_cast<int> (id_elem.size() ), pointerToDofs, M_uFESpace->map().commPtr() ) );

	M_fineVelocity.reset ( new vector_Type (*map,  Unique ) );
	M_fineVelocity->zero();

	M_finePressure.reset ( new vector_Type (*map,  Unique ) );
	M_finePressure->zero();

	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "fine scale velocity", M_uFESpace, M_fineVelocity, UInt (0),
			ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
	exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "fine scale pressure", M_pFESpace, M_finePressure, UInt (0),
				ExporterData< mesh_Type >::UnsteadyRegime, ExporterData< mesh_Type >::Cell );
}

void StabilizationSUPG_semi_implicit_ale::buildGraphs()
{
	/*
	vector_Type velocity_extrapolated_rep( M_uFESpace->map(), Repeated);

	velocity_extrapolated_rep.zero();
	
	velocity_extrapolated_rep += 1;
	
	boost::shared_ptr<SquareRoot_supg_semi_implicit_ale> squareroot(new SquareRoot_supg_semi_implicit_ale());

	{
		using namespace ExpressionAssembly;

		// Graph for block 00
		M_graph_block00.reset (new Epetra_FECrsGraph (Copy, * (M_uFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements ( M_uFESpace->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespaceUETA,
           
         TAU_M*value(M_density*M_density)*value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), phi_j )
        +TAU_M*value(M_density*M_density) * dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_j) )
        +TAU_C*div(phi_i)*div(phi_j)
        -TAU_M*value(M_density*M_viscosity)*dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), laplacian(phi_j) )
                    
                    ) >> M_graph_block00;
		
        M_graph_block00->GlobalAssemble();
		M_graph_block00->OptimizeStorage();

		// Graph for block 10
		M_graph_block10.reset (new Epetra_FECrsGraph (Copy, * (M_pFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_pFESpace->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespaceUETA,
                    
         TAU_M*value(M_density*M_alpha/M_timestep)*dot( grad(phi_i), phi_j )
        +TAU_M*value(M_density) * dot( grad(phi_i), value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_j) )
        -TAU_M*value(M_viscosity)*dot(grad(phi_i), laplacian(phi_j))

				   ) >> M_graph_block10;
        
		M_graph_block10->GlobalAssemble ( * (M_uFESpace->map().map (Unique) ), * (M_pFESpace->map().map (Unique) ) );
		M_graph_block10->OptimizeStorage();

		// Graph for block 01
		M_graph_block01.reset (new Epetra_FECrsGraph (Copy, * (M_uFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_pFESpace->mesh() ),
					 quadRuleTetra4pt,
					  M_fespaceUETA,
					  M_fespacePETA,
                    
        TAU_M*value(M_density)*dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), grad(phi_j) )
                    
				   ) >> M_graph_block01;
        
		M_graph_block01->GlobalAssemble ( * (M_pFESpace->map().map (Unique) ), * (M_uFESpace->map().map (Unique) ) );
		M_graph_block01->OptimizeStorage();

		// Graph for block 11
		M_graph_block11.reset (new Epetra_FECrsGraph (Copy, * (M_pFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_pFESpace->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespacePETA,
                    
         TAU_M*dot(  grad(phi_i), grad(phi_j) )
                    
		           ) >> M_graph_block11;
        
		M_graph_block11->GlobalAssemble ( );
		M_graph_block11->OptimizeStorage();

	}

	M_block_00.reset (new matrix_Type ( M_uFESpace->map(), *M_graph_block00 ) );
	M_block_00->zero();
	M_block_00->globalAssemble();

	M_block_01.reset (new matrix_Type ( M_uFESpace->map(), *M_graph_block01 ) );
	M_block_01->zero();
	M_block_01->globalAssemble( M_pFESpace->mapPtr(), M_uFESpace->mapPtr() );

	M_block_10.reset (new matrix_Type ( M_pFESpace->map(), *M_graph_block10 ) );
	M_block_10->zero();
	M_block_10->globalAssemble( M_uFESpace->mapPtr(), M_pFESpace->mapPtr() );

	M_block_11.reset (new matrix_Type ( M_pFESpace->map(), *M_graph_block11 ) );
	M_block_11->zero();
	M_block_11->globalAssemble( );
	*/
}

// This will be applicied to the system matrix
void StabilizationSUPG_semi_implicit_ale::apply_matrix( const vector_Type& velocityExtrapolated, const vector_Type& velocityALE )
{
    if ( !M_useGraph )
    {
        M_block_00.reset (new matrix_Type ( M_uFESpace->map() ) );
        
        M_block_01.reset (new matrix_Type ( M_uFESpace->map() ) );
        
        M_block_10.reset (new matrix_Type ( M_pFESpace->map() ) );
        
        M_block_11.reset (new matrix_Type ( M_pFESpace->map() ) );
    }
    
	// missing force

	M_block_00->zero();
	M_block_01->zero();
	M_block_10->zero();
	M_block_11->zero();

	vector_Type beta( velocityExtrapolated - velocityALE );
	vector_Type beta_rep ( beta, Repeated );

	vector_Type velocity_extrapolated_rep( velocityExtrapolated, Repeated);
	vector_Type velocity_ALE_rep( velocityALE, Repeated);

	boost::shared_ptr<SquareRoot_supg_semi_implicit_ale> squareroot(new SquareRoot_supg_semi_implicit_ale());

	using namespace ExpressionAssembly;

	if ( M_useODEfineScale )
	{
		integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA, // test  w -> phi_i
				M_fespaceUETA, // trial u^{n+1} -> phi_j

				TAU_M_TILDE*value(M_density*M_density)*value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), phi_j )
				+TAU_M_TILDE*value(M_density*M_density) * dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), value(M_fespaceUETA, beta_rep)*grad(phi_j) )
				+TAU_C_NO_DT*div(phi_i)*div(phi_j)
				-TAU_M_TILDE*value(M_density*M_viscosity)*dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), laplacian(phi_j) )

		) >> M_block_00;
	}
	else
	{
		integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA, // test  w -> phi_i
				M_fespaceUETA, // trial u^{n+1} -> phi_j

				TAU_M * (

						dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), value(M_density*M_density*M_alpha/M_timestep) * phi_j
							 											+value(M_density*M_density)*value(M_fespaceUETA, beta_rep)*grad(phi_j)
							 											-value(M_density*M_viscosity)*laplacian(phi_j)
						   )

				)

				+TAU_C*div(phi_i)*div(phi_j)

				/*
				TAU_M*value(M_density*M_density)*value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), phi_j )
				+TAU_M*value(M_density*M_density) * dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), value(M_fespaceUETA, beta_rep)*grad(phi_j) )
				+TAU_C*div(phi_i)*div(phi_j)
				-TAU_M*value(M_density*M_viscosity)*dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), laplacian(phi_j) )
				*/
		) >> M_block_00;
	}

	M_block_00->globalAssemble();

	if ( M_useODEfineScale )
	{
		integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA, // test  q -> phi_i
				M_fespaceUETA, // trial u^{n+1} -> phi_j

				TAU_M_TILDE*value(M_density*M_alpha/M_timestep)*dot( grad(phi_i), phi_j )
				+TAU_M_TILDE*value(M_density) * dot( grad(phi_i), value(M_fespaceUETA, beta_rep)*grad(phi_j) )
				-TAU_M_TILDE*value(M_viscosity)*dot(grad(phi_i), laplacian(phi_j))

		) >> M_block_10;
	}
	else
	{
		integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA, // test  q -> phi_i
				M_fespaceUETA, // trial u^{n+1} -> phi_j

				TAU_M * (
							dot( grad(phi_i), value(M_density*M_alpha/M_timestep)*phi_j
											 +value(M_density)*value(M_fespaceUETA, beta_rep)*grad(phi_j)
											 -value(M_viscosity)*laplacian(phi_j)
							   )

				)

				/*
				TAU_M*value(M_density*M_alpha/M_timestep)*dot( grad(phi_i), phi_j )
				+TAU_M*value(M_density) * dot( grad(phi_i), value(M_fespaceUETA, beta_rep)*grad(phi_j) )
				-TAU_M*value(M_viscosity)*dot(grad(phi_i), laplacian(phi_j))
				*/

		) >> M_block_10;
	}
    M_block_10->globalAssemble( M_uFESpace->mapPtr(), M_pFESpace->mapPtr() );

    if ( M_useODEfineScale )
    {
    	integrate(
    			elements(M_uFESpace->mesh()),
    			M_uFESpace->qr(),
    			M_fespaceUETA, // test  w -> phi_i
    			M_fespacePETA, // trial p^{n+1} -> phi_j

    			TAU_M_TILDE*value(M_density)*dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), grad(phi_j) )

    	) >> M_block_01;
    }
    else
    {
    	integrate(
    			elements(M_uFESpace->mesh()),
    			M_uFESpace->qr(),
    			M_fespaceUETA, // test  w -> phi_i
    			M_fespacePETA, // trial p^{n+1} -> phi_j

    			TAU_M*value(M_density)*dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), grad(phi_j) )

    	) >> M_block_01;
    }
    M_block_01->globalAssemble( M_pFESpace->mapPtr(), M_uFESpace->mapPtr() );

    if ( M_useODEfineScale )
    {
    	integrate(
    			elements(M_uFESpace->mesh()),
    			M_pFESpace->qr(),
    			M_fespacePETA, // test   q -> phi_i
    			M_fespacePETA, // trial  p^{n+1} -> phi_j

    			TAU_M_TILDE*dot(  grad(phi_j), grad(phi_i) )

    	) >> M_block_11;
    }
    else
    {
    	integrate(
    			elements(M_uFESpace->mesh()),
    			M_pFESpace->qr(),
    			M_fespacePETA, // test   q -> phi_i
    			M_fespacePETA, // trial  p^{n+1} -> phi_j

    			TAU_M*dot(  grad(phi_j), grad(phi_i) )

    	) >> M_block_11;
    }
    M_block_11->globalAssemble();
    
}

void StabilizationSUPG_semi_implicit_ale::apply_vector( vectorPtr_Type& rhs_velocity,
                                                    vectorPtr_Type& rhs_pressure,
                                                    const vector_Type& velocityExtrapolated,
                                                    const vector_Type& velocityALE,
                                                    const vector_Type& velocity_rhs)
{
	M_rhsVelocity.reset( new vector_Type ( velocity_rhs, Unique ) );

	// missing force, terms 11 and 12

	vector_Type beta( velocityExtrapolated - velocityALE );
	vector_Type beta_rep ( beta, Repeated );

    vector_Type velocity_rhs_rep( velocity_rhs, Repeated);
    vector_Type velocity_extrapolated_rep( velocityExtrapolated, Repeated);

	boost::shared_ptr<SquareRoot_supg_semi_implicit_ale> squareroot(new SquareRoot_supg_semi_implicit_ale());

    using namespace ExpressionAssembly;

    if ( M_useODEfineScale )
	{
    	integrate(
    			elements(M_uFESpace->mesh()),
    			M_uFESpace->qr(),
    			M_fespaceUETA,

    			TAU_M_TILDE*value(M_density*M_density)*dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep) )
    			+TAU_M_TILDE*value(M_density*M_density)*dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), quadpts(M_fespaceUETA, M_fineScaleVelocityRhs) )

    	) >> rhs_velocity;

    	integrate(
    			elements(M_uFESpace->mesh()),
    			M_pFESpace->qr(),
    			M_fespacePETA,

    			TAU_M_TILDE*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep))
    			+TAU_M_TILDE*value(M_density)*dot( grad(phi_i), quadpts(M_fespaceUETA, M_fineScaleVelocityRhs))

    	) >> rhs_pressure;


	}
	else
	{
		integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA,

				TAU_M*value(M_density*M_density)*dot( value(M_fespaceUETA, beta_rep)*grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep) )

		) >> rhs_velocity;

		integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA,

				TAU_M*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep))

		) >> rhs_pressure;
	}


}

} // namespace LifeV
