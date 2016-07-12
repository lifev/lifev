#include <lifev/navier_stokes_blocks/solver/StabilizationSUPGALE.hpp>

// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1)/( eval(squareroot,TAU_M_DEN) )
#define TAU_M_DEN      TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC
#define TAU_M_DEN_DT   value(M_density*M_density)*value(M_bdfOrder*M_bdfOrder)/value(M_timestep * M_timestep)
#define TAU_M_DEN_VEL  value(M_density*M_density)*dot(value(M_fespaceUETA, beta_km1), G * value(M_fespaceUETA, beta_km1))
#define TAU_M_DEN_VISC ( value(M_C_I)*value(M_viscosity*M_viscosity) *dot (G, G) )

// MACRO TO DEFINE TAU_C
#define TAU_C value(1.0)/( dot(g, TAU_M*g ) )

namespace LifeV
{

//=============================================================================
// Constructor
//=============================================================================

StabilizationSUPGALE::StabilizationSUPGALE():
		M_label("SUPG-ALE")
{
}

//=============================================================================
// Methods
//=============================================================================

void StabilizationSUPGALE::setConstant(const int & value)
{
	if ( value == 1 )
		M_C_I = 30;
	else if ( value == 2 )
		M_C_I = 60;
	else
		ASSERT(0!=0, "Please implement a suitable value for M_C_I for your velocity FE order");
}

void StabilizationSUPGALE::buildGraphs()
{
	// We assume that the fluid velocity and the mesh velocity are discretized by FE of same order
	vector_Type u_km1( M_uFESpace->map(), Repeated);
	vector_Type beta_km1( M_uFESpace->map(), Repeated);
	vector_Type p_km1( M_pFESpace->map(), Repeated);
	vector_Type u_bdf( M_uFESpace->map(), Repeated);

	u_km1.zero();
	beta_km1.zero();
	p_km1.zero();
	u_bdf.zero();

	u_km1 += 1;
	beta_km1 += 1;
	p_km1 += 1;
	u_bdf += 1;

	boost::shared_ptr<SquareRoot_SUPGALE> squareroot(new SquareRoot_SUPGALE());

	MatrixSmall<3, 3> Eye;
	Eye *= 0.0;
	Eye[0][0] = 1;
	Eye[1][1] = 1;
	Eye[2][2] = 1;

	{
		using namespace ExpressionAssembly;

		// Graph for block 00
		M_graph_block00.reset (new Epetra_FECrsGraph (Copy, * (M_uFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements ( M_uFESpace->mesh() ),
					 quadRuleTetra4pt,
					 M_fespaceUETA,
					 M_fespaceUETA,
			/*(1)*/	  TAU_M * value(M_density*M_density) * value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), phi_j )
			/*(2)*/	 +TAU_M * value(M_density*M_density) * value(M_alpha/M_timestep) * dot( phi_j*grad(phi_i), value(M_fespaceUETA, u_km1) )
			/*(3)*/	 -TAU_M * value(M_density*M_density) * dot( phi_j*grad(phi_i), value(M_fespaceUETA, u_bdf) )
			/*(4)*/	 +TAU_M * value(M_density*M_density) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), phi_j*grad(M_fespaceUETA, u_km1) )
			/*(5)*/	 +TAU_M * value(M_density*M_density) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(phi_j) )
			/*(6)*/	 +TAU_M * value(M_density*M_density) * dot( phi_j*grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(M_fespaceUETA, u_km1) )
			/*(8)*/  +TAU_M * value(M_density) * dot( phi_j*grad(phi_i), grad(M_fespacePETA, p_km1) )
			/*(9)*/  -TAU_M * value(M_density*M_viscosity) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), laplacian(phi_j) )
		   /*(10)*/  -TAU_M * value(M_density*M_viscosity) * dot( phi_j*grad(phi_i), laplacian(M_fespaceUETA, u_km1) )
		   /*(16)*/  +TAU_C * div(phi_i)*div(phi_j)
				   ) >> M_graph_block00;
		M_graph_block00->GlobalAssemble();
		M_graph_block00->OptimizeStorage();

		// Graph for block 10
		M_graph_block10.reset (new Epetra_FECrsGraph (Copy, * (M_pFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_pFESpace->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespaceUETA,
			/*(11)*/  TAU_M * value(M_alpha*M_density/M_timestep) * dot( grad(phi_i), phi_j )
			/*(12)*/ +TAU_M * value(M_density) * dot( grad(M_fespaceUETA, u_km1)*grad(phi_i), phi_j )
			/*(13)*/ +TAU_M * value(M_density) * dot( grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(phi_j) )
			/*(15)*/ -TAU_M * value(M_viscosity) * dot( grad(phi_i), laplacian(phi_j) )
				   ) >> M_graph_block10;
		M_graph_block10->GlobalAssemble ( * (M_uFESpace->map().map (Unique) ), * (M_pFESpace->map().map (Unique) ) );
		M_graph_block10->OptimizeStorage();

		// Graph for block 01
		M_graph_block01.reset (new Epetra_FECrsGraph (Copy, * (M_uFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_pFESpace->mesh() ),
					 quadRuleTetra4pt,
					  M_fespaceUETA,
					  M_fespacePETA,
			/*(7)*/	   TAU_M * value(M_density) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), grad(phi_j) )
				   ) >> M_graph_block01;
		M_graph_block01->GlobalAssemble ( * (M_pFESpace->map().map (Unique) ), * (M_uFESpace->map().map (Unique) ) );
		M_graph_block01->OptimizeStorage();

		// Graph for block 11
		M_graph_block11.reset (new Epetra_FECrsGraph (Copy, * (M_pFESpace->map().map (Unique) ), 0) );
		buildGraph ( elements (M_pFESpace->mesh() ),
					 quadRuleTetra4pt,
					 M_fespacePETA,
					 M_fespacePETA,
			/*(14)*/  TAU_M*dot(  grad(phi_i), grad(phi_j) )
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
}

void StabilizationSUPGALE::apply_matrix( const vector_Type& convective_velocity_previous_newton_step,
									 	 const vector_Type& velocity_previous_newton_step,
									 	 const vector_Type& pressure_previous_newton_step,
									 	 const vector_Type& velocity_rhs)
{
	// missing force
    if ( !M_useGraph )
    {
        M_block_00.reset (new matrix_Type ( M_uFESpace->map() ) );
        
        M_block_01.reset (new matrix_Type ( M_uFESpace->map() ) );
        
        M_block_10.reset (new matrix_Type ( M_pFESpace->map() ) );
        
        M_block_11.reset (new matrix_Type ( M_pFESpace->map() ) );
    }
    
	M_block_00->zero();
	M_block_01->zero();
	M_block_10->zero();
	M_block_11->zero();

	vector_Type beta_km1( convective_velocity_previous_newton_step, Repeated);
	vector_Type u_km1( velocity_previous_newton_step, Repeated);
	vector_Type p_km1( pressure_previous_newton_step, Repeated);
	vector_Type u_bdf( velocity_rhs, Repeated);

	boost::shared_ptr<SquareRoot_SUPGALE> squareroot(new SquareRoot_SUPGALE());

	MatrixSmall<3, 3> Eye;
	Eye *= 0.0;
	Eye[0][0] = 1;
	Eye[1][1] = 1;
	Eye[2][2] = 1;

	using namespace ExpressionAssembly;

	integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA, // test  w -> phi_i
				M_fespaceUETA, // trial \delta u -> phi_j
		/*(1)*/	  TAU_M * value(M_density*M_density) * value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), phi_j )
		/*(2)*/	 +TAU_M * value(M_density*M_density) * value(M_alpha/M_timestep) * dot( phi_j*grad(phi_i), value(M_fespaceUETA, u_km1) )
		/*(3)*/	 -TAU_M * value(M_density*M_density) * dot( phi_j*grad(phi_i), value(M_fespaceUETA, u_bdf) )
		/*(4)*/	 +TAU_M * value(M_density*M_density) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), phi_j*grad(M_fespaceUETA, u_km1) )
		/*(5)*/	 +TAU_M * value(M_density*M_density) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(phi_j) )
		/*(6)*/	 +TAU_M * value(M_density*M_density) * dot( phi_j*grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(M_fespaceUETA, u_km1) )
		/*(8)*/  +TAU_M * value(M_density) * dot( phi_j*grad(phi_i), grad(M_fespacePETA, p_km1) )
		/*(9)*/  -TAU_M * value(M_density*M_viscosity) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), laplacian(phi_j) )
	   /*(10)*/  -TAU_M * value(M_density*M_viscosity) * dot( phi_j*grad(phi_i), laplacian(M_fespaceUETA, u_km1) )
	   /*(16)*/  +TAU_C * div(phi_i)*div(phi_j)
			) >> M_block_00;
	M_block_00->globalAssemble();

	integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA, // test  q -> phi_i
				M_fespaceUETA, // trial \delta u -> phi_j
		/*(11)*/  TAU_M * value(M_alpha*M_density/M_timestep) * dot( grad(phi_i), phi_j )
		/*(12)*/ +TAU_M * value(M_density) * dot( grad(phi_i), grad(M_fespaceUETA, u_km1)*phi_j )
		/*(13)*/ +TAU_M * value(M_density) * dot( grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(phi_j) )
		/*(15)*/ -TAU_M * value(M_viscosity) * dot( grad(phi_i), laplacian(phi_j) )
	         ) >> M_block_10;
	M_block_10->globalAssemble( M_uFESpace->mapPtr(), M_pFESpace->mapPtr() );

	integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA, // test  w -> phi_i
				M_fespacePETA, // trial \delta p -> phi_j
		/*(7)*/	  TAU_M * value(M_density) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), grad(phi_j) )
		     ) >> M_block_01;
	M_block_01->globalAssemble( M_pFESpace->mapPtr(), M_uFESpace->mapPtr() );

	integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA, // test   q -> phi_i
				M_fespacePETA, // trial  \delta p -> phi_j
		/*(14)*/  TAU_M*dot(  grad(phi_i), grad(phi_j) )
		     ) >> M_block_11;
	M_block_11->globalAssemble();

}

void StabilizationSUPGALE::apply_vector( vectorPtr_Type& residual_velocity,
								     	 vectorPtr_Type& residual_pressure,
								     	 const vector_Type& convective_velocity_previous_newton_step,
								     	 const vector_Type& velocity_previous_newton_step,
								     	 const vector_Type& pressure_previous_newton_step,
								     	 const vector_Type& velocity_rhs)
{
	// missing force, terms 11 and 12

	vector_Type beta_km1( convective_velocity_previous_newton_step, Repeated);
	vector_Type u_km1( velocity_previous_newton_step, Repeated);
	vector_Type p_km1( pressure_previous_newton_step, Repeated);
	vector_Type u_bdf( velocity_rhs, Repeated);

	boost::shared_ptr<SquareRoot_SUPGALE> squareroot(new SquareRoot_SUPGALE());

	// Matrix needed to evaluate the divergence of a vector (term with TAU_C)
	MatrixSmall<3, 3> Eye;
	Eye *= 0.0;
	Eye[0][0] = 1;
	Eye[1][1] = 1;
	Eye[2][2] = 1;

    using namespace ExpressionAssembly;

	integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA,
	  /*(1)*/	  TAU_M * value(M_density*M_density)*value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), value(M_fespaceUETA, u_km1) )
	  /*(2)*/    -TAU_M * value(M_density*M_density) * dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), value(M_fespaceUETA, u_bdf) )
	  /*(3)*/    +TAU_M*value(M_density*M_density)*dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(M_fespaceUETA, u_km1) )
	  /*(4)*/	 +TAU_M*value(M_density)*dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), grad(M_fespacePETA, p_km1) )
	  /*(5)*/	 -TAU_M*value(M_density*M_viscosity)*dot( value(M_fespaceUETA, beta_km1)*grad(phi_i), laplacian(M_fespaceUETA, u_km1) )
	 /*(11)*/	 +TAU_C*div(phi_i)* trace ( grad(M_fespaceUETA, u_km1) )
			 )
			 >> residual_velocity;

	integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA,
	  /*(6)*/	  TAU_M*value(M_density*M_alpha/M_timestep)*dot( grad(phi_i), value(M_fespaceUETA, u_km1) )
	  /*(7)*/	-TAU_M*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, u_bdf))
	  /*(8)*/   +TAU_M*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, beta_km1)*grad(M_fespaceUETA, u_km1) )
	  /*(9)*/	+TAU_M*dot( grad(phi_i), grad(M_fespacePETA,p_km1) )
	 /*(10)*/   -TAU_M*value(M_viscosity)*dot(grad(phi_i), laplacian(M_fespaceUETA, u_km1))
			 )
		     >> residual_pressure;
}

} // namespace LifeV
