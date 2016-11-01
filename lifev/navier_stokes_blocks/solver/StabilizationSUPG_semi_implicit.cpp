#include <lifev/navier_stokes_blocks/solver/StabilizationSUPG_semi_implicit.hpp>

// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1.0)/( eval(squareroot,TAU_M_DEN) )
#define TAU_M_DEN      TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC
#define TAU_M_DEN_DT   value(M_density*M_density)*value(M_bdfOrder*M_bdfOrder)/value(M_timestep * M_timestep)
#define TAU_M_DEN_VEL  ( value(M_density*M_density)*dot(value(M_fespaceUETA, velocity_extrapolated_rep), G* value(M_fespaceUETA, velocity_extrapolated_rep)) )
#define TAU_M_DEN_VISC ( value(M_C_I)*value(M_viscosity*M_viscosity) *dot (G, G) )

// MACRO TO DEFINE TAU_C
#define TAU_C value(1.0)/( dot(g, TAU_M*g ) )

namespace LifeV
{

//=============================================================================
// Constructor
//=============================================================================

StabilizationSUPG_semi_implicit::StabilizationSUPG_semi_implicit():
		M_label("SUPG_semi_implicit")
{
}

//=============================================================================
// Methods
//=============================================================================

void StabilizationSUPG_semi_implicit::setConstant(const int & value)
{
	if ( value == 1 )
		M_C_I = 30;
	else if ( value == 2 )
		M_C_I = 60;
	else
		ASSERT(0!=0, "Please implement a suitable value for M_C_I for your velocity FE order");
}

// This will be applied to the system matrix
void StabilizationSUPG_semi_implicit::apply_matrix( const vector_Type& velocityExtrapolated )
{
    
	M_block_00.reset (new matrix_Type ( M_uFESpace->map() ) );

	M_block_01.reset (new matrix_Type ( M_uFESpace->map() ) );

	M_block_10.reset (new matrix_Type ( M_pFESpace->map() ) );

	M_block_11.reset (new matrix_Type ( M_pFESpace->map() ) );

	M_block_00->zero();

	M_block_01->zero();

	M_block_10->zero();

	M_block_11->zero();

	vector_Type velocity_extrapolated_rep( velocityExtrapolated, Repeated);
	boost::shared_ptr<SquareRoot_supg_semi_implicit> squareroot(new SquareRoot_supg_semi_implicit());

	using namespace ExpressionAssembly;

	integrate(
			elements(M_uFESpace->mesh()),
			M_uFESpace->qr(),
			M_fespaceUETA, // test  w -> phi_i
			M_fespaceUETA, // trial u^{n+1} -> phi_j

			TAU_M* (
					dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), value(M_density*M_density*M_alpha/M_timestep) * phi_j
							+value(M_density*M_density) * value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_j)
							-value(M_density*M_viscosity)*laplacian(phi_j)
					)

			)

			+ TAU_C*div(phi_i)*div(phi_j)

	) >> M_block_00;

	M_block_00->globalAssemble();

	integrate(
			elements(M_uFESpace->mesh()),
			M_pFESpace->qr(),
			M_fespacePETA, // test  q -> phi_i
			M_fespaceUETA, // trial u^{n+1} -> phi_j

			TAU_M *  (

					dot( grad(phi_i), value(M_density*M_alpha/M_timestep)*phi_j
							+value(M_density)*value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_j)
							-value(M_viscosity)*laplacian(phi_j)
					)

			)

	) >> M_block_10;

	M_block_10->globalAssemble( M_uFESpace->mapPtr(), M_pFESpace->mapPtr() );

	integrate(
			elements(M_uFESpace->mesh()),
			M_uFESpace->qr(),
			M_fespaceUETA, // test  w -> phi_i
			M_fespacePETA, // trial p^{n+1} -> phi_j

			TAU_M*value(M_density)*dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), grad(phi_j) )

	) >> M_block_01;

	M_block_01->globalAssemble( M_pFESpace->mapPtr(), M_uFESpace->mapPtr() );

	integrate(
			elements(M_uFESpace->mesh()),
			M_pFESpace->qr(),
			M_fespacePETA, // test   q -> phi_i
			M_fespacePETA, // trial  p^{n+1} -> phi_j

			TAU_M*dot(  grad(phi_j), grad(phi_i) )

	) >> M_block_11;

	M_block_11->globalAssemble();
    
}

void StabilizationSUPG_semi_implicit::apply_vector( vectorPtr_Type& rhs_velocity,
                                                    vectorPtr_Type& rhs_pressure,
                                                    const vector_Type& velocityExtrapolated,
                                                    const vector_Type& velocity_rhs)
{
	M_rhsVelocity.reset( new vector_Type ( velocity_rhs, Unique ) );

    vector_Type velocity_rhs_rep( velocity_rhs, Repeated);

    vector_Type velocity_extrapolated_rep( velocityExtrapolated, Repeated);

	boost::shared_ptr<SquareRoot_supg_semi_implicit> squareroot(new SquareRoot_supg_semi_implicit());

    using namespace ExpressionAssembly;

    integrate(
    		elements(M_uFESpace->mesh()),
    		M_uFESpace->qr(),
    		M_fespaceUETA,

    		TAU_M*value(M_density*M_density)*dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep) )

    ) >> rhs_velocity;

    integrate(
    		elements(M_uFESpace->mesh()),
    		M_pFESpace->qr(),
    		M_fespacePETA,

    		TAU_M*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep))

    ) >> rhs_pressure;
}

} // namespace LifeV
