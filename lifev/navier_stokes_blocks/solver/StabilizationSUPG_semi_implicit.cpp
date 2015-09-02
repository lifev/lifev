#include <lifev/navier_stokes_blocks/solver/StabilizationSUPG_semi_implicit.hpp>

// MACRO TO DEFINE TAU_M
#define TAU_M 	       value(1)/( eval(squareroot,TAU_M_DEN) )
#define TAU_M_DEN      TAU_M_DEN_DT + TAU_M_DEN_VEL + TAU_M_DEN_VISC
#define TAU_M_DEN_DT   value(M_density*M_density)*value(M_bdfOrder*M_bdfOrder)/value(M_timestep * M_timestep)
#define TAU_M_DEN_VEL  value(M_density*M_density)*dot(value(M_fespaceUETA, velocity_extrapolated_rep), value(M_fespaceUETA, velocity_extrapolated_rep))/(h_K*h_K)
#define TAU_M_DEN_VISC value(M_C_I)*value(M_viscosity*M_viscosity)/(h_K*h_K*h_K*h_K)

// MACRO TO DEFINE TAU_C
#define TAU_C (h_K*h_K)/(TAU_M)

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

void StabilizationSUPG_semi_implicit::buildGraphs()
{
	vector_Type velocity_extrapolated_rep( M_uFESpace->map(), Repeated);

	velocity_extrapolated_rep.zero();
	
	velocity_extrapolated_rep += 1;
	
	boost::shared_ptr<SquareRoot_supg_semi_implicit> squareroot(new SquareRoot_supg_semi_implicit());

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
}

// This will be applicied to the system matrix
void StabilizationSUPG_semi_implicit::apply_matrix( const vector_Type& velocityExtrapolated )
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

	vector_Type velocity_extrapolated_rep( velocityExtrapolated, Repeated);
	boost::shared_ptr<SquareRoot_supg_semi_implicit> squareroot(new SquareRoot_supg_semi_implicit());

	using namespace ExpressionAssembly;

	integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA, // test  w -> phi_i
				M_fespaceUETA, // trial u^{n+1} -> phi_j
              
		 TAU_M*value(M_density*M_density)*value(M_alpha/M_timestep) * dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), phi_j )
		+TAU_M*value(M_density*M_density) * dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_j) )
        +TAU_C*div(phi_i)*div(phi_j)
        -TAU_M*value(M_density*M_viscosity)*dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), laplacian(phi_j) )
              
			) >> M_block_00;
    
	M_block_00->globalAssemble();

	integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA, // test  q -> phi_i
				M_fespaceUETA, // trial u^{n+1} -> phi_j
              
		 TAU_M*value(M_density*M_alpha/M_timestep)*dot( grad(phi_i), phi_j )
		+TAU_M*value(M_density) * dot( grad(phi_i), value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_j) )
		-TAU_M*value(M_viscosity)*dot(grad(phi_i), laplacian(phi_j))
	         
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
	// missing force, terms 11 and 12

    vector_Type velocity_rhs_rep( velocity_rhs, Repeated);
    vector_Type velocity_extrapolated_rep( velocityExtrapolated, Repeated);

	boost::shared_ptr<SquareRoot_supg_semi_implicit> squareroot(new SquareRoot_supg_semi_implicit());

    using namespace ExpressionAssembly;

	integrate(
				elements(M_uFESpace->mesh()),
				M_uFESpace->qr(),
				M_fespaceUETA,
              
     TAU_M*value(M_density*M_density)*dot( value(M_fespaceUETA, velocity_extrapolated_rep)*grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep) )
	  			 
              )
			 >> rhs_velocity;

	integrate(
				elements(M_uFESpace->mesh()),
				M_pFESpace->qr(),
				M_fespacePETA,
    
     TAU_M*value(M_density)*dot( grad(phi_i), value(M_fespaceUETA, velocity_rhs_rep))
	  
			 )
		     >> rhs_pressure;
}

} // namespace LifeV
