#include <lifev/fsi_blocks/solver/LinearElasticity.hpp>

namespace LifeV
{

LinearElasticity::LinearElasticity( const commPtr_Type& communicator ):
		M_comm(communicator),
		M_density(0.0),
		M_young(0.0),
		M_poisson(0.0),
		M_lambda(0.0),
		M_mu(0.0),
		M_thinLayer(false),
		M_thinLayerThickness(0.0),
		M_thinLayerDensity(0.0),
		M_thinLayerLameI(0.0),
		M_thinLayerLameII(0.0),
		M_interfaceFlag(0.0)
{}

LinearElasticity::~LinearElasticity()
{}

void
LinearElasticity::setCoefficients ( const Real density, const Real young, const Real poisson)
{
	// Copying the parameters

	M_density = density;

	M_young = young;

	M_poisson = poisson;

	// Evaluation of the lame coefficients

	M_lambda = ( M_young * M_poisson ) / ( ( 1.0 + M_poisson ) * ( 1.0 - 2.0 * M_poisson ) );

	M_mu = M_young / ( 2.0 * ( 1.0 + M_poisson ) );

}

void
LinearElasticity::setCoefficientsThinLayer ( const Real density, const Real young, const Real poisson, const Real thickness, const UInt interface )
{
	M_thinLayer = true;
	M_thinLayerThickness = thickness;
	M_thinLayerDensity = density;
	M_interfaceFlag = interface;
	M_thinLayerLameI = young / ( 2 * ( 1.0 + poisson ) );
	M_thinLayerLameII = ( young * poisson ) / ( ( 1.0 + poisson ) * ( 1.0 - poisson ) );
	M_mass_no_bc_thin.reset( new matrix_Type ( M_displacementFESpace->map() ) );
	M_stiffness_no_bc_thin.reset( new matrix_Type ( M_displacementFESpace->map() ) );
}

void
LinearElasticity::setup( const meshPtr_Type& mesh, const std::string dOrder)
{
	// Create FE spaces

	M_dOrder = dOrder;

	M_displacementFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, M_dOrder, 3, M_comm) );

	M_displacementFESpace_ETA.reset( new ETFESpace_displacement ( M_displacementFESpace->mesh(), &(M_displacementFESpace->refFE()), M_comm));

	// Create matrices

	M_mass_no_bc.reset( new matrix_Type(M_displacementFESpace->map() ) );

	M_stiffness_no_bc.reset( new matrix_Type(M_displacementFESpace->map() ) );

	M_jacobian.reset( new matrix_Type(M_displacementFESpace->map() ) );

	// Initialize

	M_mass_no_bc->zero();

	M_stiffness_no_bc->zero();

	M_jacobian->zero();
}

void
LinearElasticity::assemble_matrices ( const Real timestep, const Real coeff, bcPtr_Type & bc, bool useBDF )
{
	ASSERT( M_density != 0.0, "density coefficient has not been set in LinearElasticity");
	ASSERT( M_young != 0.0,   "young coefficient has not been set in LinearElasticity");
	ASSERT( M_poisson != 0.0, "poisson coefficient has not been set in LinearElasticity");

	using namespace ExpressionAssembly;

	// mass matrix
	integrate ( elements (M_displacementFESpace_ETA->mesh() ),
				M_displacementFESpace->qr(),
				M_displacementFESpace_ETA,
				M_displacementFESpace_ETA,
				value(M_density) * dot ( phi_i, phi_j )
			  ) >> M_mass_no_bc;

	if ( M_thinLayer )
	{
		QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
		M_mass_no_bc_thin->zero();

		integrate ( boundary ( M_displacementFESpace_ETA->mesh(), M_interfaceFlag ),
				    myBDQR,
				    M_displacementFESpace_ETA,
				    M_displacementFESpace_ETA,
				    value ( M_thinLayerDensity ) *  dot ( phi_i , phi_j )
				  ) >> M_mass_no_bc_thin;

		M_mass_no_bc_thin->globalAssemble();
		*M_mass_no_bc += *M_mass_no_bc_thin;
	}

	M_mass_no_bc->globalAssemble();

	// stiffness matrix
	integrate ( elements (M_displacementFESpace_ETA->mesh() ),
				M_displacementFESpace->qr(),
				M_displacementFESpace_ETA,
				M_displacementFESpace_ETA,
				value( M_mu ) * dot ( grad(phi_i) , grad(phi_j) + transpose( grad(phi_j) ) ) +
				value( M_lambda ) * div(phi_i) * div(phi_j)
			  ) >> M_stiffness_no_bc;

	if ( M_thinLayer )
	{
		QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
		M_stiffness_no_bc_thin->zero();

		MatrixSmall<3, 3> Eye;
		Eye[0][0] = 1.0;
		Eye[1][1] = 1.0;
		Eye[2][2] = 1.0;

		using namespace ExpressionAssembly;

		integrate ( boundary ( M_displacementFESpace_ETA->mesh(), M_interfaceFlag ),
					myBDQR,
					M_displacementFESpace_ETA,
					M_displacementFESpace_ETA,
					value( M_thinLayerThickness * M_thinLayerLameII ) *
					dot (
						( grad (phi_j) + (-1) * grad (phi_j) * outerProduct ( Nface, Nface ) ) +
						transpose (grad (phi_j) + (-1) * grad (phi_j) * outerProduct ( Nface, Nface ) ),
						( grad (phi_i) + ( (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) ) )
					)
					+ value( M_thinLayerThickness * M_thinLayerLameI ) *
					dot ( value ( Eye ) , ( grad (phi_j) + (-1) * grad (phi_j) * outerProduct ( Nface, Nface ) ) ) *
					dot ( value ( Eye ) , ( grad (phi_i) + (-1) * grad (phi_i) * outerProduct ( Nface, Nface ) ) )
				) >> M_stiffness_no_bc_thin;

		M_stiffness_no_bc_thin->globalAssemble();
		*M_stiffness_no_bc += *M_stiffness_no_bc_thin;
	}

	M_stiffness_no_bc->globalAssemble();

	// jacobian matrix

	*M_jacobian += *M_mass_no_bc;

	if ( useBDF )
	{
		*M_jacobian *= ( 1.0/( timestep*timestep ) * coeff );
	}
	else
	{
		*M_jacobian *= ( 1.0/( timestep*timestep*coeff ) );
	}

	*M_jacobian += *M_stiffness_no_bc;

	// Apply BC to the jacobian

	bc->bcUpdate ( *M_displacementFESpace->mesh(), M_displacementFESpace->feBd(), M_displacementFESpace->dof() );

	bcManageMatrix( *M_jacobian, *M_displacementFESpace->mesh(), M_displacementFESpace->dof(), *bc, M_displacementFESpace->feBd(), 1.0, 0.0);

	M_jacobian->globalAssemble();

	if ( M_comm->MyPID() == 0 )
		std::cout << "\nAssembly of structure done\n\n";
}


}
