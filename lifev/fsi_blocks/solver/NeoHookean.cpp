#include <lifev/fsi_blocks/solver/NeoHookean.hpp>

namespace LifeV
{

NeoHookean::NeoHookean( const commPtr_Type& communicator ):
		M_comm(communicator),
		M_bulk(0.0),
		M_mu(0.0),
		M_offset ( 0 )
{
	M_identity (0, 0) = 1.0;
	M_identity (0, 1) = 0.0;
	M_identity (0, 2) = 0.0;
	M_identity (1, 0) = 0.0;
	M_identity (1, 1) = 1.0;
	M_identity (1, 2) = 0.0;
	M_identity (2, 0) = 0.0;
	M_identity (2, 1) = 0.0;
	M_identity (2, 2) = 1.0;
	if ( M_comm->MyPID() == 0 )
	{
		std::cout << "\nUsing NeoHookean model for the structure\n";
	}
}

NeoHookean::~NeoHookean()
{}

void
NeoHookean::setCoefficients ( const Real density, const Real young, const Real poisson)
{
	// Copying the parameters

	M_density = density;

	// Evaluation of the lame coefficients

	Real lambda = ( young * poisson ) / ( ( 1.0 + poisson ) * ( 1.0 - 2.0 * poisson ) );

	M_mu = young / ( 2.0 * ( 1.0 + poisson ) );

	M_bulk = ( 2.0 / 3.0 ) * M_mu + lambda;


}

void
NeoHookean::setup( const meshPtr_Type& mesh, const std::string dOrder)
{
	// Create FE spaces

	M_dOrder = dOrder;

	M_displacementFESpace.reset (new FESpace<mesh_Type, map_Type> (mesh, M_dOrder, 3, M_comm) );

	M_displacementFESpace_ETA.reset( new ETFESpace_displacement ( M_displacementFESpace->mesh(), &(M_displacementFESpace->refFE()), M_comm));
}

void
NeoHookean::evaluate_residual( const vectorPtr_Type& solution, const Real& coefficient, const vectorPtr_Type& csi, vectorPtr_Type& residual )
{
	// coefficient is: 1.0/ ( M_dt * M_dt * M_structureTimeAdvance->get_beta() )

	if ( M_comm->MyPID() == 0 )
	{
		std::cout << " S- Assembling residual NeoHookean structure\n";
	}

	residual.reset (new vector_Type ( M_displacementFESpace->map() ) );
	residual->zero();

	vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

	vectorPtr_Type csi_rep ( new vector_Type ( *csi, Repeated ) );

	using namespace ExpressionAssembly;

	// Definition of F
	tensorF_Type F = ExpressionDefinitions::deformationGradient( M_displacementFESpace_ETA, *solution_rep, M_offset, M_identity );

	// Definition of J
	determinantF_Type J = ExpressionDefinitions::determinantF( F );

	// Definition of tensor C
	tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

	// Definition of F^-T
	minusT_Type  F_T = ExpressionDefinitions::minusT( F );

	// Definition of tr( C )
	traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

	integrate ( elements ( M_displacementFESpace_ETA->mesh() ),
			M_displacementFESpace->qr(),
			M_displacementFESpace_ETA,
			value ( coefficient * M_density ) * dot( value( M_displacementFESpace_ETA, *solution_rep ), phi_i )  -
			value ( M_density ) * dot ( value ( M_displacementFESpace_ETA, *csi_rep ), phi_i ) +
			value ( M_mu ) * pow (J, -2.0 / 3.0) * (dot ( F - value (1.0 / 3.0) * I_C * F_T, grad (phi_i) ) )
			+ value (1.0 / 2.0) * value ( M_bulk ) * ( pow ( J , 2.0) - J + log (J) ) * dot (  F_T, grad (phi_i) )
	) >> residual;

	residual->globalAssemble();
}

void
NeoHookean::update_jacobian(const vectorPtr_Type& solution, const Real& coefficient, matrixPtr_Type& jacobian )
{
	using namespace ExpressionAssembly;

	jacobian->zero();

	vectorPtr_Type solution_rep ( new vector_Type ( *solution, Repeated ) );

	// Definition of F
	tensorF_Type F = ExpressionDefinitions::deformationGradient( M_displacementFESpace_ETA, *solution_rep, M_offset, M_identity );

	// Definition of J
	determinantF_Type J = ExpressionDefinitions::determinantF( F );

	// Definition of tensor C
	tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

	// Definition of F^-T
	minusT_Type  F_T = ExpressionDefinitions::minusT( F );

	// Definition of tr( C )
	traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

	integrate ( elements (  M_displacementFESpace_ETA->mesh() ) ,
			M_displacementFESpace->qr(),
			M_displacementFESpace_ETA,
			M_displacementFESpace_ETA,
			/* Inertia term */
			value ( coefficient * M_density ) * dot( phi_i, phi_j )
			/*Isochoric Part*/
			/*Stiffness matrix : int { -2/3 * mu * J^(-2/3) *( F^-T : \nabla \delta ) ( F : \nabla \v ) } */
			+ value (-2.0 / 3.0) * value ( M_mu ) * pow (J, - (2.0 / 3.0) )  * dot ( F_T , grad (phi_j) ) * dot ( F , grad (phi_i) )
			/*Stiffness matrix : int { 2/9 * mu * ( Ic_iso )( F^-T : \nabla \delta ) ( F^-T : \nabla \v ) } */
			+ value (2.0 / 9.0) * value ( M_mu ) * pow( J, (-2.0/3.0) ) * I_C  * dot ( F_T , grad (phi_j) ) * dot ( F_T , grad (phi_i) )
			/* Stiffness matrix : int { mu * J^(-2/3) (\nabla \delta : \nabla \v)} */
			+ value ( M_mu ) * pow (J, - (2.0 / 3.0) )  * dot ( grad (phi_j), grad (phi_i) )
			/* Stiffness matrix : int { -2/3 * mu * J^(-2/3) ( F : \nabla \delta ) ( F^-T : \nabla \v ) } */
			+ value (-2.0 / 3.0) * value ( M_mu ) * pow (J, - (2.0 / 3.0) )  * dot ( F , grad (phi_j) ) * dot ( F_T , grad (phi_i) )
			/* Stiffness matrix : int { 1/3 * mu * Ic_iso * (F^-T [\nabla \delta]^t F^-T ) : \nabla \v } */
			+ value (1.0 / 3.0) * value ( M_mu ) * pow( J, -2.0/3.0) * I_C  * dot ( ( F_T * transpose (grad (phi_j) ) * F_T ), grad (phi_i) )
			/*Volumetric Part*/
			+ value ( 1.0 / 2.0 ) * value ( M_bulk ) * ( value (2.0) *pow (J, 2.0) - J + value (1.0) ) * dot ( F_T, grad (phi_j) ) * dot ( F_T, grad (phi_i) )
			+ value ( - 1.0 / 2.0 ) * value ( M_bulk ) * ( pow (J, 2.0) - J + log (J) ) * dot ( F_T * transpose (grad (phi_j) ) * F_T,  grad (phi_i) )
	) >> jacobian;

	jacobian->globalAssemble();

}

}
