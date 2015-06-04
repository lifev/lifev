#include <lifev/fsi_blocks/solver/FSIApplyOperatorNonConforming.hpp>

#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/util/LifeChrono.hpp>


namespace LifeV
{
namespace Operators
{
//===========================================================================//
// Constructors
//===========================================================================//

FSIApplyOperatorNonConforming::FSIApplyOperatorNonConforming():
    M_label("FSIApplyOperatorNonConforming"),
    M_useTranspose(false)
{

}

FSIApplyOperatorNonConforming::~FSIApplyOperatorNonConforming()
{

}

// show information about the class
void FSIApplyOperatorNonConforming::showMe()
{

}

void
FSIApplyOperatorNonConforming::setMonolithicMap(const mapEpetraPtr_Type& monolithicMap)
{
	M_monolithicMap = monolithicMap;
}

void
FSIApplyOperatorNonConforming::setMaps( const mapEpetraPtr_Type& fluid_velocity_map,
										const mapEpetraPtr_Type& fluid_pressure_map,
										const mapEpetraPtr_Type& structure_displacement_map,
										const mapEpetraPtr_Type& lagrange_multipliers_map,
										const mapEpetraPtr_Type& ALE_map,
										const mapEpetraPtr_Type& structure_interface_map)
{
	// Setting maps
	M_u_map = fluid_velocity_map;
	M_p_map = fluid_pressure_map;
	M_ds_map = structure_displacement_map;
	M_lambda_map = lagrange_multipliers_map;
	M_ale_map = ALE_map;
	M_structure_interface_map = structure_interface_map;

	// Setting offsets
	M_fluidVelocity = M_u_map->mapSize();
	M_fluid = M_fluidVelocity + M_p_map->mapSize();
	M_structure = M_fluid + M_ds_map->mapSize();
	M_lambda = M_structure + M_lambda_map->mapSize();

	// Initialize Vectors
	M_X_velocity.reset( new VectorEpetra_Type (*M_u_map, Unique) );
	M_X_pressure.reset( new VectorEpetra_Type (*M_p_map, Unique) );
	M_X_displacement.reset( new VectorEpetra_Type (*M_ds_map, Unique) );
	M_X_lambda.reset( new VectorEpetra_Type (*M_lambda_map, Unique) );
	M_X_geometry.reset( new VectorEpetra_Type (*M_ale_map, Unique) );

	M_Y_velocity.reset( new VectorEpetra_Type (*M_u_map, Unique) );
	M_Y_pressure.reset( new VectorEpetra_Type (*M_p_map, Unique) );
	M_Y_displacement.reset( new VectorEpetra_Type (*M_ds_map, Unique) );
	M_Y_lambda.reset( new VectorEpetra_Type (*M_lambda_map, Unique) );
	M_Y_geometry.reset( new VectorEpetra_Type (*M_ale_map, Unique) );
}

void
FSIApplyOperatorNonConforming::setShapeDerivativesBlocks( const matrixEpetraPtr_Type & shapeVelocity,
   											   	   	   	  const matrixEpetraPtr_Type & shapePressure)
{
	M_shapeVelocity = shapeVelocity;
	M_shapePressure = shapePressure;
}

void
FSIApplyOperatorNonConforming::setFluidBlocks ( const matrixEpetraPtr_Type & block00,
												const matrixEpetraPtr_Type & block01,
												const matrixEpetraPtr_Type & block10)
{
	M_F_00 = block00;
	M_F_01 = block01;
	M_F_10 = block10;
	M_useStabilization = false;
}

void
FSIApplyOperatorNonConforming::setFluidBlocks ( const matrixEpetraPtr_Type & block00,
												const matrixEpetraPtr_Type & block01,
												const matrixEpetraPtr_Type & block10,
												const matrixEpetraPtr_Type & block11)
{
	M_F_00 = block00;
	M_F_01 = block01;
	M_F_10 = block10;
	M_F_11 = block11;
	M_useStabilization = true;
}

void
FSIApplyOperatorNonConforming::setInterpolants ( interpolationPtr_Type fluidToStructure,
        										 interpolationPtr_Type structureToFluid)
{
	M_FluidToStructureInterpolant = fluidToStructure;
	M_StructureToFluidInterpolant = structureToFluid;
}

void
FSIApplyOperatorNonConforming::applyInverseInterfaceFluidMass(const VectorEpetraPtr_Type&X, VectorEpetraPtr_Type&Y) const
{
	Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
	belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList_rbf3d.xml" );

	// Preconditioner
	prec_Type* precRawPtr;
	basePrecPtr_Type precPtr;
	precRawPtr = new prec_Type;
	precRawPtr->setDataFromGetPot ( M_datafile, "prec" );
	precPtr.reset ( precRawPtr );

	// Linear Solver
	LinearSolver solver;
	solver.setCommunicator ( M_comm );
	solver.setParameters ( *belosList );
	solver.setPreconditioner ( precPtr );

	// Solution
	Y->zero();

	// Solve system
	solver.setOperator ( M_fluid_interface_mass );
	solver.setRightHandSide ( X );
	solver.solve ( Y );
}

void
FSIApplyOperatorNonConforming::setInterfaceMassMatrices (  const matrixEpetraPtr_Type &  fluid_interface_mass,
    								 	 	 	 	 	   const matrixEpetraPtr_Type &  structure_interface_mass)
{
	M_fluid_interface_mass     = fluid_interface_mass;
	M_structure_interface_mass = structure_interface_mass;
}

inline int
FSIApplyOperatorNonConforming::Apply(const vector_Type & X, vector_Type & Y) const
{
	ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

//	Input vector

	const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);

	//! Extract each component of the input vector
	M_X_velocity->subset(X_vectorEpetra, *M_u_map, 0, 0);

	M_X_pressure->subset(X_vectorEpetra, *M_p_map, M_fluidVelocity, 0 );

	M_X_displacement->subset(X_vectorEpetra, *M_ds_map, M_fluid, 0 );

	M_X_lambda->subset(X_vectorEpetra, *M_lambda_map, M_structure, 0 );

	M_X_geometry->subset(X_vectorEpetra, *M_ale_map, M_lambda , 0 );

	// M_X_lambda is defined just at the fluid interface, but for coding reasons
	// I need to have it on the whole fluid domain, i.e., vector of zeros with
	// nonzero entries only at the interface
	VectorEpetra_Type lambda_omega_f(*M_u_map);
	lambda_omega_f.subset(*M_X_lambda, *M_lambda_map, 0, 0);

	//----------------------------------------------------------//
	// Applying the Jacobian to the fluid velocity and pressure //
	//----------------------------------------------------------//

	M_Y_velocity->zero();
	M_Y_pressure->zero();

	if ( M_useStabilization )
	{
		if ( M_useShapeDerivatives )
		{
			*M_Y_velocity  = (*M_F_00) * (*M_X_velocity ) + (*M_F_01) * (*M_X_pressure) + lambda_omega_f + (*M_shapeVelocity) * (*M_X_geometry);
			*M_Y_pressure  = (*M_F_10) * (*M_X_velocity ) + (*M_F_11) * (*M_X_pressure) + (*M_shapePressure) * (*M_X_geometry);
		}
		else
		{
			*M_Y_velocity  = (*M_F_00) * (*M_X_velocity ) + (*M_F_01) * (*M_X_pressure) + lambda_omega_f;
			*M_Y_pressure  = (*M_F_10) * (*M_X_velocity ) + (*M_F_11) * (*M_X_pressure);
		}
	}
	else
	{
		if ( M_useShapeDerivatives )
		{
			*M_Y_velocity  = (*M_F_00) * (*M_X_velocity ) + (*M_F_01) * (*M_X_pressure) + lambda_omega_f + (*M_shapeVelocity) * (*M_X_geometry);
			*M_Y_pressure  = (*M_F_10) * (*M_X_velocity ) + (*M_shapePressure) * (*M_X_geometry);
		}
		else
		{
			*M_Y_velocity  = (*M_F_00) * (*M_X_velocity ) + (*M_F_01) * (*M_X_pressure) + lambda_omega_f;
			*M_Y_pressure  = (*M_F_10) * (*M_X_velocity );
		}
	}

	//-----------------------------------------------------//
	// Applying the Jacobian to the structure displacement //
	//-----------------------------------------------------//

	M_Y_displacement->zero();
	*M_Y_displacement = *M_S * ( *M_X_displacement );

	// From weak to strong form of the stresses lambda,
	// need to apply the inverse of the fluid interface
	// mass
	VectorEpetraPtr_Type X_lambda_strong ( new VectorEpetra_Type ( *M_lambda_map ) );
	X_lambda_strong->zero();
	applyInverseInterfaceFluidMass(M_X_lambda, X_lambda_strong);

	// Interpolate from fluid to structure interface,
	// need to have vector defined on omega fluid
	VectorEpetraPtr_Type X_lambda_strong_omega_f ( new VectorEpetra_Type ( *M_u_map ) );
	X_lambda_strong_omega_f->subset(*X_lambda_strong, *M_lambda_map, 0, 0);
	M_FluidToStructureInterpolant->updateRhs(X_lambda_strong_omega_f);
	M_FluidToStructureInterpolant->interpolate();
	VectorEpetraPtr_Type X_lambda_strong_omega_s ( new VectorEpetra_Type ( *M_ds_map ) );
	M_FluidToStructureInterpolant->solution(X_lambda_strong_omega_s);

	// Restrict X_lambda_strong_omega_s to Gamma S
	VectorEpetraPtr_Type X_lambda_strong_gamma_s ( new VectorEpetra_Type ( *M_structure_interface_map ) );
	X_lambda_strong_gamma_s->zero();
	X_lambda_strong_gamma_s->subset(*X_lambda_strong_omega_s, *M_structure_interface_map, 0, 0);

	// From strong to weak residual on Gamma S
	VectorEpetraPtr_Type X_lambda_weak_gamma_s ( new VectorEpetra_Type ( *M_structure_interface_map ) );
	X_lambda_weak_gamma_s->zero();
	*X_lambda_weak_gamma_s = *M_structure_interface_mass * (*X_lambda_strong_gamma_s);

	// Expand vector to whole omega S
	VectorEpetraPtr_Type X_lambda_weak_omega_s ( new VectorEpetra_Type ( *M_ds_map ) );
	X_lambda_weak_omega_s->subset(*X_lambda_weak_gamma_s, *M_structure_interface_map, 0, 0);

	*M_Y_displacement -= *X_lambda_weak_omega_s;

	//---------------------------------------//
	// Applying the Jacobian to the stresses //
	//---------------------------------------//

	M_Y_lambda->zero();
	M_Y_lambda->subset(*M_X_velocity, *M_lambda_map, 0, 0);

	VectorEpetraPtr_Type structure_vel( new VectorEpetra_Type ( *M_ds_map ) );
	structure_vel->zero();
	*structure_vel += *M_X_displacement;
	*structure_vel /= M_timeStep;
	M_StructureToFluidInterpolant->updateRhs(structure_vel);
	M_StructureToFluidInterpolant->interpolate();

	VectorEpetraPtr_Type structure_vel_omega_f ( new VectorEpetra_Type ( *M_u_map ) );
	structure_vel_omega_f->zero();
	M_StructureToFluidInterpolant->solution(structure_vel_omega_f);
	VectorEpetraPtr_Type structure_vel_gamma_f ( new VectorEpetra_Type ( *M_lambda_map ) );
	structure_vel_gamma_f->zero();
	structure_vel_gamma_f->subset(*structure_vel_omega_f, *M_lambda_map, 0, 0);

	*M_Y_lambda -= *structure_vel_gamma_f;

	//----------------------------------//
	// Applying the Jacobian to the ALE //
	//----------------------------------//

	M_Y_geometry->zero();
	*M_Y_geometry = ( *M_G ) * ( *M_X_geometry );

	M_StructureToFluidInterpolant->updateRhs(M_X_displacement);
	M_StructureToFluidInterpolant->interpolate();

	VectorEpetraPtr_Type tmp ( new VectorEpetra_Type ( *M_ale_map ) );
	M_StructureToFluidInterpolant->solution(tmp);

	*M_Y_geometry -= *tmp;

//	Output vector

	VectorEpetra_Type Y_vectorEpetra(Y, M_monolithicMap, Unique);

	//! Copy the single contributions into the optput vector
	Y_vectorEpetra.subset(*M_Y_velocity, M_Y_velocity->map(), 0, 0 );
	Y_vectorEpetra.subset(*M_Y_pressure, M_Y_pressure->map(), 0, M_fluidVelocity );
	Y_vectorEpetra.subset(*M_Y_displacement, M_Y_displacement->map(), 0, M_fluid );
	Y_vectorEpetra.subset(*M_Y_lambda, M_Y_lambda->map(), 0, M_structure );
	Y_vectorEpetra.subset(*M_Y_geometry, M_Y_geometry->map(), 0, M_lambda );

	Y = dynamic_cast<Epetra_MultiVector &>( Y_vectorEpetra.epetraVector() );

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
