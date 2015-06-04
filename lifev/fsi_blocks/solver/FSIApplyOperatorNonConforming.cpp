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
										const mapEpetraPtr_Type& ALE_map)
{
	// Setting maps
	M_u_map = fluid_velocity_map;
	M_p_map = fluid_pressure_map;
	M_ds_map = structure_displacement_map;
	M_lambda_map = lagrange_multipliers_map;
	M_ale_map = ALE_map;

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

//	Output vector

	VectorEpetra_Type Y_vectorEpetra(Y, M_monolithicMap, Unique);

	//! Copy the single contributions into the optput vector
	Y_vectorEpetra.subset(*M_Y_velocity, M_Y_velocity->map(), 0, 0 );
	Y_vectorEpetra.subset(*M_Y_pressure, M_Y_pressure->map(), 0, M_fluidVelocity );
	Y_vectorEpetra.subset(*M_Y_displacement, M_Y_displacement->map(), 0, M_fluid );
	Y_vectorEpetra.subset(*M_Y_lambda, M_Y_lambda->map(), 0, M_fluid + M_structure );
	Y_vectorEpetra.subset(*M_Y_geometry, M_Y_geometry->map(), 0, M_fluid + M_structure + M_lambda );

	Y = dynamic_cast<Epetra_MultiVector &>( Y_vectorEpetra.epetraVector() );

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
