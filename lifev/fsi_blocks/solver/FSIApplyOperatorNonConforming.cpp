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

inline int
FSIApplyOperatorNonConforming::Apply(const vector_Type & X, vector_Type & Y) const
{
	ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");

//	Input vector

//	const VectorEpetra_Type X_vectorEpetra(X, M_monolithicMap, Unique);
//
//	//! Extract each component of the input vector
//	M_X_velocity->subset(X_vectorEpetra, M_F->map(), 0, 0);
//
//	M_X_pressure->subset(X_vectorEpetra, M_B->map(), M_fluidVelocity, 0 );
//
//	M_X_displacement->subset(X_vectorEpetra, M_S->map(), M_fluid, 0 );
//
//	M_X_lambda->subset(X_vectorEpetra, M_C1->map(), M_fluid + M_structure, 0 );
//
//	M_X_geometry->subset(X_vectorEpetra, M_G->map(), M_fluid + M_structure + M_lambda , 0 );

//	Output vector

//	VectorEpetra_Type Y_vectorEpetra(Y, M_monolithicMap, Unique);
//
//	//! Copy the single contributions into the optput vector
//	Y_vectorEpetra.subset(*M_Y_velocity, M_Y_velocity->map(), 0, 0 );
//	Y_vectorEpetra.subset(*M_Y_pressure, M_Y_pressure->map(), 0, M_fluidVelocity );
//	Y_vectorEpetra.subset(*M_Y_displacement, M_Y_displacement->map(), 0, M_fluid );
//	Y_vectorEpetra.subset(*M_Y_lambda, M_Y_lambda->map(), 0, M_fluid + M_structure );
//	Y_vectorEpetra.subset(*M_Y_geometry, M_Y_geometry->map(), 0, M_fluid + M_structure + M_lambda );
//
//	Y = dynamic_cast<Epetra_MultiVector &>( Y_vectorEpetra.epetraVector() );

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
