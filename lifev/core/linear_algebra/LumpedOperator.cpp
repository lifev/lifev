/*
 * LumpedOperator.cpp
 *
 *  Created on: Sep 17, 2010
 *      Author: uvilla
 */

#include <lifev/core/LifeV.hpp>
#include <lifev/operator/linear_algebra/LumpedOperator.hpp>

// Tell the compiler to ignore specific kind of warnings
LIFEV_SUPPRESS_WARNINGS

#include <Epetra_Vector.h>

// Tell the compiler to restore the warnings
LIFEV_RESTORE_WARNINGS

namespace LifeV
{
namespace Operators
{
LumpedOperator::LumpedOperator():
        M_name("LumpedOperator"),
        M_normInf(-1),
        M_useTranspose(false),
        M_alpha(1.0)
{
    // Nothing to be done here
}

void LumpedOperator::setUp(const rowMatrixPtr_Type & _matrix)
{
    ASSERT_PRE( _matrix->OperatorRangeMap().PointSameAs( _matrix->OperatorDomainMap()),
            "The matrix we want to lump must be square");
    M_matrix = _matrix;
}

int LumpedOperator::compute()
{
    M_lumped.reset(new Epetra_Vector( M_matrix->OperatorRangeMap()) );
    M_matrix->ExtractDiagonalCopy(*M_lumped);
    lumpedMatrix_Type ones(M_matrix->OperatorDomainMap()), m1(M_matrix->OperatorRangeMap());
    ones.PutScalar(1.0); m1.PutScalar(0.0);
    EPETRA_CHK_ERR(M_matrix->Multiply(false, ones, m1));
    Real totalMass, lumpMass;
    ones.Dot(m1, &totalMass);
    ones.Dot(*M_lumped, &lumpMass);
    M_lumped->Scale(totalMass/lumpMass);
    M_lumped->NormInf(&M_normInf);

    return 0;
}

int LumpedOperator::Apply(const vector_Type& X, vector_Type& Y) const
        {
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");
    ASSERT_PRE(X.Map().PointSameAs(M_lumped->Map()), "X map should be the same of M_lumped");
    ASSERT_PRE(Y.Map().PointSameAs(M_lumped->Map()), "Y map should be the same of M_lumped");

    return Y.Multiply(M_alpha, X, *M_lumped, 0.0);
        }

int LumpedOperator::ApplyInverse(const vector_Type& X, vector_Type& Y) const
        {
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "X and Y must have the same number of vectors");
    ASSERT_PRE(X.Map().PointSameAs(M_lumped->Map()), "X map should be the same of M_lumped");
    ASSERT_PRE(Y.Map().PointSameAs(M_lumped->Map()), "Y map should be the same of M_lumped");

    return Y.ReciprocalMultiply(1.0/M_alpha, *M_lumped, X, 0.0);
        }

} /* end namespace Operators */
} /*end namespace*/
