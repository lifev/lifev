/*
 * LinearOperator.hpp
 *
 *  Created on: Sep 3, 2010
 *      Author: uvilla
 */

#include <lifev/operator/linear_algebra/LinearOperator.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
namespace Operators
{

int LinearOperator::apply(const VectorEpetra & X, VectorEpetra & Y) const
{
	return Apply(X.epetraVector(), Y.epetraVector());
}

int LinearOperator::applyInverse(const VectorEpetra & X, VectorEpetra & Y)
{
    return ApplyInverse(X.epetraVector(), Y.epetraVector());
}


} /*end namespace Operators*/
} /*end namespace */
