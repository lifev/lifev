/*
 * AztecooOperator.hpp
 *
 *  Created on: Sep 3, 2010
 *      Author: uvilla
 */

#ifndef AZTECOOOPERATOR_HPP_
#define AZTECOOOPERATOR_HPP_


#include <lifev/operator/linear_algebra/InvertibleOperator.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <AztecOO.h>
#include <Teuchos_ParameterList.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
namespace Operators
{
//! @class InvertibleOperator
/*! @brief Abstract class which defines the interface of an Invertible Linear Operator.
 *
 */
class AztecooOperator : public InvertibleOperator
{
public:
	typedef AztecOO SolverType;
	typedef boost::shared_ptr<SolverType> SolverType_ptr;

	AztecooOperator();

    int numberOfIterations()
    {
        return M_linSolver->NumIters();
    }

protected:

	virtual int doApplyInverse(const vector_Type& X, vector_Type& Y) const;
	virtual void doSetOperator(){};
	virtual void doSetPreconditioner(){};
	virtual void doSetParameterList(){};

	SolverType_ptr								M_linSolver;
};

inline InvertibleOperator* createAztecooOperator() { return new AztecooOperator(); }
namespace
{
	static bool registerAztecoo = InvertibleOperatorFactory::instance().registerProduct( "AztecOO", &createAztecooOperator );
}


} // Namespace Operators

} // Namespace LifeV

#endif // LINEAROPERATOR_HPP_
