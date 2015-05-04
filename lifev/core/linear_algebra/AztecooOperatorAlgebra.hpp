/*
 * AztecooOperator.hpp
 *
 *  Created on: Sep 3, 2010
 *      Author: uvilla
 */

//@HEADER

/*!
 * \file AztecooOperator.hpp
 * \author Umberto Villa
 * \date 2010-09-03
 * Interface to the AztecOO linear solver in Trilinos. This interface requires the user to provide
 * the preconditioner in the form of a LinearOperator. If no preconditioner is provided,
 * AztecOO will use unpreconditioned Krylov methods.
 */

#ifndef AZTECOO_OPERATOR_ALGEBRA_HPP
#define AZTECOO_OPERATOR_ALGEBRA_HPP

#include <AztecOO.h>
#include <Teuchos_ParameterList.hpp>

#include <lifev/core/linear_algebra/InvertibleOperator.hpp>

namespace LifeV
{
namespace Operators
{
//! @class
/*! @brief InvertibleOperator interface to AztecOO in Trilinos.
 * AztecooOperator will use the matrix-free Krylov methods in Trilinos, therefore both the operator and
 * the preconditioner must be given in the form of LinearOperator.
 *
 * For a description of the class functionality, please refer to the parent class InvertibleOperator.
 *
 */
class AztecooOperatorAlgebra : public InvertibleOperator
{
public:
    typedef AztecOO SolverType;
    typedef boost::shared_ptr<SolverType> SolverType_ptr;

    AztecooOperatorAlgebra();

    int numberOfIterations()
    {
        return M_linSolver->NumIters();
    }

protected:

    virtual int doApplyInverse(const vector_Type& X, vector_Type& Y) const;
    virtual void doSetOperator(){};
    virtual void doSetPreconditioner(){};
    virtual void doSetParameterList(){};

    SolverType_ptr                              M_linSolver;
};

inline InvertibleOperator* createAztecooOperatorAlgebra() { return new AztecooOperatorAlgebra(); }
namespace
{
    static bool registerAztecoo = InvertibleOperatorFactory::instance().registerProduct( "AztecOO", &createAztecooOperatorAlgebra );
}


} // Namespace Operators

} // Namespace LifeV

#endif // AZTECOO_OPERATOR_HPP
