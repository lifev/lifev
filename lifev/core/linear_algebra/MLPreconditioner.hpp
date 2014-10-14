/*
 * MLPreconditioner.hpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */

/*!
 * \file MLPreconditoner.hpp
 * \author Umberto Villa
 * \date 2011-10-08
 * Interface to the ML preconditioners in Trilinos. This interface requires the user to provide
 * a matrix in Epetra_CrsFormat and the parameters of the factorization.
 */

#ifndef MLPRECONDITIONER_HPP_
#define MLPRECONDITIONER_HPP_

#include <lifev/core/linear_algebra/RowMatrixPreconditioner.hpp>

namespace LifeV
{

namespace Operators
{
//! @class
/*!
 * @brief Interface to the ML preconditioners in Trilinos.
 *
 * This class inherit from \c RowMatrixPreconditioner and it implements its abstract protected methods.
 * An \c MLPreconditioner object should be allocated by using the \c RowMatrixPreconditionerFactory factory.
 */
class MLPreconditioner : public RowMatrixPreconditioner
{
public:
	MLPreconditioner();
	virtual ~MLPreconditioner();

protected:
	int myCompute();

};

//! @fn
/* !
 * @brief helper function for the \c RowMatrixPreconditionerFactory factory.
 */
inline RowMatrixPreconditioner* createML() { return new MLPreconditioner(); }
namespace
{
	static bool registerML = RowMatrixPreconditionerFactory::instance().registerProduct( "ML", &createML );
}


}

}

#endif /* MLPRECONDITIONER_HPP_ */
