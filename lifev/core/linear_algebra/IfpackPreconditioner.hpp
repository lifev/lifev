/*
 * IfpackPreconditioner.hpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */
//@HEADER

/*!
 * \file IfpackPreconditioner.hpp
 * \author Umberto Villa
 * \date 2011-10-08
 * Interface to the Ifpack preconditioners in Trilinos. This interface requires the user to provide
 * a matrix in Epetra_CrsFormat and the parameters of the factorization.
 */

#ifndef IFPACKPRECONDITIONER_HPP_
#define IFPACKPRECONDITIONER_HPP_

#include <lifev/core/linear_algebra/RowMatrixPreconditioner.hpp>

namespace LifeV
{

namespace Operators
{
//! @class
/*!
 * @brief Interface to the Ifpack preconditioners in Trilinos.
 *
 * This class inherit from \c RowMatrixPreconditioner and it implements its abstract protected methods.
 * An \c IfpackPreconditioner object should be allocated by using the \c RowMatrixPreconditionerFactory factory.
 */
class IfpackPreconditioner : public RowMatrixPreconditioner
{

public:
	IfpackPreconditioner();
	virtual ~IfpackPreconditioner();

protected:
	int myCompute();

};

//! @fn
/* !
 * @brief helper function for the \c RowMatrixPreconditionerFactory factory.
 */
inline RowMatrixPreconditioner* createIfpack() { return new IfpackPreconditioner(); }
namespace
{
	static bool registerIF = RowMatrixPreconditionerFactory::instance().registerProduct( "Ifpack", &createIfpack );
}

} /* end Operators namespace */

} /* end LifeV namespace */



#endif /* IFPACKPRECONDITIONER_HPP_ */
