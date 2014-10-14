/*
 * TwoLevelPreconditioner.hpp
 *
 *  Created on: Oct 12, 2011
 *      Author: uvilla
 */

/*!
 * \file TwoLevelPreconditioner.hpp
 * \author Umberto Villa
 * \date 2011-10-12
 * This file contains the definition of the class \c TwoLevelPreconditioner.
 * \c TwoLevelPreconditioer require a fine level operator and extension/prolongator operator of class \c Epetra_CsrMatrix.
 * The \c ApplyInverse method of this class provides a two level approximation of the inverse of the original operator.
 */

#ifndef TWOLEVELPRECONDITIONER_HPP_
#define TWOLEVELPRECONDITIONER_HPP_

#include <lifev/core/linear_algebra/RowMatrixPreconditioner.hpp>

namespace LifeV
{

namespace Operators
{
//! @class
/*!
 * @brief An operator that provides a two level approximation of the inverse of a \c Epetra_CsrMatrix object.
 *
 * \c TwoLevelPreconditioner has exactly the same public interface of its base class \c RowMatrixPreconditioner.
 * Internally this class creates a \c TwoLevelOperator object.
 * We use the parameterList to inform a \c TwoLevelPreconditioner object of additional information as the Restriction and Exstension operator.
 * In particular the parameter list provided to an instance of this class should have:
 * <ul>
 * <li> a parameter \c EstensionMatrix which stores the extension operator (pointer to Epetra_CsrMatrix).
 * <li> a parameter \c RestrictionMatrix which stores the restriction operator (pointer to Epetra_CsrMatrix). If not provided the transpose of the extension operator will be used.
 * <li> a sublist \c CoarseLevel which contains all the information to set-up the coarse solver operator of type \c ApproximateInverseRowMatrix.
 * <li> a sublist \c FineLevel which contains all the information to set-up the fine level smoother of type \c RowMatrixPreconditioner.
 * </ul>
 *
 *
 */
class TwoLevelPreconditioner : public RowMatrixPreconditioner
{
public:
	TwoLevelPreconditioner();
	virtual ~TwoLevelPreconditioner();

protected:
    virtual int myCompute();

};

inline RowMatrixPreconditioner* createTwoLevel() { return new TwoLevelPreconditioner(); }
namespace
{
	static bool registerTL = RowMatrixPreconditionerFactory::instance().registerProduct( "TwoLevel", &createTwoLevel );
}


}

}

#endif /* TWOLEVELPRECONDITIONER_HPP_ */
