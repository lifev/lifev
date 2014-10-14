/*
 * TwoLevelOperator.hpp
 *
 *  Created on: Oct 9, 2011
 *      Author: uvilla
 */
//@HEADER

/*!
 * \file TwoLevelOperator.hpp
 * \author Umberto Villa
 * \date 2011-10-09
 * This file contains the definition of the class \c TwoLevelOperator.
 * \c TwoLevelOperator defines a two level method, and requires the fine level operator and smoother,
 * the restriction and prolongation operators, and a coarse solver operator.
 */

#ifndef TWOLEVELOPERATOR_HPP_
#define TWOLEVELOPERATOR_HPP_

#include <lifev/core/linear_algebra/LinearOperator.hpp>

namespace LifeV
{

namespace Operators
{
//! @class
/*!
 * @brief It defines a two level methods to approximately apply the inverse of a fine level operator.
 *
 * It implements the public interface of \c LinearOperator and it defines the public methods to set information
 * the fine level operator and smoother, the coarse level solver, the extension and restriction operators.
 */
class TwoLevelOperator : public LinearOperator
{
public:
	TwoLevelOperator();
	virtual ~TwoLevelOperator();

	//! @name Setters
	//@{
	//! Set the fine level operator
	void SetFineLevelOperator(const operatorPtr_Type & fineLevelOper);
	//! Set the fine level smoother (smootherOper should have a \c ApplyInverse method)
	void SetSmootherOperator(const operatorPtr_Type & smootherOper);
	//! Set the coarse level solver (coarseLevelOper should have a \c ApplyInverse method)
	void SetCoarseLevelOperator(const operatorPtr_Type & coarseLevelOper );
	//! Set the restriction operator from the fine to coarse level
	void SetRestrictionOperator(const operatorPtr_Type & restrictionOper );
	//! Set the extension operatoe from the coarse to fine level.
	void SetEstensionOperator(const operatorPtr_Type & estensionOper);
	//@}

	//! Check that the range and domains of all operators are compatible
	int checkConsistency();

	//! \warning not fully supported!
    virtual int SetUseTranspose(bool UseTranspose);
    //@}

    //! @name Mathematical functions
    //@{
    virtual int Apply(const vector_Type& X, vector_Type& Y) const;

    virtual int ApplyInverse(const vector_Type& X, vector_Type& Y) const;

    virtual double NormInf() const;
    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char * Label() const;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type & Comm() const;

    //! Returns the raw_map object associated with the domain of this operator.
    virtual const map_Type & OperatorDomainMap() const;

    //! Returns the raw_map object associated with the range of this operator.
    virtual const map_Type & OperatorRangeMap() const;
    //@}

private:
	operatorPtr_Type M_fineLevelOper;
	operatorPtr_Type M_smootherOper;
	operatorPtr_Type M_coarseLevelOper;
	operatorPtr_Type M_restrictionOper;
	operatorPtr_Type M_estensionOper;

};

} /* end namespace Operators */

} /* end namespace LifeV */

#endif /* TWOLEVELOPERATOR_HPP_ */
