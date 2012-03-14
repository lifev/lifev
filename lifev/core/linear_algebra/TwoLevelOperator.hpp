/*
 * TwoLevelOperator.hpp
 *
 *  Created on: Oct 9, 2011
 *      Author: uvilla
 */

#ifndef TWOLEVELOPERATOR_HPP_
#define TWOLEVELOPERATOR_HPP_

#include <lifev/operator/linear_algebra/LinearOperator.hpp>

namespace LifeV
{

namespace Operators
{

class TwoLevelOperator : public LinearOperator
{
public:
	TwoLevelOperator();
	virtual ~TwoLevelOperator();

	void SetFineLevelOperator(const operatorPtr_Type & fineLevelOper);
	void SetSmootherOperator(const operatorPtr_Type & smootherOper);
	void SetCoarseLevelOperator(const operatorPtr_Type & coarseLevelOper );
	void SetRestrictionOperator(const operatorPtr_Type & restrictionOper );
	void SetEstensionOperator(const operatorPtr_Type & estensionOper);

	int checkConsistency();

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
