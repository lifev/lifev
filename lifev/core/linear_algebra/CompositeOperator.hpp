/*
 * CompositeOperator.hpp
 *
 *  Created on: Jul 28, 2012
 *      Author: uvilla
 */

#ifndef COMPOSITEOPERATOR_HPP_
#define COMPOSITEOPERATOR_HPP_

#include <lifev/operator/linear_algebra/LinearOperator.hpp>

namespace LifeV
{

namespace Operators
{

class CompositeOperator : public LinearOperator
{
public:
    CompositeOperator();

    virtual ~CompositeOperator();

    //! \warning Transpose is not supported yet.
    virtual int SetUseTranspose(bool /*UseTranspose*/) { return -1; }

    //! The first operator we push is the first to be applied
    int pushBack( const operatorPtr_Type & op, const bool inverted);

    //! set whenever we want to define Apply or ApplyInverse
    void DefineAlreadyInverted(){isAlreadyInverted = true;}

    //! @name Mathematical functions
    //@{

    //! Applies all the operator (with their flags) starting from the first that was pushed
     virtual int Apply(const vector_Type& X, vector_Type& Y) const;

   //! Applies all the operator (with their flags) starting from the first that was pushed
    virtual int ApplyInverse(const vector_Type& X, vector_Type& Y) const;

    //! Returns the infinity norm of the global matrix.
    virtual double NormInf() const { return -1;}
    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char * Label() const{ return "Operators::CompositeOperator \n";}

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const { return false;}

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const {return false;}

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type & Comm() const { return (*ops.begin())->Comm(); }

    //! Returns the raw_map object associated with the domain of this operator.
    virtual const map_Type & OperatorDomainMap() const { return ops.front()->OperatorDomainMap(); }

    //! Returns the raw_map object associated with the range of this operator.
    virtual const map_Type & OperatorRangeMap() const { return ops.back()->OperatorRangeMap(); }
    //@}

private:

    int allocateTmpVects(const vector_Type& X, vector_Type& Y) const;
    void deleteTmpVects() const;

    std::vector<operatorPtr_Type> ops;
    std::vector<bool> inverted;
    mutable std::vector<Epetra_MultiVector *> vects;
    mutable int nAllocatedVectorsInMultVect;
    bool isAlreadyInverted;

};

} /* namespace Operators */
} /* namespace LifeV */
#endif /* COMPOSITEOPERATOR_HPP_ */
