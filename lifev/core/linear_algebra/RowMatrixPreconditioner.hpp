/*
 * RowMatrixPreconditioner.hpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */

#ifndef ROWMATRIXPRECONDITIONER_HPP_
#define ROWMATRIXPRECONDITIONER_HPP_

#include <lifev/operator/linear_algebra/LinearOperator.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/util/Factory.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <Epetra_CrsMatrix.h>
#include <Teuchos_ParameterList.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

namespace LifeV
{

namespace Operators
{

class RowMatrixPreconditioner : public LinearOperator
{
public:

	typedef Epetra_CrsMatrix rowMatrix_Type;
	typedef boost::shared_ptr<rowMatrix_Type> rowMatrixPtr_Type;
	typedef Teuchos::ParameterList pList_Type;

	RowMatrixPreconditioner(): LinearOperator() {};

    //! @name Destructor
    //@{
    //! Destructor
    virtual ~RowMatrixPreconditioner() {};
    //@}

    //! @name Attribute set methods
    //@{
    void SetRowMatrix(const rowMatrixPtr_Type & rowMatrix)
    {
    	ASSERT_PRE(rowMatrix.get() != 0, "RowMatrixPreconditioner::SetRowMatrix(): rowMatrix can not be a null Pointer");
    	M_rowMatrix = rowMatrix;
    }

    void SetParameterList(const pList_Type pList)
    {
    	M_pList = pList;
    }

    virtual int SetUseTranspose(bool UseTranspose)
    {
    	return M_prec->SetUseTranspose(UseTranspose);
    }
    //@}

    int Compute()
    {
    	ASSERT_PRE(M_rowMatrix.get()!= 0, "RowMatrixPreconditioner::Compute(): You need to SetRowMatrix first \n");
    	return myCompute();
    }

    //! @name Mathematical functions
    //@{
    virtual int Apply(const vector_Type& X, vector_Type& Y) const
    {
    	return M_prec->Apply(X,Y);
    }

    virtual int ApplyInverse(const vector_Type& X, vector_Type& Y) const
    {
    	return M_prec->ApplyInverse(X,Y);
    }

    virtual double NormInf() const
    {
    	return M_prec->NormInf();
    }
    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char * Label() const
    {
    	return M_prec->Label();
    }

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const
    {
    	return M_prec->UseTranspose();
    }

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const
    {
    	return M_prec->HasNormInf();
    }

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type & Comm() const
    {
    	return M_prec->Comm();
    }

    //! Returns the raw_map object associated with the domain of this operator.
    virtual const map_Type & OperatorDomainMap() const
    {
    	return M_prec->OperatorDomainMap();
    }

    //! Returns the raw_map object associated with the range of this operator.
    virtual const map_Type & OperatorRangeMap() const
    {
    	return M_prec->OperatorRangeMap();
    }
    //@}

protected:
    virtual int myCompute() = 0;

    rowMatrixPtr_Type M_rowMatrix;
    operatorPtr_Type  M_prec;
    pList_Type M_pList;


};

typedef FactorySingleton<Factory<RowMatrixPreconditioner, std::string> > RowMatrixPreconditionerFactory;

} /* end Operator namespace */
} /* end LifeV namespace */

#endif /* ROWMATRIXPRECONDITIONER_HPP_ */
