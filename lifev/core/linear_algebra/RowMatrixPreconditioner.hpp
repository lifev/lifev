/*
 * RowMatrixPreconditioner.hpp
 *
 *  Created on: Oct 8, 2011
 *      Author: uvilla
 */
/*!
 * \file RowMatrixPreconditoner.hpp
 * \author Umberto Villa
 * \date 2011-10-08
 * Abstract class to construct preconditioners from a matrix in Epetra_CsrFormat.
 */

#ifndef ROWMATRIXPRECONDITIONER_HPP_
#define ROWMATRIXPRECONDITIONER_HPP_

#include <lifev/operator/linear_algebra/LinearOperator.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/util/Factory.hpp>

#include <lifev/core/LifeV.hpp>

// Tell the compiler to ignore specific kind of warnings
LIFEV_SUPPRESS_WARNINGS

#include <Epetra_CrsMatrix.h>
#include <Teuchos_ParameterList.hpp>

// Tell the compiler to restore the warnings
LIFEV_RESTORE_WARNINGS

namespace LifeV
{

namespace Operators
{
//! @class
/*!
 * @brief Abstract class to construct preconditioners from a matrix in Epetra_CsrFormat.
 *
 * This class implementents all public methods of the parent class \c LinearOperator, and additionally introduce
 * the following public methods"
 * <ul>
 * <li> \c SetRowMatrix to specify the matrix in Epetra_CsrMatrix format
 * <li> \c SetParameterList to specify the list of parameter to be used in the preconditioner
 * <li> \c Compute do all the operator to compute the preconditioner. This method assert that the parameter are
 *         set correctly and then it calls the protected abstract method myCompute.
 * </ul>
 *
 * Concrete instances of the \c RowMatrixPreconditioner class should implement the protected method \c myCompute.
 */

class RowMatrixPreconditioner : public LinearOperator
{
public:

    //@name Typdefs
    //@{
    typedef Epetra_CrsMatrix rowMatrix_Type;
    typedef boost::shared_ptr<rowMatrix_Type> rowMatrixPtr_Type;
    typedef Teuchos::ParameterList pList_Type;
    //@}

    //! Empty constructor
    RowMatrixPreconditioner(): LinearOperator() {};


    //! Destructor
    virtual ~RowMatrixPreconditioner() {};

    //! @name Attribute set methods
    //@{
    //! Set the row matrix
    void SetRowMatrix(const rowMatrixPtr_Type & rowMatrix)
    {
        ASSERT_PRE(rowMatrix.get() != 0, "RowMatrixPreconditioner::SetRowMatrix(): rowMatrix can not be a null Pointer");
        M_rowMatrix = rowMatrix;
    }

    //! Set the list of paramenters
    void SetParameterList(const pList_Type pList)
    {
        M_pList = pList;
    }

    //! Transposition
    virtual int SetUseTranspose(bool UseTranspose)
    {
        return M_prec->SetUseTranspose(UseTranspose);
    }
    //@}

    //! Compute the preconditioner.
    /*!
     * Derived classes should implement the protected method myCompute called by this function.
     * @return: 0 success, negative number error, positive number warning
     */
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
    //! Abstract method myCompute implemented by the derived class
    virtual int myCompute() = 0;

    rowMatrixPtr_Type M_rowMatrix;
    operatorPtr_Type  M_prec;
    pList_Type M_pList;


};

//! @typedef
/*!
 * @brief \c RowMatrixPreconditionerFactory should be used for the creation of concrete instances of \c RowMatrixPreconditioner
 */
typedef FactorySingleton<Factory<RowMatrixPreconditioner, std::string> > RowMatrixPreconditionerFactory;

} /* end Operator namespace */
} /* end LifeV namespace */

#endif /* ROWMATRIXPRECONDITIONER_HPP_ */
