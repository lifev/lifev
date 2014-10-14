/*
 * ApproximatedInvertibleRowMatrix.h
 *
 *  Created on: Oct 13, 2011
 *      Author: uvilla
 */

//@HEADER

/*!
 * \file ApproximatedInvertibleRowMatrix.hpp
 * \author Umberto Villa
 * \date 2011-10-13
 * This file contains the definition of the class ApproximatedInvertibleRowMatrix.
 * ApproximatedInvertibleRowMatrix is derived by the LinearOperator class, and it can be used to
 * approximate the inverse of a Epetra_CsrMatrix, either by using the Krylov methods in Belos or AztecOO or
 * the preconditioning techniques in ML or Ifpack.
 */

#ifndef APPROXIMATEDINVERTIBLEROWMATRIX_H_
#define APPROXIMATEDINVERTIBLEROWMATRIX_H_

#include <lifev/core/linear_algebra/LinearOperator.hpp>
#include <lifev/core/linear_algebra/RowMatrixPreconditioner.hpp>
#include <lifev/core/linear_algebra/InvertibleOperator.hpp>

namespace LifeV
{

namespace Operators
{

//! @class
/*!
 * @brief Provide the \c ApplyInverse method for a Epetra_CrsMatrix object.
 *
 * ApproximatedInvertibleRowMatrix is derived by the LinearOperator class, and it can be used to
 * approximate the inverse of a Epetra_CsrMatrix, either by using the Krylov methods in Belos or AztecOO or
 * the preconditioning techniques in ML or Ifpack.
 * Internally ApproximatedInvertibleRowMatrix contains pointers to a Epetra_CsrMatrix,
 * to an InvertibleOperator object (that provides the interface to AztecOO or Belos), and RowMatrixPreconditioner
 * (the provides the interface to Ifpack or ML).
 *
 * The parameters to set up the Krylov methods or the preconditioner are given by using a Teuchos Parameter List.
 * Below how the parameter list would look like:
 *
 * \verbatim
 *         <ParameterList name="ApproximatedInvertibleRowMatrix">
 *           <Parameter name="use preconditioner as approximated inverse" type="bool" value="false"/>
 *           <Parameter name="preconditioner type" type="string" value="ML"/> <!-- Ifpack, ML, TwoLevel-->
 *           <ParameterList name="solver">
 *              <Parameter name="Linear Solver Type" type="string" value="Belos"/>
 *              <ParameterList name="AztecOO">
 *                 <Parameter name="conv" type="string" value="r0"/>
 *                 <Parameter name="max_iter" type="int" value="1500"/>
 *                 <Parameter name="output" type="string" value="warnings"/>
 *                 <Parameter name="scaling" type="string" value="none"/>
 *                 <Parameter name="solver" type="string" value="cg"/>
 *                 <Parameter name="tol" type="double" value="1e-03"/>
 *              </ParameterList>
 *              <ParameterList name="Belos">
 *                  <Parameter name="Solver Type" type="string" value="RCG"/>
 *                  <Parameter name="Preconditioner Side" type="string" value="Left"/>
 *                  <ParameterList name="options">
 *                      <Parameter name="Maximum Iterations" type="int" value="50"/>
 *                      <Parameter name="Num Blocks" type="int" value="50"/>
 *                      <Parameter name="Num Recycled Blocks" type="int" value="10"/>
 *                      <Parameter name="Convergence Tolerance" type="double" value="1e-03"/>
 *                      <Parameter name="Output Frequency" type="int" value="1"/>
 *                      <Parameter name="Verbosity" type="int" value="0"/>
 *                  </ParameterList>
 *              </ParameterList>
 *           </ParameterList>
 *           <ParameterList name="preconditioner">
 *               <ParameterList name="ML">
 *                   <ParameterList name="options">
 *                       <Parameter name="default values" type="string" value="SA"/>
 *                   </ParameterList>
 *               </ParameterList>
 *           </ParameterList>
 *       </ParameterList>
 *  \endverbatim
 *
 *
 */

class ApproximatedInvertibleRowMatrix : public LinearOperator
{
public:

	//! @name Public Typedef
	//@{
	typedef Epetra_CrsMatrix rowMatrix_Type;
	typedef boost::shared_ptr<rowMatrix_Type> rowMatrixPtr_Type;
	typedef Teuchos::ParameterList pList_Type;
	//@}

	ApproximatedInvertibleRowMatrix();

	virtual ~ApproximatedInvertibleRowMatrix();

    //! @name Attribute set methods
    //@{
    void SetRowMatrix(const rowMatrixPtr_Type & rowMatrix);

    void SetParameterList(const pList_Type pList);
    //! not fully supported!
    virtual int SetUseTranspose(bool UseTranspose);
    //@}

    int Compute();

    boost::shared_ptr<RowMatrixPreconditioner> & Preconditioner(){return M_prec;}

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

    bool usePreconditionerAsApproximatedInverse;
    rowMatrixPtr_Type M_rowMatrix;
    boost::shared_ptr<RowMatrixPreconditioner> M_prec;
    boost::shared_ptr<InvertibleOperator> M_linSolver;
    pList_Type M_pList;

};

}

}

#endif /* APPROXIMATEDINVERTIBLEROWMATRIX_H_ */
