//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Solver Operator

    @author Umberto Villa <umberto.villa@gmail.com>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-09-2010
 */

#ifndef _SOLVEROPERATOR_HPP_
#define _SOLVEROPERATOR_HPP_


#include <lifev/core/operator/LinearOperator.hpp>
#include <lifev/core/util/FactorySingleton.hpp>
#include <lifev/core/util/Factory.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Epetra_MpiComm.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

namespace LifeV
{
namespace Operators
{

//! @class SolverOperator
/*! @brief Abstract class which defines the interface of an Invertible Linear Operator.
 *
 */
class SolverOperator : public LinearOperator
{
public:
    enum SolverOperatorStatusType          { undefined, yes, no };

#ifdef HAVE_MPI
    SolverOperator ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    SolverOperator ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif
    virtual ~SolverOperator();

    //! @name Attribute set methods
    //@{

    //! If set true, transpose of this operator will be applied.
    virtual int SetUseTranspose ( bool useTranspose );

    void setOperator ( operatorPtr_Type _oper );

    void setPreconditioner ( operatorPtr_Type _prec );

    void setParameters ( const Teuchos::ParameterList& _pList );

    void setTolerance ( const Real& tolerance );

    void setUsedForPreconditioning ( const bool& enable );

    void resetCumulIterations()
    {
        M_numCumulIterations = 0;
    }

    void resetSolver();

    //@}

    //! @name Mathematical methods
    //@{

    //! Returns the result of a Epetra_Operator applied to a vector_Type X in Y.
    virtual int Apply ( const vector_Type& X, vector_Type& Y ) const;

    //! Returns the result of a Epetra_Operator inverse applied to an vector_Type X in Y.
    virtual int ApplyInverse ( const vector_Type& X, vector_Type& Y ) const;

    //! Returns the infinity norm of the global matrix.
    double NormInf() const
    {
        return M_oper->NormInf();
    }

    //! Reset the status for the state of convergence and loss of accuracy
    void resetStatus()
    {
        M_lossOfAccuracy = undefined;
        M_converged = undefined;
        M_numIterations = 0;
    }

    //@}

    //! @name Attribute access methods
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const
    {
        return M_name.c_str();
    }

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const
    {
        return M_useTranspose;
    }

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const
    {
        return M_oper->HasNormInf();
    }

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type& Comm() const
    {
        return M_oper->Comm();
    }

    //! Returns the Epetra_Map object associated with the domain of this operator.
    virtual const map_Type& OperatorDomainMap() const
    {
        return M_oper->OperatorDomainMap();
    }

    //! Returns the Epetra_Map object associated with the range of this operator.
    virtual const map_Type& OperatorRangeMap() const
    {
        return M_oper->OperatorRangeMap();
    }

    //! Returns if a loss of precision has been detected
    SolverOperatorStatusType isLossOfAccuracyDetected() const
    {
        return M_lossOfAccuracy;
    }

    //! Returns if the convergence has been achieved
    SolverOperatorStatusType hasConverged() const
    {
        return M_converged;
    }

    //! Returns the number of iterations
    /*
     * @note Negative value usually indicates problem of convergence
     */
    int numIterations() const
    {
        return M_numIterations;
    }

    //! Returns the cumul of iterations
    int numCumulIterations() const
    {
        return M_numCumulIterations;
    }

    //@}

protected:

    virtual int doApplyInverse ( const vector_Type& X, vector_Type& Y ) const = 0;
    virtual void doSetOperator() = 0;
    virtual void doSetPreconditioner() = 0;
    virtual void doSetParameterList() = 0;
    virtual void doResetSolver() = 0;

    //! The name of the Operator
    std::string M_name;

    //! The list of Parameter to feed the linear solver
    Teuchos::RCP<Teuchos::ParameterList> M_pList;

    //! The preconditioner operator
    operatorPtr_Type M_prec;

    //! The operator to be solved
    operatorPtr_Type M_oper;

    //! Whenever to use the transpose
    bool M_useTranspose;

    //! Status to see if there is a loss of accuracy
    mutable SolverOperatorStatusType M_lossOfAccuracy;

    //! Status to see if the solver has converged
    mutable SolverOperatorStatusType M_converged;

    //! Number of iterations performed by the solver
    mutable int M_numIterations;

    //! Number of cumulated iterations performed by the solver
    mutable int M_numCumulIterations;

    //! Solver tolerance
    Real M_tolerance;

    //! Print the number of iteration (used only for preconditioner LinearSolver)
    bool M_printSubiterationCount;

    //! Communicator
    boost::shared_ptr<Epetra_Comm> M_comm;

};

typedef FactorySingleton<Factory<SolverOperator, std::string> > SolverOperatorFactory;

} // Namespace Operators

} // Namespace LifeV

#endif // _SOLVEROPERATOR_HPP_
