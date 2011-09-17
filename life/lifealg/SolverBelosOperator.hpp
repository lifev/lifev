//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains the SolverBelosOperator class.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 01-09-2011

    The class provides a way to use a the SolverBelos as an operator.
    This is particularly useful to apply "exactly the inverse of a
    block in a preconditioner.
 */

#ifndef SOLVERBELOSOPERATOR_HPP
#define SOLVERBELOSOPERATOR_HPP 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Operator.h>
#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Operator.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <vector>
#include <boost/shared_ptr.hpp>
#include <life/lifecore/LifeV.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifearray/VectorEpetra.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifecore/Displayer.hpp>
#include <life/lifealg/SolverBelos.hpp>
#include <life/lifealg/Preconditioner.hpp>

class GetPot;

namespace LifeV
{

//! SolverBelosOperator - Class to wrap linear solver
/*!
  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
*/
class SolverBelosOperator
    : public Epetra_Operator
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_Operator                      prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef MatrixEpetra<Real>                   operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type> operator_type;

    typedef Displayer::comm_Type                 comm_Type;
    typedef Displayer::commPtr_Type              commPtr_Type;

    typedef Teuchos::ParameterList               list_Type;

    typedef SolverBelos                          linearSolver_type;
    typedef boost::shared_ptr<linearSolver_type> linearSolverPtr_type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Default constructor
    /*!
     * @param comm The communicator.
     */
    SolverBelosOperator( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm>() );

    //! Destructor
    virtual ~SolverBelosOperator();

    //@}

    //! @name Methods
    //@{

    //! Build a preconditioner based on the given matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
     */
    Int buildPreconditioner( operator_type& matrix, const list_Type& list );

    //! Reset the preconditioner
    void resetPreconditioner();

    //! Show informations about the preconditioner
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Epetra Operator Interface Methods
    //@{

    //! Set the matrix to be used transposed (or not)
    /*!
      @param useTranspose If true the preconditioner is transposed
     */
    Int SetUseTranspose( const bool useTranspose = false );

    //! Return true if the preconditioner is transposed
    bool UseTranspose();

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    Int Apply( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    Int ApplyInverse( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const;

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       @warning This method must not be called unless HasNormInf() returns true.
     */
    double NormInf() const;

    //! Returns a character string describing the operator
    const char* Label() const;

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    bool HasNormInf() const;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const;

    //! Return the Range map of the operator
    const Epetra_Map& OperatorRangeMap() const;

    //! Return the Domain map of the operator
    const Epetra_Map& OperatorDomainMap() const;
    //@}

private:

    linearSolverPtr_type           M_solver;
    operator_type                  M_matrix;
    boost::shared_ptr<Epetra_Comm> M_comm;
};

} // namespace LifeV

#endif /* SOLVERBELOSOPERATOR_HPP */
