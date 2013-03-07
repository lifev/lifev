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
    @brief Linear Operator

    @author Umberto Villa <umberto.villa@gmail.com>

    @date 03-09-2010
 */

#ifndef LINEAROPERATOR_HPP_
#define LINEAROPERATOR_HPP_


#include <lifev/core/LifeV.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
namespace Operators
{
//! @class LinearOperator
/*! @brief Abstract class which defines the interface of a Linear Operator.
 *
 * LinearOperator is an abstract which inherits from Epetra_Operator.
 * LinearOperator should be the base class for all the LifeV class which implements a linear operator.
 *
 * LinearOperator ensures perfect compatibility with all the Trilinos solvers,
 * plus it supports directly the LifeV::VectorEpetra data.
 *
 * Two concrete methods are implemented in LinearOperator
 * int apply(const VectorEpetra &X, VectorEpetra & Y) const ;
 * int applyInverse(const VectorEpetra &X, VectorEpetra &Y) const.
 *
 *
 * Such methods extract a raw Epetra_MultiVector from the VectorEpetra and then call the virtual methods
 * Int Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
 * or
 * Int ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
 * respectively.
 *
 *
 */
class LinearOperator : public Epetra_Operator
{
public:

    //! @name Public Types
    //@{
    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> commPtr_Type;
    typedef Epetra_Map map_Type;
    typedef boost::shared_ptr<map_Type> mapPtr_Type;
    typedef boost::shared_ptr<const map_Type> constMapPtr_Type;
    typedef Epetra_Operator operator_Type;
    typedef boost::shared_ptr<operator_Type> operatorPtr_Type;
    typedef Epetra_MultiVector vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
    //@}

    //! @name Destructor
    //@{
    //! Destructor
    virtual ~LinearOperator() {}
    //@}

    //! @name Attribute set methods
    //@{

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
    does not support transpose use, this method should return a value of -1.

    \param In
       UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
    */
    virtual int SetUseTranspose (bool UseTranspose) = 0;
    //@}

    //! @name Mathematical functions
    //@{

    //! Returns the result of a raw_operator applied to a raw_vector X in Y.
    /*!
    \param In
       X - A raw_vector of dimension NumVectors to multiply with matrix.
    \param Out
       Y -A raw_vector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    */
    virtual int Apply (const vector_Type& X, vector_Type& Y) const = 0;

    //! Returns the result of a LinearOperator applied to a VectorEpetra X in Y.
    /*!
    \param In
       X - A VectorEpetra to multiply with matrix.
    \param Out
       Y -A VectorEpetra containing result.

    \return Integer error code, set to 0 if successful.
    */
    inline int apply (const VectorEpetra& X, VectorEpetra& Y) const
    {
        return Apply (X.epetraVector(), Y.epetraVector() );
    }

    //! Returns the result of a raw_operator inverse applied to an raw_vector X in Y.
    /*!
    \param In
       X - A raw_vector of dimension NumVectors to solve for.
    \param Out
       Y -A raw_vector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must
              support the case where X and Y are the same object.
    */
    virtual int ApplyInverse (const vector_Type& X, vector_Type& Y) const = 0;

    //! Returns the result of a LinearOperator inverse applied to an VectorEpetra X in Y.
    /*!
    \param In
       X - A VectorEpetra to solve for.
    \param Out
       Y -A VectorEpetra containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method must
              support the case where X and Y are the same object.
    */

    inline int applyInverse (const VectorEpetra& X, VectorEpetra& Y)
    {
        return ApplyInverse (X.epetraVector(), Y.epetraVector() );
    }

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
    */
    virtual double NormInf() const = 0;
    //@}

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const = 0;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const = 0;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const = 0;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type& Comm() const = 0;

    //! Returns the raw_map object associated with the domain of this operator.
    virtual const map_Type& OperatorDomainMap() const = 0;

    //! Returns the raw_map object associated with the range of this operator.
    virtual const map_Type& OperatorRangeMap() const = 0;
    //@}

};

} /*end namespace Operators*/
} /*end namespace */
#endif /* LINEAROPERATOR_HPP_ */
