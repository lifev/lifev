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
    @brief File where the structures for the addition between expressions are defined.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_ADDITION_HPP
#define EXPRESSION_ADDITION_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>
#include <lifev/eta/expression/ExpressionMatrix.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionAddition  Class for representing an addition between two expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the addition in the expression tree.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the addition operation.

  <i>RExpressionType</i>: The expression on the right side of the addition operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Copiable, static display method

  <i>RExpressionType</i>: Copiable, static display method

*/
template <typename LExpressionType, typename RExpressionType>
class ExpressionAddition : public ExpressionBase< ExpressionAddition<LExpressionType, RExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // Base class (no need, just for ease of coding this class)
    typedef ExpressionBase< ExpressionAddition <LExpressionType, RExpressionType> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Complete constructor using the two expressions to be summed
    ExpressionAddition (const LExpressionType& l, const RExpressionType& r)
        : base_Type(), M_l (l), M_r (r) {}

    //! Copy constructor
    ExpressionAddition (const ExpressionAddition<LExpressionType, RExpressionType>& expression)
        : base_Type(), M_l (expression.M_l), M_r (expression.M_r) {}

    //! Destructor
    ~ExpressionAddition() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        LExpressionType::display (out);
        out << " + ";
        RExpressionType::display (out);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the expression on the left side of the addition
    const LExpressionType& left() const
    {
        return M_l;
    }

    //! Getter for the expression on the right side of the addition
    const RExpressionType& right() const
    {
        return M_r;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    ExpressionAddition();

    //@}

    // Left side of the operation
    LExpressionType M_l;

    // Right side of the operation
    RExpressionType M_r;
};


//! operator+  The generic operator for the addition between expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Operator used in construction of the expression tree. To avoid shadowing
  other operator+, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the addition operation.

  <i>RExpressionType</i>: The expression on the right side of the addition operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Same as in LifeV::ExpressionAddition

  <i>RExpressionType</i>: Same as in LifeV::ExpressionAddition

*/
template< typename LExpressionType, typename RExpressionType >
ExpressionAddition<LExpressionType, RExpressionType>
operator+ (const ExpressionBase<LExpressionType>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionAddition<LExpressionType, RExpressionType> (l.cast(), r.cast() );
}

// Specialization for the real constants
template< typename LExpressionType >
ExpressionAddition<LExpressionType, ExpressionScalar >
operator+ (const ExpressionBase<LExpressionType>& l, const Real& r)
{
    return ExpressionAddition<LExpressionType, ExpressionScalar> (l.cast(), ExpressionScalar (r) );
}

template< typename RExpressionType >
ExpressionAddition<ExpressionScalar, RExpressionType>
operator+ (const Real& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionAddition<ExpressionScalar, RExpressionType> (ExpressionScalar (l), r.cast() );
}

// Specialization for the vectorial constants
template< typename RExpressionType , UInt Vdim>
ExpressionAddition<ExpressionVector<Vdim>, RExpressionType>
operator+ (const VectorSmall<Vdim>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionAddition<ExpressionVector<Vdim>, RExpressionType> (ExpressionVector<Vdim> (l), r.cast() );
}

template< typename LExpressionType, UInt Vdim >
ExpressionAddition<LExpressionType, ExpressionVector<Vdim> >
operator+ (const ExpressionBase<LExpressionType>& l, const VectorSmall<Vdim>& r)
{
    return ExpressionAddition<LExpressionType, ExpressionVector<Vdim> > (l.cast(), ExpressionVector<Vdim> (r) );
}

// Specialization for the matrix constants
template< typename RExpressionType , UInt Dim1, UInt Dim2>
ExpressionAddition<ExpressionMatrix<Dim1, Dim2>, RExpressionType>
operator+ (const MatrixSmall<Dim1, Dim2>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionAddition<ExpressionMatrix<Dim1, Dim2>, RExpressionType> (ExpressionMatrix<Dim1, Dim2> (l), r.cast() );
}

template< typename LExpressionType, UInt Dim1, UInt Dim2 >
ExpressionAddition<LExpressionType, ExpressionMatrix<Dim1, Dim2> >
operator+ (const ExpressionBase<LExpressionType>& l, const MatrixSmall<Dim1, Dim2>& r)
{
    return ExpressionAddition<LExpressionType, ExpressionMatrix<Dim1, Dim2> > (l.cast(), ExpressionMatrix<Dim1, Dim2> (r) );
}

template<  UInt Dim1, UInt Dim2 >
ExpressionAddition< ExpressionMatrix<Dim1, Dim2>, ExpressionMatrix<Dim1, Dim2> >
operator+ (const MatrixSmall<Dim1, Dim2>& l, const MatrixSmall<Dim1, Dim2>& r)
{
    return ExpressionAddition< ExpressionMatrix<Dim1, Dim2>, ExpressionMatrix<Dim1, Dim2> > (ExpressionMatrix<Dim1, Dim2> (l), ExpressionMatrix<Dim1, Dim2> (r) );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
