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
    @brief File where the structures for the dot product between expressions are defined.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_DOT_HPP
#define EXPRESSION_DOT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionDot  Class for representing a dot product between two expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the dot product in the expression tree.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the dot product.

  <i>RExpressionType</i>: The expression on the right side of the dot product.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Copiable, static display method

  <i>RExpressionType</i>: Copiable, static display method

*/
template <typename LExpressionType, typename RExpressionType>
class ExpressionDot : public ExpressionBase< ExpressionDot<LExpressionType, RExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // No real need, just for ease of coding
    typedef ExpressionBase< ExpressionDot <LExpressionType, RExpressionType> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Full constructor, with the two expressions.
    ExpressionDot (const LExpressionType& l, const RExpressionType& r)
        : base_Type(), M_l (l), M_r (r) {}

    //! Copy constructor
    ExpressionDot (const ExpressionDot<LExpressionType, RExpressionType>& expression)
        : base_Type(), M_l (expression.M_l), M_r (expression.M_r) {}

    //! Destructor
    ~ExpressionDot() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        LExpressionType::display (out);
        out << " dot ";
        RExpressionType::display (out);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the left hand side of the dot product
    const LExpressionType& left() const
    {
        return M_l;
    }

    //! Getter for the right hand side of the dot product
    const RExpressionType& right() const
    {
        return M_r;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    ExpressionDot();

    //@}

    // Left hand side
    LExpressionType M_l;

    // Right hand side
    RExpressionType M_r;
};


//! dot  The generic function for the dot product between expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Function used in construction of the expression tree. To avoid shadowing
  other functions, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the dot product.

  <i>RExpressionType</i>: The expression on the right side of the dot product.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Same as in LifeV::ExpressionDot

  <i>RExpressionType</i>: Same as in LifeV::ExpressionDot

*/
template< typename LExpressionType, typename RExpressionType >
ExpressionDot<LExpressionType, RExpressionType>
dot (const ExpressionBase<LExpressionType>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionDot<LExpressionType, RExpressionType> (l.cast(), r.cast() );
}

// Specialization for the real constants
template< typename LExpressionType >
ExpressionDot<LExpressionType, ExpressionScalar >
dot (const ExpressionBase<LExpressionType>& l, const Real& r)
{
    return ExpressionDot<LExpressionType, ExpressionScalar> (l.cast(), ExpressionScalar (r) );
}

template< typename RExpressionType >
ExpressionDot<ExpressionScalar, RExpressionType>
dot (const Real& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionDot<ExpressionScalar, RExpressionType> (ExpressionScalar (l), r.cast() );
}

// Specialization for the vectorial constants
template< typename RExpressionType , UInt Vdim>
ExpressionDot<ExpressionVector<Vdim>, RExpressionType>
dot (const VectorSmall<Vdim>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionDot<ExpressionVector<Vdim>, RExpressionType> (ExpressionVector<Vdim> (l), r.cast() );
}

template< typename LExpressionType, UInt Vdim >
ExpressionDot<LExpressionType, ExpressionVector<Vdim> >
dot (const ExpressionBase<LExpressionType>& l, const VectorSmall<Vdim>& r)
{
    return ExpressionDot<LExpressionType, ExpressionVector<Vdim> > (l.cast(), ExpressionVector<Vdim> (r) );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
