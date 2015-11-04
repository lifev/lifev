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
    @brief File where the structures for the division between expressions are defined.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_DIVISION_HPP
#define EXPRESSION_DIVISION_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionDivision  Class for representing a quotient between two expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the division in the expression tree.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the division operation.

  <i>RExpressionType</i>: The expression on the right side of the division operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Copiable, static display method

  <i>RExpressionType</i>: Copiable, static display method

*/
template <typename LExpressionType, typename RExpressionType>
class ExpressionDivision : public ExpressionBase< ExpressionDivision<LExpressionType, RExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // Not usefull, just for ease of coding
    typedef ExpressionBase< ExpressionDivision <LExpressionType, RExpressionType> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Full constructor using the two expressions
    ExpressionDivision (const LExpressionType& l, const RExpressionType& r)
        : base_Type(), M_l (l), M_r (r) {}

    //! Copy constructor
    ExpressionDivision (const ExpressionDivision<LExpressionType, RExpressionType>& expression)
        : base_Type(), M_l (expression.M_l), M_r (expression.M_r) {}

    //! Destructor
    ~ExpressionDivision() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        LExpressionType::display (out);
        out << " / ";
        RExpressionType::display (out);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the left side of the division
    const LExpressionType& left() const
    {
        return M_l;
    }

    //! Getter for the right side of the division
    const RExpressionType& right() const
    {
        return M_r;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    ExpressionDivision();

    //@}

    // Left hand side
    LExpressionType M_l;

    // Right hand side
    RExpressionType M_r;
};

//! operator/  The generic operator for the division between expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Operator used in construction of the expression tree. To avoid shadowing
  other operator/, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the division operation.

  <i>RExpressionType</i>: The expression on the right side of the division operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Same as in LifeV::ExpressionDivision

  <i>RExpressionType</i>: Same as in LifeV::ExpressionDivision

*/
template< typename LExpressionType, typename RExpressionType >
ExpressionDivision<LExpressionType, RExpressionType>
operator/ (const ExpressionBase<LExpressionType>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionDivision<LExpressionType, RExpressionType> (l.cast(), r.cast() );
}

// Specialization for the real constants
template< typename LExpressionType >
ExpressionDivision<LExpressionType, ExpressionScalar >
operator/ (const ExpressionBase<LExpressionType>& l, const Real& r)
{
    return ExpressionDivision<LExpressionType, ExpressionScalar> (l.cast(), ExpressionScalar (r) );
}

template< typename RExpressionType >
ExpressionDivision<ExpressionScalar, RExpressionType>
operator/ (const Real& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionDivision<ExpressionScalar, RExpressionType> (ExpressionScalar (l), r.cast() );
}

// Specialization for the vectorial constants
template< typename RExpressionType , UInt Vdim>
ExpressionDivision<ExpressionVector<Vdim>, RExpressionType>
operator/ (const VectorSmall<Vdim>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionDivision<ExpressionVector<Vdim>, RExpressionType> (ExpressionVector<Vdim> (l), r.cast() );
}

template< typename LExpressionType, UInt Vdim >
ExpressionDivision<LExpressionType, ExpressionVector<Vdim> >
operator/ (const ExpressionBase<LExpressionType>& l, const VectorSmall<Vdim>& r)
{
    return ExpressionDivision<LExpressionType, ExpressionVector<Vdim> > (l.cast(), ExpressionVector<Vdim> (r) );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
