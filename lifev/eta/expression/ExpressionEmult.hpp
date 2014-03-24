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
    @brief File where the structures for the element-wise multiplication between expressions are defined.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 08-2012
 */

#ifndef EXPRESSION_EMULT_HPP
#define EXPRESSION_EMULT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>
#include <lifev/eta/expression/ExpressionMatrix.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionEmult  Class for representing a element-wise multiplication between two expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the element-wise multiplication in the expression tree.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the emult operation.

  <i>RExpressionType</i>: The expression on the right side of the emult operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Copiable, static display method

  <i>RExpressionType</i>: Copiable, static display method

*/
template <typename LExpressionType, typename RExpressionType>
class ExpressionEmult : public ExpressionBase< ExpressionEmult<LExpressionType, RExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // No real need, just for ease of coding
    typedef ExpressionBase< ExpressionEmult <LExpressionType, RExpressionType> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Full constructor, with the two expressions.
    ExpressionEmult (const LExpressionType& l, const RExpressionType& r)
        : base_Type(), M_l (l), M_r (r) {}

    //! Copy constructor
    ExpressionEmult (const ExpressionEmult<LExpressionType, RExpressionType>& expression)
        : base_Type(), M_l (expression.M_l), M_r (expression.M_r) {}

    //! Destructor
    ~ExpressionEmult() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        LExpressionType::display (out);
        out << " emult ";
        RExpressionType::display (out);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the left hand side of the emult operation
    const LExpressionType& left() const
    {
        return M_l;
    }

    //! Getter for the right hand side of the emult operation
    const RExpressionType& right() const
    {
        return M_r;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No default constructor
    ExpressionEmult();

    //@}

    // Left hand side
    LExpressionType M_l;

    // Right hand side
    RExpressionType M_r;
};


// emult  The generic function for the element-wise multiplication between expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Function used in construction of the expression tree. To avoid shadowing
  other functions, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the emult operation.

  <i>RExpressionType</i>: The expression on the right side of the emult operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Same as in LifeV::ExpressionEmult

  <i>RExpressionType</i>: Same as in LifeV::ExpressionEmult

*/
template< typename LExpressionType, typename RExpressionType >
ExpressionEmult<LExpressionType, RExpressionType>
emult (const ExpressionBase<LExpressionType>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionEmult<LExpressionType, RExpressionType> (l.cast(), r.cast() );
}


// Specialization for the matricial constants
template< typename RExpressionType , UInt Dim1, UInt Dim2>
ExpressionEmult<ExpressionMatrix<Dim1, Dim2>, RExpressionType>
emult (const MatrixSmall<Dim1, Dim2>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionEmult<ExpressionMatrix<Dim1, Dim2>, RExpressionType> (ExpressionMatrix<Dim1, Dim2> (l), r.cast() );
}

template< typename LExpressionType, UInt Dim1, UInt Dim2 >
ExpressionEmult<LExpressionType, ExpressionMatrix<Dim1, Dim2> >
emult (const ExpressionBase<LExpressionType>& l, const MatrixSmall<Dim1, Dim2>& r)
{
    return ExpressionEmult<LExpressionType, ExpressionMatrix<Dim1, Dim2> > (l.cast(), ExpressionMatrix<Dim1, Dim2> (r) );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
