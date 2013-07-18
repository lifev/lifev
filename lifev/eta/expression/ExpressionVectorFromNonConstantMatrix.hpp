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
    @brief File containing the expression to represent a vectorial constant

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_VECTORFROMNONCONSTANTMATRIX_HPP
#define EXPRESSION_VECTORFROMNONCONSTANTMATRIX_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionVector  Class representing a constant vectorial value in an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>VectorDim</i>: The dimension (size) of the vector to be represented.

*/
template <typename ExpressionType, UInt SpaceDim,  UInt FieldDim >
class ExpressionVectorFromNonConstantMatrix : public ExpressionBase<ExpressionVectorFromNonConstantMatrix<ExpressionType, SpaceDim, FieldDim> >
{
public:

    //! @name Public Types
    //@{

  typedef ExpressionBase<ExpressionVectorFromNonConstantMatrix<ExpressionType, SpaceDim, FieldDim> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the vector of values
    ExpressionVectorFromNonConstantMatrix (const ExpressionType& expression, const UInt column)
        : base_Type(), M_expr (expression), M_column(column) {}

    //! Copy constructor
    ExpressionVectorFromNonConstantMatrix (const ExpressionVectorFromNonConstantMatrix<ExpressionType, SpaceDim, FieldDim>& expr)
        : base_Type(), M_expr (expr.M_expr), M_column(expr.M_column) {}

    //! Destructor
    ~ExpressionVectorFromNonConstantMatrix() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "vector from non constant matrix ";
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the vector of values
    const ExpressionType& expr() const
    {
        return M_expr;
    }

    //! Getter for the vector of values
    const UInt& column() const
    {
        return M_column;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    ExpressionVectorFromNonConstantMatrix();

    //@}

    ExpressionType M_expr;
    UInt           M_column;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>VectorDim</i>: The dimension (size) of the vector to be represented.
*/

template< typename ExpressionType, UInt SpaceDim, UInt FieldDim >
ExpressionVectorFromNonConstantMatrix<ExpressionType, SpaceDim, FieldDim >
vectorFromMatrix (const ExpressionBase<ExpressionType>& expr, const UInt column)
{
  return ExpressionVectorFromNonConstantMatrix<ExpressionType, SpaceDim, FieldDim> (expr.cast(), column );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
