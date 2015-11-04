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
    @brief File where the structures for the product between expressions are defined.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_PRODUCT_HPP
#define EXPRESSION_PRODUCT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>
#include <lifev/eta/expression/ExpressionMatrix.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionProduct  Class for representing a product between two expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the product in the expression tree.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the product operation.

  <i>RExpressionType</i>: The expression on the right side of the product operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Copiable, static display method

  <i>RExpressionType</i>: Copiable, static display method

*/
template <typename LExpressionType, typename RExpressionType>
class ExpressionProduct : public ExpressionBase< ExpressionProduct<LExpressionType, RExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // No direct use, just ease of coding
    typedef ExpressionBase< ExpressionProduct <LExpressionType, RExpressionType> > base_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Full constructor using the two expressions
    ExpressionProduct (const LExpressionType& l, const RExpressionType& r)
        : base_Type(), M_l (l), M_r (r) {}

    //! Copy constructor
    ExpressionProduct (const ExpressionProduct<LExpressionType, RExpressionType>& expression)
        : base_Type(), M_l (expression.M_l), M_r (expression.M_r) {}

    //! Destructor
    ~ExpressionProduct() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        LExpressionType::display (out);
        out << " * ";
        RExpressionType::display (out);
    }

    //@}


    //! @name Private Methods
    //@{

    //! Getter for the left hand side
    const LExpressionType& left() const
    {
        return M_l;
    }

    //! Getter for the right hand side
    const RExpressionType& right() const
    {
        return M_r;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    ExpressionProduct();

    //@}


    // Left hand side
    LExpressionType M_l;

    // Right hand side
    RExpressionType M_r;
};


//! operator*  The generic operator for the product between expressions.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Operator used in construction of the expression tree. To avoid shadowing
  other operator*, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>LExpressionType</i>: The expression on the left side of the product operation.

  <i>RExpressionType</i>: The expression on the right side of the product operation.

  <b> Template requirements </b>

  <i>LExpressionType</i>: Same as in LifeV::ExpressionProduct

  <i>RExpressionType</i>: Same as in LifeV::ExpressionProduct

*/
template< typename LExpressionType, typename RExpressionType >
ExpressionProduct<LExpressionType, RExpressionType>
operator* (const ExpressionBase<LExpressionType>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionProduct<LExpressionType, RExpressionType> (l.cast(), r.cast() );
}

// "Specialization" for the case of a scalar
template< typename RExpressionType >
ExpressionProduct<ExpressionScalar, RExpressionType>
operator* (const Real& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionProduct<ExpressionScalar, RExpressionType> (ExpressionScalar (l), r.cast() );
}

// "Specialization" for the case of a scalar
template< typename LExpressionType >
ExpressionProduct<LExpressionType, ExpressionScalar>
operator* (const ExpressionBase<LExpressionType>& l, const Real& r)
{
    return ExpressionProduct<LExpressionType, ExpressionScalar> (l.cast(), ExpressionScalar (r) );
}

// "Specialization" for the case of a vector
template< typename RExpressionType, UInt Vdim >
ExpressionProduct<ExpressionVector<Vdim>, RExpressionType>
operator* (const VectorSmall<Vdim>& l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionProduct<ExpressionVector<Vdim>, RExpressionType> (ExpressionVector<Vdim> (l), r.cast() );
}

// "Specialization" for the case of a vector
template< typename LExpressionType, UInt Vdim>
ExpressionProduct<LExpressionType, ExpressionVector<Vdim> >
operator* (const ExpressionBase<LExpressionType>& l, const VectorSmall<Vdim>& r)
{
    return ExpressionProduct<LExpressionType, ExpressionVector<Vdim> > (l.cast(), ExpressionVector<Vdim> (r) );
}

// "Specialization" for the case of a matrix
template< typename RExpressionType, UInt Dim1 , UInt Dim2 >
ExpressionProduct<ExpressionMatrix<Dim1, Dim2>, RExpressionType>
operator* (const MatrixSmall<Dim1, Dim2>&  l, const ExpressionBase<RExpressionType>& r)
{
    return ExpressionProduct<ExpressionMatrix<Dim1, Dim2>, RExpressionType> (ExpressionMatrix<Dim1, Dim2> (l), r.cast() );
}

// "Specialization" for the case of a matrix
template< typename LExpressionType, UInt Dim1, UInt Dim2 >
ExpressionProduct<LExpressionType, ExpressionMatrix<Dim1, Dim2> >
operator* (const ExpressionBase<LExpressionType>& l, const MatrixSmall<Dim1, Dim2>& r)
{
    return ExpressionProduct<LExpressionType, ExpressionMatrix<Dim1, Dim2> > (l.cast(), ExpressionMatrix<Dim1, Dim2> (r) );
}





} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
