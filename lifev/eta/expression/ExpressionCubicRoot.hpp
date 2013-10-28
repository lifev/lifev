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

#ifndef EXPRESSION_CUBICROOT_HPP
#define EXPRESSION_CUBICROOT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionCubicRoot  Class for representing a product between two expressions.
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
template <typename BaseExpressionType>
class ExpressionCubicRoot : public ExpressionBase< ExpressionCubicRoot<BaseExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // No direct use, just ease of coding
    typedef ExpressionBase< ExpressionCubicRoot <BaseExpressionType> > base_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Full constructor using the two expressions
    ExpressionCubicRoot (const BaseExpressionType& l)
        : base_Type(), M_l (l) {}

    //! Copy constructor
    ExpressionCubicRoot (const ExpressionCubicRoot<BaseExpressionType>& expression)
        : base_Type(), M_l (expression.M_l) {}

    //! Destructor
    ~ExpressionCubicRoot() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        BaseExpressionType::display (out);
        out << " ^ ";
    }

    //@}


    //! @name Private Methods
    //@{

    //! Getter for the left hand side
    const BaseExpressionType& base() const
    {
        return M_l;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    ExpressionCubicRoot();

    //@}


    // Left hand side
    BaseExpressionType M_l;
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

template< typename  ExpressionType>
ExpressionCubicRoot<ExpressionType>
cubicroot (const ExpressionBase<ExpressionType>& l)
{
    return ExpressionCubicRoot<ExpressionType> ( l.cast() );
}


  /* 
     This expression is a specialization of the cubic root expression
     for the structure module to have a fast computation of the isochoric
     change of variable used in the structure module for the constitutive models
  */
template <typename BaseExpressionType>
class ExpressionIsochoricChangeOfVariable : 
    public ExpressionBase< ExpressionIsochoricChangeOfVariable< BaseExpressionType> >
{
public:

    //! @name Public Types
    //@{

    // No direct use, just ease of coding
    typedef ExpressionBase< ExpressionIsochoricChangeOfVariable <BaseExpressionType> > base_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Full constructor using the two expressions
    ExpressionIsochoricChangeOfVariable (const BaseExpressionType& l)
        : base_Type(), M_l (l) {}

    //! Copy constructor
    ExpressionIsochoricChangeOfVariable (const ExpressionIsochoricChangeOfVariable<BaseExpressionType>& expression)
        : base_Type(), M_l (expression.M_l) {}

    //! Destructor
    ~ExpressionIsochoricChangeOfVariable() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        BaseExpressionType::display (out);
        out << " ^ ";
    }

    //@}


    //! @name Private Methods
    //@{

    //! Getter for the left hand side
    const BaseExpressionType& base() const
    {
        return M_l;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    ExpressionIsochoricChangeOfVariable();

    //@}


    // Left hand side
    BaseExpressionType M_l;
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

template< typename  ExpressionType>
ExpressionIsochoricChangeOfVariable<ExpressionType>
isoCoV (const ExpressionBase<ExpressionType>& l)
{
    return ExpressionIsochoricChangeOfVariable<ExpressionType> ( l.cast() );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
