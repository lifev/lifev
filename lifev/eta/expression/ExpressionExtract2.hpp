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
    @brief File for the definition of the expression used for the extraction of a component from a matrix.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_EXTRACT2_HPP
#define EXPRESSION_EXTRACT2_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionExtract2  Class for representing the extraction of a component from a matrix (2 indexes specified).
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class represents the extraction of a component from a matrix in the expression tree.

  <b> Template parameters </b>

  <i>ExpressionType</i>: The expression from which we want to extract

  <b> Template requirements </b>

  <i>ExpressionType</i>: Copiable, static display method
*/

template< typename ExpressionType >
class ExpressionExtract2 : public ExpressionBase<ExpressionExtract2< ExpressionType > >
{
public:

    //! @name Public Types
    //@{

    // Base class, typedef useful to make the code cleaner
    typedef ExpressionBase<ExpressionExtract2< ExpressionType > > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Full constructor, with the expression and the specification of the position of the row to be extracted
    ExpressionExtract2 (const ExpressionType& ex, const UInt& i, const UInt& j)
        : base_Type(),
          M_ex (ex),
          M_i (i),
          M_j (j)
    {}

    //! Copy constructor
    ExpressionExtract2 (const ExpressionExtract2<ExpressionType>& expr)
        : base_Type(),
          M_ex (expr.M_ex),
          M_i (expr.M_i),
          M_j (expr.M_j)
    {}

    //! Destructor
    ~ExpressionExtract2() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "Extraction from ";
        ExpressionType::display (out);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the row index
    UInt indexI() const
    {
        return M_i;
    }

    //! Getter for the column index
    UInt indexJ() const
    {
        return M_j;
    }

    //! Getter for the expression from which we extract
    const ExpressionType& exprEx() const
    {
        return M_ex;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    ExpressionExtract2();

    //@}

    // Storage for the row index
    UInt M_i;

    // Storage for the column index
    UInt M_j;

    // Storage for the expression from which we extract
    ExpressionType M_ex;
};

//! extract The generic function for the extraction
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  Function used in construction of the expression tree. To avoid shadowing
  other functions, it uses the ExpressionBase type to distinguish expressions
  from other types.

  <b> Template parameters </b>

  <i>ExpressionType</i>: The expression from which we extract

  <b> Template requirements </b>

  <i>ExpressionType</i>: Same as LifeV::ExpressionExtract2
*/

template< typename ExpressionType >
inline ExpressionExtract2<ExpressionType>
extract (const ExpressionBase<ExpressionType>& ex, const UInt& i, const UInt& j)
{
    return ExpressionExtract2<ExpressionType> (ex.cast(), i, j);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
