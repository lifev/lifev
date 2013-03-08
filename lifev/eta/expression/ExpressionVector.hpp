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

#ifndef EXPRESSION_VECTOR_HPP
#define EXPRESSION_VECTOR_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>

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
template < UInt VectorDim >
class ExpressionVector : public ExpressionBase<ExpressionVector<VectorDim> >
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionVector<VectorDim> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the vector of values
    ExpressionVector (const VectorSmall<VectorDim>& myValue)
        : base_Type(), M_value (myValue) {}

    //! Copy constructor
    ExpressionVector (const ExpressionVector<VectorDim>& expr)
        : base_Type(), M_value (expr.M_value) {}

    //! Destructor
    ~ExpressionVector() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "vector[" << VectorDim << "] ";
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the vector of values
    const VectorSmall<VectorDim>& value() const
    {
        return M_value;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    ExpressionVector();

    //@}

    VectorSmall<VectorDim> M_value;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>VectorDim</i>: The dimension (size) of the vector to be represented.
*/
template<UInt VectorDim>
inline ExpressionVector<VectorDim>
value (const VectorSmall<VectorDim>& myValue)
{
    return ExpressionVector<VectorDim> (myValue);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
