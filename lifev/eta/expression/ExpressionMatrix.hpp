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
    @brief File containing the expression to represent a matricial constant

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_MATRIX_HPP
#define EXPRESSION_MATRIX_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionMatrix  Class representing a matricial constant value in an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>iDim</i>: The i dimension (size) of the matrix to be represented.
  <i>jDim</i>: The j dimension (size) of the matrix to be represented.

*/
template < UInt iDim, UInt jDim >
class ExpressionMatrix : public ExpressionBase<ExpressionMatrix<iDim, jDim> >
{
public:

    //! @name Public Types
    //@{

    typedef ExpressionBase<ExpressionMatrix<iDim, jDim> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the vector of values
    ExpressionMatrix (const MatrixSmall<iDim, jDim>& myValue)
        : base_Type(), M_value (myValue) {}

    //! Copy constructor
    ExpressionMatrix (const ExpressionMatrix<iDim, jDim>& expr)
        : base_Type(), M_value (expr.M_value) {}

    //! Destructor
    ~ExpressionMatrix() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "matrix[" << iDim << "][" << jDim << "] ";
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the vector of values
    const MatrixSmall<iDim, jDim>& value() const
    {
        return M_value;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    ExpressionMatrix();

    //@}

    MatrixSmall<iDim, jDim> M_value;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>VectorDim</i>: The dimension (size) of the vector to be represented.
*/
template<UInt iDim, UInt jDim>
inline ExpressionMatrix<iDim, jDim>
value (const MatrixSmall<iDim, jDim>& myValue)
{
    return ExpressionMatrix<iDim, jDim> (myValue);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
