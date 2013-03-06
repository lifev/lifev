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

#ifndef EXPRESSION_MATRIX_HPP
#define EXPRESSION_MATRIX_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/array/MatrixSmall.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionVector  Class representing a constant matrix value in an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>MatrixDim1, MatrixDim2</i>: The dimensions (size) of the matrix to be represented.

*/
template < UInt MatrixDim1, UInt MatrixDim2 >
class ExpressionMatrix : public ExpressionBase<ExpressionMatrix<MatrixDim1,MatrixDim2> >
{
public:

    //! @name Public Types
    //@{

	typedef ExpressionBase<ExpressionMatrix<MatrixDim1,MatrixDim2> > base_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor using the vector of values
	ExpressionMatrix(const MatrixSmall<MatrixDim1,MatrixDim2>& myValue)
	: base_Type(), M_value(myValue) {}

    //! Copy constructor
	ExpressionMatrix(const ExpressionMatrix<MatrixDim1,MatrixDim2>& expr)
	: base_Type(), M_value(expr.M_value) {}

    //! Destructor
    ~ExpressionMatrix(){}

    //@}


    //! @name Methods
    //@{

    //! Display method
	static void display(std::ostream& out= std::cout)
	{ out << "matrix[" << MatrixDim1 << "][" << MatrixDim2 << "] ";}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the vector of values
	const MatrixSmall<MatrixDim1,MatrixDim2>& value() const { return M_value; }

    //@}

private:

    //! @name Private Methods
    //@{

    ExpressionMatrix();

    //@}

	MatrixSmall<MatrixDim1,MatrixDim2> M_value;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  <b> Template parameters </b>

  <i>MatrixDim1, MatrixDim2</i>: The dimensions (size) of the matrix to be represented.
*/
template<UInt MatrixDim1, UInt MatrixDim2>
inline ExpressionMatrix<MatrixDim1,MatrixDim2>
value(const MatrixSmall<MatrixDim1,MatrixDim2>& myValue)
{
	return ExpressionMatrix<MatrixDim1,MatrixDim2>(myValue);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
