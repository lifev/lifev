//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the definition of the EvaluationVector class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */
#ifndef EVALUATION_VECTORFROMNONCONSTANTMATRIX_HPP
#define EVALUATION_VECTORFROMNONCONSTANTMATRIX_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/eta/expression/ExpressionVectorFromNonConstantMatrix.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for a vectorial constant
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing a vectorial constant in the assembly

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <typename EvaluationType, UInt SpaceDim, UInt FieldDim>
class EvaluationVectorFromNonConstantMatrix
{
public:

    //! @name Public Types
    //@{

    //! Type of the value returned by this class
    typedef VectorSmall<FieldDim>                 return_Type;
    typedef std::vector<return_Type>              vector_Type;
    typedef MatrixSmall<FieldDim, SpaceDim>       matrix_Type;

    //@}


    //! @name Static constants
    //@{

    //! Flag for the global current FE
    const static flag_Type S_globalUpdateFlag;

    //! Flag for the test current FE
    const static flag_Type S_testUpdateFlag;

    //! Flag for the solution current FE
    const static flag_Type S_solutionUpdateFlag;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Copy constructor
    EvaluationVectorFromNonConstantMatrix (const EvaluationVectorFromNonConstantMatrix<EvaluationType, SpaceDim, FieldDim >& evaluation)
        : M_evaluation( evaluation.M_evaluation ), M_value (evaluation.M_value), M_column( evaluation.M_column)
    {
        if (evaluation.M_quadrature != 0)
        {
            M_quadrature = new QuadratureRule (* (evaluation.M_quadrature) );
        }

    }

    //! Expression-based constructor
    template< typename Expression>
    explicit EvaluationVectorFromNonConstantMatrix (const ExpressionVectorFromNonConstantMatrix<Expression, SpaceDim, FieldDim>& expression)
        : M_evaluation( expression.expr() ), M_value(0) , M_column( expression.column() ), M_quadrature(0)
    {}

    //! Destructor
    ~EvaluationVectorFromNonConstantMatrix()
    {
        // Deleting quadrature rule object
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }

    }

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& iElement)
    {
        zero();

        M_evaluation.update ( iElement );

        // Loop on each quad point
        matrix_Type matrixDOF;

	for( UInt q(0); q < M_quadrature->nbQuadPt(); ++q )
	{
	    // Getting the small matrix from the inner expression
	    matrixDOF = M_evaluation.value_qij( q, 0, 0 );

	    // Extracting the wanted column
	    M_value[q] = matrixDOF.extractColumn( M_column );
	}
    }

    //! Re-initiliaze method
    void zero ( )
    {
      for( UInt q(0); q < M_quadrature->nbQuadPt(); ++q )
       {
	 for( UInt i(0); i < FieldDim; ++i )
	   {
	       M_value[ q ][ i ] = 0.0;
	   }
       }
    }

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "vector from a non constant matrix[" << FieldDim << "]";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* globalCFE)
    {
        M_evaluation.setGlobalCFE (globalCFE);
    }

    //! Do nothing setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        M_evaluation.setTestCFE (testCFE);
    }

    //! Do nothing setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_evaluation.setSolutionCFE (solutionCFE);
    }

    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
        M_quadrature = new QuadratureRule (qr);	
        M_value.resize ( qr.nbQuadPt() );

        M_evaluation.setQuadrature (qr);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        return M_value[q];
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& /*i*/) const
    {
        return M_value[q];
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_value[q];
    }

    //@}

private:

    // Storage
    EvaluationType        M_evaluation;
    vector_Type           M_value;
    UInt                  M_column;
    QuadratureRule*       M_quadrature;
};


template<typename EvaluationType,UInt SpaceDim , UInt FieldDim>
const flag_Type EvaluationVectorFromNonConstantMatrix<EvaluationType,SpaceDim, FieldDim>::S_globalUpdateFlag = EvaluationType::S_globalUpdateFlag;

template<typename EvaluationType,UInt SpaceDim , UInt FieldDim>
const flag_Type EvaluationVectorFromNonConstantMatrix<EvaluationType,SpaceDim, FieldDim>::S_testUpdateFlag = EvaluationType::S_testUpdateFlag;

template<typename EvaluationType,UInt SpaceDim , UInt FieldDim>
const flag_Type EvaluationVectorFromNonConstantMatrix<EvaluationType,SpaceDim, FieldDim>::S_solutionUpdateFlag = EvaluationType::S_solutionUpdateFlag;


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
