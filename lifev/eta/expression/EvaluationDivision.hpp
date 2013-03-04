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
     @brief This file contains the definition of the EvaluationDivision class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef EVALUTATION_DIVISION_HPP
#define EVALUTATION_DIVISION_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/array/OperationSmallDivision.hpp>

#include <lifev/eta/expression/ExpressionDivision.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for the quotient of two other Evaluations
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing a division operation during the assembly

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <typename EvaluationLType, typename EvaluationRType>
class EvaluationDivision
{
public:

    //! @name Public Types
    //@{

    //! Type of the value returned by the left operand
    typedef typename EvaluationLType::return_Type Lreturn_Type;

    //! Type of the value returned by the right operand
    typedef typename EvaluationRType::return_Type Rreturn_Type;

    //! Type of the value returned by this class
    typedef typename OperationSmallDivision<Lreturn_Type, Rreturn_Type>::result_Type return_Type;

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
    EvaluationDivision (const EvaluationDivision& eval)
        : M_evaluationL (eval.M_evaluationL),
          M_evaluationR (eval.M_evaluationR)
    {}

    //! Constructor from the corresponding expression
    template <typename L, typename R>
    explicit EvaluationDivision (const ExpressionDivision<L, R>& expression)
        : M_evaluationL (expression.left() ),
          M_evaluationR (expression.right() )
    {}

    //! Destructor
    ~EvaluationDivision() {}

    //@}


    //! @name Methods
    //@{

    //! Internal update method
    void update (const UInt& iElement)
    {
        M_evaluationL.update (iElement);
        M_evaluationR.update (iElement);
    }

    //! Display method
    static void display ( std::ostream& out = std::cout)
    {
        EvaluationLType::display (out);
        out << " / ";
        EvaluationRType::display (out);
    }

    //@}


    //! @name Set Methods
    //@{

    //! Setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* globalCFE)
    {
        M_evaluationL.setGlobalCFE (globalCFE);
        M_evaluationR.setGlobalCFE (globalCFE);
    }

    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        M_evaluationL.setTestCFE (testCFE);
        M_evaluationR.setTestCFE (testCFE);
    }

    //! Setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_evaluationL.setSolutionCFE (solutionCFE);
        M_evaluationR.setSolutionCFE (solutionCFE);
    }

    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        M_evaluationL.setQuadrature (qr);
        M_evaluationR.setQuadrature (qr);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        return M_evaluationL.value_q (q) / M_evaluationR.value_q (q);
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& i) const
    {
        return M_evaluationL.value_qi (q, i) / M_evaluationR.value_qi (q, i);
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& i, const UInt& j) const
    {
        return M_evaluationL.value_qij (q, i, j) / M_evaluationR.value_qij (q, i, j);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationDivision();

    //@}

    // Internal storage
    EvaluationLType M_evaluationL;
    EvaluationRType M_evaluationR;
};


template< typename EvaluationLType, typename EvaluationRType>
const flag_Type EvaluationDivision<EvaluationLType, EvaluationRType>::S_globalUpdateFlag
    = EvaluationLType::S_globalUpdateFlag | EvaluationRType::S_globalUpdateFlag;

template< typename EvaluationLType, typename EvaluationRType>
const flag_Type EvaluationDivision<EvaluationLType, EvaluationRType>::S_testUpdateFlag
    = EvaluationLType::S_testUpdateFlag | EvaluationRType::S_testUpdateFlag;

template< typename EvaluationLType, typename EvaluationRType>
const flag_Type EvaluationDivision<EvaluationLType, EvaluationRType>::S_solutionUpdateFlag
    = EvaluationLType::S_solutionUpdateFlag | EvaluationRType::S_solutionUpdateFlag;


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
