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
     @brief This file contains the definition of the EvaluationExtract1 class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUTATION_EXTRACT1_HPP
#define EVALUTATION_EXTRACT1_HPP


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/eta/expression/ExpressionExtract1.hpp>
#include <lifev/eta/expression/EvaluationMatrix.hpp>
#include <lifev/eta/array/OperationSmallExtract.hpp>


#include <lifev/core/fem/QuadratureRule.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for the extraction of a row resp. component from a matrix resp. a vector.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the extraction of a row resp. component from a matrix resp. a vector in the assembly.

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.

  <b> Template requirement </b>

  There is one template argument for this class: The EvaluationType.

  <i> EvaluationType </i> It has to be an Evaluation class.
 */
template <typename EvaluationType>
class EvaluationExtract1
{
public:

    //! @name Public Types
    //@{

    //! Type of the value returned by this class
    typedef typename OperationSmallExtract1< typename EvaluationType::return_Type >::result_Type return_Type;

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
    EvaluationExtract1 (const EvaluationExtract1< EvaluationType >& eval)
        : M_i (eval.M_i),
          M_evaluation (eval.M_evaluation)
    {}

    //! Constructor from the corresponding expression
    template< typename ExpressionType >
    explicit EvaluationExtract1 (const ExpressionExtract1< ExpressionType >& expression)
        : M_i (expression.indexI() ),
          M_evaluation (expression.exprEx() )
    {}


    //! Destructor
    ~EvaluationExtract1()
    {}

    //@}


    //! @name Methods
    //@{

    //! Internal update
    void update (const UInt& iElement)
    {
        M_evaluation.update (iElement);
    }

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "Extraction from ";
        EvaluationType::display (out);
    }

    //@}


    //! @name Set Methods
    //@{

    //! Setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* globalCFE)
    {
        M_evaluation.setGlobalCFE (globalCFE);
    }

    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        M_evaluation.setTestCFE (testCFE);
    }

    //! Setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_evaluation.setSolutionCFE (solutionCFE);
    }

    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        M_evaluation.setQuadrature (qr);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        return M_evaluation.value_q (q).extract (M_i);;
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& i) const
    {
        return M_evaluation.value_qi (q, i).extract (M_i);
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& i, const UInt& j) const
    {
        return M_evaluation.value_qij (q, i, j).extract (M_i);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationExtract1();

    //@}

    // Internal storage
    UInt M_i;
    EvaluationType M_evaluation;
};



template<typename EvaluationType>
const flag_Type EvaluationExtract1<EvaluationType>::S_globalUpdateFlag
    = EvaluationType::S_globalUpdateFlag;

template<typename EvaluationType>
const flag_Type EvaluationExtract1<EvaluationType>::S_testUpdateFlag
    = EvaluationType::S_testUpdateFlag;

template<typename EvaluationType>
const flag_Type EvaluationExtract1<EvaluationType>::S_solutionUpdateFlag
    = EvaluationType::S_solutionUpdateFlag;


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
