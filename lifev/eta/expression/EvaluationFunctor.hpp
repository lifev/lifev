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
     @brief This file contains the definition of the EvaluationFunctor class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUTATION_FUNCTOR_HPP
#define EVALUTATION_FUNCTOR_HPP


#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionFunctor.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for a generic functor with 1 argument
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing a functor with 1 arguement during the assembly

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.

  <b> Template requirement </b>

  There are two template arguments for this class: The functor and the ArgumentEvaluationType.

  <i> FunctorType </i> The functor has to be copiable, must implement the operator() const, has
  a type named ReturnType, the operator() must return a ReturnType.

  <i> Argument Evaluation </i> It has to be an Evaluation class.
 */
template <typename FunctorType, typename ArgumentEvaluationType>
class EvaluationFunctor1
{
public:

    //! @name Public Types
    //@{

    //! Type of the values returned by this class
    typedef typename FunctorType::return_Type return_Type;

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
    EvaluationFunctor1 (const EvaluationFunctor1<FunctorType, ArgumentEvaluationType>& eval)
        : M_functor (eval.M_functor),
          M_evaluation (eval.M_evaluation)
    {}

    //! Constructor from the corresponding expression
    template<typename Argument>
    explicit EvaluationFunctor1 (const ExpressionFunctor1<FunctorType, Argument>& expression)
        : M_functor (expression.functor() ),
          M_evaluation (expression.argument() )
    {}

    //! Destructor
    ~EvaluationFunctor1()
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
        out << "f( ";
        ArgumentEvaluationType::display (out);
        out << " )";
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
        return (*M_functor) (M_evaluation.value_q (q) );
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& i) const
    {
        return (*M_functor) (M_evaluation.value_qi (q, i) );
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& i, const UInt& j) const
    {
        return (*M_functor) (M_evaluation.value_qij (q, i, j) );
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationFunctor1();

    //@}

    // Internal storage
    std::shared_ptr<FunctorType> M_functor;
    ArgumentEvaluationType M_evaluation;
};


template<typename FunctorType, typename ArgumentEvaluationType>
const flag_Type EvaluationFunctor1<FunctorType, ArgumentEvaluationType>::S_globalUpdateFlag
    = ArgumentEvaluationType::S_globalUpdateFlag;

template<typename FunctorType, typename ArgumentEvaluationType>
const flag_Type EvaluationFunctor1<FunctorType, ArgumentEvaluationType>::S_testUpdateFlag
    = ArgumentEvaluationType::S_testUpdateFlag;

template<typename FunctorType, typename ArgumentEvaluationType>
const flag_Type EvaluationFunctor1<FunctorType, ArgumentEvaluationType>::S_solutionUpdateFlag
    = ArgumentEvaluationType::S_solutionUpdateFlag;



//! Evaluation for a generic functor with 2 arguments
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing a functor with 2 arguements during the assembly

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.

  <b> Template requirement </b>

  There are two template arguments for this class: The functor and the ArgumentEvaluationType.

  <i> FunctorType </i> The functor has to be copiable, must implement the operator() const, has
  a type named return_Type, the operator() must return a return_Type.

  <i> Argument1EvaluationType </i> It has to be an Evaluation class.

  <i> Argument2EvaluationType </i> It has to be an Evaluation class.

 */
template <typename FunctorType, typename Argument1EvaluationType, typename Argument2EvaluationType>
class EvaluationFunctor2
{
public:

    //! @name Public Types
    //@{

    //! Type returned by this class
    typedef typename FunctorType::return_Type return_Type;

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
    EvaluationFunctor2 (const EvaluationFunctor2<FunctorType, Argument1EvaluationType, Argument2EvaluationType>& eval)
        : M_functor (eval.M_functor),
          M_evaluation1 (eval.M_evaluation1),
          M_evaluation2 (eval.M_evaluation2)
    {}

    //! Constructor from the corresponding expression
    template<typename Argument1, typename Argument2>
    explicit EvaluationFunctor2 (const ExpressionFunctor2<FunctorType, Argument1, Argument2>& expression)
        : M_functor (expression.functor() ),
          M_evaluation1 (expression.argument1() ),
          M_evaluation2 (expression.argument2() ) {}

    //! Destructor
    ~EvaluationFunctor2()
    {}

    //@}


    //! @name Methods
    //@{

    //! Internal update
    void update (const UInt& iElement)
    {
        M_evaluation1.update (iElement);
        M_evaluation2.update (iElement);
    }

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "f( ";
        Argument1EvaluationType::display (out);
        out << " , ";
        Argument2EvaluationType::display (out);
        out << " )";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* globalCFE)
    {
        M_evaluation1.setGlobalCFE (globalCFE);
        M_evaluation2.setGlobalCFE (globalCFE);
    }

    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        M_evaluation1.setTestCFE (testCFE);
        M_evaluation2.setTestCFE (testCFE);
    }

    //! Setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_evaluation1.setSolutionCFE (solutionCFE);
        M_evaluation2.setSolutionCFE (solutionCFE);
    }

    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        M_evaluation1.setQuadrature (qr);
        M_evaluation2.setQuadrature (qr);
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        return (*M_functor) (M_evaluation1.value_q (q), M_evaluation2.value_q (q) );
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& i) const
    {
        return (*M_functor) (M_evaluation1.value_qi (q, i), M_evaluation2.value_qi (q, i) );
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& i, const UInt& j) const
    {
        return (*M_functor) (M_evaluation1.value_qij (q, i, j), M_evaluation2.value_qij (q, i, j) );
    }

    //@}

private:


    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationFunctor2();

    //@}

    // Internal storage
    std::shared_ptr<FunctorType> M_functor;
    Argument1EvaluationType M_evaluation1;
    Argument2EvaluationType M_evaluation2;
};


template<typename FunctorType, typename Argument1EvaluationType, typename Argument2EvaluationType>
const flag_Type EvaluationFunctor2<FunctorType, Argument1EvaluationType, Argument2EvaluationType>::S_globalUpdateFlag
    = Argument1EvaluationType::S_globalUpdateFlag | Argument2EvaluationType::S_globalUpdateFlag;

template<typename FunctorType, typename Argument1EvaluationType, typename Argument2EvaluationType>
const flag_Type EvaluationFunctor2<FunctorType, Argument1EvaluationType, Argument2EvaluationType>::S_testUpdateFlag
    = Argument1EvaluationType::S_testUpdateFlag | Argument2EvaluationType::S_testUpdateFlag;

template<typename FunctorType, typename Argument1EvaluationType, typename Argument2EvaluationType>
const flag_Type EvaluationFunctor2<FunctorType, Argument1EvaluationType, Argument2EvaluationType>::S_solutionUpdateFlag
    = Argument1EvaluationType::S_solutionUpdateFlag | Argument2EvaluationType::S_solutionUpdateFlag;

} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
