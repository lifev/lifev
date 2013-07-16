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
    @brief File for the definition of the expression used for a general functor.

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */

#ifndef EXPRESSION_SAVEVECTORACCORDINGTOFUNCTOR_HPP
#define EXPRESSION_SAVEVECTORACCORDINGTOFUNCTOR_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/expression/ExpressionBase.hpp>
#include <lifev/eta/expression/ExpressionScalar.hpp>
#include <lifev/eta/expression/ExpressionVector.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>

namespace LifeV
{

namespace ExpressionAssembly
{

//! class ExpressionSaveVectorAccordingToFunctor1  Class representing a functor with 1 expression as arguement.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This is an experimental expression.

  This expression is usefull to build expressions that are "function-like":
  In case one wants to have an expression that is a function, the functor
  passed to this expression will have no members and its operator() will
  operate the function evaluation.

  A functor, however, is more general in the sense that it can contain
  data. These data can be very simple (e.g. just a Real to represent
  the current time for a function that would depend on time as well)
  or arbitrarily complicated.

  It is required to pass the functor as a shared pointer in order to
  avoid the copy of heavy data that could be stored in the functor.

  <b> Template parameters </b>

  <i>FunctorType</i>: The type of the functor

  <i>ArgumentType</i>: The type of the argument, that is an expression
  (it is usually different from the input type expected by the functor)

  <b> Template requirements </b>

  <i>FunctorType</i>: None

  <i>ArgumentType</i>: Copiable, has a static method for display
*/

template< typename FunctorType, typename ArgumentType >
class ExpressionSaveVectorAccordingToFunctor1 : public ExpressionBase<ExpressionSaveVectorAccordingToFunctor1<FunctorType, ArgumentType> >
{
public:

    //! @name Public Types
    //@{

    // Base class, typedef usefull to make the code cleaner
    typedef ExpressionBase<ExpressionSaveVectorAccordingToFunctor1<FunctorType, ArgumentType> > base_Type;
    typedef VectorEpetra                                                                        vector_Type;
    typedef boost::Shared_ptr<vector_Type>                                                      vectorPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Argument with the required data in argument
    ExpressionSaveVectorAccordingToFunctor1 (boost::shared_ptr<FunctorType> fct,
                                             const ArgumentType& arg,
                                             const vectorPtr_Type originVector,
                                             const vectorPtr_Type saveVector)
        : base_Type(),
          M_functor (fct),
          M_argument (arg),
          M_originVector (originVector)
          M_saveVector (saveVector)
    {}

    //! Copy constructor
    ExpressionSaveVectorAccordingToFunctor1 (const ExpressionSaveVectorAccordingToFunctor1<FunctorType, ArgumentType>& expr)
        : base_Type(),
          M_functor (expr.M_functor),
          M_argument (expr.M_argument),
          M_originVector (expr.M_originVector)
          M_saveVector (expr.M_saveVector)
    {}

    //! Destructor
    ~ExpressionSaveVectorAccordingToFunctor1() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "fct( ";
        ArgumentType::display (out);
        out << " )";
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the functor
    boost::shared_ptr<FunctorType> functor() const
    {
        return M_functor;
    }

    //! Getter for the expression to be placed as argument
    const ArgumentType& argument() const
    {
        return M_argument;
    }

    //! Getter for the origin vector
    const vectorPtr_Type originVector() const
    {
        return M_originVector;
    }

    //! Getter for the save vector
    const vectorPtr_Type saveVector() const
    {
        return M_saveVector;
    }


    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    ExpressionSaveVectorAccordingToFunctor1();

    //@}

    // Storage for the functor
    boost::shared_ptr<FunctorType> M_functor;

    // Storage for the argument
    ArgumentType M_argument;

    vectorPtr_Type M_originVector;
    vectorPtr_Type M_saveVector;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This helper function builds the expression associated with the
  action of a functor on an expression.

  <b> Template parameters </b>

  <i>FunctorType</i>: The type of the functor

  <i>ArgumentType</i>: The type of the argument, that is an expression
  (it is usually different from the input type expected by the functor)

  <b> Template requirements </b>

  <i>FunctorType</i>: Same as LifeV::ExpressionSaveVectorAccordingToFunctor1

  <i>ArgumentType</i>: Same as LifeV::ExpressionSaveVectorAccordingToFunctor1
*/

template< typename FunctorType, typename ArgumentType>
inline ExpressionSaveVectorAccordingToFunctor1<FunctorType, ArgumentType>
eval (boost::shared_ptr<FunctorType> fct,
      const ExpressionBase<ArgumentType>& argument,
      const boost::shared_ptr<VectorEpetra> originVector,
      const boost::shared_ptr<VectorEpetra> saveVector)
{
    return ExpressionSaveVectorAccordingToFunctor1<FunctorType, ArgumentType> (fct, argument.cast(), originVector, saveVector );
}

///
///
///  2 ARGUMENTS
///
///

//! class ExpressionSaveVectorAccordingToFunctor2  Class representing a functor with 2 expressions as arguement.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  See the documentation of the class LifeV::ExpressionSaveVectorAccordingToFunctor1

  <b> Template parameters </b>

  <i>FunctorType</i>: The type of the functor

  <i>ArgumentType1</i>: The type of the first argument, that is an expression
  (it is usually different from the input type expected by the functor)

  <i>ArgumentType2</i>: The type of the second argument, that is an expression
  (it is usually different from the input type expected by the functor)

  <b> Template requirements </b>

  <i>FunctorType</i>: None

  <i>ArgumentType1</i>: Copiable, has a static method for display

  <i>ArgumentType2</i>: Copiable, has a static method for display
*/

template< typename FunctorType, typename ArgumentType1, typename ArgumentType2 >
class ExpressionSaveVectorAccordingToFunctor2 : public ExpressionBase<ExpressionSaveVectorAccordingToFunctor2<FunctorType, ArgumentType1, ArgumentType2> >
{
public:

    //! @name Public Types
    //@{

    // Base class, typedef usefull to make the code cleaner
    typedef ExpressionBase<ExpressionSaveVectorAccordingToFunctor2<FunctorType, ArgumentType1, ArgumentType2> > base_Type;
    typedef VectorEpetra                                                                        vector_Type;
    typedef boost::Shared_ptr<vector_Type>                                                      vectorPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Argument with the required data in argument
    ExpressionSaveVectorAccordingToFunctor2 (boost::shared_ptr<FunctorType> fct, const ArgumentType1& arg1, const ArgumentType2& arg2,
                                             const vectorPtr_Type originVector, const vectorPtr_Type saveVector)
        : base_Type(),
          M_functor (fct),
          M_argument1 (arg1),
          M_argument2 (arg2),
          M_originVector (originVector),
          M_saveVector (saveVector)
    {}

    //! Copy constructor
    ExpressionSaveVectorAccordingToFunctor2 (const ExpressionSaveVectorAccordingToFunctor2<FunctorType, ArgumentType1, ArgumentType2>& expr)
        : base_Type(),
          M_functor (expr.M_functor),
          M_argument1 (expr.M_argument1),
          M_argument2 (expr.M_argument2),
          M_originVector (expr.M_originVector),
          M_saveVector (expr.M_saveVector)
    {}

    //! Destructor
    ~ExpressionSaveVectorAccordingToFunctor2() {}

    //@}


    //! @name Methods
    //@{

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "fct( ";
        ArgumentType1::display (out);
        out << " , " << ArgumentType2::display (out);
        out << " )";
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the functor
    boost::shared_ptr<FunctorType> functor() const
    {
        return M_functor;
    }

    //! Getter for the expression to be placed as first argument
    const ArgumentType1& argument1() const
    {
        return M_argument1;
    }

    //! Getter for the expression to be placed as second argument
    const ArgumentType2& argument2() const
    {
        return M_argument2;
    }

    const vectorPtr_Type originVector() const
    {
        return M_originVector;
    }

    const vectorPtr_Type saveVector() const
    {
        return M_saveVector;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    ExpressionSaveVectorAccordingToFunctor2();

    //@}

    // Storage for the functor
    boost::shared_ptr<FunctorType> M_functor;

    // Storage for the arguments
    ArgumentType1 M_argument1;
    ArgumentType2 M_argument2;

    vectorPtr_Type M_originVector;
    vectorPtr_Type M_saveVector;
};

//! Simple function to be used in the construction of an expression
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This helper function builds the expression associated with the
  action of a functor on an expression.

  <b> Template parameters </b>

  <i>FunctorType</i>: The type of the functor

  <i>ArgumentType1</i>: The type of the first argument, that is an expression
  (it is usually different from the input type expected by the functor)

  <i>ArgumentType2</i>: The type of the second argument, that is an expression
  (it is usually different from the input type expected by the functor)

  <b> Template requirements </b>

  <i>FunctorType</i>: Same as LifeV::ExpressionSaveVectorAccordingToFunctor2

  <i>ArgumentType1</i>: Same as LifeV::ExpressionSaveVectorAccordingToFunctor2

  <i>ArgumentType2</i>: Same as LifeV::ExpressionSaveVectorAccordingToFunctor2
*/

template< typename FunctorType, typename ArgumentType1, typename ArgumentType2>
inline ExpressionSaveVectorAccordingToFunctor2<FunctorType, ArgumentType1, ArgumentType2>
eval (boost::shared_ptr<FunctorType> fct, const ArgumentType1& arg1, const ArgumentType2& arg2,
      const boost::shared_ptr<VectorEpetra> originVector,  const boost::shared_ptr<VectorEpetra> saveVector)
{
    return ExpressionSaveVectorAccordingToFunctor2<FunctorType, ArgumentType1, ArgumentType2> (fct, arg1, arg2,
                                                                                               originVector, saveVector);
}

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
