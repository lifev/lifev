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
     @brief This file contains the definition of the EvaluationDivJ class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef EVALUATION_DIV_J_HPP
#define EVALUATION_DIV_J_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionDivJ.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

template <UInt fieldDim, UInt spaceDim>
class EvaluationDivJ
{

private:

    //! @name Constructors, destructor
    //@{

    //! Empty constructor
    EvaluationDivJ();

    //! Copy constructor
    EvaluationDivJ (const EvaluationDivJ& provider);

    //! Destructor
    ~EvaluationDivJ();

    //@}

};

//! Evaluation of the basis function div(phi_j) in the case of a vectorial FE.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the divergence of the solution in the assembly.

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <UInt spaceDim>
class EvaluationDivJ<3, spaceDim>
{
public:

    //! @name Public Types
    //@{

    //! The type of the values returned by this class
    typedef Real return_Type;

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

    //! Empty constructor
    EvaluationDivJ() {}

    //! Copy constructor
    EvaluationDivJ (const EvaluationDivJ& provider)
        : M_valuesPtr (provider.M_valuesPtr)
    {}

    //! Expression-based constructor
    explicit EvaluationDivJ (const ExpressionDivJ& /*expression*/) {}

    //! Destructor
    ~EvaluationDivJ() {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& /*iElement*/) {}

    //! Display method
    static void display ( std::ostream& out = std::cout)
    {
        out << "div_j";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* /*globalCFE*/) {}

    //! Do nothing setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* /*testCFE*/) {}

    //! Setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_valuesPtr = & (solutionCFE->M_divergence);
    }

    //! Do nothing setter for the quadrature rule
    void setQuadrature (const QuadratureRule&) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value for a matrix
    const return_Type& value_qij (const UInt& q, const UInt& /*i*/, const UInt& j) const
    {
        return (*M_valuesPtr) [q][j];
    }

    //@}

private:

    //! Pointer to the data
    std::vector< std::vector < return_Type > > const* M_valuesPtr;

};


template<UInt spaceDim>
const flag_Type EvaluationDivJ<3, spaceDim>::S_globalUpdateFlag = ET_UPDATE_NONE;

template<UInt spaceDim>
const flag_Type EvaluationDivJ<3, spaceDim>::S_testUpdateFlag = ET_UPDATE_NONE;

template<UInt spaceDim>
const flag_Type EvaluationDivJ<3, spaceDim>::S_solutionUpdateFlag = ET_UPDATE_DIVERGENCE;


//! Evaluation of the basis function div(phi_j) in the case of a 2D vectorial FE.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the divergence of the solution in the assembly.

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <UInt spaceDim>
class EvaluationDivJ<2, spaceDim>
{
public:

    //! @name Public Types
    //@{

    //! The type of the values returned by this class
    typedef Real return_Type;

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

    //! Empty constructor
    EvaluationDivJ() {}

    //! Copy constructor
    EvaluationDivJ (const EvaluationDivJ& provider)
        : M_valuesPtr (provider.M_valuesPtr)
    {}

    //! Expression-based constructor
    explicit EvaluationDivJ (const ExpressionDivJ& /*expression*/) {}

    //! Destructor
    ~EvaluationDivJ() {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& /*iElement*/) {}

    //! Display method
    static void display ( std::ostream& out = std::cout)
    {
        out << "div_j";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* /*globalCFE*/) {}

    //! Do nothing setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* /*testCFE*/) {}

    //! Setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* solutionCFE)
    {
        M_valuesPtr = & (solutionCFE->M_divergence);
    }

    //! Do nothing setter for the quadrature rule
    void setQuadrature (const QuadratureRule&) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value for a matrix
    const return_Type& value_qij (const UInt& q, const UInt& /*i*/, const UInt& j) const
    {
        return (*M_valuesPtr) [q][j];
    }

    //@}

private:

    //! Pointer to the data
    std::vector< std::vector < return_Type > > const* M_valuesPtr;

};


template<UInt spaceDim>
const flag_Type EvaluationDivJ<2, spaceDim>::S_globalUpdateFlag = ET_UPDATE_NONE;

template<UInt spaceDim>
const flag_Type EvaluationDivJ<2, spaceDim>::S_testUpdateFlag = ET_UPDATE_NONE;

template<UInt spaceDim>
const flag_Type EvaluationDivJ<2, spaceDim>::S_solutionUpdateFlag = ET_UPDATE_DIVERGENCE;


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
