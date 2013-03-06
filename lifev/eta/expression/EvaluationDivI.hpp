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
     @brief This file contains the definition of the EvaluationDivI class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUATION_DIV_I_HPP
#define EVALUATION_DIV_I_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/eta/expression/ExpressionDivI.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

template <UInt fieldDim, UInt spaceDim>
class EvaluationDivI
{
private:

    //! @name Constructors, destructor
    //@{

    //! Empty constructor
    EvaluationDivI();

    //! Copy constructor
    EvaluationDivI (const EvaluationDivI& provider);

    //! Destructor
    ~EvaluationDivI();

    //@}
};


//! Evaluation of the basis function divergence in the case of a vectorial FE.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the divergence of a test function in the assembly.

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <UInt spaceDim>
class EvaluationDivI<3, spaceDim>
{
public:

    //! @name Public Types
    //@{

    //! Type of the values returned by this class
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
    EvaluationDivI() {}

    //! Copy constructor
    EvaluationDivI (const EvaluationDivI<3, spaceDim>& provider) : M_valuesPtr (provider.M_valuesPtr) {}

    //! Expression-based constructor
    explicit EvaluationDivI (const ExpressionDivI& /*expression*/) {}

    //! Destructor
    ~EvaluationDivI() {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& /*iElement*/) {}

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "div_i";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* /*globalCFE*/) {}

    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        ASSERT (testCFE != 0, "Nul pointer to the testCFE cannot be set");
        M_valuesPtr = & (testCFE->M_divergence);
    }

    //! Do nothing setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* /*solutionCFE*/) {}

    //! Do nothing setter for the quadrature rule
    void setQuadrature (const QuadratureRule&) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value for a vector
    const return_Type& value_qi (const UInt& q, const UInt& i) const
    {
        ASSERT ( q < M_valuesPtr->size(), "Quadrature point index invalid");
        ASSERT ( i < (*M_valuesPtr) [q].size(), "Dof index invalid");
        return (*M_valuesPtr) [q][i];
    }

    //! Getter for the value for a matrix
    const return_Type& value_qij (const UInt& q, const UInt& i, const UInt& /*j*/) const
    {
        ASSERT ( q < M_valuesPtr->size(), "Quadrature point index invalid");
        ASSERT ( i < (*M_valuesPtr) [q].size(), "Dof index invalid");
        return (*M_valuesPtr) [q][i];
    }

    //@}

private:

    //! Storage of the pointer to the data
    std::vector< std::vector < return_Type > > const* M_valuesPtr;

};


template<UInt spaceDim>
const flag_Type EvaluationDivI<3, spaceDim>::S_globalUpdateFlag = ET_UPDATE_NONE;

template<UInt spaceDim>
const flag_Type EvaluationDivI<3, spaceDim>::S_testUpdateFlag = ET_UPDATE_DIVERGENCE;

template<UInt spaceDim>
const flag_Type EvaluationDivI<3, spaceDim>::S_solutionUpdateFlag = ET_UPDATE_NONE;


//! Evaluation of the basis function divergence in the case of a 2D vectorial FE.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the divergence of a test function in the assembly.

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <UInt spaceDim>
class EvaluationDivI<2, spaceDim>
{
public:

    //! @name Public Types
    //@{

    //! Type of the values returned by this class
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
    EvaluationDivI() {}

    //! Copy constructor
    EvaluationDivI (const EvaluationDivI<2, spaceDim>& provider) : M_valuesPtr (provider.M_valuesPtr) {}

    //! Expression-based constructor
    explicit EvaluationDivI (const ExpressionDivI& /*expression*/) {}

    //! Destructor
    ~EvaluationDivI() {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& /*iElement*/) {}

    //! Display method
    static void display (std::ostream& out = std::cout)
    {
        out << "div_i";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* /*globalCFE*/) {}

    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* testCFE)
    {
        ASSERT (testCFE != 0, "Nul pointer to the testCFE cannot be set");
        M_valuesPtr = & (testCFE->M_divergence);
    }

    //! Do nothing setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* /*solutionCFE*/) {}

    //! Do nothing setter for the quadrature rule
    void setQuadrature (const QuadratureRule&) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value for a vector
    const return_Type& value_qi (const UInt& q, const UInt& i) const
    {
        ASSERT ( q < M_valuesPtr->size(), "Quadrature point index invalid");
        ASSERT ( i < (*M_valuesPtr) [q].size(), "Dof index invalid");
        return (*M_valuesPtr) [q][i];
    }

    //! Getter for the value for a matrix
    const return_Type& value_qij (const UInt& q, const UInt& i, const UInt& /*j*/) const
    {
        ASSERT ( q < M_valuesPtr->size(), "Quadrature point index invalid");
        ASSERT ( i < (*M_valuesPtr) [q].size(), "Dof index invalid");
        return (*M_valuesPtr) [q][i];
    }

    //@}

private:

    //! Storage of the pointer to the data
    std::vector< std::vector < return_Type > > const* M_valuesPtr;

};


template<UInt spaceDim>
const flag_Type EvaluationDivI<2, spaceDim>::S_globalUpdateFlag = ET_UPDATE_NONE;

template<UInt spaceDim>
const flag_Type EvaluationDivI<2, spaceDim>::S_testUpdateFlag = ET_UPDATE_DIVERGENCE;

template<UInt spaceDim>
const flag_Type EvaluationDivI<2, spaceDim>::S_solutionUpdateFlag = ET_UPDATE_NONE;


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
