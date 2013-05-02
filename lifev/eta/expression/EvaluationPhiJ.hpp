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
     @brief This file contains the definition of the EvaluationPhiJ class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef EVALUATION_PHI_J_HPP
#define EVALUATION_PHI_J_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionPhiJ.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation of the basis function phi_j in the case of a scalar FE.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the solution in the assembly.

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template<UInt solutionDim>
class EvaluationPhiJ
{
public:

    //! @name Public Types
    //@{

    //! Type of the values returned by this class
    typedef VectorSmall<solutionDim> return_Type;

    //@}


    //! @name Static constants
    //@{

    //! Flag for the global current FE
    const static flag_Type S_globalUpdateFlag;

    //! Flag for the test current FE
    const static flag_Type S_testUpdateFlag;

    //! Flag fot the solution current FE
    const static flag_Type S_solutionUpdateFlag;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Empty constructor
    EvaluationPhiJ() {}

    //! Copy constructor
    EvaluationPhiJ (const EvaluationPhiJ& provider)
        : M_valuesPtr (provider.M_valuesPtr)
    {}

    //! Expression-based constructor
    explicit EvaluationPhiJ (const ExpressionPhiJ& /*expression*/) {}

    //! Destructor
    ~EvaluationPhiJ() {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& /*iElement*/) {}

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "phi_j";
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
        ASSERT (solutionCFE != 0, "Nul pointer to the solutionCFE cannot be set");
        M_valuesPtr = & (solutionCFE->M_phi);
    }

    //! Do nothing setter for the quadrature
    void setQuadrature (const QuadratureRule&) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value for a matrix
    const return_Type& value_qij (const UInt& q, const UInt& /*i*/, const UInt& j) const
    {
        ASSERT ( q < M_valuesPtr->size(), "Quadrature point index invalid");
        ASSERT ( j < (*M_valuesPtr) [q].size(), "Dof index invalid");
        return (*M_valuesPtr) [q][j];
    }

    //@}

private:

    //! Storage for the pointer to the data
    std::vector< std::vector < return_Type > > const* M_valuesPtr;

};


template<UInt solutionDim>
const flag_Type EvaluationPhiJ<solutionDim>::S_globalUpdateFlag = ET_UPDATE_NONE;

template<UInt solutionDim>
const flag_Type EvaluationPhiJ<solutionDim>::S_testUpdateFlag = ET_UPDATE_NONE;

template<UInt solutionDim>
const flag_Type EvaluationPhiJ<solutionDim>::S_solutionUpdateFlag = ET_UPDATE_PHI;


//! Evaluation of the basis function phi_j in the case of a scalar FE.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing the solution in the assembly.

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template <>
class EvaluationPhiJ<1>
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

    //! Flag fot the solution current FE
    const static flag_Type S_solutionUpdateFlag;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Empty constructor
    EvaluationPhiJ() {}

    //! Copy constructor
    EvaluationPhiJ (const EvaluationPhiJ& provider)
        : M_valuesPtr (provider.M_valuesPtr)
    {}

    //! Expression-based constructor
    explicit EvaluationPhiJ (const ExpressionPhiJ& /*expression*/) {}

    //! Destructor
    ~EvaluationPhiJ() {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing internal update
    void update (const UInt& /*iElement*/) {}

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "phi_j";
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
        ASSERT (solutionCFE != 0, "Nul pointer to the solutionCFE cannot be set");
        M_valuesPtr = & (solutionCFE->M_phi);
    }

    //! Do nothing setter for the quadrature
    void setQuadrature (const QuadratureRule&) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value for a matrix
    const return_Type& value_qij (const UInt& q, const UInt& /*i*/, const UInt& j) const
    {
        ASSERT ( q < M_valuesPtr->size(), "Quadrature point index invalid");
        ASSERT ( j < (*M_valuesPtr) [q].size(), "Dof index invalid");
        return (*M_valuesPtr) [q][j];
    }

    //@}

private:

    //! Storage for the pointer to the data
    std::vector< std::vector < Real > > const* M_valuesPtr;

};


const flag_Type EvaluationPhiJ<1>::S_globalUpdateFlag = ET_UPDATE_NONE;

const flag_Type EvaluationPhiJ<1>::S_testUpdateFlag = ET_UPDATE_NONE;

const flag_Type EvaluationPhiJ<1>::S_solutionUpdateFlag = ET_UPDATE_PHI;


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
