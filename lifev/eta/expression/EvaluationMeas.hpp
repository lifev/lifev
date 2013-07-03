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
    @brief This file contains the definition of the evaluation for the measure

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 03 Jan 2012
 */

#ifndef EVALUATION_MEASURE_HPP
#define EVALUATION_MEASURE_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionMeas.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! EvaluationMeas - Evaluation of the measure of the cell
/*!
    @author Samuel Quinodoz
 */
template< UInt spaceDim>
class EvaluationMeas
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
    EvaluationMeas() {}

    //! Copy constructor
    EvaluationMeas (const EvaluationMeas<spaceDim>& evaluation)
        : M_valuePtr (evaluation.M_valuePtr)
    {}

    //! Expression-based constructor
    explicit EvaluationMeas (const ExpressionMeas& /*expression*/) {}

    //! Destructor
    ~EvaluationMeas() {}

    //@}


    //! @name Methods
    //@{

    //! Do nothing update
    void update (const UInt& /*iElement*/) {}

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "meas_K";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* globalCFE)
    {
        ASSERT (globalCFE != 0, "Nul pointer to the globalCFE cannot be set");
        M_valuePtr = & (globalCFE->M_measure);

        std::cout << "ciao ciao:" << *M_valuePtr << std::endl;
        int n;
        std::cin >> n;
    }

    //! Setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* /*testCFE*/) {}

    //! Do nothing setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* /*solutionCFE*/) {}

    //! Do nothing setter for the quadrature rule
    void setQuadrature (const QuadratureRule&) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the value for a value
    return_Type value_q (const UInt& /*q*/) const
    {
        return *M_valuePtr;

        std::cout << "ciao ciao:" << *M_valuePtr << std::endl;
        int n;
        std::cin >> n;

    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& /*q*/, const UInt& /*i*/) const
    {
        return *M_valuePtr;

        std::cout << "ciao ciao:" << *M_valuePtr << std::endl;
        int n;
        std::cin >> n;

    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& /*q*/, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return *M_valuePtr;

        std::cout << "ciao ciao:" << *M_valuePtr << std::endl;
        int n;
        std::cin >> n;

    }

    //@}

private:

    //! Storage for the pointer to the data
    Real const* M_valuePtr;

};

template<UInt spaceDim>
const flag_Type EvaluationMeas<spaceDim>::S_globalUpdateFlag = ET_UPDATE_MEASURE;

template<UInt spaceDim>
const flag_Type EvaluationMeas<spaceDim>::S_testUpdateFlag = ET_UPDATE_NONE;

template<UInt spaceDim>
const flag_Type EvaluationMeas<spaceDim>::S_solutionUpdateFlag = ET_UPDATE_NONE;


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif /* EVALUATION_HK_HPP */
