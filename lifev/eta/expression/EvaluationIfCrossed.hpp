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
     @brief This file contains the definition of the EvaluationInterpolateValue class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUATION_IF_CROSSED_HPP
#define EVALUATION_IF_CROSSED_HPP

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionIfCrossed.hpp>

#include <boost/shared_ptr.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for the interpolation of a FE function
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */
template<typename MeshType, typename MapType, UInt SpaceDim>
class EvaluationIfCrossed
{
public:

    //! @name Public Types
    //@{

    //! Type of the value returned by this class
    typedef Real return_Type;

    //! Type of the FESpace that has to be used with this class
    typedef ETFESpace<MeshType, MapType, SpaceDim, 1> fespace_Type;

    //! Pointer on the FESpace
    typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;

    //! Vector of the values
    typedef VectorEpetra vector_Type;

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
    EvaluationIfCrossed (const EvaluationIfCrossed<MeshType, MapType, SpaceDim>& evaluation)
        :
        M_fespace ( evaluation.M_fespace),
        M_vector ( evaluation.M_vector, Repeated),
        M_value ( evaluation.M_value)
    {}

    //! Expression-based constructor
    explicit EvaluationIfCrossed (const ExpressionIfCrossed<MeshType, MapType, SpaceDim>& expression)
        :
        M_fespace ( expression.fespace() ),
        M_vector ( expression.vector(), Repeated ),
        M_value (0.0)
    {}

    //! Destructor
    ~EvaluationIfCrossed()
    {}

    //@}


    //! @name Methods
    //@{

    //! Interal update, computes the interpolated values.
    void update (const UInt& iElement)
    {
        bool existPositive (false);
        bool existNegative (false);

        for (UInt i (0); i < 4; ++i)
        {
            UInt globalID (M_fespace->dof().localToGlobalMap (iElement, i) );

            Real interpolatedValue = M_vector[globalID];

            if (interpolatedValue >= 0)
            {
                existPositive = true;
            }
            else
            {
                existNegative = true;
            }
        }

        if (existPositive && existNegative)
        {
            M_value = 1.0;
        }
        else
        {
            M_value = 0.0;
        }
    }

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "ifCrossed";
    }

    //@}


    //! @name Set Methods
    //@{

    //! Do nothing setter for the global current FE
    template< typename CFEType >
    void setGlobalCFE (const CFEType* /*globalCFE*/)
    {}

    //! Do nothing setter for the test current FE
    template< typename CFEType >
    void setTestCFE (const CFEType* /*testCFE*/)
    {}

    //! Do nothing setter for the solution current FE
    template< typename CFEType >
    void setSolutionCFE (const CFEType* /*solutionCFE*/)
    {}

    //! Setter for the quadrature rule (deep copy)
    void setQuadrature (const QuadratureRule& )
    {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& /*q*/) const
    {
        return M_value;
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& /*q*/, const UInt& /*i*/) const
    {
        return M_value;
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& /*q*/, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_value;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationIfCrossed();

    //@}

    fespacePtr_Type M_fespace;
    vector_Type M_vector;

    Real M_value;
};

template<typename MeshType, typename MapType, UInt SpaceDim>
const flag_Type
EvaluationIfCrossed<MeshType, MapType, SpaceDim>::
S_globalUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim>
const flag_Type
EvaluationIfCrossed<MeshType, MapType, SpaceDim>::
S_testUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim>
const flag_Type
EvaluationIfCrossed<MeshType, MapType, SpaceDim>::
S_solutionUpdateFlag = ET_UPDATE_NONE;


} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
