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

#ifndef EVALUATION_PATCHAREA_HPP
#define EVALUATION_PATCHAREA_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionInterpolateValue.hpp>

#include <boost/shared_ptr.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! Evaluation for the interpolation of a FE function
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing an interpolated FE value.

  This is the generic implementation, so representing a vectorial FE

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
class EvaluationPatchArea
{
public:

    //! @name Public Types
    //@{

    //! Type of the value returned by this class
    typedef VectorSmall<FieldDim>   return_Type;
    typedef ETFESpace<MeshType, MapType, SpaceDim, FieldDim> ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type>                ETFESpacePtr_Type;
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
    EvaluationPatchArea (const EvaluationPatchArea& evaluation)
        :
        M_fespace ( evaluation.M_fespace),
        M_quadrature (0),
        //M_currentFE(M_fespace->refFE(),M_fespace->geoMap()),
        M_currentFE (evaluation.M_currentFE),
        M_patchArea (evaluation.M_patchArea)
    {
        if (evaluation.M_quadrature != 0)
        {
            M_quadrature = new QuadratureRule (* (evaluation.M_quadrature) );
            //M_currentFE.setQuadratureRule(*(evaluation.M_quadrature));
        }
    }

    //! Expression-based constructor
    explicit EvaluationPatchArea (const ExpressionPatchArea<MeshType, MapType, SpaceDim, FieldDim>& expression)
        :
        M_fespace ( expression.fespace() ),
        M_quadrature (0),
        M_currentFE (M_fespace->refFE(), M_fespace->geoMap() ),
        M_patchArea (0)
    {}

    //! Destructor
    ~EvaluationPatchArea()
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
    }

    //@}


    //! @name Methods
    //@{

    //! Interal update, computes the interpolated values.
    void update (const UInt& iElement)
    {
        zero();

        M_currentFE.update (M_fespace->mesh()->element (iElement), ET_UPDATE_PHI | ET_UPDATE_WDET );

        for (UInt q (0); q < M_quadrature->nbQuadPt(); ++q)
        {
            for (UInt iDim (0); iDim < FieldDim; ++iDim)
            {
                M_patchArea[q][iDim] += M_currentFE.measure();
            }
        }
    }


    //! Erase all the interpolated values stored internally
    void zero()
    {
        for (UInt q (0); q < M_quadrature->nbQuadPt(); ++q)
        {
            for (UInt iDim (0); iDim < FieldDim; ++iDim)
            {
                M_patchArea[q][iDim] = 0.0;
            }
        }
    }

    //! Show the values interpolated
    void showValues() const
    {
        std::cout << " Values : " << std::endl;

        for (UInt i (0); i < M_quadrature->nbQuadPt(); ++i)
        {
            std::cout << M_patchArea[i] << std::endl;
        }
    }

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "patchAread" << "\n";
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
    void setQuadrature (const QuadratureRule& qr)
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
        M_quadrature = new QuadratureRule (qr);
        M_currentFE.setQuadratureRule (qr);
        M_patchArea.resize (qr.nbQuadPt() );
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        return M_patchArea[q];
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& /*i*/) const
    {
        return M_patchArea[q];
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_patchArea[q];
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationPatchArea();

    //@}

    ETFESpacePtr_Type M_fespace;

    QuadratureRule* M_quadrature;
    ETCurrentFE<SpaceDim, 1> M_currentFE;

    std::vector<return_Type> M_patchArea;
};


template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
const flag_Type
EvaluationPatchArea<MeshType, MapType, SpaceDim, FieldDim>::
S_globalUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
const flag_Type
EvaluationPatchArea<MeshType, MapType, SpaceDim, FieldDim>::
S_testUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
const flag_Type
EvaluationPatchArea<MeshType, MapType, SpaceDim, FieldDim>::
S_solutionUpdateFlag = ET_UPDATE_NONE;

} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
