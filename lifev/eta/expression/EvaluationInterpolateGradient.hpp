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
     @brief This file contains the definition of the EvaluationInterpolateGradient class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef EVALUATION_INTERPOLATE_GRADIENT_HPP
#define EVALUATION_INTERPOLATE_GRADIENT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorSmall.hpp>

#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentFlag.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <lifev/eta/expression/ExpressionInterpolateGradient.hpp>

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

  !! NOT DEFINED YET !! (Reason: miss SimpleTensor)
 */
template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
class EvaluationInterpolateGradient
{
public:

    //! @name Public Types
    //@{

    //! Type of the value returned by this class
    typedef MatrixSmall<FieldDim, SpaceDim> return_Type;

    //! Type of the FESpace to be used in this class
    typedef ETFESpace<MeshType, MapType, SpaceDim, FieldDim> fespace_Type;

    //! Type of the pointer on the FESpace
    typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;

    //! Type of the vector to be used
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
    EvaluationInterpolateGradient (const EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, FieldDim>& evaluation)
        :
        M_fespace ( evaluation.M_fespace),
        M_vector ( evaluation.M_vector, Repeated),
        M_quadrature (0),
        M_currentFE (evaluation.M_currentFE),
        M_interpolatedGradients (evaluation.M_interpolatedGradients)
    {
        if (evaluation.M_quadrature != 0)
        {
            M_quadrature = new QuadratureRule (* (evaluation.M_quadrature) );
        }
    }

    //! Expression-based constructor
    explicit EvaluationInterpolateGradient (const ExpressionInterpolateGradient<MeshType, MapType, SpaceDim, FieldDim>& expression)
        :
        M_fespace ( expression.fespace() ),
        M_vector ( expression.vector(), Repeated ),
        M_quadrature (0),
        M_currentFE (M_fespace->refFE(), M_fespace->geoMap() ),
        M_interpolatedGradients (0)
    {}

    //! Destructor
    ~EvaluationInterpolateGradient()
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
    }

    //@}


    //! @name Methods
    //@{
    //! Internal update: computes the interpolated gradients
    void update (const UInt& iElement)
    {
        zero();

        M_currentFE.update (M_fespace->mesh()->element (iElement), ET_UPDATE_DPHI);

        for (UInt i (0); i < M_fespace->refFE().nbDof(); ++i)
        {
            for (UInt q (0); q < M_quadrature->nbQuadPt(); ++q)
            {
                for (UInt iDim (0); iDim < SpaceDim; ++iDim)
                {
                    for (UInt jDim (0); jDim < FieldDim; ++jDim)
                    {
                        UInt globalID (M_fespace->dof().localToGlobalMap (iElement, i)
                                       + jDim * M_fespace->dof().numTotalDof() );

                        M_interpolatedGradients[q][jDim][iDim] +=
                            M_currentFE.dphi (jDim * M_currentFE.nbFEDof() + i, jDim, iDim, q)
                            * M_vector[globalID];
                    }
                }
            }
        }
    }

    //! Erase the interpolated gradients stored internally
    void zero()
    {
        for (UInt q (0); q < M_quadrature->nbQuadPt(); ++q)
        {
            for (UInt i (0); i < SpaceDim; ++i)
            {
                for (UInt j (0); j < FieldDim; ++j)
                {
                    M_interpolatedGradients[q][j][i] = 0.0;
                }
            }
        }
    }

    //! Show the values
    void showValues() const
    {
        std::cout << " Gradients : " << std::endl;

        for (UInt i (0); i < M_quadrature->nbQuadPt(); ++i)
        {
            std::cout << M_interpolatedGradients[i] << std::endl;
        }
    }

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "interpolated[ " << FieldDim << " ][ " << SpaceDim << " ]";
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

    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
        M_quadrature = new QuadratureRule (qr);
        M_currentFE.setQuadratureRule (qr);
        M_interpolatedGradients.resize (qr.nbQuadPt() );
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        return M_interpolatedGradients[q];
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& /*i*/) const
    {
        return M_interpolatedGradients[q];
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_interpolatedGradients[q];
    }

    //@}


private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationInterpolateGradient();

    //@}

    //! Data storage
    fespacePtr_Type M_fespace;
    vector_Type M_vector;
    QuadratureRule* M_quadrature;

    //! Structure for the computations
    ETCurrentFE<SpaceDim, FieldDim> M_currentFE;

    //! Storage for the temporary values
    std::vector<return_Type> M_interpolatedGradients;
};

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
const flag_Type
EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, FieldDim>::
S_globalUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
const flag_Type
EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, FieldDim>::
S_testUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
const flag_Type
EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, FieldDim>::
S_solutionUpdateFlag = ET_UPDATE_NONE;


//! Evaluation for the interpolation of the gradient of a FE function
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class aims at representing an interpolated FE gradient.

  This is the specialized (partially) implementation representing a scalar FE

  This class is an Evaluation class, and therefore, has all the methods
  required to work within the Evaluation trees.
 */
template<typename MeshType, typename MapType, UInt SpaceDim>
class EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, 1>
{

public:

    //! @name Public Types
    //@{

    //! Type of the value returned by this class
    typedef VectorSmall<SpaceDim> return_Type;

    //! Type of the FESpace to be used in this class
    typedef ETFESpace<MeshType, MapType, SpaceDim, 1> fespace_Type;

    //! Type of the pointer on the FESpace
    typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;

    //! Type of the vector to be used
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
    EvaluationInterpolateGradient (const EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, 1>& evaluation)
        :
        M_fespace ( evaluation.M_fespace),
        M_vector ( evaluation.M_vector, Repeated),
        M_quadrature (0),
        M_currentFE (evaluation.M_currentFE),
        M_interpolatedGradients (evaluation.M_interpolatedGradients)
    {
        if (evaluation.M_quadrature != 0)
        {
            M_quadrature = new QuadratureRule (* (evaluation.M_quadrature) );
        }
    }

    //! Expression-based constructor
    explicit EvaluationInterpolateGradient (const ExpressionInterpolateGradient<MeshType, MapType, SpaceDim, 1>& expression)
        :
        M_fespace ( expression.fespace() ),
        M_vector ( expression.vector(), Repeated ),
        M_quadrature (0),
        M_currentFE (M_fespace->refFE(), M_fespace->geoMap() ),
        M_interpolatedGradients (0)
    {}

    //! Destructor
    ~EvaluationInterpolateGradient()
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
    }

    //@}


    //! @name Methods
    //@{

    //! Internal update: computes the interpolated gradients
    void update (const UInt& iElement)
    {
        zero();

        M_currentFE.update (M_fespace->mesh()->element (iElement), ET_UPDATE_DPHI);

        for (UInt i (0); i < M_fespace->refFE().nbDof(); ++i)
        {
            for (UInt q (0); q < M_quadrature->nbQuadPt(); ++q)
            {
                UInt globalID (M_fespace->dof().localToGlobalMap (iElement, i) );

                for (UInt iDim (0); iDim < SpaceDim; ++iDim)
                {
                    M_interpolatedGradients[q][iDim] +=
                        M_currentFE.dphi (i, iDim, q)
                        * M_vector[globalID];
                }
            }
        }
    }

    //! Erase the interpolated gradients stored internally
    void zero()
    {
        for (UInt q (0); q < M_quadrature->nbQuadPt(); ++q)
        {
            for (UInt i (0); i < SpaceDim; ++i)
            {
                M_interpolatedGradients[q][i] = 0.0;
            }
        }
    }

    //! Show the values
    void showValues() const
    {
        std::cout << " Gradients : " << std::endl;

        for (UInt i (0); i < M_quadrature->nbQuadPt(); ++i)
        {
            std::cout << M_interpolatedGradients[i] << std::endl;
        }
    }

    //! Display method
    static void display (ostream& out = std::cout)
    {
        out << "interpolated[ " << SpaceDim << " ]";
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

    //! Setter for the quadrature rule
    void setQuadrature (const QuadratureRule& qr)
    {
        if (M_quadrature != 0)
        {
            delete M_quadrature;
        }
        M_quadrature = new QuadratureRule (qr);
        M_currentFE.setQuadratureRule (qr);
        M_interpolatedGradients.resize (qr.nbQuadPt() );
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for a value
    return_Type value_q (const UInt& q) const
    {
        return M_interpolatedGradients[q];
    }

    //! Getter for the value for a vector
    return_Type value_qi (const UInt& q, const UInt& /*i*/) const
    {
        return M_interpolatedGradients[q];
    }

    //! Getter for the value for a matrix
    return_Type value_qij (const UInt& q, const UInt& /*i*/, const UInt& /*j*/) const
    {
        return M_interpolatedGradients[q];
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    EvaluationInterpolateGradient();

    //@}

    //! Data storage
    fespacePtr_Type M_fespace;
    vector_Type M_vector;
    QuadratureRule* M_quadrature;

    //! Structure for the computations
    ETCurrentFE<SpaceDim, 1> M_currentFE;

    //! Storage for the temporary values
    std::vector<return_Type> M_interpolatedGradients;
};


template<typename MeshType, typename MapType, UInt SpaceDim>
const flag_Type
EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, 1>::
S_globalUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim>
const flag_Type
EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, 1>::
S_testUpdateFlag = ET_UPDATE_NONE;

template<typename MeshType, typename MapType, UInt SpaceDim>
const flag_Type
EvaluationInterpolateGradient<MeshType, MapType, SpaceDim, 1>::
S_solutionUpdateFlag = ET_UPDATE_NONE;

} // Namespace ExpressionAssembly

} // Namespace LifeV
#endif
