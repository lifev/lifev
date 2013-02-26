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
     @brief This file contains the definition of the IntegrateValueElement class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_VALUE_ELEMENT_HPP
#define INTEGRATE_VALUE_ELEMENT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <boost/shared_ptr.hpp>



namespace LifeV
{

namespace ExpressionAssembly
{

//! The class to actually perform the loop over the elements to compute a value
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is used to store the data required for the integration of a value and
  perform that assembly with a loop over the elements, and then, for each elements,
  using the Evaluation corresponding to the Expression (This convertion is done
  within a typedef).
 */
template < typename MeshType, typename ExpressionType>
class IntegrateValueElement
{
public:

    //! @name Public Types
    //@{

    //! Evaluation type
    typedef typename ExpressionToEvaluation < ExpressionType,
            0,
            0,
            3 >::evaluation_Type evaluation_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateValueElement (const boost::shared_ptr<MeshType>& mesh,
                           const QuadratureRule& quadrature,
                           const ExpressionType& expression);

    //! Copy constructor
    IntegrateValueElement ( const IntegrateValueElement < MeshType, ExpressionType>& integrator);

    //! Destructor
    ~IntegrateValueElement();

    //@}


    //! @name Operators
    //@{

    //! Operator wrapping the addTo method
    inline void operator>> (Real& value)
    {
        addTo (value);
    }

    //@}


    //! @name Methods
    //@{

    //! Ouput method
    void check (std::ostream& out = std::cout);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, sum over the quadrature nodes,
      sum into the global value.
     */
    void addTo (Real& value);

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    IntegrateValueElement();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Quadrature to be used
    QuadratureRule M_quadrature;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    ETCurrentFE<3, 1>* M_globalCFE;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename ExpressionType>
IntegrateValueElement < MeshType, ExpressionType>::
IntegrateValueElement (const boost::shared_ptr<MeshType>& mesh,
                       const QuadratureRule& quadrature,
                       const ExpressionType& expression)
    :   M_mesh (mesh),
        M_quadrature (quadrature),
        M_evaluation (expression),

        M_globalCFE (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), quadrature) )

{
    M_evaluation.setQuadrature (quadrature);
    M_evaluation.setGlobalCFE (M_globalCFE);
}


template < typename MeshType, typename ExpressionType>
IntegrateValueElement < MeshType, ExpressionType>::
IntegrateValueElement ( const IntegrateValueElement < MeshType, ExpressionType>& integrator)
    :   M_mesh (integrator.M_mesh),
        M_quadrature (integrator.M_quadrature),
        M_evaluation (integrator.M_evaluation),

        M_globalCFE (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_quadrature) )

{
    M_evaluation.setQuadrature (M_quadrature);
    M_evaluation.setGlobalCFE (M_globalCFE);
}


template < typename MeshType, typename ExpressionType>
IntegrateValueElement < MeshType, ExpressionType>::
~IntegrateValueElement()
{
    delete M_globalCFE;
}


// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename ExpressionType>
void
IntegrateValueElement < MeshType, ExpressionType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
}


template < typename MeshType, typename ExpressionType>
void
IntegrateValueElement < MeshType, ExpressionType>::
addTo (Real& value)
{
    UInt nbElements (M_mesh->numElements() );
    UInt nbQuadPt (M_quadrature.nbQuadPt() );

    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Update the currentFEs
        M_globalCFE->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);

        // Update the evaluation
        M_evaluation.update (iElement);


        // Make the assembly
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {
            value += M_evaluation.value_q (iQuadPt)
                     * M_globalCFE->wDet (iQuadPt);
        }
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
