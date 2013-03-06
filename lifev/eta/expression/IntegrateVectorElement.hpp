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
     @brief This file contains the definition of the IntegrateVectorElement class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_VECTOR_ELEMENT_HPP
#define INTEGRATE_VECTOR_ELEMENT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETVectorElemental.hpp>

#include <boost/shared_ptr.hpp>



namespace LifeV
{

namespace ExpressionAssembly
{

//! The class to actually perform the loop over the elements to assemble a vector
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is used to store the data required for the assembly of a vector and
  perform that assembly with a loop over the elements, and then, for each elements,
  using the Evaluation corresponding to the Expression (this convertion is done
  within a typedef).
 */
template < typename MeshType, typename TestSpaceType, typename ExpressionType>
class IntegrateVectorElement
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
            TestSpaceType::S_fieldDim,
            0,
            MeshType::S_geoDimensions >::evaluation_Type evaluation_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateVectorElement (const boost::shared_ptr<MeshType>& mesh,
                            const QuadratureRule& quadrature,
                            const boost::shared_ptr<TestSpaceType>& testSpace,
                            const ExpressionType& expression);

    //! Copy constructor
    IntegrateVectorElement ( const IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType>& integrator);

    //! Destructor
    ~IntegrateVectorElement();

    //@}


    //! @name Operator
    //@{

    //! Operator wrapping the addTo method
    template <typename VectorType>
    inline void operator>> (VectorType& vector)
    {
        addTo (vector);
    }

    //! Operator wrapping the addTo method (for shared_ptr)
    template <typename VectorType>
    inline void operator>> (boost::shared_ptr<VectorType> vector)
    {
        addTo (*vector);
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
      performed: update the values, update the local vector,
      sum over the quadrature nodes, assemble in the global
      vector.
     */
    template <typename VectorType>
    void addTo (VectorType& vec);

    //@}

private:

    //! @name Private Methods
    //@{

    // No default constructor
    IntegrateVectorElement();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Quadrature to be used
    QuadratureRule M_quadrature;

    // Shared pointer on the Space
    boost::shared_ptr<TestSpaceType> M_testSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    ETCurrentFE<MeshType::S_geoDimensions, 1>* M_globalCFE;
    ETCurrentFE<TestSpaceType::S_spaceDim, TestSpaceType::S_fieldDim>* M_testCFE;

    //ETVectorElemental<1> M_elementalVector;
    ETVectorElemental M_elementalVector;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType>::
IntegrateVectorElement (const boost::shared_ptr<MeshType>& mesh,
                        const QuadratureRule& quadrature,
                        const boost::shared_ptr<TestSpaceType>& testSpace,
                        const ExpressionType& expression)
    :   M_mesh (mesh),
        M_quadrature (quadrature),
        M_testSpace (testSpace),
        M_evaluation (expression),

        //M_globalCFE(new ETCurrentFE<3,1>(feTetraP0,geometricMapFromMesh<MeshType>(),quadrature)),
        M_testCFE (new ETCurrentFE<TestSpaceType::S_spaceDim, TestSpaceType::S_fieldDim> (testSpace->refFE(), testSpace->geoMap(), quadrature) ),

        //M_elementalVector(testSpace->refFE().nbDof())
        M_elementalVector (TestSpaceType::S_fieldDim * testSpace->refFE().nbDof() )
{
    switch (MeshType::geoShape_Type::BasRefSha::S_shape)
    {
        case LINE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case TRIANGLE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case QUAD:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case TETRA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        case HEXA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), quadrature);
            break;
        default:
            ERROR_MSG ("Unrecognized element shape");
    }
    M_evaluation.setQuadrature (quadrature);
    M_evaluation.setGlobalCFE (M_globalCFE);
    M_evaluation.setTestCFE (M_testCFE);
}


template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType>::
IntegrateVectorElement ( const IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType>& integrator)
    :   M_mesh (integrator.M_mesh),
        M_quadrature (integrator.M_quadrature),
        M_testSpace (integrator.M_testSpace),
        M_evaluation (integrator.M_evaluation),

        //M_globalCFE(new ETCurrentFE<3,1>(feTetraP0,geometricMapFromMesh<MeshType>(),M_quadrature)),
        M_testCFE (new ETCurrentFE<TestSpaceType::S_spaceDim, TestSpaceType::S_fieldDim> (M_testSpace->refFE(), M_testSpace->geoMap(), M_quadrature) ),

        M_elementalVector (integrator.M_elementalVector)
{
    switch (MeshType::geoShape_Type::BasRefSha::S_shape)
    {
        case LINE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case TRIANGLE:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case QUAD:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case TETRA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        case HEXA:
            M_globalCFE = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), M_quadrature);
            break;
        default:
            ERROR_MSG ("Unrecognized element shape");
    }
    M_evaluation.setQuadrature (M_quadrature);
    M_evaluation.setGlobalCFE (M_globalCFE);
    M_evaluation.setTestCFE (M_testCFE);
}


template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType>::
~IntegrateVectorElement()
{
    delete M_globalCFE;
    delete M_testCFE;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType>
void
IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
    out << " Elemental vector : " << std::endl;
    M_elementalVector.showMe (out);
    out << std::endl;
}

template < typename MeshType, typename TestSpaceType, typename ExpressionType>
template <typename VectorType>
void
IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType>::
addTo (VectorType& vec)
{
    UInt nbElements (M_mesh->numElements() );
    UInt nbQuadPt (M_quadrature.nbQuadPt() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );

    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Zeros out the elemental vector
        M_elementalVector.zero();

        // Update the currentFEs
        M_globalCFE->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
        M_testCFE->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);

        // Update the evaluation
        M_evaluation.update (iElement);

        // Loop on the blocks

        for (UInt iblock (0); iblock < TestSpaceType::S_fieldDim; ++iblock)
        {
            // Set the row global indices in the local vector
            for (UInt i (0); i < nbTestDof; ++i)
            {
                /*M_elementalVector.block(0).setRowIndex
                    (i + iblock*nbTestDof,
                    M_testSpace->dof().localToGlobalMap(iElement,i)+ iblock*M_testSpace->dof().numTotalDof());*/
                M_elementalVector.setRowIndex
                (i + iblock * nbTestDof,
                 M_testSpace->dof().localToGlobalMap (iElement, i) + iblock * M_testSpace->dof().numTotalDof() );
            }

            // Make the assembly
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    /*M_elementalVector.block(0).element(i+iblock*nbTestDof) +=
                        M_evaluation.value_qi(iQuadPt,i+iblock*nbTestDof)
                        * M_globalCFE->wDet(iQuadPt);*/
                    M_elementalVector.element (i + iblock * nbTestDof) +=
                        M_evaluation.value_qi (iQuadPt, i + iblock * nbTestDof)
                        * M_globalCFE->wDet (iQuadPt);

                }
            }
        }

        M_elementalVector.pushToGlobal (vec);
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
