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
     @brief This file contains the definition of the IntegrateVectorElementLSAdapted class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_VECTOR_ELEMENT_LS_ADAPTED_HPP
#define INTEGRATE_VECTOR_ELEMENT_LS_ADAPTED_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/level_set/fem/LevelSetQRAdapter.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETVectorElemental.hpp>

#include <boost/shared_ptr.hpp>



namespace LifeV
{

namespace ExpressionAssembly
{

//! The class to actually perform the loop over the elements to assemble a vector using special QR.
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is used to store the data required for the assembly of a vector and
  perform that assembly with a loop over the elements, and then, for each elements,
  using the Evaluation corresponding to the Expression (this convertion is done
  within a typedef).

  The quadrature used is given by the structure adapting the quadrature to the
  situation.
 */
template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
class IntegrateVectorElementLSAdapted
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
            TestSpaceType::S_fieldDim,
            0,
            3 >::evaluation_Type evaluation_Type;

    typedef LevelSetQRAdapter<LSFESpaceType, LSVectorType> QRAdapter_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateVectorElementLSAdapted (const boost::shared_ptr<MeshType>& mesh,
                                     const QRAdapter_Type& QRAdapter,
                                     const boost::shared_ptr<TestSpaceType>& testSpace,
                                     const ExpressionType& expression);

    //! Copy constructor
    IntegrateVectorElementLSAdapted ( const IntegrateVectorElementLSAdapted
                                      < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>& integrator);

    //! Destructor
    ~IntegrateVectorElementLSAdapted();

    //@}


    //! @name Operator
    //@{

    //! Operator wrapping the addTo method
    template <typename VectorType>
    inline void operator>> (VectorType& vector)
    {
        addTo (vector);
    }

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

    template <typename VectorType>
    void addTo (boost::shared_ptr<VectorType> vec)
    {
        addTo (*vec);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    // No default constructor
    IntegrateVectorElementLSAdapted();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Quadrature adapter
    QRAdapter_Type M_QRAdapter;

    // Shared pointer on the Space
    boost::shared_ptr<TestSpaceType> M_testSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    // CurrentFEs are duplicated to consider the
    // adapted and unadapted quadrature rules.
    ETCurrentFE<3, 1>* M_globalCFE_unadapted;
    ETCurrentFE<3, 1>* M_globalCFE_adapted;

    ETCurrentFE<3, TestSpaceType::S_fieldDim>* M_testCFE_unadapted;
    ETCurrentFE<3, TestSpaceType::S_fieldDim>* M_testCFE_adapted;

    //ETVectorElemental<1> M_elementalVector;
    ETVectorElemental M_elementalVector;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateVectorElementLSAdapted < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
IntegrateVectorElementLSAdapted (const boost::shared_ptr<MeshType>& mesh,
                                 const QRAdapter_Type& QRAdapter,
                                 const boost::shared_ptr<TestSpaceType>& testSpace,
                                 const ExpressionType& expression)
    :   M_mesh (mesh),
        M_QRAdapter (QRAdapter),
        M_testSpace (testSpace),
        M_evaluation (expression),

        M_globalCFE_unadapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) ),
        M_globalCFE_adapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) ),

        M_testCFE_unadapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (testSpace->refFE(), testSpace->geoMap(), M_QRAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (testSpace->refFE(), testSpace->geoMap(), M_QRAdapter.standardQR() ) ),

        M_elementalVector (TestSpaceType::S_fieldDim * testSpace->refFE().nbDof() )
{
    M_evaluation.setQuadrature (M_QRAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_unadapted);
    M_evaluation.setTestCFE (M_testCFE_unadapted);
}


template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateVectorElementLSAdapted < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
IntegrateVectorElementLSAdapted ( const IntegrateVectorElementLSAdapted
                                  < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>& integrator)
    :   M_mesh (integrator.M_mesh),
        M_QRAdapter (integrator.M_QRAdapter),
        M_testSpace (integrator.M_testSpace),
        M_evaluation (integrator.M_evaluation),

        M_globalCFE_unadapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) ),
        M_globalCFE_adapted (new ETCurrentFE<3, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) ),

        M_testCFE_unadapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (M_testSpace->refFE(), M_testSpace->geoMap(), M_QRAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<3, TestSpaceType::S_fieldDim> (M_testSpace->refFE(), M_testSpace->geoMap(), M_QRAdapter.standardQR() ) ),

        M_elementalVector (integrator.M_elementalVector)
{
    M_evaluation.setQuadrature (M_QRAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_unadapted);
    M_evaluation.setTestCFE (M_testCFE_unadapted);
}


template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateVectorElementLSAdapted < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
~IntegrateVectorElementLSAdapted()
{
    delete M_globalCFE_unadapted;
    delete M_globalCFE_adapted;
    delete M_testCFE_unadapted;
    delete M_testCFE_adapted;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
void
IntegrateVectorElementLSAdapted < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
    out << " Elemental vector : " << std::endl;
    M_elementalVector.showMe (out);
    out << std::endl;
}

template < typename MeshType,
         typename TestSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
template <typename VectorType>
void
IntegrateVectorElementLSAdapted < MeshType, TestSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
addTo (VectorType& vec)
{
    const UInt nbElements (M_mesh->numElements() );
    const UInt nbQuadPt_unadapted (M_QRAdapter.standardQR().nbQuadPt() );
    const UInt nbTestDof (M_testSpace->refFE().nbDof() );

    // Even if this is not the case, this
    // prevents many potentially buggy behaviours
    // and does not harm.
    bool isPreviousAdapted (true);

    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Zeros out the elemental vector
        M_elementalVector.zero();

        // Update the quadrature rule
        M_QRAdapter.update (iElement);


        if (M_QRAdapter.isAdaptedElement() )
        {
            // Change the QR
            M_evaluation.setQuadrature ( M_QRAdapter.adaptedQR() );
            M_globalCFE_adapted->setQuadratureRule ( M_QRAdapter.adaptedQR() );
            M_testCFE_adapted->setQuadratureRule ( M_QRAdapter.adaptedQR() );

            // Reset the CFEs in the evaluation
            M_evaluation.setGlobalCFE ( M_globalCFE_adapted );
            M_evaluation.setTestCFE ( M_testCFE_adapted );

            // Update the currentFEs
            M_globalCFE_adapted->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
            M_testCFE_adapted->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);

            // Update the evaluation
            M_evaluation.update (iElement);

            // Loop on the blocks

            for (UInt iblock (0); iblock < TestSpaceType::S_fieldDim; ++iblock)
            {
                // Set the row global indices in the local vector
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    M_elementalVector.setRowIndex
                    (i + iblock * nbTestDof,
                     M_testSpace->dof().localToGlobalMap (iElement, i) + iblock * M_testSpace->dof().numTotalDof() );
                }

                // Make the assembly
                for (UInt iQuadPt (0); iQuadPt < M_QRAdapter.adaptedQR().nbQuadPt(); ++iQuadPt)
                {
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        M_elementalVector.element (i + iblock * nbTestDof) +=
                            M_evaluation.value_qi (iQuadPt, i + iblock * nbTestDof)
                            * M_globalCFE_adapted->wDet (iQuadPt);
                    }
                }
            }

            M_elementalVector.pushToGlobal (vec);

            isPreviousAdapted = true;

        }
        else
        {
            if (isPreviousAdapted)
            {
                M_evaluation.setQuadrature ( M_QRAdapter.standardQR() );
                M_evaluation.setGlobalCFE ( M_globalCFE_unadapted );
                M_evaluation.setTestCFE ( M_testCFE_unadapted );
            }

            // Update the currentFEs
            M_globalCFE_unadapted->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
            M_testCFE_unadapted->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);

            // Update the evaluation
            M_evaluation.update (iElement);

            // Loop on the blocks

            for (UInt iblock (0); iblock < TestSpaceType::S_fieldDim; ++iblock)
            {
                // Set the row global indices in the local vector
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    M_elementalVector.setRowIndex
                    (i + iblock * nbTestDof,
                     M_testSpace->dof().localToGlobalMap (iElement, i) + iblock * M_testSpace->dof().numTotalDof() );
                }

                // Make the assembly
                for (UInt iQuadPt (0); iQuadPt < nbQuadPt_unadapted; ++iQuadPt)
                {
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        M_elementalVector.element (i + iblock * nbTestDof) +=
                            M_evaluation.value_qi (iQuadPt, i + iblock * nbTestDof)
                            * M_globalCFE_unadapted->wDet (iQuadPt);
                    }
                }
            }

            M_elementalVector.pushToGlobal (vec);

            isPreviousAdapted = false;

        }
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif // INTEGRATE_VECTOR_ELEMENT_LS_ADAPTED_HPP
