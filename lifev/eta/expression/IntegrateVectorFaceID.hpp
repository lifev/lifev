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

#ifndef INTEGRATE_VECTOR_FACE_ID_HPP
#define INTEGRATE_VECTOR_FACE_ID_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/QuadratureBoundary.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentBDFE.hpp>
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
class IntegrateVectorFaceID
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
            TestSpaceType::field_dim,
            0,
            3 >::evaluation_Type evaluation_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateVectorFaceID (const boost::shared_ptr<MeshType>& mesh,
                           const UInt boundaryID,
                           const QuadratureBoundary& quadratureBD,
                           const boost::shared_ptr<TestSpaceType>& testSpace,
                           const ExpressionType& expression);

    //! Copy constructor
    IntegrateVectorFaceID ( const IntegrateVectorFaceID < MeshType, TestSpaceType, ExpressionType>& integrator);

    //! Destructor
    ~IntegrateVectorFaceID();

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
    IntegrateVectorFaceID();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Identifier for the boundary
    UInt M_boundaryId;

    // Quadrature to be used
    QuadratureBoundary M_quadratureBoundary;

    // Shared pointer on the Space
    boost::shared_ptr<TestSpaceType> M_testSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    std::vector<ETCurrentBDFE<3>*> M_globalCFE;
    std::vector<ETCurrentFE<3, TestSpaceType::field_dim>*> M_testCFE;

    ETVectorElemental M_elementalVector;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorFaceID < MeshType, TestSpaceType, ExpressionType>::
IntegrateVectorFaceID (const boost::shared_ptr<MeshType>& mesh,
                       const UInt boundaryID,
                       const QuadratureBoundary& quadratureBD,
                       const boost::shared_ptr<TestSpaceType>& testSpace,
                       const ExpressionType& expression)
    :   M_mesh (mesh),
        M_boundaryId (boundaryID),
        M_quadratureBoundary (quadratureBD),
        M_testSpace (testSpace),
        M_evaluation (expression),

        M_globalCFE (4),
        M_testCFE (4),

        M_elementalVector (TestSpaceType::field_dim * testSpace->refFE().nbDof() )
{
    for (UInt i (0); i < 4; ++i)
    {
        M_globalCFE[i] = new ETCurrentBDFE<3> (geometricMapFromMesh<MeshType>()
                                               , M_quadratureBoundary.qr (i) );
        M_testCFE[i] = new ETCurrentFE<3, TestSpaceType::field_dim> (testSpace->refFE()
                                                                     , testSpace->geoMap()
                                                                     , M_quadratureBoundary.qr (i) );
    }

    // Set the tangent on the different faces
    std::vector< VectorSmall<3> > t0 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    t0[0][0] = 1;
    t0[0][1] = 0;
    t0[0][2] = 0;
    t0[1][0] = 0;
    t0[1][1] = 1;
    t0[1][2] = 0;
    std::vector< VectorSmall<3> > t1 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    t1[0][0] = 0;
    t1[0][1] = 0;
    t1[0][2] = 1;
    t1[1][0] = 1;
    t1[1][1] = 0;
    t1[1][2] = 0;
    std::vector< VectorSmall<3> > t2 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    //t2[0][0]=-1/std::sqrt(6);    t2[0][1]=-1/std::sqrt(6);    t2[0][2]=2/std::sqrt(6);
    //t2[1][0]=-1/std::sqrt(2);    t2[1][1]=1/std::sqrt(2);    t2[1][2]=0;
    t2[0][0] = -1;
    t2[0][1] = 0;
    t2[0][2] = 1;
    t2[1][0] = -1;
    t2[1][1] = 1;
    t2[1][2] = 0;

    std::vector< VectorSmall<3> > t3 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    t3[0][0] = 0;
    t3[0][1] = 1;
    t3[0][2] = 0;
    t3[1][0] = 0;
    t3[1][1] = 0;
    t3[1][2] = 1;

    M_globalCFE[0]->setRefTangents (t0);
    M_globalCFE[1]->setRefTangents (t1);
    M_globalCFE[2]->setRefTangents (t2);
    M_globalCFE[3]->setRefTangents (t3);


    M_evaluation.setQuadrature (M_quadratureBoundary.qr (0) );
    M_evaluation.setGlobalCFE (M_globalCFE[0]);
    M_evaluation.setTestCFE (M_testCFE[0]);
}


template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorFaceID < MeshType, TestSpaceType, ExpressionType>::
IntegrateVectorFaceID ( const IntegrateVectorFaceID < MeshType, TestSpaceType, ExpressionType>& integrator)
    :   M_mesh (integrator.M_mesh),
        M_boundaryId (integrator.M_boundaryId),
        M_quadratureBoundary (integrator.M_quadratureBoundary),
        M_testSpace (integrator.M_testSpace),
        M_evaluation (integrator.M_evaluation),

        M_globalCFE (4),
        M_testCFE (4),

        M_elementalVector (integrator.M_elementalVector)
{
    for (UInt i (0); i < 4; ++i)
    {
        M_globalCFE[i] = new ETCurrentBDFE<3> (geometricMapFromMesh<MeshType>()
                                               , M_quadratureBoundary.qr (i) );
        M_testCFE[i] = new ETCurrentFE<3, TestSpaceType::field_dim> (M_testSpace->refFE()
                                                                     , M_testSpace->geoMap()
                                                                     , M_quadratureBoundary.qr (i) );
    }

    // Set the tangent on the different faces
    std::vector< VectorSmall<3> > t0 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    t0[0][0] = 1;
    t0[0][1] = 0;
    t0[0][2] = 0;
    t0[1][0] = 0;
    t0[1][1] = 1;
    t0[1][2] = 0;
    std::vector< VectorSmall<3> > t1 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    t1[0][0] = 0;
    t1[0][1] = 0;
    t1[0][2] = 1;
    t1[1][0] = 1;
    t1[1][1] = 0;
    t1[1][2] = 0;
    std::vector< VectorSmall<3> > t2 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    //t2[0][0]=-1/std::sqrt(6);    t2[0][1]=-1/std::sqrt(6);    t2[0][2]=2/std::sqrt(6);
    //t2[1][0]=-1/std::sqrt(2);    t2[1][1]=1/std::sqrt(2);    t2[1][2]=0;
    t2[0][0] = -1;
    t2[0][1] = 0;
    t2[0][2] = 1;
    t2[1][0] = -1;
    t2[1][1] = 1;
    t2[1][2] = 0;

    std::vector< VectorSmall<3> > t3 (2, VectorSmall<3> (0.0, 0.0, 0.0) );
    t3[0][0] = 0;
    t3[0][1] = 1;
    t3[0][2] = 0;
    t3[1][0] = 0;
    t3[1][1] = 0;
    t3[1][2] = 1;

    M_globalCFE[0]->setRefTangents (t0);
    M_globalCFE[1]->setRefTangents (t1);
    M_globalCFE[2]->setRefTangents (t2);
    M_globalCFE[3]->setRefTangents (t3);


    M_evaluation.setQuadrature (M_quadratureBoundary.qr (0) );
    M_evaluation.setGlobalCFE (M_globalCFE[0]);
    M_evaluation.setTestCFE (M_testCFE[0]);
}


template < typename MeshType, typename TestSpaceType, typename ExpressionType>
IntegrateVectorFaceID < MeshType, TestSpaceType, ExpressionType>::
~IntegrateVectorFaceID()
{
    for (UInt i (0); i < 4; ++i)
    {
        delete M_globalCFE[i];
        delete M_testCFE[i];
    }
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType>
void
IntegrateVectorFaceID < MeshType, TestSpaceType, ExpressionType>::
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
IntegrateVectorFaceID < MeshType, TestSpaceType, ExpressionType>::
addTo (VectorType& vec)
{
    UInt nbBoundaryFaces (M_mesh->numBFaces() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );

    /*CurrentBoundaryFE origFE(M_testSpace->refFE().boundaryFE(),
                             M_testSpace->geoMap().boundaryMap(),
                             quadRuleTria3pt);*/

    for (UInt iFace (0); iFace < nbBoundaryFaces; ++iFace)
    {
        // Check the identifier
        if ( M_mesh->face (iFace).markerID() != M_boundaryId )
        {
            continue;
        }

        //origFE.updateMeasNormalQuadPt(M_mesh->face(iFace));

        // Zeros out the elemental vector
        M_elementalVector.zero();

        // Get the number of the face in the adjacent element
        UInt faceIDinAdjacentElement (M_mesh->face (iFace).firstAdjacentElementPosition() );

        // Get the ID of the adjacent element
        UInt adjacentElementID (M_mesh->face (iFace).firstAdjacentElementIdentity() );

        // Update the currentFEs
        M_globalCFE[faceIDinAdjacentElement]
        ->update (M_mesh->element (adjacentElementID) );
        M_testCFE[faceIDinAdjacentElement]
        ->update (M_mesh->element (adjacentElementID), evaluation_Type::S_testUpdateFlag);

        // Update the evaluation
        M_evaluation.setQuadrature (M_quadratureBoundary.qr (faceIDinAdjacentElement) );
        M_evaluation.setGlobalCFE (M_globalCFE[faceIDinAdjacentElement]);
        M_evaluation.setTestCFE (M_testCFE[faceIDinAdjacentElement]);

        M_evaluation.update (adjacentElementID);

        // Loop on the blocks
        for (UInt iblock (0); iblock < TestSpaceType::field_dim; ++iblock)
        {
            // Set the row global indices in the local vector
            for (UInt i (0); i < nbTestDof; ++i)
            {
                M_elementalVector.setRowIndex
                (i + iblock * nbTestDof,
                 M_testSpace->dof().localToGlobalMap (adjacentElementID, i) + iblock * M_testSpace->dof().numTotalDof() );
            }

            // Make the assembly
            for (UInt iQuadPt (0); iQuadPt < M_quadratureBoundary.qr (faceIDinAdjacentElement).nbQuadPt(); ++iQuadPt)
            {
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    M_elementalVector.element (i + iblock * nbTestDof) +=
                        M_evaluation.value_qi (iQuadPt, i + iblock * nbTestDof)
                        * M_globalCFE[faceIDinAdjacentElement]->M_wMeas[iQuadPt];

                }
            }
        }

        M_elementalVector.pushToGlobal (vec);
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
