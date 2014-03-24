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
     @brief This file contains the definition of the IntegrateMatrixFaceID class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_MATRIX_FACE_ID_LSADAPTED_HPP
#define INTEGRATE_MATRIX_FACE_ID_LSADAPTED_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/eta/fem/QuadratureBoundary.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/ETCurrentBDFE.hpp>

#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETMatrixElemental.hpp>

#include <lifev/eta/fem/LevelSetBDQRAdapter.hpp>

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
template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
class IntegrateMatrixFaceIDLSAdapted
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
					      TestSpaceType::field_dim,
					      SolutionSpaceType::field_dim,
					      3 >::evaluation_Type evaluation_Type;

    typedef LevelSetBDQRAdapter<LSFESpaceType, LSVectorType> BDQRAdapter_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateMatrixFaceIDLSAdapted (const boost::shared_ptr<MeshType>& mesh,
                                    const UInt boundaryID,
                                    const BDQRAdapter_Type& quadratureBD,
                                    const boost::shared_ptr<TestSpaceType> testSpace,
                                    const boost::shared_ptr<SolutionSpaceType> solutionSpace,
                                    const ExpressionType& expression);

    //! Copy constructor
    IntegrateMatrixFaceIDLSAdapted ( const IntegrateMatrixFaceIDLSAdapted
                                     < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>& integrator);

    //! Destructor
    ~IntegrateMatrixFaceIDLSAdapted();

    //@}


    //! @name Operator
    //@{

    //! Operator wrapping the addTo method
    template <typename MatrixType>
    inline void operator>> (MatrixType& mat)
    {
        addTo (mat);
    }

    //! Operator wrapping the addTo method (for shared_ptr)
    template <typename MatrixType>
    inline void operator>> (boost::shared_ptr<MatrixType> mat)
    {
        addTo (mat);
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
    template <typename MatrixType>
    void addTo (MatrixType& mat);

    template <typename MatrixType>
    inline void addTo (boost::shared_ptr<MatrixType> mat)
    {
        ASSERT (mat != 0, " Cannot assemble with an empty matrix");
        addTo (*mat);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    // No default constructor
    IntegrateMatrixFaceIDLSAdapted();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Identifier for the boundary
    UInt M_boundaryId;

    // Quadrature to be used
    BDQRAdapter_Type M_qrAdapter;

    // Shared pointer on the Space
    boost::shared_ptr<TestSpaceType> M_testSpace;
    boost::shared_ptr<SolutionSpaceType> M_solutionSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    std::vector<ETCurrentBDFE<3>*> M_globalCFE;

    std::vector<ETCurrentFE<3, TestSpaceType::field_dim>*> M_testCFE;
    std::vector<ETCurrentFE<3, SolutionSpaceType::field_dim>*> M_solutionCFE;

    ETMatrixElemental M_elementalMatrix;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
IntegrateMatrixFaceIDLSAdapted (const boost::shared_ptr<MeshType>& mesh,
                                const UInt boundaryID,
                                const BDQRAdapter_Type& quadratureBD,
                                const boost::shared_ptr<TestSpaceType> testSpace,
                                const boost::shared_ptr<SolutionSpaceType> solutionSpace,
                                const ExpressionType& expression)
    :   M_mesh (mesh),
        M_boundaryId (boundaryID),
        M_qrAdapter (quadratureBD),
        M_testSpace (testSpace),
        M_solutionSpace (solutionSpace),
        M_evaluation (expression),

        M_globalCFE (4),
        M_testCFE (4),
        M_solutionCFE (4),

        M_elementalMatrix (TestSpaceType::field_dim * testSpace->refFE().nbDof(), SolutionSpaceType::field_dim * solutionSpace->refFE().nbDof() )
{
    for (UInt i (0); i < 4; ++i)
    {
        M_globalCFE[i] = new ETCurrentBDFE<3> (geometricMapFromMesh<MeshType>()
                                               , M_qrAdapter.adaptedBdQR (i) );
        M_testCFE[i] = new ETCurrentFE<3, TestSpaceType::field_dim> (testSpace->refFE()
                                                                     , testSpace->geoMap()
                                                                     , M_qrAdapter.adaptedBdQR (i) );
        M_solutionCFE[i] = new ETCurrentFE<3, SolutionSpaceType::field_dim> (solutionSpace->refFE()
                                                                             , solutionSpace->geoMap()
                                                                             , M_qrAdapter.adaptedBdQR (i) );

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


    M_evaluation.setQuadrature (M_qrAdapter.adaptedBdQR (0) );
    M_evaluation.setGlobalCFE (M_globalCFE[0]);
    M_evaluation.setTestCFE (M_testCFE[0]);
    M_evaluation.setSolutionCFE (M_solutionCFE[0]);
}

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
IntegrateMatrixFaceIDLSAdapted ( const IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>& integrator)
    :   M_mesh (integrator.M_mesh),
        M_boundaryId (integrator.M_boundaryId),
        M_qrAdapter (integrator.M_qrAdapter),
        M_testSpace (integrator.M_testSpace),
        M_solutionSpace (integrator.M_solutionSpace),
        M_evaluation (integrator.M_evaluation),

        M_globalCFE (4),
        M_testCFE (4),
        M_solutionCFE (4),

        M_elementalMatrix (integrator.M_elementalMatrix)
{
    for (UInt i (0); i < 4; ++i)
    {
        M_globalCFE[i] = new ETCurrentBDFE<3> (geometricMapFromMesh<MeshType>()
                                               , M_qrAdapter.adaptedBdQR (i) );
        M_testCFE[i] = new ETCurrentFE<3, TestSpaceType::field_dim> (M_testSpace->refFE()
                                                                     , M_testSpace->geoMap()
                                                                     , M_qrAdapter.adaptedBdQR (i) );
        M_solutionCFE[i] = new ETCurrentFE<3, SolutionSpaceType::field_dim> (M_solutionSpace->refFE()
                                                                             , M_solutionSpace->geoMap()
                                                                             , M_qrAdapter.adaptedBdQR (i) );
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


    M_evaluation.setQuadrature (M_qrAdapter.adaptedBdQR (0) );
    M_evaluation.setGlobalCFE (M_globalCFE[0]);
    M_evaluation.setTestCFE (M_testCFE[0]);
    M_evaluation.setSolutionCFE (M_solutionCFE[0]);
}


template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
~IntegrateMatrixFaceIDLSAdapted()
{
    for (UInt i (0); i < 4; ++i)
    {
        delete M_globalCFE[i];
        delete M_testCFE[i];
        delete M_solutionCFE[i];
    }
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
void
IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
    out << " Elemental matrix : " << std::endl;
    M_elementalMatrix.showMe (out);
    out << std::endl;
}


template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename LSFESpaceType,
         typename LSVectorType >
template <typename MatrixType>
void
IntegrateMatrixFaceIDLSAdapted < MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, LSFESpaceType, LSVectorType>::
addTo (MatrixType& mat)
{
    UInt nbBoundaryFaces (M_mesh->numBFaces() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );
    UInt nbSolutionDof (M_solutionSpace->refFE().nbDof() );

    for (UInt iFace (0); iFace < nbBoundaryFaces; ++iFace)
    {
        // Check the identifier
        if ( M_mesh->face (iFace).markerID() != M_boundaryId )
        {
            continue;
        }

        // Zeros out the elemental vector
        M_elementalMatrix.zero();

        // Get the number of the face in the adjacent element
        UInt faceIDinAdjacentElement (M_mesh->face (iFace).firstAdjacentElementPosition() );

        // Get the ID of the adjacent element
        UInt adjacentElementID (M_mesh->face (iFace).firstAdjacentElementIdentity() );

        // Update the QR
        M_qrAdapter.update (adjacentElementID, faceIDinAdjacentElement);

        // Change the qr
        M_globalCFE[faceIDinAdjacentElement]->setQuadratureRule (M_qrAdapter.adaptedBdQR (faceIDinAdjacentElement) );
        M_testCFE[faceIDinAdjacentElement]->setQuadratureRule (M_qrAdapter.adaptedBdQR (faceIDinAdjacentElement) );
        M_solutionCFE[faceIDinAdjacentElement]->setQuadratureRule (M_qrAdapter.adaptedBdQR (faceIDinAdjacentElement) );

        // Update the currentFEs
        M_globalCFE[faceIDinAdjacentElement]
        ->update (M_mesh->element (adjacentElementID) );
        M_testCFE[faceIDinAdjacentElement]
        ->update (M_mesh->element (adjacentElementID), evaluation_Type::S_testUpdateFlag);
        M_solutionCFE[faceIDinAdjacentElement]
        ->update (M_mesh->element (adjacentElementID), evaluation_Type::S_solutionUpdateFlag);


        // Update the evaluation
        M_evaluation.setQuadrature (M_qrAdapter.adaptedBdQR (faceIDinAdjacentElement) );
        M_evaluation.setGlobalCFE (M_globalCFE[faceIDinAdjacentElement]);
        M_evaluation.setTestCFE (M_testCFE[faceIDinAdjacentElement]);
        M_evaluation.setSolutionCFE (M_solutionCFE[faceIDinAdjacentElement]);

        M_evaluation.update (adjacentElementID);

        // Loop on the blocks
        for (UInt iblock (0); iblock < TestSpaceType::field_dim; ++iblock)
        {
            for (UInt jblock (0); jblock < SolutionSpaceType::field_dim; ++jblock)
            {

                // Set the row global indices in the local matrix
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    M_elementalMatrix.setRowIndex
                    (i + iblock * nbTestDof,
                     M_testSpace->dof().localToGlobalMap (adjacentElementID, i) + iblock * M_testSpace->dof().numTotalDof() );
                }

                for (UInt j (0); j < nbSolutionDof; ++j)
                {
                    M_elementalMatrix.setColumnIndex
                    (j + jblock * nbSolutionDof,
                     M_solutionSpace->dof().localToGlobalMap (adjacentElementID, j) + jblock * M_solutionSpace->dof().numTotalDof() );
                }


                // Make the assembly
                for (UInt iQuadPt (0); iQuadPt < M_qrAdapter.adaptedBdQR (faceIDinAdjacentElement).nbQuadPt(); ++iQuadPt)
                {
                    for (UInt i (0); i < nbTestDof; ++i)
                    {
                        for (UInt j (0); j < nbSolutionDof; ++j)
                        {
                            M_elementalMatrix.element (i + iblock * nbTestDof, j + jblock * nbSolutionDof) +=
                                M_evaluation.value_qij (iQuadPt, i + iblock * nbTestDof, j + jblock * nbSolutionDof)
                                * M_globalCFE[faceIDinAdjacentElement]->M_wMeas[iQuadPt];
                        }
                    }
                }
            }
        }

        M_elementalMatrix.pushToGlobal (mat);
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
