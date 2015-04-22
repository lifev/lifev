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
     @brief This file contains the definition of the IntegrateMatrixElement class.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef INTEGRATE_MATRIX_ELEMENT_HPP
#define INTEGRATE_MATRIX_ELEMENT_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include <lifev/core/LifeV.hpp>

#include <lifev/core/util/OpenMPParameters.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>
#include <lifev/eta/fem/QRAdapterBase.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETMatrixElemental.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>



namespace LifeV
{

namespace ExpressionAssembly
{


//! The class to actually perform the loop over the elements to assemble a matrix
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  This class is used to store the data required for the assembly of a matrix and
  perform that assembly with a loop over the elements, and then, for each elements,
  using the Evaluation corresponding to the Expression (This convertion is done
  within a typedef).
 */
template < typename MeshType,
         typename TestSpaceType,
         typename SolutionSpaceType,
         typename ExpressionType,
         typename QRAdapterType >

class IntegrateMatrixElement
{
public:

    //! @name Public Types
    //@{

    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
            TestSpaceType::field_dim,
            SolutionSpaceType::field_dim,
            MeshType::S_geoDimensions >::evaluation_Type  evaluation_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
    IntegrateMatrixElement (const boost::shared_ptr<MeshType>& mesh,
                            const QRAdapterType& qrAdapter,
                            const boost::shared_ptr<TestSpaceType>& testSpace,
                            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                            const ExpressionType& expression,
                            const UInt offsetUp = 0,
                            const UInt offsetLeft = 0,
                            const UInt regionFlag = 0,
                            const UInt * const volumeElements = NULL );

    //! Full data constructor
    IntegrateMatrixElement (const boost::shared_ptr<MeshType>& mesh,
                            const QRAdapterType& qrAdapter,
                            const boost::shared_ptr<TestSpaceType>& testSpace,
                            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                            const ExpressionType& expression,
                            const OpenMPParameters& ompParams,
                            const UInt offsetUp = 0,
                            const UInt offsetLeft = 0,
                            const UInt regionFlag = 0,
                            const UInt * const volumeElements = NULL );

    //! Copy constructor
    IntegrateMatrixElement (const IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>& integrator);

    //! Destructor
    ~IntegrateMatrixElement();

    //@}


    //! @name Operators
    //@{

    //! Operator wrapping the addTo method
    template <typename MatrixType>
    inline void operator>> (MatrixType& mat)
    {
        if ( mat.filled() )
        {
            addToClosed (mat);
        }
        else
        {
            addTo (mat);
        }
    }

    //! Operator wrapping the addTo method (for shared_ptr)
    template <typename MatrixType>
    inline void operator>> (boost::shared_ptr<MatrixType> mat)
    {
        if (mat->filled() )
        {
            addToClosed (mat);
        }
        else
        {
            addTo (mat);
        }
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
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.
     */
    template <typename MatrixType>
    void addTo (MatrixType& mat);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.
      The method is used for closed matrices
     */
    template <typename MatrixType>
    void addToClosed (MatrixType& mat);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.

      Specialized for the case where the matrix is passed as a shared_ptr
     */
    template <typename MatrixType>
    inline void addTo (boost::shared_ptr<MatrixType> mat)
    {
        ASSERT (mat != 0, " Cannot assemble with an empty matrix");
        addTo (*mat);
    }

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local matrix,
      sum over the quadrature nodes, assemble in the global
      matrix.
      This method is used with closed matrices.

      Specialized for the case where the matrix is passed as a shared_ptr
     */
    template <typename MatrixType>
    inline void addToClosed (boost::shared_ptr<MatrixType> mat)
    {
        ASSERT (mat != 0, " Cannot assemble with an empty matrix");
        addToClosed (*mat);
    }
    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    IntegrateMatrixElement();

    //! Perform the computations for a single element
    /*!
     * This method computes the elemental matrix for a given element
     * index
     */
    void integrateElement  (const UInt iElement,
                            const UInt nbQuadPt,
                            const UInt nbTestDof,
                            const UInt nbSolutionDof,
                            ETMatrixElemental& elementalMatrix,
                            evaluation_Type& evaluation,
                            ETCurrentFE<MeshType::S_geoDimensions, 1>& globalCFE,
                            ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>& testCFE,
                            ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim>& solutionCFE);
    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;

    // Quadrature to be used
    QRAdapterType M_qrAdapter;

    // Shared pointer on the Spaces
    boost::shared_ptr<TestSpaceType> M_testSpace;
    boost::shared_ptr<SolutionSpaceType> M_solutionSpace;

    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;

    ETCurrentFE<MeshType::S_geoDimensions, 1>* M_globalCFE_std;
    ETCurrentFE<MeshType::S_geoDimensions, 1>* M_globalCFE_adapted;

    ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>* M_testCFE_std;
    ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>* M_testCFE_adapted;

    ETCurrentFE<TestSpaceType::space_dim, SolutionSpaceType::field_dim>* M_solutionCFE_std;
    ETCurrentFE<TestSpaceType::space_dim, SolutionSpaceType::field_dim>* M_solutionCFE_adapted;

    UInt M_offsetUp;
    UInt M_offsetLeft;

    // Data for multi-threaded assembly
    OpenMPParameters M_ompParams;

    // Data for integration on one subRegion, flag and elements on which perform the integration
    const UInt M_regionFlag;
    const UInt * const M_volumeElements;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
IntegrateMatrixElement (const boost::shared_ptr<MeshType>& mesh,
                        const QRAdapterType& qrAdapter,
                        const boost::shared_ptr<TestSpaceType>& testSpace,
                        const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                        const ExpressionType& expression,
                        const UInt offsetUp,
                        const UInt offsetLeft,
                        const UInt regionFlag,
                        const UInt * const volumeElements )
    :   M_mesh (mesh),
        M_qrAdapter (qrAdapter),
        M_testSpace (testSpace),
        M_solutionSpace (solutionSpace),
        M_evaluation (expression),

        M_testCFE_std (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),

        M_solutionCFE_std (new ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim> (solutionSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),
        M_solutionCFE_adapted (new ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim> (solutionSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),

        M_offsetUp (offsetUp),
        M_offsetLeft (offsetLeft),
        M_ompParams(),
        M_regionFlag( regionFlag ),
        M_volumeElements(volumeElements)
{
    switch (MeshType::geoShape_Type::BasRefSha::S_shape)
    {
        case LINE:
            M_globalCFE_std  = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted  = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) ;
            break;
        case TRIANGLE:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        case QUAD:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        case TETRA:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        case HEXA:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        default:
            ERROR_MSG ("Unrecognized element shape");
    }
    M_evaluation.setQuadrature (qrAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_std);
    M_evaluation.setTestCFE (M_testCFE_std);
    M_evaluation.setSolutionCFE (M_solutionCFE_std);
}


template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
IntegrateMatrixElement (const boost::shared_ptr<MeshType>& mesh,
                        const QRAdapterType& qrAdapter,
                        const boost::shared_ptr<TestSpaceType>& testSpace,
                        const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                        const ExpressionType& expression,
                        const OpenMPParameters& ompParams,
                        const UInt offsetUp,
                        const UInt offsetLeft,
                        const UInt regionFlag,
                        const UInt * const volumeElements )
    :   M_mesh (mesh),
        M_qrAdapter (qrAdapter),
        M_testSpace (testSpace),
        M_solutionSpace (solutionSpace),
        M_evaluation (expression),
        M_testCFE_std (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),

        M_solutionCFE_std (new ETCurrentFE<TestSpaceType::space_dim, SolutionSpaceType::field_dim> (solutionSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),
        M_solutionCFE_adapted (new ETCurrentFE<TestSpaceType::space_dim, SolutionSpaceType::field_dim> (solutionSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),

        M_offsetUp (offsetUp),
        M_offsetLeft (offsetLeft),
        M_ompParams (ompParams),
        M_regionFlag( regionFlag ),
        M_volumeElements(volumeElements)
{
    switch (MeshType::geoShape_Type::BasRefSha::S_shape)
    {
        case LINE:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        case TRIANGLE:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        case QUAD:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        case TETRA:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        case HEXA:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() );
            break;
        default:
            ERROR_MSG ("Unrecognized element shape");
    }
    M_evaluation.setQuadrature (qrAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_std);
    M_evaluation.setTestCFE (M_testCFE_std);
    M_evaluation.setSolutionCFE (M_solutionCFE_std);
}



template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
IntegrateMatrixElement (const IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>& integrator)
    :   M_mesh (integrator.M_mesh),
        M_qrAdapter (integrator.M_qrAdapter),
        M_testSpace (integrator.M_testSpace),
        M_solutionSpace (integrator.M_solutionSpace),
        M_evaluation (integrator.M_evaluation),

        M_testCFE_std (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (M_testSpace->refFE(), M_testSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),
        M_testCFE_adapted (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (M_testSpace->refFE(), M_testSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),

        M_solutionCFE_std (new ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim> (M_solutionSpace->refFE(), M_solutionSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),
        M_solutionCFE_adapted (new ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim> (M_solutionSpace->refFE(), M_solutionSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),

        M_offsetUp (integrator.M_offsetUp),
        M_offsetLeft (integrator.M_offsetLeft),

        M_ompParams (integrator.M_ompParams),
        M_regionFlag( integrator.M_regionFlag ),
        M_volumeElements( integrator.M_volumeElements )
{
    switch (MeshType::geoShape_Type::BasRefSha::S_shape)
    {
        case LINE:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feSegP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            break;
        case TRIANGLE:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            break;
        case QUAD:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            break;
        case TETRA:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            break;
        case HEXA:
            M_globalCFE_std = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            M_globalCFE_adapted = new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), integrator.M_qrAdapter.standardQR() );
            break;
        default:
            ERROR_MSG ("Unrecognized element shape");
    }
    M_evaluation.setQuadrature (integrator.M_qrAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_std);
    M_evaluation.setTestCFE (M_testCFE_std);
    M_evaluation.setSolutionCFE (M_solutionCFE_std);
}

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
~IntegrateMatrixElement()
{
    delete M_globalCFE_std;
    delete M_globalCFE_adapted;
    delete M_testCFE_std;
    delete M_testCFE_adapted;
    delete M_solutionCFE_std;
    delete M_solutionCFE_adapted;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
void
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
}

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
void
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
integrateElement (const UInt iElement, const UInt nbQuadPt,
                  const UInt nbTestDof, const UInt nbSolutionDof,
                  ETMatrixElemental& elementalMatrix,
                  evaluation_Type& evaluation,
                  ETCurrentFE<MeshType::S_geoDimensions, 1>& globalCFE,
                  ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>& testCFE,
                  ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim>& solutionCFE)
{
    // Zeros out the matrix
    elementalMatrix.zero();

    evaluation.update (iElement);

    // Update the currentFEs
    globalCFE.update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
    testCFE.update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);
    solutionCFE.update (M_mesh->element (iElement), evaluation_Type::S_solutionUpdateFlag);


    // Loop on the blocks

    for (UInt iblock (0); iblock < TestSpaceType::field_dim; ++iblock)
    {
        for (UInt jblock (0); jblock < SolutionSpaceType::field_dim; ++jblock)
        {

            // Set the row global indices in the local matrix
            for (UInt i (0); i < nbTestDof; ++i)
            {
                elementalMatrix.setRowIndex
                (i + iblock * nbTestDof,
                 M_testSpace->dof().localToGlobalMap (iElement, i) + iblock * M_testSpace->dof().numTotalDof() + M_offsetUp);
            }

            // Set the column global indices in the local matrix
            for (UInt j (0); j < nbSolutionDof; ++j)
            {
                elementalMatrix.setColumnIndex
                (j + jblock * nbSolutionDof,
                 M_solutionSpace->dof().localToGlobalMap (iElement, j) + jblock * M_solutionSpace->dof().numTotalDof() + M_offsetLeft);
            }

            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                for (UInt i (0); i < nbTestDof; ++i)
                {
                    for (UInt j (0); j < nbSolutionDof; ++j)
                    {
                        elementalMatrix.element (i + iblock * nbTestDof, j + jblock * nbSolutionDof) +=
                            evaluation.value_qij (iQuadPt, i + iblock * nbTestDof, j + jblock * nbSolutionDof)
                            * globalCFE.wDet (iQuadPt);

                    }
                }
            }
        }
    }
}


template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
template <typename MatrixType>
void
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
addTo (MatrixType& mat)
{
    UInt nbElements (M_mesh->numElements() );
    //UInt nbQuadPt_std (M_qrAdapter.standardQR().nbQuadPt() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );
    UInt nbSolutionDof (M_solutionSpace->refFE().nbDof() );

    // TODO: Shall these be members or not? I think it depends on OMP
    // globalCFE => if yes: move the above cases here
    // elementalMatrix => if no, add back the member and remove the variable
    // Including a copy of evaluation: if no, use the member

    ETMatrixElemental elementalMatrix (TestSpaceType::field_dim * M_testSpace->refFE().nbDof(),
                                       SolutionSpaceType::field_dim * M_solutionSpace->refFE().nbDof() );

    evaluation_Type evaluation (M_evaluation);

    // Defaulted to true for security
    bool isPreviousAdapted (true);

    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Extracting the marker
        UInt markerID = M_testSpace->mesh()->element ( iElement ).markerID( );

//        std::cout << "M_regionFlag is " << M_regionFlag << " markerID is " << markerID << " ID " << M_mesh->comm()->MyPID() << " element " << iElement << std::endl;


        if ( M_regionFlag == 0 )
        {
            markerID = M_regionFlag;
        }

        elementalMatrix.zero();

        if ( markerID == M_regionFlag )
        {

            // Update the quadrature rule adapter
            M_qrAdapter.update (iElement);

            if (M_qrAdapter.isAdaptedElement() )
            {
                // Set the quadrature rule everywhere
                evaluation.setQuadrature ( M_qrAdapter.adaptedQR() );
                M_globalCFE_adapted -> setQuadratureRule ( M_qrAdapter.adaptedQR() );
                M_testCFE_adapted -> setQuadratureRule ( M_qrAdapter.adaptedQR() );
                M_solutionCFE_adapted -> setQuadratureRule ( M_qrAdapter.adaptedQR() );

                // Reset the CurrentFEs in the evaluation
                evaluation.setGlobalCFE ( M_globalCFE_adapted );
                evaluation.setTestCFE ( M_testCFE_adapted );
                evaluation.setSolutionCFE ( M_solutionCFE_adapted );

                integrateElement (iElement, M_qrAdapter.adaptedQR().nbQuadPt(), nbTestDof, nbSolutionDof,
                                  elementalMatrix, evaluation, *M_globalCFE_adapted , //*globalCFE,
                                  *M_testCFE_adapted, *M_solutionCFE_adapted);

                isPreviousAdapted = true;

            }
            else
            {
                // Change in the evaluation if needed
                if (isPreviousAdapted)
                {
                    M_evaluation.setQuadrature ( M_qrAdapter.standardQR() );
                    M_evaluation.setGlobalCFE ( M_globalCFE_std );
                    M_evaluation.setTestCFE ( M_testCFE_std );
                    M_evaluation.setSolutionCFE ( M_solutionCFE_std );

                    isPreviousAdapted = false;
                }

                integrateElement (iElement, M_qrAdapter.standardQR().nbQuadPt(), nbTestDof, nbSolutionDof,
                                  elementalMatrix, evaluation, *M_globalCFE_std , //*globalCFE,
                                  *M_testCFE_std, *M_solutionCFE_std);

            }
        }

        elementalMatrix.pushToGlobal (mat);
    }
}

template < typename MeshType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType, typename QRAdapterType>
template <typename MatrixType>
void
IntegrateMatrixElement<MeshType, TestSpaceType, SolutionSpaceType, ExpressionType, QRAdapterType>::
addToClosed (MatrixType& mat)
{
    UInt nbElements (M_mesh->numElements() );
    //UInt nbQuadPt (M_qrAdapter.standardQR().nbQuadPt() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );
    UInt nbSolutionDof (M_solutionSpace->refFE().nbDof() );

    // OpenMP setup and pragmas around the loop
    M_ompParams.apply();

    #pragma omp parallel
    {
        QRAdapterType qrAdapter (M_qrAdapter);


        // Update the currentFEs
        boost::scoped_ptr<ETCurrentFE<MeshType::S_geoDimensions, 1> > globalCFE_std;
        boost::scoped_ptr<ETCurrentFE<MeshType::S_geoDimensions, 1> > globalCFE_adapted;

        switch (MeshType::geoShape_Type::BasRefSha::S_shape)
        {
            case LINE:
                globalCFE_std.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                globalCFE_adapted.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                break;
            case TRIANGLE:
                globalCFE_std.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                globalCFE_adapted.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTriaP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                break;
            case QUAD:
                globalCFE_std.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                globalCFE_adapted.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feQuadQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                break;
            case TETRA:
                globalCFE_std.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                globalCFE_adapted.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feTetraP0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                break;
            case HEXA:
                globalCFE_std.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                globalCFE_adapted.reset ( new ETCurrentFE<MeshType::S_geoDimensions, 1> (feHexaQ0, geometricMapFromMesh<MeshType>(), qrAdapter.standardQR() ) );
                break;
            default:
                ERROR_MSG ("Unrecognized element shape");
        }

        ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>
        testCFE_std (M_testSpace->refFE(), M_testSpace->geoMap(), M_qrAdapter.standardQR() );

        ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>
        testCFE_adapted (M_testSpace->refFE(), M_testSpace->geoMap(), M_qrAdapter.standardQR() );


        ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim>
        solutionCFE_std (M_solutionSpace->refFE(), M_testSpace->geoMap(),
                         M_qrAdapter.standardQR() );

        ETCurrentFE<SolutionSpaceType::space_dim, SolutionSpaceType::field_dim>
        solutionCFE_adapted (M_solutionSpace->refFE(), M_testSpace->geoMap(),
                             M_qrAdapter.standardQR() );

        evaluation_Type evaluation (M_evaluation);
        // Update the evaluation is done within the if statement
        /*
        evaluation.setQuadrature (M_qrAdapter.standardQR());
        evaluation.setGlobalCFE (&(*globalCFE));
        evaluation.setTestCFE (&testCFE);
        evaluation.setSolutionCFE (&solutionCFE);
        */

        ETMatrixElemental elementalMatrix (TestSpaceType::field_dim * M_testSpace->refFE().nbDof(),
                                           SolutionSpaceType::field_dim * M_solutionSpace->refFE().nbDof() );

        // Defaulted to true for security
        bool isPreviousAdapted (true);

        #pragma omp for schedule(runtime)
        for (UInt iElement (0); iElement < nbElements; ++iElement)
        {
            // Update the quadrature rule adapter
            qrAdapter.update (iElement);

            // TODO: move QRule choice inside a common method for AddTo and AddToClosed
            // TODO: Remove the members repeated here
            // TODO: use a policy to say if: 1) matrix open/closed (with graph) 2) with or without QR adapter

            if (qrAdapter.isAdaptedElement() )
            {
                // Set the quadrature rule everywhere
                evaluation.setQuadrature ( qrAdapter.adaptedQR() );
                globalCFE_adapted -> setQuadratureRule ( qrAdapter.adaptedQR() );
                testCFE_adapted.setQuadratureRule ( qrAdapter.adaptedQR() );
                solutionCFE_adapted. setQuadratureRule ( qrAdapter.adaptedQR() );

                // Reset the CurrentFEs in the evaluation
                evaluation.setGlobalCFE ( globalCFE_adapted.get() );
                evaluation.setTestCFE ( &testCFE_adapted );
                evaluation.setSolutionCFE ( &solutionCFE_adapted );

                integrateElement (iElement, qrAdapter.adaptedQR().nbQuadPt(), nbTestDof, nbSolutionDof,
                                  elementalMatrix, evaluation, *globalCFE_adapted ,
                                  testCFE_adapted, solutionCFE_adapted);

                isPreviousAdapted = true;

            }
            else
            {
                // Change in the evaluation if needed
                if (isPreviousAdapted)
                {
                    evaluation.setQuadrature ( qrAdapter.standardQR() );
                    evaluation.setGlobalCFE ( globalCFE_std.get() );
                    evaluation.setTestCFE ( &testCFE_std );
                    evaluation.setSolutionCFE ( &solutionCFE_std );

                    isPreviousAdapted = false;
                }

                integrateElement (iElement, M_qrAdapter.standardQR().nbQuadPt(), nbTestDof, nbSolutionDof,
                                  elementalMatrix, evaluation, *globalCFE_std ,
                                  testCFE_std, solutionCFE_std);

            }

            elementalMatrix.pushToGlobal (mat);
        }

        M_ompParams.restorePreviousNumThreads();
    }
}

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
