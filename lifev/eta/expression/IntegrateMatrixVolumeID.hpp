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
     @brief This file contains the definition of the IntegrateMatrixVolumeID class.

     @date 06/2011
     @author Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef INTEGRATE_MATRIX_VOLUME_ID_HPP
#define INTEGRATE_MATRIX_VOLUME_ID_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>

#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <lifev/eta/array/ETMatrixElemental.hpp>

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
template < typename MeshType, typename FunctorType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
class IntegrateMatrixVolumeID
{
public:

	//! @name Public Types
    //@{

    //! Type of the Evaluation
	typedef typename ExpressionToEvaluation< ExpressionType,
                                             TestSpaceType::field_dim,
                                             SolutionSpaceType::field_dim,
                                             3>::evaluation_Type evaluation_Type;

    typedef typename MeshType::element_Type       elementMesh_Type;
    typedef std::vector<elementMesh_Type const* > vectorVolumes_Type;
    typedef boost::shared_ptr<vectorVolumes_Type> vectorVolumesPtr_Type;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full data constructor
	IntegrateMatrixVolumeID(const boost::shared_ptr<MeshType>& mesh,
                            const boost::shared_ptr<FunctorType>& functorSelection,
                            const QuadratureRule& quadrature,
                            const boost::shared_ptr<TestSpaceType>& testSpace,
                            const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                            const ExpressionType& expression);

    //! Copy constructor
	IntegrateMatrixVolumeID( const IntegrateMatrixVolumeID < MeshType, FunctorType, TestSpaceType, SolutionSpaceType, ExpressionType> & integrator);

    //! Destructor
    ~IntegrateMatrixVolumeID();

    //@}


    //! @name Operator
    //@{

    //! Operator wrapping the addTo method
    template <typename MatrixType>
    inline void operator>>(MatrixType& mat)
    {
        addTo(mat);
    }

    //! Operator wrapping the addTo method (for shared_ptr)
    template <typename MatrixType>
    inline void operator>>(boost::shared_ptr<MatrixType> mat)
    {
        addTo(mat);
    }

    //@}


    //! @name Methods
    //@{

    //! Ouput method
	void check(std::ostream& out = std::cout);

    //! Method that performs the assembly
    /*!
      The loop over the elements is located right
      in this method. Everything for the assembly is then
      performed: update the values, update the local vector,
      sum over the quadrature nodes, assemble in the global
      vector.
     */
	template <typename MatrixType>
	void addTo(MatrixType& mat);

	template <typename MatrixType>
	inline void addTo(boost::shared_ptr<MatrixType> mat)
    {
        ASSERT(mat!=0, " Cannot assemble with an empty matrix");
        addTo(*mat);
    }

    //@}

private:

    //! @name Private Methods
    //@{

    // No default constructor
	IntegrateMatrixVolumeID();

    //@}

    // Pointer on the mesh
	boost::shared_ptr<MeshType> M_mesh;

    // Pointer to the selectionObject
    boost::shared_ptr<FunctorType> M_functorSelection;
    // Quadrature to be used
	QuadratureRule M_quadrature;

    // Shared pointer on the Space
	boost::shared_ptr<TestSpaceType> M_testSpace;
	boost::shared_ptr<SolutionSpaceType> M_solutionSpace;

    // Tree to compute the values for the assembly
	evaluation_Type M_evaluation;

	ETCurrentFE<3,1>* M_globalCFE;
	ETCurrentFE<3,TestSpaceType::field_dim>* M_testCFE;
	ETCurrentFE<3,SolutionSpaceType::field_dim>* M_solutionCFE;

    ETMatrixElemental M_elementalMatrix;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename FunctorType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixVolumeID<MeshType,FunctorType,TestSpaceType,SolutionSpaceType,ExpressionType>::
IntegrateMatrixVolumeID(const boost::shared_ptr<MeshType>& mesh,
                        const boost::shared_ptr<FunctorType>& functorSelection,
                        const QuadratureRule& quadrature,
                        const boost::shared_ptr<TestSpaceType>& testSpace,
                        const boost::shared_ptr<SolutionSpaceType>& solutionSpace,
                        const ExpressionType& expression)
    :	M_mesh( mesh ),
        M_functorSelection( functorSelection ),
        M_quadrature(quadrature),
        M_testSpace(testSpace),
        M_solutionSpace(solutionSpace),
        M_evaluation(expression),

        M_globalCFE(new ETCurrentFE<3,1>(feTetraP0,geometricMapFromMesh<MeshType>(),quadrature)),
        M_testCFE(new ETCurrentFE<3,TestSpaceType::field_dim>(testSpace->refFE(),testSpace->geoMap(),quadrature)),
        M_solutionCFE(new ETCurrentFE<3,SolutionSpaceType::field_dim>(solutionSpace->refFE(),testSpace->geoMap(),quadrature)),

        M_elementalMatrix(TestSpaceType::field_dim*testSpace->refFE().nbDof(),
                          SolutionSpaceType::field_dim*solutionSpace->refFE().nbDof())
{
    M_evaluation.setQuadrature(quadrature);
    M_evaluation.setGlobalCFE(M_globalCFE);
    M_evaluation.setTestCFE(M_testCFE);
    M_evaluation.setSolutionCFE(M_solutionCFE);
}

template < typename MeshType, typename FunctorType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixVolumeID<MeshType,FunctorType, TestSpaceType,SolutionSpaceType,ExpressionType>::
IntegrateMatrixVolumeID(const IntegrateMatrixVolumeID<MeshType, FunctorType,TestSpaceType,SolutionSpaceType,ExpressionType>& integrator)
    :	M_mesh(integrator.M_mesh),
        M_functorSelection(integrator.M_functorSelection),
        M_quadrature(integrator.M_quadrature),
        M_testSpace(integrator.M_testSpace),
        M_solutionSpace(integrator.M_solutionSpace),
        M_evaluation(integrator.M_evaluation),

        M_globalCFE(new ETCurrentFE<3,1>(feTetraP0,geometricMapFromMesh<MeshType>(),M_quadrature)),
        M_testCFE(new ETCurrentFE<3,TestSpaceType::field_dim>(M_testSpace->refFE(), M_testSpace->geoMap(),M_quadrature)),
        M_solutionCFE(new ETCurrentFE<3,SolutionSpaceType::field_dim>(M_solutionSpace->refFE(), M_solutionSpace->geoMap(),M_quadrature)),

        M_elementalMatrix(integrator.M_elementalMatrix)
{
    M_evaluation.setQuadrature(M_quadrature);
    M_evaluation.setGlobalCFE(M_globalCFE);
    M_evaluation.setTestCFE(M_testCFE);
    M_evaluation.setSolutionCFE(M_solutionCFE);
}

template < typename MeshType, typename FunctorType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
IntegrateMatrixVolumeID<MeshType,FunctorType,TestSpaceType,SolutionSpaceType,ExpressionType>::
~IntegrateMatrixVolumeID()
{
    delete M_globalCFE;
    delete M_testCFE;
    delete M_solutionCFE;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename FunctorType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
void
IntegrateMatrixVolumeID<MeshType,FunctorType, TestSpaceType,SolutionSpaceType,ExpressionType>::
check(std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display(out);
    out << std::endl;
    out << " Elemental matrix : " << std::endl;
    M_elementalMatrix.showMe(out);
    out << std::endl;
}


template < typename MeshType, typename FunctorType, typename TestSpaceType, typename SolutionSpaceType, typename ExpressionType>
template <typename MatrixType>
void
IntegrateMatrixVolumeID<MeshType, FunctorType, TestSpaceType,SolutionSpaceType,ExpressionType>::
addTo(MatrixType& mat)
{

    //First: extract the list of volumes over which we integrate according to the
    //       functor that is received

    UInt numExtractedVolumes = this->M_mesh->elementList().countAccordingToPredicate( *(this->M_functorSelection) );
    vectorVolumes_Type extractedVolumes( numExtractedVolumes );

    //Extracting the volumes
    extractedVolumes = this->M_mesh->elementList().extractAccordingToPredicate( *(this->M_functorSelection) );

    //number of volumes
    UInt nbElements( extractedVolumes.size() );
    UInt nbQuadPt(M_quadrature.nbQuadPt());
    UInt nbTestDof(M_testSpace->refFE().nbDof());
    UInt nbSolutionDof(M_solutionSpace->refFE().nbDof());

    for (UInt iElement(0); iElement< nbElements; ++iElement)
    {
        // Zeros out the matrix
        M_elementalMatrix.zero();

        // Update the currentFEs need to access to the current element of the loop
        M_globalCFE->update( *extractedVolumes[iElement] ,evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
        M_testCFE->update( *extractedVolumes[iElement] ,evaluation_Type::S_testUpdateFlag);
        M_solutionCFE->update( *extractedVolumes[iElement] ,evaluation_Type::S_solutionUpdateFlag);

        // Update the evaluation
        M_evaluation.update(iElement);

        // Loop on the blocks

        for (UInt iblock(0); iblock < TestSpaceType::field_dim; ++iblock)
        {
            for (UInt jblock(0); jblock < SolutionSpaceType::field_dim; ++jblock)
            {

                // Set the row global indices in the local matrix
                for (UInt i(0); i<nbTestDof; ++i)
                {
                    M_elementalMatrix.setRowIndex
						(i+iblock*nbTestDof,
                         M_testSpace->dof().localToGlobalMap(iElement,i)+ iblock*M_testSpace->dof().numTotalDof());
                }

                // Set the column global indices in the local matrix
                for (UInt j(0); j<nbSolutionDof; ++j)
                {
                    M_elementalMatrix.setColumnIndex
						(j+jblock*nbSolutionDof,
                         M_solutionSpace->dof().localToGlobalMap(iElement,j)+ jblock*M_solutionSpace->dof().numTotalDof());
                }

                for (UInt iQuadPt(0); iQuadPt< nbQuadPt; ++iQuadPt)
                {
                    for (UInt i(0); i<nbTestDof; ++i)
                    {
                        for (UInt j(0); j<nbSolutionDof; ++j)
                        {
                            M_elementalMatrix.element(i+iblock*nbTestDof,j+jblock*nbSolutionDof) +=
								M_evaluation.value_qij(iQuadPt,i+iblock*nbTestDof,j+jblock*nbSolutionDof)
								* M_globalCFE->wDet(iQuadPt);

                        }
                    }
                }
            }
        }

        M_elementalMatrix.pushToGlobal(mat);
    }
}

} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
