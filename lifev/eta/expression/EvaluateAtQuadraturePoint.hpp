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
 @brief This file contains the evaluation of vectors at quadrature points.
 
 @date 06/2011
 @author Davide Forti <davide.forti@epfl.ch>
 */

#ifndef EVALUATE_QUAD_PT_HPP
#define EVALUATE_QUAD_PT_HPP

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>
#include <lifev/eta/fem/QRAdapterBase.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>
#include <lifev/eta/expression/EvaluationPhiI.hpp>

#include <lifev/eta/array/ETVectorElemental.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{
    
namespace ExpressionAssembly
{

//! The class to actually perform the loop over the elements to compute the stresses at the center of each element
/*!
 @author Davide Forti <davide.forti@epfl.ch>
 
 */
template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
class EvaluateAtQuadraturePoint
{
public:
    
    //! @name Public Types
    //@{
    
    //! Type of the Evaluation
    typedef typename ExpressionToEvaluation < ExpressionType,
    TestSpaceType::field_dim,
    0,
    MeshType::S_geoDimensions >::evaluation_Type evaluation_Type;
    
    //@}
    
    
    //! @name Constructors, destructor
    //@{
    
    //! Full data constructor
    EvaluateAtQuadraturePoint (const boost::shared_ptr<MeshType>& mesh,
                            const QRAdapterType& qrAdapter,
                            const boost::shared_ptr<TestSpaceType>& testSpace,
                            const ExpressionType& expression,
                            const UInt offset = 0);
    
    //! Copy constructor
    EvaluateAtQuadraturePoint ( const IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType, QRAdapterType>& integrator);
    
    //! Destructor
    ~EvaluateAtQuadraturePoint();
    
    //@}
    
    
    //! @name Operator
    //@{
    
    //! Operator wrapping the addTo method
    inline void operator>> (std::vector<std::vector<Real>>& vector)
    {
        evaluate (vector);
    }

    //! Operator wrapping the addTo method
    inline void operator>> (std::vector<std::vector<VectorSmall<3>>>& vector)
    {
        evaluate (vector);
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
    void evaluate (std::vector<std::vector<Real>>& vec); // for scalar field
    
    void evaluate (std::vector<std::vector<VectorSmall<3>>>& vec); // for 3D field
    
    //@}
    
private:
    
    //! @name Private Methods
    //@{
    
    // No default constructor
    EvaluateAtQuadraturePoint();
    
    //@}
    
    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;
    
    // Quadrature to be used
    QRAdapterType M_qrAdapter;
    
    // Shared pointer on the Space
    boost::shared_ptr<TestSpaceType> M_testSpace;
    
    // Tree to compute the values for the assembly
    evaluation_Type M_evaluation;
    
    ETCurrentFE<MeshType::S_geoDimensions, 1>* M_globalCFE_std;
    ETCurrentFE<MeshType::S_geoDimensions, 1>* M_globalCFE_adapted;
    
    ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>* M_testCFE_std;
    ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim>* M_testCFE_adapted;
    
    ETVectorElemental M_elementalVector;
    
    // Offset
    UInt M_offset;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
EvaluateAtQuadraturePoint < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
EvaluateAtQuadraturePoint (const boost::shared_ptr<MeshType>& mesh,
                        const QRAdapterType& qrAdapter,
                        const boost::shared_ptr<TestSpaceType>& testSpace,
                        const ExpressionType& expression,
                        const UInt offset)
:   M_mesh (mesh),
M_qrAdapter (qrAdapter),
M_testSpace (testSpace),
M_evaluation (expression),

M_testCFE_std (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),
M_testCFE_adapted (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (testSpace->refFE(), testSpace->geoMap(), qrAdapter.standardQR() ) ),

M_elementalVector (TestSpaceType::field_dim * testSpace->refFE().nbDof() ),

M_offset (offset)
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
    M_evaluation.setQuadrature ( qrAdapter.standardQR() );
    
    M_evaluation.setGlobalCFE (M_globalCFE_std);
    M_evaluation.setTestCFE (M_testCFE_std);
}


template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
EvaluateAtQuadraturePoint < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
EvaluateAtQuadraturePoint ( const IntegrateVectorElement < MeshType, TestSpaceType, ExpressionType, QRAdapterType>& integrator)
:   M_mesh (integrator.M_mesh),
M_qrAdapter (integrator.M_qrAdapter),
M_testSpace (integrator.M_testSpace),
M_evaluation (integrator.M_evaluation),

M_testCFE_std (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (M_testSpace->refFE(), M_testSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),
M_testCFE_adapted (new ETCurrentFE<TestSpaceType::space_dim, TestSpaceType::field_dim> (M_testSpace->refFE(), M_testSpace->geoMap(), integrator.M_qrAdapter.standardQR() ) ),

M_elementalVector (integrator.M_elementalVector),
M_offset (integrator.M_offset)
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
}


template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
EvaluateAtQuadraturePoint < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
~EvaluateAtQuadraturePoint()
{
    delete M_globalCFE_std;
    delete M_globalCFE_adapted;
    delete M_testCFE_std;
    delete M_testCFE_adapted;
}

// ===================================================
// Methods
// ===================================================

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
void
EvaluateAtQuadraturePoint < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
check (std::ostream& out)
{
    out << " Checking the integration : " << std::endl;
    M_evaluation.display (out);
    out << std::endl;
    out << " Elemental vector : " << std::endl;
    M_elementalVector.showMe (out);
    out << std::endl;
}

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
void
EvaluateAtQuadraturePoint < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
evaluate (std::vector<std::vector<Real>>& vec)
{
    UInt nbElements (M_mesh->numElements() );
    UInt nbQuadPt_std (M_qrAdapter.standardQR().nbQuadPt() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );
    
    // Defaulted to true for security
    bool isPreviousAdapted (true);
    
    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        
        // Update the quadrature rule adapter
        M_qrAdapter.update (iElement);
        
        // Check if the last one was adapted
        if (isPreviousAdapted)
        {
            M_evaluation.setQuadrature ( M_qrAdapter.standardQR() );
            M_evaluation.setGlobalCFE ( M_globalCFE_std );
            M_evaluation.setTestCFE ( M_testCFE_std );
            
            isPreviousAdapted = false;
        }
        
        // Update the currentFEs
        M_globalCFE_std->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
        M_testCFE_std->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);
        
        // Update the evaluation
        M_evaluation.update (iElement);
        
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt_std; ++iQuadPt)
        {
            vec[iElement][iQuadPt] = M_evaluation.value_qi (iQuadPt, 0);
        }
    }
    
}

template < typename MeshType, typename TestSpaceType, typename ExpressionType, typename QRAdapterType>
void
EvaluateAtQuadraturePoint < MeshType, TestSpaceType, ExpressionType, QRAdapterType>::
evaluate (std::vector<std::vector<VectorSmall<3>>>& vec)
{
    UInt nbElements (M_mesh->numElements() );
    UInt nbQuadPt_std (M_qrAdapter.standardQR().nbQuadPt() );
    UInt nbTestDof (M_testSpace->refFE().nbDof() );
    
    // Defaulted to true for security
    bool isPreviousAdapted (true);
    
    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        // Update the quadrature rule adapter
        M_qrAdapter.update (iElement);
        
        // Check if the last one was adapted
        if (isPreviousAdapted)
        {
            M_evaluation.setQuadrature ( M_qrAdapter.standardQR() );
            M_evaluation.setGlobalCFE ( M_globalCFE_std );
            M_evaluation.setTestCFE ( M_testCFE_std );
            
            isPreviousAdapted = false;
        }
        
        // Update the currentFEs
        M_globalCFE_std->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);
        M_testCFE_std->update (M_mesh->element (iElement), evaluation_Type::S_testUpdateFlag);
        
        // Update the evaluation
        M_evaluation.update (iElement);
        
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt_std; ++iQuadPt)
        {
            vec[iElement][iQuadPt] = M_evaluation.value_qi (iQuadPt , 1 );
        }
    }
    
}
    
    
} // Namespace ExpressionAssembly
    
} // Namespace LifeV

#endif
