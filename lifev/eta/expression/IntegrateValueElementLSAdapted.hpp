//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief File containing the class IntegrateValueElementLSAdapted

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 21 Dec 2011


 */

#ifndef INTEGRATEVALUEELEMENTLSADAPTED_H
#define INTEGRATEVALUEELEMENTLSADAPTED_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/level_set/fem/LevelSetQRAdapter.hpp>
#include <lifev/eta/fem/ETCurrentFE.hpp>
#include <lifev/eta/fem/MeshGeometricMap.hpp>

#include <lifev/eta/expression/ExpressionToEvaluation.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

namespace ExpressionAssembly
{

//! IntegrateValueElementLSAdapted - Class to integrate a value using a quadrature conforming to a level set
/*!
    @author Samuel Quinodoz
    @see Reference to papers (if available)

 */

template< typename MeshType, typename ExpressionType, typename LSFESpaceType, typename VectorType>
class IntegrateValueElementLSAdapted
{
public:

    //! @name Public Types
    //@{

    typedef typename ExpressionToEvaluation < ExpressionType,
            0,
            0,
            3 >::evaluation_Type evaluation_Type;

    typedef LevelSetQRAdapter<LSFESpaceType, VectorType> QRAdapter_Type;


    //@}


    //! @name Constructor & Destructor
    //@{

    IntegrateValueElementLSAdapted ( const std::shared_ptr<MeshType>& mesh,
                                     const QRAdapter_Type& QRAdapter,
                                     const ExpressionType& expression);

    IntegrateValueElementLSAdapted ( const IntegrateValueElementLSAdapted
                                     <MeshType, ExpressionType, LSFESpaceType, VectorType>& integrator);

    ~IntegrateValueElementLSAdapted();


    //@}


    //! @name Operators
    //@{

    inline void operator>> (Real& value)
    {
        addTo (value);
    }

    //@}


    //! @name Methods
    //@{

    void addTo (Real& value);

    //@}

private:

    //! @name Private Methods
    //@{

    IntegrateValueElementLSAdapted();

    //@}

    // Mesh
    std::shared_ptr<MeshType> M_mesh;

    // Quadrature Adapter
    QRAdapter_Type M_QRAdapter;

    // Evaluation tree
    evaluation_Type M_evaluation;

    // CurrentFE
    ETCurrentFE<3, 1>* M_globalCFE_unadapted;
    ETCurrentFE<3, 1>* M_globalCFE_adapted;

};


template< typename MeshType, typename ExpressionType, typename LSFESpaceType, typename VectorType>
IntegrateValueElementLSAdapted<MeshType, ExpressionType, LSFESpaceType, VectorType>::
IntegrateValueElementLSAdapted ( const std::shared_ptr<MeshType>& mesh,
                                 const QRAdapter_Type& QRAdapter,
                                 const ExpressionType& expression)
    :
    M_mesh (mesh),
    M_QRAdapter (QRAdapter),
    M_evaluation (expression),

    M_globalCFE_unadapted ( new ETCurrentFE<3, 1> ( feTetraP0, geometricMapFromMesh<MeshType>(), QRAdapter.standardQR() ) ),
    M_globalCFE_adapted ( new ETCurrentFE<3, 1> ( feTetraP0, geometricMapFromMesh<MeshType>(), QRAdapter.standardQR() ) )
{
    M_evaluation.setQuadrature ( QRAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_unadapted);
}

template< typename MeshType, typename ExpressionType, typename LSFESpaceType, typename VectorType>
IntegrateValueElementLSAdapted<MeshType, ExpressionType, LSFESpaceType, VectorType>::
IntegrateValueElementLSAdapted ( const IntegrateValueElementLSAdapted
                                 <MeshType, ExpressionType, LSFESpaceType, VectorType>& integrator)
    :
    M_mesh ( integrator.M_mesh ),
    M_QRAdapter ( integrator.M_QRAdapter ),
    M_evaluation ( integrator.M_evaluation ),

    M_globalCFE_unadapted ( new ETCurrentFE<3, 1> ( feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) ),
    M_globalCFE_adapted ( new ETCurrentFE<3, 1> ( feTetraP0, geometricMapFromMesh<MeshType>(), M_QRAdapter.standardQR() ) )
{
    M_evaluation.setQuadrature (M_QRAdapter.standardQR() );
    M_evaluation.setGlobalCFE (M_globalCFE_unadapted);
}


template< typename MeshType, typename ExpressionType, typename LSFESpaceType, typename VectorType>
IntegrateValueElementLSAdapted<MeshType, ExpressionType, LSFESpaceType, VectorType>::
~IntegrateValueElementLSAdapted()
{
    delete M_globalCFE_unadapted;
    delete M_globalCFE_adapted;
}

template< typename MeshType, typename ExpressionType, typename LSFESpaceType, typename VectorType>
void
IntegrateValueElementLSAdapted<MeshType, ExpressionType, LSFESpaceType, VectorType>::
addTo (Real& value)
{
    const UInt nbElements (M_mesh->numElements() );
    const UInt nbQuadPt_unadapted (M_QRAdapter.standardQR().nbQuadPt() );

    // Even if this is not the case, this
    // prevents many potentially buggy behaviours
    // and does not harm.
    bool isPreviousAdapted (true);


    for (UInt iElement (0); iElement < nbElements; ++iElement)
    {
        M_QRAdapter.update (iElement);

        if ( M_QRAdapter.isAdaptedElement() )
        {
            // Set the new QR everywhere
            M_evaluation.setQuadrature ( M_QRAdapter.adaptedQR() );
            M_globalCFE_adapted->setQuadratureRule ( M_QRAdapter.adaptedQR() );

            M_evaluation.setGlobalCFE ( M_globalCFE_adapted );

            // Update the currentFEs
            M_globalCFE_adapted->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);

            // Update the evaluation
            M_evaluation.update (iElement);

            // Make the assembly
            for (UInt iQuadPt (0); iQuadPt < M_QRAdapter.adaptedQR().nbQuadPt(); ++iQuadPt)
            {
                value += M_evaluation.value_q (iQuadPt)
                         * M_globalCFE_adapted->wDet (iQuadPt);
            }

            isPreviousAdapted = true;
        }
        else
        {
            // if the previous element was adapted, reset the QR
            if (isPreviousAdapted)
            {
                M_evaluation.setQuadrature ( M_QRAdapter.standardQR() );
                M_evaluation.setGlobalCFE ( M_globalCFE_unadapted );
            }

            // Update the currentFEs
            M_globalCFE_unadapted->update (M_mesh->element (iElement), evaluation_Type::S_globalUpdateFlag | ET_UPDATE_WDET);

            // Update the evaluation
            M_evaluation.update (iElement);


            // Make the assembly
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt_unadapted; ++iQuadPt)
            {
                value += M_evaluation.value_q (iQuadPt)
                         * M_globalCFE_unadapted->wDet (iQuadPt);
            }

            isPreviousAdapted = false;
        }
    }
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif /* INTEGRATEVALUEELEMENTLSADAPTED_H */
