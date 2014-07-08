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
     @brief This file contains the definition of the ETCurrentFE.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef ETCURRENTFE_HPP
#define ETCURRENTFE_HPP


#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/eta/fem/ETCurrentFlag.hpp>

#include <lifev/core/fem/GeometricMap.hpp>

#include <lifev/core/fem/ReferenceFEScalar.hpp>
#include <lifev/core/fem/ReferenceFEHdiv.hpp>
#include <lifev/core/fem/ReferenceFEHybrid.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>


#include <vector>

namespace LifeV
{

namespace ExpressionAssembly
{
/*
   Predeclaration of the classes that are friend to this class.
   This is not mandatory with the newest standards (C++03 and later?)
   but is required for some earlier standards. So, better to
   have them.
*/

template <UInt dim>
class EvaluationPhiI;

template <UInt dim>
class EvaluationPhiJ;

template <UInt dim, UInt FSpaceDim>
class EvaluationDphiI;

template <UInt dim, UInt FSpaceDim>
class EvaluationDphiJ;

template <UInt dim, UInt FSpaceDim>
class EvaluationLaplacianPhiJ;
    
template <UInt dim, UInt FSpaceDim>
class EvaluationDivI;

template <UInt dim, UInt FSpaceDim>
class EvaluationDivJ;

template <UInt dim>
class EvaluationHK;

template <UInt dim>
class EvaluationPosition;

template <UInt dim>
class EvaluationMeas;

} // Namespace Expression Assembly

/*!
  ETCurrenteFE is a template class. If fieldDim the general
  case is treated as representing a vectorial FE (only the case
  where fieldDim represents a scalar FE, but this is a partial
  specialization of this class).

*/
// This header contains the non-specialized version of the class
#include "ETCurrentFE_FD3.hpp"

/*!
  The ETCurrentFE is the class to be used to compute quantities related
  to the quadrature and the basis functions in a real element (current element).

  These quantities are computed using the update method.
  The update is managed through a system of flags that have to be passed to that method
  to require some of the quantities to be computed.

  The access to the data, getters have to be used (excepted for some specific classes
  which have a direct access to the stored data through friendship).

  This is the partial specialization of the generic ETCurrentFE for scalar finite element spaces.

*/
template <UInt spaceDim>
class ETCurrentFE<spaceDim, 1>
{

    //! @name Friends
    //@{

    //!Friend to allow direct access to the raw data
    template <UInt dim>
    friend class ExpressionAssembly::EvaluationPhiI;

    //!Friend to allow direct access to the raw data
    template <UInt dim>
    friend class ExpressionAssembly::EvaluationPhiJ;

    //!Friend to allow direct access to the raw data
    template <UInt dim, UInt FSpaceDim>
    friend class ExpressionAssembly::EvaluationDphiI;

    //!Friend to allow direct access to the raw data
    template <UInt dim, UInt FSpaceDim>
    friend class ExpressionAssembly::EvaluationDphiJ;
    
    //!Friend to allow direct access to the raw data
    template <UInt dim, UInt FSpaceDim>
    friend class ExpressionAssembly::EvaluationLaplacianPhiJ;

    //!Friend to allow direct access to the raw data
    template <UInt dim>
    friend class ExpressionAssembly::EvaluationHK;

    //!Friend to allow direct access to the raw data
    template <UInt dim>
    friend class ExpressionAssembly::EvaluationPosition;

    //!Friend to allow direct access to the raw data
    template <UInt dim>
    friend class ExpressionAssembly::EvaluationMeas;


    //@}

private:
    
    // Vector return type for phi
    typedef VectorSmall< spaceDim > array1D_Return_Type;
    
public:

    //! @name Static constants
    //@{

    //! Static value for the space dimension
    static const UInt S_spaceDimension; // = spaceDim;

    //! Static value for the dimension of the field (1 here as it is a scalar FE)
    static const UInt S_fieldDimension; // = 1;

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full constructor
    /*!
      @param refFE The reference element for the FE
      @param geoMap The geometric map from the reference element to the current element
      @param qr The quadrature rule
     */
    ETCurrentFE (const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr);

    //! Constructor without quadrature rule
    /*!
      @param refFE The reference element for the FE
      @param geoMap The geometric map from the reference element to the current element
     */
    ETCurrentFE (const ReferenceFE& refFE, const GeometricMap& geoMap);

    //! Copy constructor
    /*!
      @param otherFE The currentFE to be copied
     */
    ETCurrentFE (const ETCurrentFE<spaceDim, 1>& otherFE);

    //! Destructor
    virtual ~ETCurrentFE();

    //@}


    //! @name Methods
    //@{

    //! Update method
    /*!
      The update method computes the required quantities (indicated by the flag,
      see in \ref Update_flag) on the element.

      @param element The current element were to compute the values
      @param flag The flag indicating the quantities to update

      <b>Template parameters</b>

      <i>elementType</i>: The type of the element

      <b>Template requirements</b>

      <i>elementType</i>: has a id() method returning the identifier (global); has a localId() method returning
      the identifier (local); has a method point(UInt i) that returns the ith point.
     */
    template<typename elementType>
    void update (const elementType& element, const flag_Type& flag);

    //! ShowMe method
    /*!
      @param out Output stream were to print the informations
     */
    void showMe (std::ostream& out = std::cout) const;

    //@}


    //! @name Set Methods
    //@{

    //! Setter for the quadrature rule
    /*!
      Beware that when this setter is called, all the internal
      structures have to be reshaped (according to the number
      of quadrature nodes).

      @param qr The new quadrature to use.
     */
    void setQuadratureRule (const QuadratureRule& qr);

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the number of degrees of freedom of this element
    /*!
      @return The number of number of degrees of freedom of this element
     */
    UInt nbFEDof() const
    {
        return M_nbFEDof;
    }

    //! Getter for the values of the basis functions in the quadrature nodes in the current element
    /*!
      @param i The index of the basis function
      @param q The index of the quadrature node
      @return The value of the ith basis function in the qth quadrature node
     */
    Real phi (const UInt& i, const UInt& q) const
    {
        ASSERT ( i < M_nbFEDof, "No basis function with this index");
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index");
        return M_phi[q][i];
    }

    //! Getter for the coordinates of the quadrature nodes in the current element
    /*!
      @param q The index of the quadrature node
      @param coord The coordinate required
      @return The coordth component of the coordinate of the qth quadrature point
     */
    Real quadNode (const UInt& q, const UInt& coord) const
    {
        ASSERT ( M_isQuadNodeUpdated, "Quadrature nodes have not been updated");
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index");
        ASSERT ( coord < spaceDim, "No such coordinate index");
        return M_quadNode[q][coord];
    }

    //! Getter for the weighted jacobian determinant (current element)
    /*!
      @param q The index of the quadrature point
      @return The weighted determinant of the jacobian transform in the qth quadrature node
     */
    Real wDet (const UInt& q) const
    {
        ASSERT ( M_isWDetUpdated, "Weighted determinant has not been updated");
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index");
        return M_wDet[q];
    }

    //! Getter for the derivatives of the basis function in the quadrature nodes (current element)
    /*!
      @param i The index of the basis function
      @param dxi The direction of the derivative required (0 for d/dx, 1 for d/dy...)
      @param q The index of the quadrature node
      @return The value of the ith basis function derived w.r. to dxi, in the qth quadrature node.
     */
    Real dphi (const UInt& i, const UInt& dxi, const UInt& q) const
    {
        ASSERT ( M_isDphiUpdated, "Derivative of the basis functions have not been updated");
        ASSERT ( i < M_nbFEDof, "No basis function with this index");
        ASSERT ( dxi < spaceDim, "No such coordinate index");
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index");
        return M_dphi[q][i][dxi];
    }

    //! Getter for the derivatives of the basis function in the quadrature nodes (current element)
    /*!
      @param i The index of the basis function
      @param q The index of the quadrature node
      @return The local vector of the basis functions derived w.r. to dxi, in the qth quadrature node.
     */
    VectorSmall<spaceDim> const& dphi (const UInt& i, const UInt& q) const
    {
        ASSERT ( M_isDphiUpdated, "Derivative of the basis functions have not been updated");
        ASSERT ( i < M_nbFEDof, "No basis function with this index");
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index");
        return M_dphi[q][i];
    }

    //! Getter for the laplacian of the basis functions in the quadrature nodes in the current element
    /*!
     @param i The index of the basis function
     @param q The index of the quadrature node
     @return The laplacian of the ith basis function in the qth quadrature node
     */
    
    const VectorSmall<spaceDim> laplacian (const UInt& i, const UInt& q) const
    {
        ASSERT ( M_isLaplacianUpdated, "Divergence of the basis functions have not been updated");
        ASSERT ( i < M_nbFEDof, "No basis function with this index" );
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index" );
        
        return ( M_laplacian[q][i] );
    }
    
    //! Getter for the identifier of the current element
    /*!
      @return The (global) identifier of the current element.
     */
    UInt currentId() const
    {
        ASSERT (M_isCellNodeUpdated, "Cell has not been updated");
        return M_currentId;
    }

    //! Getter for the diameter of the current element
    /*!
      @return The diameter of the current element
    */
    Real diameter() const
    {
        ASSERT (M_isDiameterUpdated, "Diameter has not been updated");
        return M_diameter;
    }

    //! Getter for the measure of the current element
    /*!
      @return The measure of the current element
    */
    Real measure() const
    {
        //ASSERT(M_isMeasure, "Measure has not been updated");

        return M_measure;
    }


    //@}

private:
    
    //Private typedefs for the 1D array
    typedef std::vector< Real > array1D_Type;

    //Private typedefs for the 2D array (array of 1D array)
    typedef std::vector< array1D_Type > array2D_Type;

    //Private typedefs for the 3D array (array of 2D array)
    typedef std::vector< array2D_Type > array3D_Type;
    
    //Private typedefs for the 4D array (array of 3D array)
    typedef std::vector< array3D_Type > array4D_Type;
    
    //Private typedefs for the 1D array of vector
    typedef std::vector< VectorSmall<spaceDim> > array1D_vector_Type;

    //Private typedefs for the 2D array of vector
    typedef std::vector< std::vector< VectorSmall<spaceDim> > > array2D_vector_Type;
    
    //To contain the second derivatives
    typedef std::vector< std::vector< MatrixSmall<spaceDim,spaceDim> > > array2D_matrix_Type;

    //Private typedefs for the 3D array of vector
    typedef std::vector< std::vector< std::vector< VectorSmall<spaceDim> > > > array3D_vector_Type;
    
    //! @name Private Methods
    //@{

    //! No default constructor
    ETCurrentFE();

    //! No assignement
    void operator= (const ETCurrentFE<spaceDim, 1>&);

    //! Resize all the internal containers w.r. to the stored data and compute the constant values
    void setupInternalConstants();

    //! Update the cell nodes
    template< typename ElementType >
    void updateCellNode (const ElementType& element);

    //! Update the cell nodes with a std::vector of points
    /*!
      Note here that this method is NOT a template specialization (which
      is impossible in this case, full specialization being forbidden in
      non specialized template classes), but an overload of the template
      method (which represents a better match in case the right arguments
      are passed).
     */
    void updateCellNode (const std::vector<VectorSmall<spaceDim> >& ptsCoordinates);

    //! Update the diameter of the cell
    void updateDiameter();

    //! Update the measure of the cell
    void updateMeasure();

    //! Update the quadrature nodes
    void updateQuadNode (const UInt& iQuadPt);

    //! Update Jacobian
    void updateJacobian (const UInt& iQuadPt);

    //! Update DetJacobian
    void updateDetJacobian (const UInt& iQuadPt);

    //! Update InverseJacobian
    void updateInverseJacobian (const UInt& iQuadPt);

    //! Update WeightedJacobian
    void updateWDet (const UInt& iQuadPt);

    //! Update Dphi
    void updateDphi (const UInt& iQuadPt);
    
    //! Update D2phi
    void updateD2phi (const UInt& iQuadPt);
    
    //! Update Laplacian
    void updateLaplacian (const UInt& iQuadPt);

    //@}


    // Pointer on the reference FE
    const ReferenceFE* M_referenceFE;

    // Pointer on the geometric map
    const GeometricMap* M_geometricMap;

    // Pointer on the quadrature rule
    const QuadratureRule* M_quadratureRule;

    // Number of FE dof (current element)
    UInt M_nbFEDof;

    // Number of dof for the geometric map
    UInt M_nbMapDof;

    // Number of quadrature point
    UInt M_nbQuadPt;

    // Global identifier
    UInt M_currentId;

    // Local identifier
    UInt M_currentLocalId;

    // Diameter of the element
    Real M_diameter;

    // Measure of the element
    Real M_measure;

    // Storage for the values of the basis functions
    array2D_Type M_phi;

    // Storage for the values of the geometric map
    array2D_Type M_phiMap;

    // Storage for the derivatives of the basis functions
    array3D_Type M_dphiReferenceFE;

    // Storage for the derivatives of the geometric map
    array3D_Type M_dphiGeometricMap;
    
    // Storage for the second derivatives of the basis functions (note: the geometric map used for the 2nd derivative is M_dphiGeometricMap, no support for isoparametric elements)
    array4D_Type M_d2phiReferenceFE;

    // Storage for the coordinates of the nodes of the current element
    array2D_Type M_cellNode;

    // Storage for the position the quadrature nodes (current element)
    array1D_vector_Type M_quadNode;

    // Storage for the jacobian of the transformation
    array3D_Type M_jacobian;

    // Storage for the determinant of the jacobian of the transformation
    array1D_Type M_detJacobian;

    // Storage for the weighted determinant
    array1D_Type M_wDet;

    // Storage for the inverse of the jacobian
    array3D_Type M_tInverseJacobian;

    // Storage for the derivative of the basis functions
    array2D_vector_Type M_dphi;
    
    // Storage for the second derivative of the basis functions
    array2D_matrix_Type M_d2phi;
    
    // Storage for the laplacian
    array2D_Type M_laplacian;

#ifdef HAVE_LIFEV_DEBUG
    // Debug informations, defined only if the code
    // is compiled in debug mode. These booleans store the
    // information about what the last call to "update"
    // has actually computed (names are self explanatory)

    bool M_isCellNodeUpdated;
    bool M_isDiameterUpdated;
    bool M_isMeasureUpdated;
    bool M_isQuadNodeUpdated;
    bool M_isJacobianUpdated;
    bool M_isDetJacobianUpdated;
    bool M_isInverseJacobianUpdated;
    bool M_isWDetUpdated;
    bool M_isPhiUpdated;
    bool M_isDphiUpdated;
    bool M_isD2phiUpdated;
    bool M_isLaplacianUpdated;
#endif

};


// ===================================================
// IMPLEMENTATION
// ===================================================


template <UInt spaceDim>
const UInt ETCurrentFE<spaceDim, 1>::
S_spaceDimension = spaceDim;

template <UInt spaceDim>
const UInt ETCurrentFE<spaceDim, 1>::
S_fieldDimension = 1;


// ===================================================
// Constructors & Destructor
// ===================================================


template< UInt spaceDim>
ETCurrentFE<spaceDim, 1>::
ETCurrentFE (const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr)
    :
    M_referenceFE (&refFE),
    M_geometricMap (&geoMap),
    M_quadratureRule (new QuadratureRule (qr) ),

    M_nbFEDof (M_referenceFE->nbDof() ),
    M_nbMapDof (M_geometricMap->nbDof() ),
    M_nbQuadPt (M_quadratureRule->nbQuadPt() ),

    M_currentId(),
    M_currentLocalId(),

    M_diameter(),
    M_measure(),
    M_phi(),
    M_phiMap(),
    M_dphiReferenceFE(),
    M_dphiGeometricMap(),

    M_cellNode(),
    M_quadNode(),
    M_jacobian(),
    M_detJacobian(),
    M_wDet(),
    M_tInverseJacobian(),
    M_dphi(),
    M_d2phi(),
    M_laplacian()

#ifdef HAVE_LIFEV_DEBUG
    , M_isCellNodeUpdated (false),
    M_isDiameterUpdated (false),
    M_isMeasureUpdated (false),
    M_isQuadNodeUpdated (false),
    M_isJacobianUpdated (false),
    M_isDetJacobianUpdated (false),
    M_isInverseJacobianUpdated (false),
    M_isWDetUpdated (false),
    M_isDphiUpdated (false),
    M_isD2phiUpdated (false),
    M_isLaplacianUpdated (false)
#endif


{
    // Everything's there, so reshape the arrays
    setupInternalConstants();
}


template< UInt spaceDim>
ETCurrentFE<spaceDim, 1>::
ETCurrentFE (const ReferenceFE& refFE, const GeometricMap& geoMap)
    :
    M_referenceFE (&refFE),
    M_geometricMap (&geoMap),
    M_quadratureRule (0),

    M_nbFEDof (M_referenceFE->nbDof() ),
    M_nbMapDof (M_geometricMap->nbDof() ),
    M_nbQuadPt (0),

    M_currentId(),
    M_currentLocalId(),

    M_diameter(),
    M_measure(),
    M_phi(),
    M_phiMap(),
    M_dphiReferenceFE(),
    M_dphiGeometricMap(),

    M_cellNode(),
    M_quadNode(),
    M_jacobian(),
    M_detJacobian(),
    M_wDet(),
    M_tInverseJacobian(),
    M_dphi(),
    M_d2phi(),
    M_laplacian()
#ifdef HAVE_LIFEV_DEBUG
    , M_isCellNodeUpdated (false),
    M_isDiameterUpdated (false),
    M_isMeasureUpdated (false),
    M_isQuadNodeUpdated (false),
    M_isJacobianUpdated (false),
    M_isDetJacobianUpdated (false),
    M_isInverseJacobianUpdated (false),
    M_isWDetUpdated (false),
    M_isDphiUpdated (false),
    M_isD2phiUpdated (false),
    M_isLaplacianUpdated (false)
#endif

{
    // Miss the QR, so no reshape
}

template< UInt spaceDim>
ETCurrentFE<spaceDim, 1>::
ETCurrentFE (const ETCurrentFE<spaceDim, 1>& otherFE)
    :
    M_referenceFE (otherFE.M_referenceFE),
    M_geometricMap (otherFE.M_geometricMap),
    M_quadratureRule (new QuadratureRule (*otherFE.M_quadratureRule) ),

    M_nbFEDof (otherFE.M_nbFEDof),
    M_nbMapDof (otherFE.M_nbMapDof),
    M_nbQuadPt (otherFE.M_nbQuadPt),

    M_currentId (otherFE.M_currentId),
    M_currentLocalId (otherFE.M_currentLocalId),

    M_diameter (otherFE.M_diameter),
    M_measure (otherFE.M_measure),
    M_phi (otherFE.M_phi),
    M_phiMap (otherFE.M_phiMap),
    M_dphiReferenceFE (otherFE.M_dphiReferenceFE),
    M_dphiGeometricMap (otherFE.M_dphiGeometricMap),

    M_cellNode (otherFE.M_cellNode),
    M_quadNode (otherFE.M_quadNode),
    M_jacobian (otherFE.M_jacobian),
    M_detJacobian (otherFE.M_detJacobian),
    M_wDet (otherFE.M_wDet),
    M_tInverseJacobian (otherFE.M_tInverseJacobian),
    M_dphi (otherFE.M_dphi),
    M_d2phi(otherFE.M_d2phi),
    M_laplacian(otherFE.M_laplacian)
#ifdef HAVE_LIFEV_DEBUG
    //Beware for the comma at the begining of this line!
    , M_isCellNodeUpdated ( otherFE.M_isCellNodeUpdated ),
    M_isDiameterUpdated ( otherFE.M_isDiameterUpdated ),
    M_isMeasureUpdated ( otherFE.M_isMeasureUpdated ),
    M_isQuadNodeUpdated ( otherFE.M_isQuadNodeUpdated ),
    M_isJacobianUpdated ( otherFE.M_isJacobianUpdated ),
    M_isDetJacobianUpdated ( otherFE.M_isDetJacobianUpdated ),
    M_isInverseJacobianUpdated ( otherFE.M_isInverseJacobianUpdated ),
    M_isWDetUpdated ( otherFE.M_isWDetUpdated ),
    M_isDphiUpdated ( otherFE.M_isDphiUpdated ),
    M_isD2phiUpdated ( otherFE.M_isD2phiUpdated ),
    M_isLaplacianUpdated (otherFE.M_isLaplacianUpdated )
#endif

{}


template< UInt spaceDim>
ETCurrentFE<spaceDim, 1>::
~ETCurrentFE()
{
    if (M_quadratureRule != 0)
    {
        delete M_quadratureRule;
    }
}


// ===================================================
// Methods
// ===================================================


template< UInt spaceDim>
template<typename elementType>
void
ETCurrentFE<spaceDim, 1>::
update (const elementType& element, const flag_Type& flag)
{
    ASSERT (M_referenceFE != 0, "No reference FE for the update");
    ASSERT (M_geometricMap != 0, "No geometric mapping for the update");
    ASSERT (M_quadratureRule != 0, "No quadrature rule for the update");

#ifdef HAVE_LIFEV_DEBUG
    // Reset all the flags to false
    M_isCellNodeUpdated = false;
    M_isDiameterUpdated = false;
    M_isMeasureUpdated = false;
    M_isQuadNodeUpdated = false;
    M_isJacobianUpdated = false;
    M_isDetJacobianUpdated = false;
    M_isInverseJacobianUpdated = false;
    M_isWDetUpdated = false;
    M_isDphiUpdated = false;
    M_isD2phiUpdated = false;
    M_isLaplacianUpdated = false;
#endif

    // update the cell informations if required
    if ( flag & ET_UPDATE_ONLY_CELL_NODE )
    {
        updateCellNode (element);
    }
    if ( flag & ET_UPDATE_ONLY_DIAMETER )
    {
        updateDiameter();
    }

    // Loop over the quadrature nodes
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        // and update the required quantities
        if ( flag & ET_UPDATE_ONLY_QUAD_NODE )
        {
            updateQuadNode (i);
        }
        if ( flag & ET_UPDATE_ONLY_JACOBIAN )
        {
            updateJacobian (i);
        }
        if ( flag & ET_UPDATE_ONLY_DET_JACOBIAN )
        {
            updateDetJacobian (i);
        }
        if ( flag & ET_UPDATE_ONLY_T_INVERSE_JACOBIAN )
        {
            updateInverseJacobian (i);
        }
        if ( flag & ET_UPDATE_ONLY_W_DET_JACOBIAN )
        {
            updateWDet (i);
        }
        if ( flag & ET_UPDATE_ONLY_DPHI )
        {
            updateDphi (i);
        }
        if ( flag & ET_UPDATE_ONLY_D2PHI )
        {
            updateD2phi (i);
        }
        if ( flag & ET_UPDATE_ONLY_LAPLACIAN )
        {
            updateLaplacian (i);
        }
    }

    if ( flag & ET_UPDATE_ONLY_MEASURE )
    {
        updateMeasure();
    }
}


template<UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
showMe (std::ostream& out) const
{
    out << " Number of FE Dof   : " << M_nbFEDof << std::endl;
    out << " Number of Map Dof  : " << M_nbMapDof << std::endl;
    out << " Number of QR point : " << M_nbQuadPt << std::endl;
    out << std::endl;

    out << " Cell Nodes : " << std::endl;
    for ( UInt i (0); i < M_nbMapDof; ++i )
    {
        for ( UInt icoor (0); icoor < S_spaceDimension; ++icoor)
        {
            out << M_cellNode[i][icoor] << " ";
        }
        out << std::endl;
    }
    out << std::endl;

    out << " Jacobian : " << std::endl;
    for (UInt iQuad (0); iQuad < M_nbQuadPt; ++iQuad)
    {
        for (UInt iDim (0); iDim < S_spaceDimension; ++iDim)
        {
            for (UInt jDim (0); jDim < S_spaceDimension; ++jDim)
            {
                out << M_jacobian[iQuad][iDim][jDim] << " ";
            }
            out << std::endl;
        }
        out << std::endl;
    }

    out << " Det jacobian : " << std::endl;
    for (UInt iQuad (0); iQuad < M_nbQuadPt; ++iQuad)
    {
        out << M_detJacobian[iQuad] << " ";
    }
    out << std::endl;

    out << " T inverse Jacobian : " << std::endl;
    for (UInt iQuad (0); iQuad < M_nbQuadPt; ++iQuad)
    {
        for (UInt iDim (0); iDim < S_spaceDimension; ++iDim)
        {
            for (UInt jDim (0); jDim < S_spaceDimension; ++jDim)
            {
                out << M_tInverseJacobian[iQuad][iDim][jDim] << " ";
            }
            out << std::endl;
        }
        out << std::endl;
    }

    out << " DPhi : " << std::endl;
    for (UInt iQuad (0); iQuad < M_nbQuadPt; ++iQuad)
    {
        for (UInt iDof (0); iDof < M_nbFEDof; ++iDof)
        {
            for (UInt iCoor (0); iCoor < S_spaceDimension; ++iCoor)
            {
                out << M_dphi[iQuad][iDof][iCoor] << " ";
            }
            out << std::endl;
        }
        out << std::endl;
    }
    
    out << " D2Phi : " << std::endl;
    for (UInt iQuad (0); iQuad < M_nbQuadPt; ++iQuad)
    {
        for (UInt iDof (0); iDof < M_nbFEDof; ++iDof)
        {
            for (UInt iCoor (0); iCoor < S_spaceDimension; ++iCoor)
            {
                for (UInt jCoor (0); jCoor < S_spaceDimension; ++jCoor)
                {
                    out << M_d2phi[iQuad][iDof][iCoor][jCoor] << " ";
                }
                out << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
    }
    
    out << " Laplacian : " << std::endl;
    for (UInt iQuad (0); iQuad < M_nbQuadPt; ++iQuad)
    {
        for (UInt iDof (0); iDof < M_nbFEDof; ++iDof)
        {
            out << M_laplacian[iQuad][iDof] << " ";
        }
        out << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================

template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
setQuadratureRule (const QuadratureRule& qr)
{
    if (M_quadratureRule != 0)
    {
        delete M_quadratureRule;
    }
    M_quadratureRule = new QuadratureRule (qr);
    M_nbQuadPt = qr.nbQuadPt();
    setupInternalConstants();
}


// ===================================================
// Private Methods
// ===================================================


template< UInt spaceDim >
void
ETCurrentFE<spaceDim, 1>::
setupInternalConstants()
{
    // The first group of values can be computed as it
    // it does not depend on the current element

    // PHI
    M_phi.resize (M_nbQuadPt);
    for (UInt q (0); q < M_nbQuadPt; ++q)
    {
        M_phi[q].resize (M_nbFEDof);
        for (UInt j (0); j < M_nbFEDof; ++j)
        {
            M_phi[q][j] = M_referenceFE->phi (j, M_quadratureRule->quadPointCoor (q) );
        }
    }

#ifdef HAVE_LIFEV_DEBUG
    M_isDphiUpdated = true;
#endif

    // PHI MAP
    M_phiMap.resize (M_nbQuadPt);
    for (UInt q (0); q < M_nbQuadPt; ++q)
    {
        M_phiMap[q].resize (M_nbMapDof);
        for (UInt i (0); i < M_nbMapDof; ++i)
        {
            M_phiMap[q][i] = M_geometricMap->phi (i, M_quadratureRule->quadPointCoor (q) );
        }
    }

    // DPHIREFERENCEFE
    M_dphiReferenceFE.resize (M_nbQuadPt);
    for (UInt q (0); q < M_nbQuadPt; ++q)
    {
        M_dphiReferenceFE[q].resize (M_nbFEDof);
        for (UInt i (0); i < M_nbFEDof; ++i)
        {
            M_dphiReferenceFE[q][i].resize (spaceDim);
            for (UInt j (0); j < spaceDim; ++j)
            {
                M_dphiReferenceFE[q][i][j] = M_referenceFE->dPhi (i, j, M_quadratureRule->quadPointCoor (q) );
            }
        }
    }

    // DPHIGEOMETRICMAP
    M_dphiGeometricMap.resize (M_nbQuadPt);
    for (UInt q (0); q < M_nbQuadPt; ++q)
    {
        M_dphiGeometricMap[q].resize (M_nbMapDof);
        for (UInt i (0); i < M_nbMapDof; ++i)
        {
            M_dphiGeometricMap[q][i].resize (spaceDim);
            for (UInt j (0); j < spaceDim; ++j)
            {
                M_dphiGeometricMap[q][i][j] = M_geometricMap->dPhi (i, j, M_quadratureRule->quadPointCoor (q) );
            }
        }
    }
    
    // DPHI2REFERENCEFE
    M_d2phiReferenceFE.resize (M_nbQuadPt);
    for (UInt q (0); q < M_nbQuadPt; ++q)
    {
        M_d2phiReferenceFE[q].resize (M_nbFEDof);
        for (UInt i (0); i < M_nbFEDof; ++i)
        {
            M_d2phiReferenceFE[q][i].resize (spaceDim);
            for (UInt j (0); j < spaceDim; ++j)
            {
                M_d2phiReferenceFE[q][i][j].resize (spaceDim);
                for (UInt k (0); k < spaceDim; ++k)
                {
                    M_d2phiReferenceFE[q][i][j][k] = M_referenceFE->d2Phi (i, j, k, M_quadratureRule->quadPointCoor (q) );
                }
            }
        }
    }

    // The second group of values cannot be computed
    // now because it depends on the current element.
    // So, we just make space for it.

    // Cell nodes
    M_cellNode.resize (M_nbMapDof);
    for (UInt i (0); i < M_nbMapDof; ++i)
    {
        M_cellNode[i].resize (spaceDim);
    }

    // Quad nodes
    M_quadNode.resize (M_nbQuadPt);
    /*for (UInt i(0); i<M_nbQuadPt; ++i)
    {
        M_quadNode[i].resize(spaceDim);
        }*/

    // Jacobian
    M_jacobian.resize (M_nbQuadPt);
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        M_jacobian[i].resize (spaceDim);
        for (UInt j (0); j < spaceDim; ++j)
        {
            M_jacobian[i][j].resize (spaceDim);
        }
    }

    // Det jacobian
    M_detJacobian.resize (M_nbQuadPt);

    // wDet
    M_wDet.resize (M_nbQuadPt);

    // tInverseJacobian
    M_tInverseJacobian.resize (M_nbQuadPt);
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        M_tInverseJacobian[i].resize (spaceDim);
        for (UInt j (0); j < spaceDim; ++j)
        {
            M_tInverseJacobian[i][j].resize (spaceDim);
        }
    }

    // dphi
    M_dphi.resize (M_nbQuadPt);
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        M_dphi[i].resize (M_nbFEDof);
    }
    
    // d2phi
    M_d2phi.resize (M_nbQuadPt);
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        M_d2phi[i].resize (M_nbFEDof);
    }

    // laplacian
    M_laplacian.resize (M_nbQuadPt);
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        M_laplacian[i].resize (M_nbFEDof);
    }

    
}


template <UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateQuadNode (const UInt& iQuadPt)
{
    // Check the requirements
    ASSERT (M_isCellNodeUpdated, "Cell must be updated to compute the quadrature node position");

    // Set the check boolean
#ifdef HAVE_LIFEV_DEBUG
    M_isQuadNodeUpdated = true;
#endif

    for (UInt iDim (0); iDim < S_spaceDimension; ++iDim)
    {
        M_quadNode[iQuadPt][iDim] = 0.0;

        for (UInt iDof (0); iDof < M_nbMapDof; ++iDof)
        {
            M_quadNode[iQuadPt][iDim] += M_cellNode[iDof][iDim] * M_phiMap[iQuadPt][iDof];
        }
    }
}


template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateJacobian (const UInt& iQuadPt)
{
    // Check the requirements
    ASSERT (M_isCellNodeUpdated, "Cell must be updated to compute the jacobian");

    // Set the check boolean
#ifdef HAVE_LIFEV_DEBUG
    M_isJacobianUpdated = true;
#endif

    Real partialSum (0.0);
    for (UInt iDim (0); iDim < S_spaceDimension; ++iDim)
    {
        for (UInt jDim (0); jDim < S_spaceDimension; ++jDim)
        {
            partialSum = 0.0;
            for (UInt iterNode (0); iterNode < M_nbMapDof; ++iterNode)
            {
                partialSum += M_cellNode[iterNode][iDim] * M_dphiGeometricMap[iQuadPt][iterNode][jDim];
            }

            M_jacobian[iQuadPt][iDim][jDim] = partialSum;
        }
    }
}

template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateWDet (const UInt& iQuadPt)
{
    ASSERT (M_isDetJacobianUpdated, "Determinant of the jacobian must be updated to compute WDet");

#ifdef HAVE_LIFEV_DEBUG
    M_isWDetUpdated = true;
#endif

    M_wDet[iQuadPt] = M_detJacobian[iQuadPt] * M_quadratureRule->weight (iQuadPt);
}


template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateDphi (const UInt& iQuadPt)
{
    ASSERT (M_isInverseJacobianUpdated,
            "Inverse jacobian must be updated to compute the derivative of the basis functions");

#ifdef HAVE_LIFEV_DEBUG
    M_isDphiUpdated = true;
#endif

    Real partialSum (0.0);

    for (UInt iDof (0); iDof < M_nbFEDof; ++iDof)
    {
        for (UInt iCoor (0); iCoor < S_spaceDimension; ++iCoor)
        {
            partialSum = 0.0;
            for (UInt jCoor (0); jCoor < S_spaceDimension; ++jCoor)
            {
                partialSum += M_tInverseJacobian[iQuadPt][iCoor][jCoor] * M_dphiReferenceFE[iQuadPt][iDof][jCoor];
            }
            M_dphi[iQuadPt][iDof][iCoor] = partialSum;
        }
    }
}
    
template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateD2phi (const UInt& iQuadPt)
{
    ASSERT (M_isInverseJacobianUpdated,
            "Inverse jacobian must be updated to compute the derivative of the basis functions");
    
#ifdef HAVE_LIFEV_DEBUG
    M_isD2phiUpdated = true;
#endif
    
    Real partialSum (0.0);
    
    std::cout << "\n---------- BEGIN D2PHI ---------\n";
    
    for (UInt iDof (0); iDof < M_nbFEDof; ++iDof)
    {
        std::cout << "\n begin Dof " << iDof << ":\n";
        
        for (UInt iCoor (0); iCoor < S_spaceDimension; ++iCoor)
        {
            for (UInt jCoor (0); jCoor < S_spaceDimension; ++jCoor)
            {
                partialSum = 0.0;
                for ( UInt k1 (0); k1 < S_spaceDimension; ++k1 )
                {
                    for ( UInt k2 (0) ; k2 < S_spaceDimension; ++k2 )
                    {
                        partialSum += M_tInverseJacobian[iQuadPt][iCoor][k1] * M_d2phiReferenceFE[iQuadPt][iDof][k1][k2] * M_tInverseJacobian[iQuadPt][jCoor][k2];
                    }
                }
                M_d2phi[iQuadPt][iDof][iCoor][jCoor] = partialSum;
            }
        }
        std::cout << M_d2phi[iQuadPt][iDof];
        std::cout << "\n\n";
        std::cout << "End dof ---------------------\n";
    }
    
    std::cout << "\n---------- END D2PHI ---------\n";
}
    
template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateLaplacian (const UInt& iQuadPt)
{
    ASSERT ( M_isD2phiUpdated,
            "Basis function second derivatives must be updated to compute the laplacian" );
    
#ifdef HAVE_LIFEV_DEBUG
    M_isLaplacianUpdated = true;
#endif
    
    Real partialSum ( 0.0 );
    
    std::cout << "\n-------- LAPLACIAN BEGIN -------\n";
    
    for ( UInt iDof ( 0 ); iDof < M_nbFEDof; ++iDof )
    {
        partialSum = 0.0;
        for ( UInt iCoor ( 0 ); iCoor < S_spaceDimension; ++iCoor ) 
        {
            partialSum += M_d2phi[iQuadPt][iDof][iCoor][iCoor];
        }
        M_laplacian[iQuadPt][iDof] = partialSum;
        std::cout << "Laplacian of dof " << iDof << " equal to ";
        std::cout << M_laplacian[iQuadPt][iDof];
        std::cout << "\n---------------\n";
    }
}
    
    
template< UInt spaceDim>
template< typename ElementType >
void
ETCurrentFE<spaceDim, 1>::
updateCellNode (const ElementType& element)
{

#ifdef HAVE_LIFEV_DEBUG
    M_isCellNodeUpdated = true;
#endif

    M_currentId      = element.id();
    M_currentLocalId = element.localId();

    for ( UInt i (0); i < M_nbMapDof; ++i )
    {
        for ( UInt icoor (0); icoor < S_spaceDimension; ++icoor)
        {
            M_cellNode[i][icoor] = element.point (i).coordinate (icoor);
        }
    }
}

template< UInt spaceDim >
void
ETCurrentFE<spaceDim, 1>::
updateCellNode (const std::vector<VectorSmall<spaceDim> >& ptsCoordinates)
{
#ifdef HAVE_LIFEV_DEBUG
    M_isCellNodeUpdated = true;
#endif

    ASSERT ( ptsCoordinates.size() == M_nbMapDof, "Number of points does not define the right geometric shape.");

    for ( UInt i (0); i < M_nbMapDof; ++i )
    {
        for ( UInt icoor (0); icoor < S_spaceDimension; ++icoor)
        {
            M_cellNode[i][icoor] = ptsCoordinates[i][icoor];
        }
    }
}

template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateDiameter()
{
    ASSERT (M_isCellNodeUpdated, "Cell must be updated to compute the diameter");

#ifdef HAVE_LIFEV_DEBUG
    M_isDiameterUpdated = true;
#endif

    M_diameter = 0.0;
    Real dist;

    for ( UInt i (0); i < M_nbMapDof; ++i )
    {
        for ( UInt j (i + 1); j < M_nbMapDof; ++j)
        {
            dist = 0.0;
            for ( UInt icoor (0); icoor < S_spaceDimension; ++icoor)
            {
                dist += (M_cellNode[i][icoor] - M_cellNode[j][icoor]) * (M_cellNode[i][icoor] - M_cellNode[j][icoor]);
            }
            M_diameter = std::max (M_diameter, dist);
        }
    }

    M_diameter = std::sqrt (M_diameter);
}

template< UInt spaceDim>
void
ETCurrentFE<spaceDim, 1>::
updateMeasure()
{
    ASSERT (M_isWDetUpdated, "Wdet must be updated to compute the measure");

#ifdef HAVE_LIFEV_DEBUG
    M_isMeasureUpdated = true;
#endif

    M_measure = 0.0;

    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
    {
        M_measure += M_wDet[iq];
    }
}


} // Namespace LifeV

#endif /* ETCURRENTFE_HPP */
