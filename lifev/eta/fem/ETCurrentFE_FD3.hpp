///*!
// ETCurrenteFE is a template class.
// This is the non-specialized class used with fieldDim = 3 or 2 (using faster VectorSmall return objects)
//
// */

template< UInt spaceDim, UInt fieldDim >
class ETCurrentFE
{

    //! @name Friends
    //@{

    //!Friend to allow direct access to the raw data
    template< UInt dim >
    friend class ExpressionAssembly::EvaluationPhiI;

    //!Friend to allow direct access to the raw data
    template< UInt dim >
    friend class ExpressionAssembly::EvaluationPhiJ;

    //!Friend to allow direct access to the raw data
    template< UInt dim, UInt FSpaceDim >
    friend class ExpressionAssembly::EvaluationDphiI;

    //!Friend to allow direct access to the raw data
    template< UInt dim, UInt FSpaceDim >
    friend class ExpressionAssembly::EvaluationDphiJ;

    //!Friend to allow direct access to the raw data
    template< UInt dim, UInt FSpaceDim >
    friend class ExpressionAssembly::EvaluationDivI;

    //!Friend to allow direct access to the raw data
    template< UInt dim, UInt FSpaceDim >
    friend class ExpressionAssembly::EvaluationDivJ;

    //@}

private:

    // Vector return type for phi
    typedef VectorSmall< fieldDim > array1D_Return_Type;

    // Matrix return type for dphi
    typedef MatrixSmall< fieldDim, spaceDim > matrix_Return_Type;

    //Private typedefs for the 2D array of vector
    typedef std::vector< std::vector< array1D_Return_Type > > array2D_vector_Type;

    //Private typedefs for the 3D array of vector
    typedef std::vector< std::vector< matrix_Return_Type > > array2D_matrix_Type;

public:

    //! @name Static constants
    //@{

    //! Static value for the space dimension
    static const UInt S_spaceDimension; // = spaceDim;

    //! Static value for the dimension of the field (fieldDim here as it is a vectorial FE)
    static const UInt S_fieldDimension; // = fieldDim;

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
    ETCurrentFE (const ETCurrentFE<spaceDim, fieldDim>& otherFE);

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
      @return The vector<3> of the ith basis function in the qth quadrature node
     */
    const array1D_Return_Type& phi (const UInt& i, const UInt& q) const
    {
        ASSERT ( i < fieldDim * M_nbFEDof, "No basis function with this index" );
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index" );

        return ( M_phi[q][i] );
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
      @param iCoor The component of the basis function to be derived
      @param dxi The direction of the derivative required (0 for d/dx, 1 for d/dy...)
      @param q The index of the quadrature node
      @return The iCoor component of the ith basis function derived w.r. to dxi, in the qth quadrature node.
     */
    const Real& dphi (const UInt& i, const UInt& iCoor, const UInt& dxi, const UInt& q) const
    {
        ASSERT ( M_isDphiUpdated, "Derivative of the basis functions have not been updated");
        ASSERT ( i < fieldDim * M_nbFEDof, "No basis function with this index");
        ASSERT ( iCoor < fieldDim, "No such coordinate index");
        ASSERT ( dxi < spaceDim, "No such coordinate index");
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index");

        return M_dphi[q][i][iCoor][dxi];
    }

    //! Getter for the derivatives of the basis function in the quadrature nodes (current element)
    /*!
      @param i The index of the basis function
      @param q The index of the quadrature node
      @return The local vector of the basis functions derived w.r. to dxi, in the qth quadrature node.
     */
    MatrixSmall<spaceDim, 3> const& dphi (const UInt& i, const UInt& q) const
    {
        ASSERT ( M_isDphiUpdated, "Derivative of the basis functions have not been updated");
        ASSERT ( i < 3 * M_nbFEDof, "No basis function with this index");
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index");
        return M_dphi[q][i];
    }


    //! Getter for the divergence of the basis functions in the quadrature nodes in the current element
    /*!
      @param i The index of the basis function
      @param q The index of the quadrature node
      @return The divergence of the ith basis function in the qth quadrature node
     */
    const array1D_Return_Type& divergence (const UInt& i, const UInt& q) const
    {
        ASSERT ( M_isDivergenceUpdated, "Divergence of the basis functions have not been updated");
        ASSERT ( i < fieldDim * M_nbFEDof, "No basis function with this index" );
        ASSERT ( q < M_nbQuadPt, "No quadrature point with this index" );

        return ( M_divergence[q][i] );
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

    //@}

private:

    //Private typedefs for the 1D array
    typedef std::vector< Real > array1D_Type;

    //Private typedefs for the 2D array (array of 1D array)
    typedef std::vector< array1D_Type > array2D_Type;

    //Private typedefs for the 3D array (array of 2D array)
    typedef std::vector< array2D_Type > array3D_Type;

    //! @name Private Methods
    //@{

    //! No default constructor
    ETCurrentFE();

    //! No assignement
    void operator= ( const ETCurrentFE< spaceDim, fieldDim >& );

    //! Resize all the internal containers w.r. to the stored data and compute the constant values
    void setupInternalConstants();

    //! Update the cell nodes
    template< typename ElementType >
    void updateCellNode (const ElementType& element);

    //! Update the quadrature nodes
    void updateQuadNode (const UInt& iQuadPt);

    //! Update Jacobian
    void updateJacobian (const UInt& iQuadPt);

    //! Update DetJacobian
    void updateDetJacobian ( const UInt& iQuadPt );

    //! Update InverseJacobian
    void updateInverseJacobian ( const UInt& iQuadPt );

    //! Update WeightedJacobian
    void updateWDet ( const UInt& iQuadPt );

    //! Update Dphi
    void updateDphi ( const UInt& iQuadPt );

    //! Update Divergence
    void updateDivergence (const UInt& iQuadPt);

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

    // Storage for the values of the basis functions
    array2D_vector_Type M_phi;

    // Storage for the values of the geometric map
    array2D_Type M_phiMap;
    // Storage for the derivatives of the basis functions
    array3D_Type M_dphiReferenceFE;
    // Storage for the derivatives of the geometric map
    array3D_Type M_dphiGeometricMap;

    // Storage for the coordinates of the nodes of the current element
    array2D_Type M_cellNode;
    // Storage for the position the quadrature nodes (current element)
    array2D_Type M_quadNode;
    // Storage for the jacobian of the transformation
    array3D_Type M_jacobian;
    // Storage for the determinant of the jacobian of the transformation
    array1D_Type M_detJacobian;
    // Storage for the weighted determinant
    array1D_Type M_wDet;
    // Storage for the inverse of the jacobian
    array3D_Type M_tInverseJacobian;

    // Storage for the derivative of the basis functions
    array2D_matrix_Type M_dphi;

    // Storage for the divergence of the basis functions
    array2D_Type M_divergence;

#ifdef HAVE_LIFEV_DEBUG
    // Debug informations, defined only if the code
    // is compiled in debug mode. These booleans store the
    // information about what the last call to "update"
    // has actually computed (names are self explanatory)
    bool M_isCellNodeUpdated;
    bool M_isQuadNodeUpdated;
    bool M_isJacobianUpdated;
    bool M_isDetJacobianUpdated;
    bool M_isInverseJacobianUpdated;
    bool M_isWDetUpdated;
    bool M_isPhiUpdated;
    bool M_isDphiUpdated;
    bool M_isDivergenceUpdated;
#endif

};

// ===================================================
// IMPLEMENTATION
// ===================================================

template< UInt spaceDim, UInt fieldDim >
const UInt ETCurrentFE< spaceDim, fieldDim >::S_spaceDimension = spaceDim;

template< UInt spaceDim, UInt fieldDim >
const UInt ETCurrentFE< spaceDim, fieldDim >::S_fieldDimension = fieldDim;

// ===================================================
// Constructors & Destructor
// ===================================================

template< UInt spaceDim, UInt fieldDim >
ETCurrentFE<spaceDim, fieldDim>::
ETCurrentFE (const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr)
    :
    M_referenceFE (&refFE),
    M_geometricMap (&geoMap),
    M_quadratureRule (&qr),

    M_nbFEDof (M_referenceFE->nbDof() ),
    M_nbMapDof (M_geometricMap->nbDof() ),
    M_nbQuadPt (M_quadratureRule->nbQuadPt() ),

    M_currentId(),
    M_currentLocalId(),

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
    M_divergence()

#ifdef HAVE_LIFEV_DEBUG
    , M_isCellNodeUpdated (false),
    M_isQuadNodeUpdated (false),
    M_isJacobianUpdated (false),
    M_isDetJacobianUpdated (false),
    M_isInverseJacobianUpdated (false),
    M_isWDetUpdated (false),
    M_isPhiUpdated (false),
    M_isDphiUpdated (false),
    M_isDivergenceUpdated (false)
#endif

{
    // Everything's there, so reshape the arrays
    setupInternalConstants();
}

template< UInt spaceDim, UInt fieldDim >
ETCurrentFE<spaceDim, fieldDim>::
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
    M_divergence()

#ifdef HAVE_LIFEV_DEBUG
    , M_isCellNodeUpdated (false),
    M_isQuadNodeUpdated (false),
    M_isJacobianUpdated (false),
    M_isDetJacobianUpdated (false),
    M_isInverseJacobianUpdated (false),
    M_isWDetUpdated (false),
    M_isPhiUpdated (false),
    M_isDphiUpdated (false),
    M_isDivergenceUpdated (false)
#endif

{
    // Miss the QR, so no reshape
}

template< UInt spaceDim, UInt fieldDim >
ETCurrentFE<spaceDim, fieldDim>::
ETCurrentFE (const ETCurrentFE<spaceDim, fieldDim>& otherFE)
    :
    M_referenceFE (otherFE.M_referenceFE),
    M_geometricMap (otherFE.M_geometricMap),
    M_quadratureRule (otherFE.M_quadratureRule),

    M_nbFEDof (otherFE.M_nbFEDof),
    M_nbMapDof (otherFE.M_nbMapDof),
    M_nbQuadPt (otherFE.M_nbQuadPt),

    M_currentId (otherFE.M_currentId),
    M_currentLocalId (otherFE.M_currentLocalId),

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
    M_divergence (otherFE.M_divergence)

#ifdef HAVE_LIFEV_DEBUG
    //Beware for the comma at the begining of this line!
    , M_isCellNodeUpdated ( otherFE.M_isCellNodeUpdated ),
    M_isQuadNodeUpdated ( otherFE.M_isQuadNodeUpdated ),
    M_isJacobianUpdated ( otherFE.M_isJacobianUpdated ),
    M_isDetJacobianUpdated ( otherFE.M_isDetJacobianUpdated ),
    M_isInverseJacobianUpdated ( otherFE.M_isInverseJacobianUpdated ),
    M_isWDetUpdated ( otherFE.M_isWDetUpdated ),
    M_isPhiUpdated ( otherFE.M_isPhiUpdated ),
    M_isDphiUpdated ( otherFE.M_isDphiUpdated ),
    M_isDivergenceUpdated ( otherFE.M_isDivergenceUpdated )
#endif

{}

template< UInt spaceDim, UInt fieldDim >
ETCurrentFE<spaceDim, fieldDim>::
~ETCurrentFE()
{}

// ===================================================
// Methods
// ===================================================

template< UInt spaceDim, UInt fieldDim >
template<typename elementType>
void
ETCurrentFE<spaceDim, fieldDim>::
update (const elementType& element, const flag_Type& flag)
{
    ASSERT (M_referenceFE != 0, "No reference FE for the update");
    ASSERT (M_geometricMap != 0, "No geometric mapping for the update");
    ASSERT (M_quadratureRule != 0, "No quadrature rule for the update");

#ifdef HAVE_LIFEV_DEBUG
    // Reset all the flags to false
    M_isCellNodeUpdated = false;
    M_isQuadNodeUpdated = false;
    M_isJacobianUpdated = false;
    M_isDetJacobianUpdated = false;
    M_isInverseJacobianUpdated = false;
    M_isWDetUpdated = false;
    M_isDphiUpdated = false;
    M_isDivergenceUpdated = false;
#endif

    // update the cell informations if required
    if (flag & ET_UPDATE_ONLY_CELL_NODE)
    {
        updateCellNode (element);
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
        if ( flag & ET_UPDATE_ONLY_DIVERGENCE )
        {
            updateDivergence (i);
        }
    }
}

template< UInt spaceDim, UInt fieldDim >
void
ETCurrentFE<spaceDim, fieldDim>::
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
                for (UInt jCoor (0); jCoor < S_spaceDimension; ++jCoor)
                {
                    out << M_dphi[iQuad][iDof][iCoor][jCoor] << " ";
                }
            }
            out << std::endl;
        }
        out << std::endl;
    }

    out << " Divergence : " << std::endl;
    for (UInt iQuad (0); iQuad < M_nbQuadPt; ++iQuad)
    {
        for (UInt iDof (0); iDof < fieldDim * M_nbFEDof; ++iDof)
        {
            out << M_divergence[iQuad][iDof] << " ";
        }
        out << std::endl;
    }

}

// ===================================================
// Set Methods
// ===================================================

template< UInt spaceDim, UInt fieldDim >
void
ETCurrentFE<spaceDim, fieldDim>::
setQuadratureRule (const QuadratureRule& qr)
{
    M_quadratureRule = &qr;
    M_nbQuadPt = qr.nbQuadPt();
    setupInternalConstants();
}

// ===================================================
// Private Methods
// ===================================================

template< UInt spaceDim, UInt fieldDim >
void
ETCurrentFE<spaceDim, fieldDim>::
setupInternalConstants()
{
    // The first group of values can be computed as it
    // it does not depend on the current element

    // PHI
    M_phi.resize ( M_nbQuadPt );
    for ( UInt q ( 0 ); q < M_nbQuadPt; ++q )
    {
        // we have M_nbFEDof * 3 basis functions
        M_phi[q].resize ( M_nbFEDof * 3 );

        // set only appropriate values, other are initialized to 0 by default constructor (of VectorSmall)
        for ( UInt j ( 0 ); j < M_nbFEDof; ++j )
        {
            M_phi[q][j][0] = M_referenceFE->phi ( j, M_quadratureRule->quadPointCoor ( q ) );

            // copy other values according to the vectorial basis functions
            for ( UInt k ( 1 ); k < S_fieldDimension; ++k )
            {
                M_phi[q][k * M_nbFEDof + j][k] = M_phi[q][j][0];
            }
        }
    }
#ifdef HAVE_LIFEV_DEBUG
    M_isPhiUpdated = true;
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
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        M_quadNode[i].resize (spaceDim);
    }

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
        // we have fieldDim * DoF basis functions
        M_dphi[i].resize ( fieldDim * M_nbFEDof );
    }

    // divergence
    M_divergence.resize (M_nbQuadPt);
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        // we have fieldDim * DoF basis functions
        M_divergence[i].resize ( fieldDim * M_nbFEDof );
    }
}

template <UInt spaceDim, UInt fieldDim >
void
ETCurrentFE<spaceDim, fieldDim>::
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

template< UInt spaceDim, UInt fieldDim >
void
ETCurrentFE<spaceDim, fieldDim>::
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

template< UInt spaceDim, UInt fieldDim >
void ETCurrentFE< spaceDim, fieldDim >::updateWDet ( const UInt& iQuadPt )
{
    ASSERT ( M_isDetJacobianUpdated,
             "Determinant of the jacobian must be updated to compute WDet" );

#ifdef HAVE_LIFEV_DEBUG
    M_isWDetUpdated = true;
#endif

    M_wDet[iQuadPt] = M_detJacobian[iQuadPt] * M_quadratureRule->weight ( iQuadPt );
}

template< UInt spaceDim, UInt fieldDim >
void ETCurrentFE< spaceDim, fieldDim >::updateDphi ( const UInt& iQuadPt )
{
    ASSERT ( M_isInverseJacobianUpdated,
             "Inverse jacobian must be updated to compute the derivative of the basis functions" );

#ifdef HAVE_LIFEV_DEBUG
    M_isDphiUpdated = true;
#endif

    Real partialSum ( 0.0 );

    for ( UInt iDof ( 0 ); iDof < M_nbFEDof; ++iDof )
    {
        for ( UInt iCoor ( 0 ); iCoor < S_spaceDimension; ++iCoor )
        {
            partialSum = 0.0;
            for ( UInt jCoor ( 0 ); jCoor < S_spaceDimension; ++jCoor )
            {
                partialSum += M_tInverseJacobian[iQuadPt][iCoor][jCoor] * M_dphiReferenceFE[iQuadPt][iDof][jCoor];
            }

            // set only appropriate values, other are initialized to 0 by default constructor (of VectorSmall)
            M_dphi[iQuadPt][iDof][0][iCoor] = partialSum;

            // copy other values according to the vectorial basis functions
            for ( UInt k ( 1 ); k < fieldDim; ++k)
            {
                M_dphi[iQuadPt][k * M_nbFEDof + iDof][k][iCoor] = partialSum;
            }
        }
    }
}

template< UInt spaceDim, UInt fieldDim >
void ETCurrentFE< spaceDim, fieldDim >::updateDivergence ( const UInt& iQuadPt )
{
    ASSERT ( M_isDphiUpdated,
             "Basis function derivatives must be updated to compute the divergence" );

#ifdef HAVE_LIFEV_DEBUG
    M_isDivergenceUpdated = true;
#endif

    Real partialSum ( 0.0 );

    for ( UInt iDof ( 0 ); iDof < fieldDim * M_nbFEDof; ++iDof )
    {
        partialSum = 0.0;
        for ( UInt iCoor ( 0 ); iCoor < S_spaceDimension; ++iCoor )
        {
            partialSum += M_dphi[iQuadPt][iDof][iCoor][iCoor];
        }
        M_divergence[iQuadPt][iDof] = partialSum;
    }
}


template< UInt spaceDim, UInt fieldDim >
template< typename ElementType >
void
ETCurrentFE<spaceDim, fieldDim>::
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
