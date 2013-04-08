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
    @brief File containing a class with methods to manage fields defined on the domain boundary

    @author Alessandro Veneziani <ale@mathcs.emory.edu>
    @author Tiziano Passerini <tiziano@mathcs.emory.edu>

    @maintainer Tiziano Passerini <tiziano@mathcs.emory.edu>

    @date 09-2008
 */

#ifndef POSTPROCESSINGBOUNDARY_H
#define POSTPROCESSINGBOUNDARY_H 1

//#include <string>
//#include <iostream>
//#include <sstream>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! Class with methods to manage fields defined on the domain boundary
/*!
    The data structures are based on vector M_vectorNumberingPerFacetVector:
    <ol>
        <li> the size is equal to the number of boundary facets in the mesh
        <li> each component is a vector of identifiers
        <li> in each vector of identifiers, the i-th component is the position in vector
            M_dofGlobalIdVector where to find the global Id for i-th dof in the facet
    </ol>

    Vector M_dofGlobalIdVector has this usage:
    <ol>
        <li> the size is equal to the total number of boundary dof's
        <li> the k-th component is the global numbering of the k-th boundary dof
    </ol>

    How to use M_vectorNumberingPerFacetVector and M_dofGlobalIdVector together?
    <ol>
        <li> M_vectorNumberingPerFacetVector[ j ] is a vector containing the vector indices
             for boundary face j
        <li> M_dofGlobalIdVector[ M_vectorNumberingPerFacetVector[ j ][ i ] ] is the global ID
             associated to the i-th dof on j-th face
    </ol>

    NB: here "global" is intended wrt the global mesh (prior to partitioning, if this applies)
 */
template <typename MeshType>
class PostProcessingBoundary
{

    //! @name Public Types
    //@{
    typedef typename MeshType::elementShape_Type             elementGeometricShape_Type;
    typedef typename elementGeometricShape_Type::GeoBShape   facetGeometricShape_Type;
    typedef MeshType                                         mesh_Type;
    typedef boost::shared_ptr<MeshType>                      meshPtr_Type;
    typedef CurrentFEManifold*                               currentBdFEPtr_Type;
    typedef DOF*                                             dofPtr_Type;
    //@}

public:

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    PostProcessingBoundary<MeshType>()
    {}

    //! Copy constructor
    /*!
        @note all pointer members are copied (no new allocation is done in the class)
     */
    PostProcessingBoundary ( const PostProcessingBoundary& /*postProc*/ )
    {}

    //! Constructor
    /*!
        In this general case we allow to pass an arbitrary number (nFESpaces) of dof/fe classes

        NOTE: in most parts of the code, only currentBdFEVector[0] and dofVector[0] will be considered

        \param mesh the mesh
        \param currentBdFEVector a vector of finite element prototypes for boundary elements
        \param dofVector a vector of classes for the description of degrees of freedom
        \param epetraMap the map describing the partition over different processors
        \param nFESpaces number of elements in vectors currentBdFEVector, dofVector
     */
    PostProcessingBoundary<MeshType> ( meshPtr_Type mesh,
                                       std::vector<currentBdFEPtr_Type > currentBdFEVector,
                                       std::vector<dofPtr_Type > dofVector,
                                       const MapEpetra& epetraMap, UInt nFESpaces = 1 );

    /*!
     \brief Constructor for the case in which we have only one currentBdFE and dof
     */
    PostProcessingBoundary<MeshType> ( meshPtr_Type mesh,
                                       currentBdFEPtr_Type currentBdFE, dofPtr_Type dof,
                                       const MapEpetra& epetraMap );

    /*!
     \brief Constructor for the case in which we have two currentBdFE's and dof's

     This is the case of NS problem, in which we will want to consider both
     velocity and pressure dof/fe classes
     */
    PostProcessingBoundary<MeshType> ( meshPtr_Type mesh,
                                       currentBdFEPtr_Type feBdu, dofPtr_Type dofu,
                                       currentBdFEPtr_Type feBdp, dofPtr_Type dofp,
                                       const MapEpetra& epetraMap );

    //@}

    /*! @name Methods */
    //@{

    /*! @defgroup boundary_methods
     These methods compute quantities on boundary sections associated to a given flag
     */

    /*!
       This method computes the measure of boundary section "flag"
       @ingroup boundary_methods

       \param flag is the marker of the considered boundary section
     */
    Real measure ( const markerID_Type& flag );

    /*!
       This method computes the flux of vectorField across boundary section "flag"
       @ingroup boundary_methods

      \tparam VectorType Vector type. Basic policy for type VectorType: operator[] available
      \param vectorField is intended to be a vector field
      \param flag is the marker of the considered boundary section
      \param feSpace is the identifier of the desired FE Space in M_feSpaceVector
      \param nDim is the dimension of vectorField
     */
    template< typename VectorType >
    Real flux ( const VectorType& vectorField, const markerID_Type& flag, UInt feSpace = 0, UInt nDim = nDimensions );

    /*! Compute the kinetic normal stress (i.e., the normal stress due to the kinetic energy) on a boundary face.
     *
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *  @ingroup boundary_methods
     *
     *  This method computes the following quantity:
     *
     *  \f[
     *  \mathcal{K} = \frac{1}{2}\rho_\textrm{F}\frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}\left({\mathbf{u}}_\textrm{F} \mathbf{\cdot} {\mathbf{n}}_\textrm{F}\right)^2  \textrm{d} \Gamma
     *  \f]
     *
     *  @param velocity velocity
     *  @param density density of the fluid
     *  @param flag the flag of the boundary face
     *  @param feSpace the FE space
     *  @param nDim the dimension size
     *  @return the kinetic normal stress
     */
    template< typename VectorType >
    Real kineticNormalStress ( const VectorType& velocity, const Real& density, const markerID_Type& flag, UInt feSpace = 0, UInt nDim = nDimensions );

    /*! Compute the derivative of the kinetic normal stress (i.e., the derivative of the normal stress due to the kinetic energy) on a boundary face.
     *
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *  @ingroup boundary_methods
     *
     *  This method computes the following quantity:
     *
     *  \f[
     *  \begin{array}{r@{\,\,}c@{\,\,}l@{\qquad}l}
     *  \textrm{D}\mathcal{K} &=&\displaystyle\frac{1}{2}\rho_\textrm{F} \displaystyle\frac{1}{{\left|\Gamma^t_{\textrm{F},j_1}\right|}^2}
     *  \left(\displaystyle\int_{\Gamma^t_{\textrm{F},j_1}}{\mathbf\nabla}_\Gamma \mathbf \cdot \delta {\mathbf d}_\textrm{F}  \textrm{d} \Gamma\right)
     *  \left(\displaystyle\int_{\Gamma^t_{\textrm{F},j_1}}{\left({\mathbf u}_\textrm{F} \mathbf \cdot {\mathbf n}_\textrm{F}\right)}^2  \textrm{d} \Gamma \right) \\[4ex]
     *  &-&\displaystyle\frac{1}{2}\rho_\textrm{F}\displaystyle\frac{1}{\left|\Gamma^t_{\textrm{F},j_1}\right|}
     *  \left(\displaystyle\int_{\Gamma^t_{\textrm{F},j_1}}2({\mathbf u}_\textrm{F}\mathbf \cdot {\mathbf n}_\textrm{F})\left(\delta {\mathbf u}_\textrm{F} \mathbf \cdot {\mathbf n}_\textrm{F}\right)  \textrm{d} \Gamma
     *  +\displaystyle\int_{\Gamma^t_{\textrm{F},j_1}}\left({\mathbf\nabla}_\Gamma \mathbf \cdot \delta {\mathbf d}_\textrm{F}\right)\left({\mathbf u}_\textrm{F} \mathbf \cdot {\mathbf n}_\textrm{F}\right)^2  \textrm{d} \Gamma \right)
     *  \end{array}
     *  \f]
     *
     *  @param velocity velocity
     *  @param velocityDerivative velocity derivative
     *  @param density density of the fluid
     *  @param flag the flag of the boundary face
     *  @param feSpace the FE space
     *  @param nDim the dimension size
     *  @return the kinetic normal stress derivative
     */
    template< typename VectorType >
    Real kineticNormalStressDerivative ( const VectorType& velocity, const VectorType& velocityDerivative, const Real& density, const markerID_Type& flag, UInt feSpace = 0, UInt nDim = nDimensions );

    /*!
       This method computes the average value of a field on the boundary section "flag"
       @ingroup boundary_methods

      \tparam VectorType Vector type. Basic policy for type VectorType: operator[] available

      \param field is intended to be a vector or a scalar (pressure in NS problem),
       this method computes the average value of field on section "flag"

       \return the averaged vector
     */
    template< typename VectorType >
    Vector average ( const VectorType& field, const markerID_Type& flag, UInt feSpace = 0, UInt nDim = 1 );

    /*!
       This method computes an approximate normal vector on the boundary section "flag"
       @ingroup boundary_methods

      \return the approximate normal vector
     */
    Vector normal ( const markerID_Type& flag, UInt feSpace = 0, UInt nDim = nDimensions );

    /*! Compute the geometric center of a boundary face.
     *
     *  This method computes the geometric center of a boundary section "flag"
     *
     *  @param flag the flag of the boundary face
     *  @param feSpace the FE space
     *  @param nDim the dimension size
     *  @return the vector containing the x-y-z coordinates of the geometric center
     *
     */
    Vector geometricCenter ( const markerID_Type& flag, UInt feSpace = 0, UInt nDim = nDimensions );

    // NOT READY!
#if 0
    /*!
      This procedure computes the tangential stresses on the boundary of the domain
      starting from an estimate of the viscous stresses on the boundary.

      The basic approach resorts to weak (residual-based) computing of the stresses.

      \param r contains the estimate of the viscous stresses
      \param residual switches two possibilities:
        If (residual == true), vector r contains the residual of the fluid problem,
        i. e. the integral over the domain of the product of the stress and the test functions
          r_i = \int_Omega tau \cdot \phi_i domega

        If (residual == false), vector r represents an estimate of the stress
        in each dof of the _boundary_
          r_i = tau_i
     */
    template<typename VectorType>
    VectorType compute_sstress ( const VectorType& r, UInt ncomp, bool residual ) const;

    /*!
      Here we compute the stress over the boundary, _given_ the velocity gradient

      s(n) = nu ( grad(u) + grad(u)^T ) \cdot n
      component-wise:
      s_i = nu ( d u_i / d x_j + d u_j / d x_i ) n_j

      \param grad the velocity gradient
        grad[ dofGlobalId-1 + ncomp*dim + scomp*nDimensions*dim ] =
         = ( d u_{scomp} / d x_ncomp +  d u_{ucomp} / d x_scomp )
         (at node dofGlobalId)
     */
    template<typename VectorType>
    VectorType compute_stress ( const VectorType& grad, const UInt& dim, const Real& visc  ) const;
#endif
    //@}

    /*! @name Methods */
    //@{
    /*!
     These methods print to screen information about the class data structures
     */
    void showDOFIndexMap ( std::ostream& output = std::cout ) const;

    void showPatchesMeasure ( std::ostream& output = std::cout ) const;

    void showPatchesNormal ( std::ostream& output = std::cout ) const;

    void showPatchesPhi ( std::ostream& output = std::cout ) const;

    /*!
     Access to private members
     */
    const Vector& patchMeasureInFESpace ( const UInt& feSpace = 0 ) const
    {
        return M_patchMeasureVector[feSpace];
    }

    const Vector& patchNormalInFESpace ( const UInt& feSpace = 0 ) const
    {
        return M_patchNormalVector[feSpace];
    }

    const std::vector< ID >& dofGlobalIdInFESpace ( const UInt& feSpace = 0 ) const
    {
        return M_dofGlobalIdVector[feSpace];
    }

    //! Number of boundary DOF for the mesh at hand
    const UInt& numBoundaryDofInFESpace ( const UInt& feSpace = 0 ) const
    {
        return M_numBoundaryDofVector[feSpace];
    }
    //@}


private:
    /*! @name Private Methods */
    //@{
    /*!
     These methods build the data structures describing patches of boundary elements
     (patch measure, patch normal and integral of the test function over the patch)
     */
    void                                         computePatchesMeasure();
    void                                         computePatchesNormal();
    void                                         computePatchesPhi();
    void                                         buildVectors();
    //@{

    UInt                                         M_numFESpaces;

    std::vector<UInt>                            M_numBoundaryDofVector;
    std::vector<UInt>                            M_numDofPerPeakVector, M_numDofPerRidgeVector, M_numDofPerFacetVector;
    UInt                                         M_numRidgesPerFacet, M_numFacetPerFacet;
    UInt                                         M_numPeaksPerElement, M_numRidgesPerElement;
    std::vector<UInt>                            M_numPeakDofPerFacetVector, M_numRidgeDofPerFacetVector;
    std::vector<UInt>                            M_numTotalDofPerFacetVector;
    std::vector<UInt>                            M_numPeakDofPerElement, M_numRidgeDofPerElementVector;
    UInt                                         M_numBoundaryFacets;
    std::vector<UInt>                            M_numTotalDofVector;
    // vector whose ith-component is the measure of the patch of the ith-node on the boundary
    std::vector<Vector >                         M_patchMeasureVector;
    // vector with the components of the average normal vector of each boundary node
    std::vector<Vector >                         M_patchNormalVector;
    // vector with \int_{Patch(i)} \phi_i dx where \phi_i is the basis function.
    std::vector<Vector >                         M_patchIntegratedPhiVector;

    // store once for all a map, with key={boundary flag}, value={ID list}
    std::map< markerID_Type, std::list<ID> >   M_boundaryMarkerToFacetIdMap;

    // for each boundary face, it contains the numbering of the dof of the face
    std::vector< std::vector< std::vector<ID> > > M_vectorNumberingPerFacetVector;
    // it converts from a local numbering over the boundary faces on the global numbering of the mesh
    std::vector< std::vector< ID > >             M_dofGlobalIdVector;

    // reference to the current boundary FE
    std::vector<currentBdFEPtr_Type >            M_currentBdFEPtrVector;
    // reference to a DOF class
    std::vector<dofPtr_Type >                    M_dofPtrVector;
    // pointer to the mesh
    meshPtr_Type                                 M_meshPtr;
    // pointer to the processor mapping
    boost::shared_ptr<MapEpetra>                 M_epetraMapPtr;

    const Int                                    M_geoDimension;

};

//
// IMPLEMENTATIONS
//
template <typename MeshType>
PostProcessingBoundary<MeshType>::PostProcessingBoundary ( meshPtr_Type meshPtr,
                                                           std::vector<currentBdFEPtr_Type > currentBdFEVector,
                                                           std::vector<dofPtr_Type > dofVector,
                                                           const MapEpetra& epetraMap, UInt nvar ) :
    M_numFESpaces (nvar),
    M_numBoundaryDofVector (M_numFESpaces),
    M_numDofPerPeakVector (M_numFESpaces), M_numDofPerRidgeVector (M_numFESpaces), M_numDofPerFacetVector (M_numFESpaces),
    M_numPeakDofPerFacetVector (M_numFESpaces), M_numRidgeDofPerFacetVector (M_numFESpaces),
    M_numTotalDofPerFacetVector (M_numFESpaces),
    M_numPeakDofPerElement (M_numFESpaces), M_numRidgeDofPerElementVector (M_numFESpaces),
    M_numTotalDofVector (M_numFESpaces),
    M_patchMeasureVector (M_numFESpaces), M_patchNormalVector (M_numFESpaces),
    M_patchIntegratedPhiVector (M_numFESpaces),
    M_vectorNumberingPerFacetVector (M_numFESpaces), M_dofGlobalIdVector (M_numFESpaces),
    M_currentBdFEPtrVector (currentBdFEVector), M_dofPtrVector (dofVector),
    M_meshPtr ( meshPtr ), M_epetraMapPtr ( new MapEpetra (epetraMap) ),
    M_geoDimension (MeshType::S_geoDimensions)
{
    for (UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace)
    {
        M_vectorNumberingPerFacetVector[M_numFESpaces].clear();
        M_dofGlobalIdVector[iFESpace].clear();
    }
    // Some useful local variables, to save some typing
    M_numRidgesPerFacet = facetGeometricShape_Type::S_numRidges; // Number of face's vertices
    M_numFacetPerFacet = facetGeometricShape_Type::S_numFacets;    // Number of face's ridges

    M_numPeaksPerElement = elementGeometricShape_Type::S_numPeaks; // Number of element's vertices
    M_numRidgesPerElement = elementGeometricShape_Type::S_numRidges;    // Number of element's ridges

    M_numBoundaryFacets = M_meshPtr->numBoundaryFacets();    // number of faces on boundary

    // Construction of M_vectorNumberingPerFacetVector & other data structures
    buildVectors();

    computePatchesMeasure();
    computePatchesNormal();
    computePatchesPhi();
}

template <typename MeshType>
PostProcessingBoundary<MeshType>::PostProcessingBoundary ( meshPtr_Type mesh,
                                                           currentBdFEPtr_Type currentBdFE, dofPtr_Type dof,
                                                           const MapEpetra& epetraMap ) :
    M_numFESpaces (1),
    M_numBoundaryDofVector (M_numFESpaces),
    M_numDofPerPeakVector (M_numFESpaces), M_numDofPerRidgeVector (M_numFESpaces), M_numDofPerFacetVector (M_numFESpaces),
    M_numPeakDofPerFacetVector (M_numFESpaces), M_numRidgeDofPerFacetVector (M_numFESpaces),
    M_numTotalDofPerFacetVector (M_numFESpaces),
    M_numPeakDofPerElement (M_numFESpaces), M_numRidgeDofPerElementVector (M_numFESpaces),
    M_numTotalDofVector (M_numFESpaces),
    M_patchMeasureVector (M_numFESpaces), M_patchNormalVector (M_numFESpaces), M_patchIntegratedPhiVector (M_numFESpaces),
    M_vectorNumberingPerFacetVector (M_numFESpaces), M_dofGlobalIdVector (M_numFESpaces),
    M_currentBdFEPtrVector (M_numFESpaces), M_dofPtrVector (M_numFESpaces),
    M_meshPtr ( mesh ), M_epetraMapPtr ( new MapEpetra (epetraMap) ),
    M_geoDimension (MeshType::S_geoDimensions)
{
    M_currentBdFEPtrVector[0] = currentBdFE;
    M_dofPtrVector[0] = dof;

    // Some useful local variables, to save some typing
    M_numRidgesPerFacet = facetGeometricShape_Type::S_numRidges; // Number of facet's ridges
    M_numFacetPerFacet = facetGeometricShape_Type::S_numFacets;    // Number of facet's facets

    M_numPeaksPerElement = elementGeometricShape_Type::S_numPeaks; // Number of element's vertices
    M_numRidgesPerElement = elementGeometricShape_Type::S_numRidges;    // Number of element's ridges

    M_numBoundaryFacets = M_meshPtr->numBoundaryFacets();    // number of faces on boundary

    // Construction of M_vectorNumberingPerFacetVector & other data structures
    buildVectors();

    computePatchesMeasure();
    computePatchesNormal();
    computePatchesPhi();
}

template <typename MeshType>
PostProcessingBoundary<MeshType>::PostProcessingBoundary ( meshPtr_Type mesh,
                                                           currentBdFEPtr_Type feBdu, dofPtr_Type dofu,
                                                           currentBdFEPtr_Type feBdp, dofPtr_Type dofp,
                                                           const MapEpetra& epetraMap ) :
    M_numFESpaces (2),
    M_numBoundaryDofVector (M_numFESpaces),
    M_numDofPerPeakVector (M_numFESpaces), M_numDofPerRidgeVector (M_numFESpaces), M_numDofPerFacetVector (M_numFESpaces),
    M_numPeakDofPerFacetVector (M_numFESpaces), M_numRidgeDofPerFacetVector (M_numFESpaces),
    M_numTotalDofPerFacetVector (M_numFESpaces),
    M_numPeakDofPerElement (M_numFESpaces), M_numRidgeDofPerElementVector (M_numFESpaces),
    M_numTotalDofVector (M_numFESpaces),
    M_patchMeasureVector (M_numFESpaces), M_patchNormalVector (M_numFESpaces), M_patchIntegratedPhiVector (M_numFESpaces),
    M_vectorNumberingPerFacetVector (M_numFESpaces), M_dofGlobalIdVector (M_numFESpaces),
    M_currentBdFEPtrVector (M_numFESpaces), M_dofPtrVector (M_numFESpaces),
    M_meshPtr ( mesh ), M_epetraMapPtr ( new MapEpetra (epetraMap) ),
    M_geoDimension (MeshType::S_geoDimensions)
{
    M_currentBdFEPtrVector[0] = feBdu;
    M_dofPtrVector[0] = dofu;
    M_currentBdFEPtrVector[1] = feBdp;
    M_dofPtrVector[1] = dofp;

    // Some useful local variables, to save some typing
    M_numRidgesPerFacet = facetGeometricShape_Type::S_numRidges; // Number of face's vertices
    M_numFacetPerFacet = facetGeometricShape_Type::S_numFacets;    // Number of face's ridges

    M_numPeaksPerElement = elementGeometricShape_Type::S_numPeaks; // Number of element's vertices
    M_numRidgesPerElement = elementGeometricShape_Type::S_numRidges;    // Number of element's ridges

    M_numBoundaryFacets = M_meshPtr->numBoundaryFacets();    // number of faces on boundary

    // Construction of M_vectorNumberingPerFacetVector & other data structures
    buildVectors();

    computePatchesMeasure();
    computePatchesNormal();
    computePatchesPhi();
}

// Measure of faces with a certain marker
template<typename MeshType>
void PostProcessingBoundary<MeshType>::buildVectors()
{
    std::vector<ID>            boundaryDofGlobalIdVector;

    UInt                      iFirstAdjacentElement, iPeakLocalId, iFacetLocalId, iRidgeLocalId;
    ID                        dofLocalId, dofGlobalId, dofAuxiliaryId;
    std::vector<ID>           numBoundaryDofVector (M_numFESpaces);
    std::vector<ID>::iterator dofGlobalIdVectorIterator;
    markerID_Type           boundaryFlag;

    for (UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace)
    {

        numBoundaryDofVector[iFESpace] = 0;

        // Some useful local variables, to save some typing
        M_numDofPerPeakVector[iFESpace] = M_dofPtrVector[iFESpace]->localDofPattern().nbDofPerPeak(); // number of DOF per vertices
        M_numDofPerRidgeVector[iFESpace] = M_dofPtrVector[iFESpace]->localDofPattern().nbDofPerRidge();   // number of DOF per ridges
        M_numDofPerFacetVector[iFESpace] = M_dofPtrVector[iFESpace]->localDofPattern().nbDofPerFacet();   // number of DOF per faces

        // number of peak's DOF on a face
        M_numPeakDofPerFacetVector[iFESpace] = M_numDofPerPeakVector[iFESpace] * M_numRidgesPerFacet;
        // number of ridge's DOF on a face
        M_numRidgeDofPerFacetVector[iFESpace] = M_numDofPerRidgeVector[iFESpace] * M_numFacetPerFacet;

        // number of total DOF on a face
        M_numTotalDofPerFacetVector[iFESpace] =
            M_numPeakDofPerFacetVector[iFESpace] + M_numRidgeDofPerFacetVector[iFESpace] + M_numDofPerFacetVector[iFESpace];
        // number of peak's DOF on a Element
        M_numPeakDofPerElement[iFESpace] = M_numPeaksPerElement * M_numDofPerPeakVector[iFESpace];
        // number of ridge's DOF on a Element
        M_numRidgeDofPerElementVector[iFESpace] = M_numRidgesPerElement * M_numDofPerRidgeVector[iFESpace];

        M_numTotalDofVector[iFESpace] = M_dofPtrVector[iFESpace]->numTotalDof();
    }

    // ===================================================
    // Loop on boundary faces
    // ===================================================
    for ( ID iboundaryFacet = 0 ; iboundaryFacet < M_numBoundaryFacets; ++iboundaryFacet )
    {

        iFirstAdjacentElement = M_meshPtr->boundaryFacet ( iboundaryFacet ).firstAdjacentElementIdentity(); // id of the element adjacent to the face
        iFacetLocalId = M_meshPtr->boundaryFacet ( iboundaryFacet ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

        boundaryFlag = M_meshPtr->boundaryFacet (iboundaryFacet ).markerID();
        M_boundaryMarkerToFacetIdMap[boundaryFlag].push_back ( iboundaryFacet ); // fill the flag-to-faceIdList map

        for (UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace)
        {
            clearVector (boundaryDofGlobalIdVector);

            //OLD: boundaryDofGlobalIdVector.clearVector();
            boundaryDofGlobalIdVector.resize ( M_numTotalDofPerFacetVector[iFESpace] );

            // updating finite element information
            M_currentBdFEPtrVector[iFESpace]->update ( M_meshPtr->boundaryFacet ( iboundaryFacet ), UPDATE_W_ROOT_DET_METRIC );

            // ===================================================
            // Peak based Dof
            // ===================================================
            if ( M_numDofPerPeakVector[iFESpace] )
            {

                // loop on face vertices
                for ( ID iRidgePerFacet = 0; iRidgePerFacet < M_numRidgesPerFacet; ++iRidgePerFacet )
                {
                    // local peak number (in element)
                    iPeakLocalId = elementGeometricShape_Type::faceToPoint ( iFacetLocalId, iRidgePerFacet );

                    // Loop number of DOF per peak
                    for ( ID localDof = 0; localDof < M_numDofPerPeakVector[iFESpace]; ++localDof )
                    {
                        dofLocalId = iRidgePerFacet * M_numDofPerPeakVector[iFESpace] + localDof ; // local Dof
                        dofGlobalId = M_dofPtrVector[iFESpace]->localToGlobalMap (
                                          iFirstAdjacentElement, ( iPeakLocalId ) * M_numDofPerPeakVector[iFESpace] + localDof ); // global Dof
                        dofGlobalIdVectorIterator = find (
                                                        M_dofGlobalIdVector[iFESpace].begin(), M_dofGlobalIdVector[iFESpace].end(), dofGlobalId );
                        if ( dofGlobalIdVectorIterator == M_dofGlobalIdVector[iFESpace].end() )
                        {
                            // the dofGlobalId has been encountered for the first time
                            boundaryDofGlobalIdVector[ dofLocalId ] = numBoundaryDofVector[iFESpace];
                            M_dofGlobalIdVector[iFESpace].push_back ( dofGlobalId ); // local to boundary global on this face
                            numBoundaryDofVector[iFESpace]++;
                        }
                        else
                        {
                            // the dofGlobalId has been already inserted in the M_dofGlobalIdVector vector
                            dofAuxiliaryId = ( ID ) ( ( dofGlobalIdVectorIterator - M_dofGlobalIdVector[iFESpace].begin() ) );
                            boundaryDofGlobalIdVector[ dofLocalId ] = dofAuxiliaryId; // local to boundary global on this face
                        }
                    }
                }
            }
            // ===================================================
            // Ridge based Dof
            // ===================================================
            if ( M_numDofPerRidgeVector[iFESpace] )
            {
                // loop on face ridges
                for ( ID iFacetPerFacet = 0; iFacetPerFacet < M_numFacetPerFacet; ++iFacetPerFacet )
                {
                    // local ridge number (in element)
                    iRidgeLocalId = elementGeometricShape_Type::facetToRidge ( iFacetLocalId, iFacetPerFacet );
                    // Loop number of DOF per ridge
                    for ( ID localDof = 0; localDof < M_numDofPerRidgeVector[iFESpace]; ++localDof )
                    {
                        dofLocalId = M_numPeakDofPerFacetVector[iFESpace] +
                                     iFacetPerFacet * M_numDofPerRidgeVector[iFESpace] + localDof ; // local Dof
                        dofGlobalId = M_dofPtrVector[iFESpace]->localToGlobalMap (
                                          iFirstAdjacentElement, M_numPeakDofPerElement[iFESpace] + iRidgeLocalId *
                                          M_numDofPerRidgeVector[iFESpace] + localDof ); // global Dof
                        dofGlobalIdVectorIterator = find (
                                                        M_dofGlobalIdVector[iFESpace].begin(), M_dofGlobalIdVector[iFESpace].end(), dofGlobalId );
                        if ( dofGlobalIdVectorIterator == M_dofGlobalIdVector[iFESpace].end() )
                        {
                            // the dofGlobalId has been encountered for the first time
                            boundaryDofGlobalIdVector[ dofLocalId ] = numBoundaryDofVector[iFESpace];
                            M_dofGlobalIdVector[iFESpace].push_back ( dofGlobalId ); // local to boundary global on this facet
                            numBoundaryDofVector[iFESpace]++;
                        }
                        else
                        {
                            // the dofGlobalId has been already inserted in the M_dofGlobalIdVector vector
                            dofAuxiliaryId = ( ID ) ( dofGlobalIdVectorIterator - M_dofGlobalIdVector[iFESpace].begin() );
                            boundaryDofGlobalIdVector[ dofLocalId ] = dofAuxiliaryId; // local to boundary global on this facet
                        }
                    }
                }
            }
            // ===================================================
            // Facet based Dof
            // ===================================================
            // Loop on number of DOF per facet
            for ( ID localDof = 0; localDof < M_numDofPerFacetVector[iFESpace]; ++localDof )
            {
                // local Dof
                dofLocalId = M_numRidgeDofPerFacetVector[iFESpace] + M_numPeakDofPerFacetVector[iFESpace] + localDof;
                dofGlobalId = M_dofPtrVector[iFESpace]->localToGlobalMap (
                                  iFirstAdjacentElement, M_numRidgeDofPerElementVector[iFESpace] + M_numPeakDofPerElement[iFESpace] +
                                  iFacetLocalId * M_numDofPerFacetVector[iFESpace] + localDof ); // global Dof
                dofGlobalIdVectorIterator = find (
                                                M_dofGlobalIdVector[iFESpace].begin(), M_dofGlobalIdVector[iFESpace].end(), dofGlobalId );
                if ( dofGlobalIdVectorIterator == M_dofGlobalIdVector[iFESpace].end() )
                {
                    // the dofGlobalId has been encountered for the first time
                    boundaryDofGlobalIdVector[ dofLocalId ] = numBoundaryDofVector[iFESpace];
                    M_dofGlobalIdVector[iFESpace].push_back ( dofGlobalId ); // local to boundary global on this facet
                    numBoundaryDofVector[iFESpace]++;
                }
                else
                {
                    // the dofGlobalId has been already inserted in the M_dofGlobalIdVector vector
                    dofAuxiliaryId = ( ID ) ( dofGlobalIdVectorIterator - M_dofGlobalIdVector[iFESpace].begin() ) ;
                    boundaryDofGlobalIdVector[ dofLocalId ] = dofAuxiliaryId; // local to boundary global on this facet
                }
            }

            M_vectorNumberingPerFacetVector[iFESpace].push_back ( boundaryDofGlobalIdVector );
        }
    }

    for (UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace)
    {
        // each processor holds information on HIS OWN patches
        M_numBoundaryDofVector[iFESpace] = M_dofGlobalIdVector[iFESpace].size();

        M_patchMeasureVector[iFESpace].resize ( M_numBoundaryDofVector[iFESpace] );
        for ( Vector::iterator it = M_patchMeasureVector[iFESpace].begin(); it < M_patchMeasureVector[iFESpace].end(); it++ )
        {
            *it = 0.0;
        }

        M_patchNormalVector[iFESpace].resize ( M_numBoundaryDofVector[iFESpace] * M_geoDimension );
        for ( Vector::iterator it = M_patchNormalVector[iFESpace].begin(); it < M_patchNormalVector[iFESpace].end(); it++ )
        {
            *it = 0.0;
        }

        M_patchIntegratedPhiVector[iFESpace].resize ( M_numBoundaryDofVector[iFESpace] );
        for ( Vector::iterator it = M_patchIntegratedPhiVector[iFESpace].begin();
                it < M_patchIntegratedPhiVector[iFESpace].end(); it++ )
        {
            *it = 0.0;
        }
    }
}

///////////////////////////////////////////////


// Measure of facets with a certain marker
template<typename MeshType>
Real PostProcessingBoundary<MeshType>::measure ( const markerID_Type& flag )
{
    // Each processor computes the measure across his own flagged facets --> measureScatter
    // At the end I'll reduce the process measures --> measure
    Real measureScatter (0.0), measure (0.);

    std::list<ID> facetList ( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    //
    // Loop on flagged processor facets
    //
    for ( Iterator j = facetList.begin(); j != facetList.end(); ++j )
    {

        M_currentBdFEPtrVector[0]->update ( M_meshPtr->boundaryFacet ( *j ), UPDATE_W_ROOT_DET_METRIC ); // updating finite element information

        measureScatter += M_currentBdFEPtrVector[0]->measure();

    }

    // reducing per-processor information
    M_epetraMapPtr->comm().SumAll ( &measureScatter, &measure, 1 );

    return measure;
}


// flux of vector field "field" through facets with a certain marker
template<typename MeshType>
template<typename VectorType>
Real PostProcessingBoundary<MeshType>::flux ( const VectorType& field, const markerID_Type& flag, UInt feSpace, UInt nDim )
{
    // Each processor computes the flux across his own flagged facets --> fluxScatter
    // At the end I'll reduce the process fluxes --> flux
    Real fluxScatter (0.0), flux (0.);

    // I need the global DOF ID to query the vector
    // dofVectorIndex is the index of the dof in the data structure of PostProcessingBoundary class
    // dofGlobalId is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt dofVectorIndex, dofGlobalId;

    // list of flagged facets on current processor
    std::list<ID> facetList ( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of field in the current facet
    Vector localFieldVector (nDim * M_numTotalDofPerFacetVector[feSpace]);

    // Loop on facetList
    for (Iterator j = facetList.begin(); j != facetList.end(); ++j)
    {

        // Updating quadrature data on the current facet
        M_currentBdFEPtrVector[feSpace]->update (M_meshPtr->boundaryFacet (*j), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS );

        // Quadrature formula
        // Loop on quadrature points
        for (UInt iq = 0; iq < M_currentBdFEPtrVector[feSpace]->nbQuadPt(); ++iq)
        {

            // Dot product
            // Loop on components
            for (UInt iComponent = 0; iComponent < nDim; ++iComponent)
            {

                // Interpolation
                // Loop on local dof
                for (ID iDof = 0; iDof < M_numTotalDofPerFacetVector[feSpace]; ++iDof)
                {

                    // Extracting nodal values of field in the current facet
                    dofVectorIndex = M_vectorNumberingPerFacetVector[feSpace][ ( UInt ) * j ][ iDof ];
                    dofGlobalId = M_dofGlobalIdVector[feSpace][dofVectorIndex]; // this is in the GLOBAL mesh

                    localFieldVector[iComponent * M_numTotalDofPerFacetVector[feSpace] + iDof] =
                        field[iComponent * M_numTotalDofVector[feSpace] + dofGlobalId];

                    fluxScatter += M_currentBdFEPtrVector[feSpace]->wRootDetMetric (iq)
                                   * localFieldVector[iComponent * M_numTotalDofPerFacetVector[feSpace] + iDof]
                                   * M_currentBdFEPtrVector[feSpace]->phi (iDof, iq)
                                   * M_currentBdFEPtrVector[feSpace]->normal (iComponent, iq);
                }
            }
        }
    }
    // Reducing per-processor values
    M_epetraMapPtr->comm().SumAll ( &fluxScatter, &flux, 1 );

    return flux;
}

template<typename MeshType>
template<typename VectorType>
Real PostProcessingBoundary<MeshType>::kineticNormalStress ( const VectorType& velocity, const Real& density, const markerID_Type& flag, UInt feSpace, UInt /*nDim*/ )
{
    // Each processor computes the quantities across his own flagged facets
    Real kineticNormalStressScatter (0.0), kineticNormalStress (0.0);
    Real areaScatter (0.0), area (0.0);
    Real temp (0.0);

    // Compute the normal
    Vector faceNormal = normal ( flag );

    // I need the global DOF ID to query the vector
    // dofVectorIndex is the index of the dof in the data structure of PostProcessingBoundary class
    // dofGlobalId is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt dofVectorIndex, dofGlobalId;

    // List of flagged facets on current processor
    std::list<ID> facetList ( M_boundaryMarkerToFacetIdMap[flag] );

    // Loop on facetList
    for ( std::list<ID>::iterator j (facetList.begin() ); j != facetList.end(); ++j )
    {
        // Updating quadrature data on the current facet
        M_currentBdFEPtrVector[feSpace]->update ( M_meshPtr->boundaryFacet ( *j ), UPDATE_NORMALS | UPDATE_W_ROOT_DET_METRIC  );

        // Computing the area
        areaScatter += M_currentBdFEPtrVector[0]->measure();

        // Quadrature formula (loop on quadrature points)
        for ( UInt iq (0); iq < M_currentBdFEPtrVector[feSpace]->nbQuadPt(); ++iq )
        {
            // Loop on local dof
            for ( ID iDof (0); iDof < M_numTotalDofPerFacetVector[feSpace]; ++iDof )
            {
                // Extracting nodal values of field in the current facet
                dofVectorIndex = M_vectorNumberingPerFacetVector[feSpace][ ( UInt ) * j ][ iDof ];
                dofGlobalId    = M_dofGlobalIdVector[feSpace][dofVectorIndex]; // this is in the GLOBAL mesh

                temp = velocity[0 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[0] // u_x * n_x
                       + velocity[1 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[1] // u_y * n_y
                       + velocity[2 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[2]; // u_z * n_z

                kineticNormalStressScatter += M_currentBdFEPtrVector[feSpace]->wRootDetMetric (iq)
                                              * M_currentBdFEPtrVector[feSpace]->phi (Int (iDof), iq)
                                              * temp * temp;
            }
        }
    }

    // Reducing per-processor values
    M_epetraMapPtr->comm().SumAll ( &areaScatter, &area, 1 );
    M_epetraMapPtr->comm().SumAll ( &kineticNormalStressScatter, &kineticNormalStress, 1 );

    return 0.5 * density * kineticNormalStress / area;
}

template<typename MeshType>
template<typename VectorType>
Real PostProcessingBoundary<MeshType>::kineticNormalStressDerivative ( const VectorType& velocity, const VectorType& velocityDerivative,
                                                                       const Real& density, const markerID_Type& flag, UInt feSpace, UInt /*nDim*/ )
{
    //TODO The two terms which depend on the displacement of the fluid have still to be coded

    // Each processor computes the quantities across his own flagged facets
    Real kineticNormalStressScatter (0.0), kineticNormalStress (0.0);
    Real areaScatter (0.0), area (0.0);
    Real temp (0.0);

    // Compute the normal
    Vector faceNormal = normal ( flag );

    // I need the global DOF ID to query the vector
    // dofVectorIndex is the index of the dof in the data structure of PostProcessingBoundary class
    // dofGlobalId is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt dofVectorIndex, dofGlobalId;

    // List of flagged facets on current processor
    std::list<ID> facetList ( M_boundaryMarkerToFacetIdMap[flag] );

    // Loop on facetList
    for ( std::list<ID>::iterator j (facetList.begin() ); j != facetList.end(); ++j )
    {
        // Updating quadrature data on the current facet
        M_currentBdFEPtrVector[feSpace]->update ( M_meshPtr->boundaryFacet ( *j ), UPDATE_NORMALS | UPDATE_W_ROOT_DET_METRIC );

        // Computing the area
        areaScatter += M_currentBdFEPtrVector[0]->measure();

        // Quadrature formula (loop on quadrature points)
        for ( UInt iq (0); iq < M_currentBdFEPtrVector[feSpace]->nbQuadPt(); ++iq )
        {
            // Loop on local dof
            for ( ID iDof (0); iDof < M_numTotalDofPerFacetVector[feSpace]; ++iDof )
            {
                // Extracting nodal values of field in the current facet
                dofVectorIndex = M_vectorNumberingPerFacetVector[feSpace][ ( UInt ) * j ][ iDof ];
                dofGlobalId    = M_dofGlobalIdVector[feSpace][dofVectorIndex]; // this is in the GLOBAL mesh

                temp = (
                           velocity[0 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[0] // u_x * n_x
                           + velocity[1 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[1] // u_y * n_y
                           + velocity[2 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[2] // u_z * n_z
                       )
                       * (
                           velocityDerivative[0 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[0] // du_x * n_x
                           + velocityDerivative[1 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[1] // du_y * n_y
                           + velocityDerivative[2 * M_numTotalDofVector[feSpace] + dofGlobalId] * faceNormal[2] // du_z * n_z
                       );

                kineticNormalStressScatter += M_currentBdFEPtrVector[feSpace]->wRootDetMetric (iq)
                                              * M_currentBdFEPtrVector[feSpace]->phi (Int (iDof), iq)
                                              * temp;
            }
        }
    }

    // Reducing per-processor values
    M_epetraMapPtr->comm().SumAll ( &areaScatter, &area, 1 );
    M_epetraMapPtr->comm().SumAll ( &kineticNormalStressScatter, &kineticNormalStress, 1 );

    return density * kineticNormalStress / area;
}

// Average value of field on facets with a certain marker
template<typename MeshType>
template<typename VectorType>
Vector PostProcessingBoundary<MeshType>::average ( const VectorType& field, const markerID_Type& flag, UInt feSpace, UInt nDim )
{
    // Each processor computes the average value on his own flagged facets --> fieldAverageScatter
    // At the end I'll reduce the process values --> fieldAverage
    Vector fieldAverageScatter (nDim), fieldAverage (nDim), localField (nDim);
    // basic policy for type VectorType: operator[] available
    for ( UInt iComponent = 0; iComponent < nDim; ++iComponent )
    {
        fieldAverageScatter[iComponent] = 0.;
        fieldAverage[iComponent] = 0.;
        localField[iComponent] = 0.;
    }

    // The total measure of the considered facets
    Real measureScatter (0.), measure;

    // I need the global Dof ID to query the Oseen solution vector
    // dofVectorIndex is the id of the DOF in the data structure of PostProcessingBoundary class
    // dofGlobalId is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt dofVectorIndex, dofGlobalId;

    // list of flagged facets on current processor
    std::list<ID> facetList ( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of field in the current facet
    Vector localFieldVector (M_numTotalDofPerFacetVector[feSpace]);

    // Loop on facets
    for (Iterator j = facetList.begin(); j != facetList.end(); ++j)
    {

        // basic policy for type VectorType: operator[] available
        for ( UInt iComponent = 0; iComponent < nDim; ++iComponent )
        {
            localField[iComponent] = 0.;
        }

        // Updating quadrature data on the current facet
        M_currentBdFEPtrVector[feSpace]->update (M_meshPtr->boundaryFacet (*j), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS );

        // Loop on components
        for (UInt iComponent = 0; iComponent < nDim; ++iComponent)
        {

            // Quadrature formula
            // Loop on quadrature points
            for (UInt iq = 0; iq < M_currentBdFEPtrVector[feSpace]->nbQuadPt(); ++iq)
            {

                // Interpolation
                // Loop on local dof
                for (ID iDof = 0; iDof < M_numTotalDofPerFacetVector[feSpace]; ++iDof)
                {

                    // Extracting nodal values of field in the current facet
                    dofVectorIndex = M_vectorNumberingPerFacetVector[feSpace][ ( UInt ) * j ][ iDof];
                    dofGlobalId = M_dofGlobalIdVector[feSpace][dofVectorIndex]; // this is in the GLOBAL mesh

                    // basic policy for type VectorType: operator[] available
                    localFieldVector[iDof] = field[iComponent * M_numTotalDofVector[feSpace] + dofGlobalId];

                    localField[iComponent] += M_currentBdFEPtrVector[feSpace]->wRootDetMetric (iq)
                                              * localFieldVector[iDof] * M_currentBdFEPtrVector[feSpace]->phi (Int (iDof), iq);
                }
            }
            // Computing the field integral over the boundary facets
            fieldAverageScatter[iComponent] += localField[iComponent];
        }

        // Computing the measure
        measureScatter += M_currentBdFEPtrVector[feSpace]->measure();
    }

    M_epetraMapPtr->comm().SumAll ( &measureScatter, &measure, 1 );

    // Reducing per-processor values
    for ( UInt iComponent (0); iComponent < nDim; ++iComponent )
    {
        M_epetraMapPtr->comm().SumAll ( &fieldAverageScatter[iComponent], &fieldAverage[iComponent], 1 );
    }

    return fieldAverage / measure;
}

// approximate normal for a certain marker
template<typename MeshType>
Vector PostProcessingBoundary<MeshType>::normal ( const markerID_Type& flag, UInt feSpace, UInt nDim )
{
    // Each processor computes the normal on his own flagged facets --> normalScatter
    // At the end I'll reduce the process normals --> normal
    Vector normalScatter (3) , normal (3);

    // Initialize vectors to zero
    normal[0] = 0.0;
    normal[1] = 0.0;
    normal[2] = 0.0;

    normalScatter[0] = 0.0;
    normalScatter[1] = 0.0;
    normalScatter[2] = 0.0;

    // list of flagged facets on current processor
    std::list<ID> facetList ( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of field in the current facet
    Vector localFieldVector ( nDim * M_numTotalDofPerFacetVector[feSpace] );

    // Loop on facetList
    for (Iterator j = facetList.begin(); j != facetList.end(); ++j)
    {
        // Updating quadrature data on the current facet
        M_currentBdFEPtrVector[feSpace]->update ( M_meshPtr->boundaryFacet (*j), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS );

        // Quadrature formula (loop on quadrature points)
        for ( UInt iq (0); iq < M_currentBdFEPtrVector[feSpace]->nbQuadPt(); ++iq )
        {
            // Dot product (loop on components)
            for ( UInt iComponent (0); iComponent < nDim; ++iComponent )
            {
                // Interpolation (loop on local dof)
                for (ID iDof (0); iDof < M_numTotalDofPerFacetVector[ feSpace ]; ++iDof)
                {
                    normalScatter (iComponent) += M_currentBdFEPtrVector[feSpace]->wRootDetMetric ( iq )
                                                  * M_currentBdFEPtrVector[feSpace]->phi ( Int (iDof), iq )
                                                  * M_currentBdFEPtrVector[feSpace]->normal ( Int (iComponent), iq );
                }
            }
        }
    }

    // Reducing per-processor values
    M_epetraMapPtr->comm().SumAll ( &normalScatter (0), &normal (0), 1 );
    M_epetraMapPtr->comm().SumAll ( &normalScatter (1), &normal (1), 1 );
    M_epetraMapPtr->comm().SumAll ( &normalScatter (2), &normal (2), 1 );

    // Scale normal to unity length
    Real nn = std::sqrt ( normal (0) * normal (0) + normal (1) * normal (1) + normal (2) * normal (2) );

#ifdef DEBUG
    if ( std::fabs ( nn ) < 1e-6 )
    {
        debugStream ( 5000 ) << "Approximate surface normal could not be reliably computed.\n";
        debugStream ( 5000 ) << "Modulus of the integrated normal vector was: " << nn  << "\n";
    }
#endif

    return ( normal / nn );
}

// approximate geometric center for a certain marker
template<typename MeshType>
Vector PostProcessingBoundary<MeshType>::geometricCenter ( const markerID_Type& flag, UInt feSpace, UInt nDim )
{
    // Each processor computes the geometric center on his own flagged facets --> geometricCenterScatter
    // At the end I'll reduce the process geometricCenterScatter --> geometricCenter
    Real areaScatter (0.0), area (0.0);
    Vector geometricCenterScatter (3), geometricCenter (3);

    geometricCenter[0] = 0.0;
    geometricCenter[1] = 0.0;
    geometricCenter[2] = 0.0;

    geometricCenterScatter[0] = 0.0;
    geometricCenterScatter[1] = 0.0;
    geometricCenterScatter[2] = 0.0;

    // list of flagged facets on current processor
    std::list<ID> facetList ( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of field in the current facet
    Vector localFieldVector (nDim * M_numTotalDofPerFacetVector[feSpace]);

    // Loop on facetList
    for (Iterator j = facetList.begin(); j != facetList.end(); ++j)
    {
        // Updating quadrature data on the current facet
        M_currentBdFEPtrVector[feSpace]->update (M_meshPtr->boundaryFacet (*j), UPDATE_QUAD_NODES | UPDATE_NORMALS | UPDATE_W_ROOT_DET_METRIC );

        // Compute the area of the facet
        areaScatter += M_currentBdFEPtrVector[0]->measure();

        // Quadrature formula (loop on quadrature points)
        for (UInt iq (0); iq < M_currentBdFEPtrVector[feSpace]->nbQuadPt(); ++iq)
        {
            // Dot product (loop on components)
            for (UInt iComponent (0); iComponent < nDim; ++iComponent)
            {
                // Interpolation (loop on local dof)
                for (ID iDof (0); iDof < M_numTotalDofPerFacetVector[feSpace]; ++iDof)
                {
                    geometricCenterScatter (iComponent) += M_currentBdFEPtrVector[feSpace]->wRootDetMetric (iq)
                                                           * M_currentBdFEPtrVector[feSpace]->phi (Int (iDof), iq)
                                                           * M_currentBdFEPtrVector[feSpace]->quadPt (iq, Int (iComponent) );
                }
            }
        }
    }

    // Reducing per-processor values
    M_epetraMapPtr->comm().SumAll ( &geometricCenterScatter (0), &geometricCenter (0), 1 );
    M_epetraMapPtr->comm().SumAll ( &geometricCenterScatter (1), &geometricCenter (1), 1 );
    M_epetraMapPtr->comm().SumAll ( &geometricCenterScatter (2), &geometricCenter (2), 1 );

    M_epetraMapPtr->comm().SumAll ( &areaScatter, &area, 1 );

    // Didive by the area
    geometricCenter /= area;

    return geometricCenter;
}


// Measure of patches on the boundary
template<typename MeshType>
void PostProcessingBoundary<MeshType>::computePatchesMeasure()
{
    for ( UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace )
    {

        // measure of the mesh facet
        Real localMeasure;

        // index of the considered dof in this class' vectors
        ID dofVectorIndex;

        // ===================================================
        // Loop on boundary facets
        // ===================================================
        for ( ID iboundaryFacet = 0 ; iboundaryFacet < M_numBoundaryFacets; ++iboundaryFacet )
        {
            // updating finite element information
            M_currentBdFEPtrVector[iFESpace]->update ( M_meshPtr->boundaryFacet ( iboundaryFacet ), UPDATE_W_ROOT_DET_METRIC );

            localMeasure = M_currentBdFEPtrVector[iFESpace]->measure();
            // Loop on the total DOF per Facet
            for ( ID iDof = 0; iDof < M_numTotalDofPerFacetVector[iFESpace]; ++iDof )
            {
                // Extracting local ID of iDof
                dofVectorIndex = M_vectorNumberingPerFacetVector[iFESpace][ iboundaryFacet ][ iDof];
                M_patchMeasureVector[iFESpace][dofVectorIndex] += localMeasure;
            }
        }
    }
}

template <typename MeshType>
void PostProcessingBoundary<MeshType>::showPatchesMeasure ( std::ostream& output ) const
{
    for ( UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace )
    {

        output << "\n***** Post Proc: Measure of the patches *****"
               << "\n\tfor variable " << iFESpace << std::endl;
        ID counter = 0;

        for ( Vector::const_iterator it = M_patchMeasureVector[iFESpace].begin();
                it != M_patchMeasureVector[iFESpace].end(); it++ )
        {
            output << "Boundary DOF: " << counter
                   << ", corresponding to Global DOF: " << M_dofGlobalIdVector[iFESpace][ counter ]
                   << " has patch measure: " << *it << std::endl;
            counter++;
        }
    }
}

template<typename MeshType>
void PostProcessingBoundary<MeshType>::showDOFIndexMap ( std::ostream& output ) const
{
    for ( UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace )
    {
        output << "\n***** Data for variable " << iFESpace << "*****" << std::endl;

        Int counter = 0;
        output << "\n***** Post Proc: Vector Indexes per Facet *****" << std::endl;
        output << M_vectorNumberingPerFacetVector[iFESpace].size() << std::endl;
        for ( std::vector<std::vector<ID> >::iterator it1 = M_vectorNumberingPerFacetVector[iFESpace].begin();
                it1 < M_vectorNumberingPerFacetVector[iFESpace].end(); it1++ )
        {
            counter++;
            output << "Boundary Facet " << counter << ", indexes: " << std::endl;
            for ( std::vector<ID>::iterator it2 = it1->begin(); it2 < it1->end(); it2++ )
            {
                output << *it2 << ",";
            }
            output << std::endl;
        }

        output << "***** Post Proc: From Vector Index to Global DOF *****" << std::endl;
        output << M_dofGlobalIdVector[iFESpace].size() << std::endl;

        for ( std::vector<ID>::iterator it3 = M_dofGlobalIdVector[iFESpace].begin();
                it3 < M_dofGlobalIdVector[iFESpace].end(); it3++ )
        {
            output << "Index :" << it3 - M_dofGlobalIdVector[iFESpace].begin()
                   << ", Global DOF: " << *it3 << std::endl;
        }
    }
}

/////////////////////////////////////////////////

///////////////////////////////////////////////


// Normal vectors of patches on the boundary
template<typename MeshType>
void PostProcessingBoundary<MeshType>::computePatchesNormal()
{
    for ( UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace )
    {

        // for each patch, the average of each component of the normal vector
        Real sum;

        // index of the considered dof in this class' vectors
        ID dofVectorIndex;

        // ===================================================
        // Loop on boundary facets
        // ===================================================
        for ( ID iboundaryFacet = 0 ; iboundaryFacet < M_numBoundaryFacets; ++iboundaryFacet )
        {
            // updating finite element information
            M_currentBdFEPtrVector[iFESpace]->update ( M_meshPtr->boundaryFacet ( iboundaryFacet ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS );

            // Loop on the components
            for ( Int iComponent = 0; iComponent < M_geoDimension; iComponent++ )
            {
                sum = 0.;
                // Loop on the quadrature points
                for ( UInt iQuadraturePoint = 0; iQuadraturePoint < M_currentBdFEPtrVector[iFESpace]->nbQuadPt();
                        ++iQuadraturePoint )
                {
                    sum += M_currentBdFEPtrVector[iFESpace]->normal ( iComponent, iQuadraturePoint ) *
                           M_currentBdFEPtrVector[iFESpace]->wRootDetMetric ( iQuadraturePoint );
                }
                for ( ID iDof = 0; iDof < M_numTotalDofPerFacetVector[iFESpace]; ++iDof )
                {
                    // Extracting local ID of iDof
                    dofVectorIndex = M_vectorNumberingPerFacetVector[iFESpace][ iboundaryFacet ][ iDof ];
                    M_patchNormalVector[iFESpace][ iComponent * M_numBoundaryDofVector[iFESpace] + dofVectorIndex ] += sum;
                }
            }
        }
        // Normalization of the averaged normals with the patch measure
        for ( UInt iBoundaryDof = 0; iBoundaryDof < M_numBoundaryDofVector[iFESpace]; ++iBoundaryDof )
        {
            Real localMeasure = M_patchMeasureVector[iFESpace][ iBoundaryDof ];
            for ( Int icc = 0; icc < M_geoDimension; icc++ )
            {
                M_patchNormalVector[iFESpace][ icc * M_numBoundaryDofVector[iFESpace] + iBoundaryDof ] *= 1. / localMeasure;
            }
        }
    }
}


template <typename MeshType>
void PostProcessingBoundary<MeshType>::showPatchesNormal ( std::ostream& output ) const
{
    for ( UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace )
    {

        output << "\n***** Post Proc: Normal vector on the patches *****"
               << "\n\tfor variable " << iFESpace << std::endl;

        ID counter = 0;

        for ( Vector::const_iterator it = M_patchMeasureVector[iFESpace].begin();
                it != M_patchMeasureVector[iFESpace].end(); it++ )
        {
            output << "Boundary DOF: " << counter
                   << ", corresponding to Global DOF: " << M_dofGlobalIdVector[iFESpace][ counter ]
                   << " has patch measure: " << *it << std::endl;
            output << "and normal components " ;
            for ( Int iComponent = 0; iComponent < M_geoDimension; iComponent++ )
                output <<
                       M_patchNormalVector[iFESpace][ iComponent * M_numBoundaryDofVector[iFESpace] + counter ] << " ";

            output << std::endl;
            counter++;
        }
    }

    output << "End SHOW NORMAL" << std::endl;
}

//////////////////////////////////////////////////
//
///////////////////////////////////////////////////

// Vector with the integral of the shape functions on the patches on the boundary
template<typename MeshType>
void PostProcessingBoundary<MeshType>::computePatchesPhi()
{
    for ( UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace )
    {

        // sum contributions from each facet of the patch
        Real sum;

        // ID of the considered dof in this class' vectors
        ID dofVectorIndex;

        // ===================================================
        // Loop on boundary facets
        // ===================================================
        for ( ID iboundaryFacet = 0 ; iboundaryFacet < M_numBoundaryFacets; ++iboundaryFacet )
        {
            // updating finite element information
            M_currentBdFEPtrVector[iFESpace]->update ( M_meshPtr->boundaryFacet ( iboundaryFacet ), UPDATE_W_ROOT_DET_METRIC );

            for ( ID iDof = 0; iDof < M_numTotalDofPerFacetVector[iFESpace]; ++iDof )
            {
                sum = 0.0;
                // global dof
                dofVectorIndex = M_vectorNumberingPerFacetVector[iFESpace][ iboundaryFacet ][ iDof ];

                // Loop on the quadrature points
                for ( UInt iQuadraturePoint = 0; iQuadraturePoint < M_currentBdFEPtrVector[iFESpace]->nbQuadPt();
                        ++iQuadraturePoint )
                {
                    sum += M_currentBdFEPtrVector[iFESpace]->phi ( ( Int ) ( iDof), iQuadraturePoint )
                           * M_currentBdFEPtrVector[iFESpace]->wRootDetMetric ( iQuadraturePoint );
                }
                M_patchIntegratedPhiVector[iFESpace][ dofVectorIndex ] += sum;
            }
        }
    }
}


template <typename MeshType>
void PostProcessingBoundary<MeshType>::showPatchesPhi ( std::ostream& output ) const
{
    for ( UInt iFESpace = 0; iFESpace < M_numFESpaces; ++iFESpace )
    {

        output << "\n***** Post Proc: Average phi on the patches *****"
               << "\n\tfor variable " << iFESpace << std::endl;

        ID counter = 0;

        for ( Vector::const_iterator it = M_patchMeasureVector[iFESpace].begin();
                it != M_patchMeasureVector[iFESpace].end(); it++ )
        {
            output << "Boundary DOF: " << counter
                   << ", corresponding to Global DOF: " << M_dofGlobalIdVector[iFESpace][ counter ]
                   << " has patch measure: " << *it << std::endl;
            output << "and average phi  " << M_patchIntegratedPhiVector[iFESpace][ counter ] << std::endl ;
            counter++;
        }
    }
    output << "End SHOW PHI" << std::endl;
}

// the following part is definitely not ready yet TP 09/2008
#if 0
////////////////////////////////
////////////////////////////////
////////////////////////////////
template<typename MeshType>
template<typename VectorType>
VectorType PostProcessingBoundary<MeshType>::compute_sstress ( const VectorType& r, UInt ncomp, bool residual = true ) const
{
    ASSERT ( ncomp = M_geoDimension, "Error: Shear stress computation possible only for vector unknowns" );

    // prepare the vectors
    VectorType stress ( M_numBoundaryDofVector[0] * M_geoDimension );
    stress.clear();
    VectorType nstress ( M_numBoundaryDofVector[0] * M_geoDimension );
    nstress.clear();
    VectorType sstress ( M_numBoundaryDofVector[0] * M_geoDimension );
    sstress.clear();
    ID counter, dofGlobalId;

    // number of DOFs for each component
    UInt dim = r.size() / ncomp;

    // helper structures to avoid "if" statement
    // if residual==true, vector r has ncomp*M_numBoundaryDofVector components
    // if residual==false, vector r has ncomp*M_numTotalDofVector components
    Vector coef (2);
    std::vector<ID> index (2);

    // loop on locally stored DOFs
    for ( counter = 0; counter < M_numBoundaryDofVector[0]; ++counter )
    {
        // find global ID for boundary dof
        dofGlobalId = M_dofGlobalIdVector[0][counter]; // this is in the GLOBAL mesh

        coef[0] = 1.;
        index[0] = counter;
        coef[1] = 1. / M_patchIntegratedPhiVector[0][counter];
        index[1] = dofGlobalId;

        for ( UInt ind_comp = 0; ind_comp < ncomp; ++ind_comp )
        {
            stress[ counter + ind_comp * M_numBoundaryDofVector[0] ] = //damned conventions trouble : 0 or 1
                coef[residual] * r[ index[residual] + ind_comp * dim ];
        }
    }

    Real sn = 0.;
    ///// Normal stress
    for ( counter = 0; counter < M_numBoundaryDofVector[0]; counter++ )
    {
        sn = 0.;
        for ( UInt ind_comp = 0; ind_comp < ncomp; ind_comp++ )
            sn += stress[ counter + ind_comp * M_numBoundaryDofVector[0] ] *
                  M_patchNormalVector[0][ counter + ind_comp * M_numBoundaryDofVector[0] ];

        for ( UInt ind_comp = 0; ind_comp < ncomp; ind_comp++ )
            nstress[ counter + ind_comp * M_numBoundaryDofVector[0] ] = sn *
                                                                        M_patchNormalVector[0][ counter + ind_comp * M_numBoundaryDofVector[0] ];
    }

    // Shear Stress: this vector lives on the patches (ncomp*M_numBoundaryDofVector components)
    sstress = stress - nstress;

    return sstress;
} // compute_sstress



////////////////////////////////
////////////////////////////////
///////////////////////////////
template<typename MeshType>
template<typename VectorType>
VectorType PostProcessingBoundary<MeshType>::compute_stress ( const VectorType& grad, const UInt& dim,
                                                              const Real& visc ) const
{
    VectorType stress ( M_numBoundaryDofVector[0] * M_geoDimension );
    stress.clear();

    ID dofGlobalId;

    // cycle over boundary dof
    for ( ID bound_dof = 0; bound_dof < M_numBoundaryDofVector[0]; ++bound_dof )
    {
        // find global ID for boundary dof
        dofGlobalId = M_dofGlobalIdVector[0][bound_dof]; // this is in the GLOBAL mesh

        // cycle over stress components
        for ( UInt scomp = 0; scomp < M_geoDimension; ++scomp )

            // cycle over normal components
            for ( UInt ncomp = 0; ncomp < M_geoDimension; ++ncomp )
            {
                stress[ bound_dof + scomp * M_numBoundaryDofVector[0] ] += visc *
                                                                           // grad!
                                                                           ( grad[ dofGlobalId + ncomp * dim + scomp * nDimensions * dim ] +
                                                                             // transpose grad!
                                                                             grad[ dofGlobalId + scomp * dim + ncomp * nDimensions * dim ] )
                                                                           * M_patchNormalVector[ bound_dof + ncomp * M_numBoundaryDofVector[0] ];
            }
    }


    return stress;
} // compute_stress
#endif
} // namespace LifeV

#endif /* POSTPROCESSINGBOUNDARY_H */

