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

#ifndef POSTPROC_H
#define POSTPROC_H 1

#include <string>
#include <iostream>
#include <sstream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

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
class PostProc
{

    //! @name Public Types
    //@{
    typedef typename MeshType::VolumeShape                   elementGeometricShape_Type;
    typedef typename elementGeometricShape_Type::GeoBShape   facetGeometricShape_Type;
    typedef MeshType                                         mesh_Type;
    typedef boost::shared_ptr<MeshType>                      meshPtr_Type;
    typedef CurrentBdFE*                                     currentBdFEPtr_Type;
    typedef Dof*                                             dofPtr_Type;
    //@}

public:

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    PostProc<MeshType>()
    {}

    //! Copy constructor
    /*!
        @note all pointer members are copied (no new allocation is done in the class)
     */
    PostProc( const PostProc& postProc )
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
    PostProc<MeshType>( meshPtr_Type mesh,
                    std::vector<currentBdFEPtr_Type > currentBdFEVector,
                    std::vector<dofPtr_Type > dofVector,
                    const EpetraMap& epetraMap, UInt nFESpaces = 1 );

    /*!
     \brief Constructor for the case in which we have only one currentBdFE and dof
     */
    PostProc<MeshType>( meshPtr_Type mesh,
                    currentBdFEPtr_Type currentBdFE, dofPtr_Type dof,
                    const EpetraMap& epetraMap );

    /*!
     \brief Constructor for the case in which we have two currentBdFE's and dof's

     This is the case of NS problem, in which we will want to consider both
     velocity and pressure dof/fe classes
     */
    PostProc<MeshType>( meshPtr_Type mesh,
                    currentBdFEPtr_Type feBdu, dofPtr_Type dofu,
                    currentBdFEPtr_Type feBdp, dofPtr_Type dofp,
                    const EpetraMap& epetraMap );

    //@}

    /*! @name Methods */
    //@{
    /*!
     These methods build the data structures describing patches of boundary elements
     (patch measure, patch normal and integral of the test function over the patch)
     */
    void set_area();
    void set_normal();
    void set_phi();

    /*! @defgroup boundary_methods
     These methods compute quantities on boundary sections associated to a given flag
     */

    /*!
       This method computes the measure of boundary section "flag"
       @ingroup boundary_methods

       \param flag is the marker of the considered boundary section
     */
    Real area( const entityFlag_Type& flag );

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
    Real flux( const VectorType& vectorField, const entityFlag_Type& flag, UInt feSpace = 0, UInt nDim = nDimensions);

    /*!
       This method computes the average value of a field on the boundary section "flag"
       @ingroup boundary_methods

      \tparam VectorType Vector type. Basic policy for type VectorType: operator[] available

      \param field is intended to be a vector or a scalar (pressure in NS problem),
       this method computes the average value of field on section "flag"

       \return the averaged vector
     */
    template< typename VectorType >
    Vector average( const VectorType& field, const entityFlag_Type& flag, UInt feSpace = 0, UInt nDim = 1);

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
    VectorType compute_sstress( const VectorType& r, UInt ncomp, bool residual ) const;

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
    VectorType compute_stress( const VectorType& grad, const UInt& dim, const Real& visc  ) const;
#endif
    //@}

    /*! @name Methods */
    //@{
    /*!
     These methods print to screen information about the class data structures
     */
    void show_bdLtoG();

    void show_area();

    void show_normal();

    void show_phi();
    //@}

    /*! @name Getters */
    //@{
    /*!
     Access to private members
     */
    Vector patch_area() const { return M_patchMeasureVector; }

    Vector patch_normal() const { return M_patchNormalVector; }

    std::vector< ID > fBdToIn() const { return M_dofGlobalIdVector; }

    //! Number of boundary DOF for the mesh at hand
    UInt nBdDof() const { return M_numBoundaryDofVector; }
    //@}


private:
    void                                         build_vectors();

    UInt                                         M_numFESpaces;

    std::vector<UInt>                            M_numBoundaryDofVector;
    std::vector<UInt>                            M_numDofPerVertexVector, M_numDofPerEdgeVector, M_numDofPerFaceVector;
    UInt                                         M_numVerticesPerFace, M_numEdgesPerFace;
    UInt                                         M_numVerticesPerElement, M_numEdgesPerElement;
    std::vector<UInt>                            M_numVertexDofPerFaceVector, M_numEdgeDofPerFaceVector;
    std::vector<UInt>                            M_numTotalDofPerFaceVector;
    std::vector<UInt>                            M_numVertexDofPerElement, M_numEdgeDofPerElementVector;
    UInt                                         M_numBoundaryFaces;
    std::vector<UInt>                            M_numTotalDofVector;
    // vector whose ith-component is the area of the patch of the ith-node on the boundary
    std::vector<Vector >              M_patchMeasureVector;
    // vector with the components of the average normal vector of each boundary node
    std::vector<Vector >              M_patchNormalVector;
    // vector with \int_{Patch(i)} \phi_i dx where \phi_i is the basis function.
    std::vector<Vector >              M_patchIntegratedPhiVector;

    // store once for all a map, with key={boundary flag}, value={ID list}
    std::map< entityFlag_Type, std::list<ID> >        M_boundaryMarkerToFacetIdMap;

    // for each boundary face, it contains the numbering of the dof of the face
    std::vector< std::vector< SimpleVect<ID> > > M_vectorNumberingPerFacetVector;
    // it converts from a local numbering over the boundary faces on the global numbering of the mesh
    std::vector< std::vector< ID > >             M_dofGlobalIdVector;

    // reference to the current boundary FE
    std::vector<currentBdFEPtr_Type >            M_currentBdFEPtrVector;
    // reference to a Dof class
    std::vector<dofPtr_Type >                    M_dofPtrVector;
    // pointer to the mesh
    meshPtr_Type                                 M_meshPtr;
    // pointer to the processor mapping
    boost::shared_ptr<EpetraMap>                 M_epetraMapPtr;

};

//
// IMPLEMENTATIONS
//
template <typename MeshType>
PostProc<MeshType>::PostProc( meshPtr_Type meshPtr,
                          std::vector<currentBdFEPtr_Type > currentBdFEVector,
                          std::vector<dofPtr_Type > dofVector,
                          const EpetraMap& epetraMap, UInt nvar ) :
        M_numFESpaces(nvar),
        M_numBoundaryDofVector(M_numFESpaces),
        M_numDofPerVertexVector(M_numFESpaces), M_numDofPerEdgeVector(M_numFESpaces), M_numDofPerFaceVector(M_numFESpaces),
        M_numVertexDofPerFaceVector(M_numFESpaces), M_numEdgeDofPerFaceVector(M_numFESpaces),
        M_numTotalDofPerFaceVector(M_numFESpaces),
        M_numVertexDofPerElement(M_numFESpaces), M_numEdgeDofPerElementVector(M_numFESpaces),
        M_numTotalDofVector(M_numFESpaces),
        M_patchMeasureVector(M_numFESpaces), M_patchNormalVector(M_numFESpaces),
        M_patchIntegratedPhiVector(M_numFESpaces),
        M_vectorNumberingPerFacetVector(M_numFESpaces), M_dofGlobalIdVector(M_numFESpaces),
        M_currentBdFEPtrVector(currentBdFEVector), M_dofPtrVector(dofVector),
        M_meshPtr( meshPtr ), M_epetraMapPtr( new EpetraMap(epetraMap) )
{
    for (UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace)
    {
        M_vectorNumberingPerFacetVector[M_numFESpaces].clear();
        M_dofGlobalIdVector[iFESpace].clear();
    }
    // Some useful local variables, to save some typing
    M_numVerticesPerFace = facetGeometricShape_Type::numVertices; // Number of face's vertices
    M_numEdgesPerFace = facetGeometricShape_Type::numEdges;    // Number of face's edges

    M_numVerticesPerElement = elementGeometricShape_Type::numVertices; // Number of element's vertices
    M_numEdgesPerElement = elementGeometricShape_Type::numEdges;    // Number of element's edges

    M_numBoundaryFaces = M_meshPtr->numBFaces();    // number of faces on boundary

    // Construction of M_vectorNumberingPerFacetVector & other data structures
    build_vectors();

    set_area();
    set_normal();
    set_phi();
}

template <typename MeshType>
PostProc<MeshType>::PostProc( meshPtr_Type mesh,
                          currentBdFEPtr_Type currentBdFE, dofPtr_Type dof,
                          const EpetraMap& epetraMap ) :
        M_numFESpaces(1),
        M_numBoundaryDofVector(M_numFESpaces),
        M_numDofPerVertexVector(M_numFESpaces), M_numDofPerEdgeVector(M_numFESpaces), M_numDofPerFaceVector(M_numFESpaces),
        M_numVertexDofPerFaceVector(M_numFESpaces), M_numEdgeDofPerFaceVector(M_numFESpaces),
        M_numTotalDofPerFaceVector(M_numFESpaces),
        M_numVertexDofPerElement(M_numFESpaces), M_numEdgeDofPerElementVector(M_numFESpaces),
        M_numTotalDofVector(M_numFESpaces),
        M_patchMeasureVector(M_numFESpaces), M_patchNormalVector(M_numFESpaces), M_patchIntegratedPhiVector(M_numFESpaces),
        M_vectorNumberingPerFacetVector(M_numFESpaces), M_dofGlobalIdVector(M_numFESpaces),
        M_currentBdFEPtrVector(M_numFESpaces), M_dofPtrVector(M_numFESpaces),
        M_meshPtr( mesh ), M_epetraMapPtr( new EpetraMap(epetraMap) )
{
    M_currentBdFEPtrVector[0]=currentBdFE;
    M_dofPtrVector[0]=dof;

    // Some useful local variables, to save some typing
    M_numVerticesPerFace = facetGeometricShape_Type::numVertices; // Number of face's vertices
    M_numEdgesPerFace = facetGeometricShape_Type::numEdges;    // Number of face's edges

    M_numVerticesPerElement = elementGeometricShape_Type::numVertices; // Number of element's vertices
    M_numEdgesPerElement = elementGeometricShape_Type::numEdges;    // Number of element's edges

    M_numBoundaryFaces = M_meshPtr->numBFaces();    // number of faces on boundary

    // Construction of M_vectorNumberingPerFacetVector & other data structures
    build_vectors();

    set_area();
    set_normal();
    set_phi();
}

template <typename MeshType>
PostProc<MeshType>::PostProc( meshPtr_Type mesh,
                          currentBdFEPtr_Type feBdu, dofPtr_Type dofu,
                          currentBdFEPtr_Type feBdp, dofPtr_Type dofp,
                          const EpetraMap& epetraMap ) :
        M_numFESpaces(2),
        M_numBoundaryDofVector(M_numFESpaces),
        M_numDofPerVertexVector(M_numFESpaces), M_numDofPerEdgeVector(M_numFESpaces), M_numDofPerFaceVector(M_numFESpaces),
        M_numVertexDofPerFaceVector(M_numFESpaces), M_numEdgeDofPerFaceVector(M_numFESpaces),
        M_numTotalDofPerFaceVector(M_numFESpaces),
        M_numVertexDofPerElement(M_numFESpaces), M_numEdgeDofPerElementVector(M_numFESpaces),
        M_numTotalDofVector(M_numFESpaces),
        M_patchMeasureVector(M_numFESpaces), M_patchNormalVector(M_numFESpaces), M_patchIntegratedPhiVector(M_numFESpaces),
        M_vectorNumberingPerFacetVector(M_numFESpaces), M_dofGlobalIdVector(M_numFESpaces),
        M_currentBdFEPtrVector(M_numFESpaces), M_dofPtrVector(M_numFESpaces),
        M_meshPtr( mesh ), M_epetraMapPtr( new EpetraMap(epetraMap) )
{
    M_currentBdFEPtrVector[0] = feBdu;
    M_dofPtrVector[0] = dofu;
    M_currentBdFEPtrVector[1] = feBdp;
    M_dofPtrVector[1] = dofp;

    // Some useful local variables, to save some typing
    M_numVerticesPerFace = facetGeometricShape_Type::numVertices; // Number of face's vertices
    M_numEdgesPerFace = facetGeometricShape_Type::numEdges;    // Number of face's edges

    M_numVerticesPerElement = elementGeometricShape_Type::numVertices; // Number of element's vertices
    M_numEdgesPerElement = elementGeometricShape_Type::numEdges;    // Number of element's edges

    M_numBoundaryFaces = M_meshPtr->numBFaces();    // number of faces on boundary

    // Construction of M_vectorNumberingPerFacetVector & other data structures
    build_vectors();

    set_area();
    set_normal();
    set_phi();
}

// Area of faces with a certain marker
template<typename MeshType>
void PostProc<MeshType>::build_vectors()
{
    SimpleVect<ID>            boundaryDofGlobalIdVector;

    UInt                      iFirstAdjacentElement, iVertexLocalId, iFaceLocalId, iEdgeLocalId;
    ID                        dofLocalId, dofGlobalId, dofAuxiliaryId;
    std::vector<ID>           numBoundaryDofVector(M_numFESpaces);
    std::vector<ID>::iterator dofGlobalIdVectorIterator;
    entityFlag_Type                boundaryFlag;

    for (UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace)
    {

        numBoundaryDofVector[iFESpace] = 1;

        // Some useful local variables, to save some typing
        M_numDofPerVertexVector[iFESpace] = M_dofPtrVector[iFESpace]->localDofPattern().nbDofPerVertex(); // number of Dof per vertices
        M_numDofPerEdgeVector[iFESpace] = M_dofPtrVector[iFESpace]->localDofPattern().nbDofPerEdge();   // number of Dof per edges
        M_numDofPerFaceVector[iFESpace] = M_dofPtrVector[iFESpace]->localDofPattern().nbDofPerFace();   // number of Dof per faces

        // number of vertex's Dof on a face
        M_numVertexDofPerFaceVector[iFESpace] = M_numDofPerVertexVector[iFESpace] * M_numVerticesPerFace;
        // number of edge's Dof on a face
        M_numEdgeDofPerFaceVector[iFESpace] = M_numDofPerEdgeVector[iFESpace] * M_numEdgesPerFace;

        // number of total Dof on a face
        M_numTotalDofPerFaceVector[iFESpace] =
        		M_numVertexDofPerFaceVector[iFESpace] + M_numEdgeDofPerFaceVector[iFESpace] + M_numDofPerFaceVector[iFESpace];
        // number of vertex's Dof on a Element
        M_numVertexDofPerElement[iFESpace] = M_numVerticesPerElement * M_numDofPerVertexVector[iFESpace];
        // number of edge's Dof on a Element
        M_numEdgeDofPerElementVector[iFESpace] = M_numEdgesPerElement * M_numDofPerEdgeVector[iFESpace];

        M_numTotalDofVector[iFESpace] = M_dofPtrVector[iFESpace]->numTotalDof();
    }

    // ===================================================
    // Loop on boundary faces
    // ===================================================
    for ( ID iboundaryFace = 1 ; iboundaryFace <= M_numBoundaryFaces; ++iboundaryFace )
    {

        iFirstAdjacentElement = M_meshPtr->bElement( iboundaryFace ).ad_first();  // id of the element adjacent to the face
        iFaceLocalId = M_meshPtr->bElement( iboundaryFace ).pos_first(); // local id of the face in its adjacent element

        boundaryFlag = M_meshPtr->bElement(iboundaryFace).marker();
        M_boundaryMarkerToFacetIdMap[boundaryFlag].push_back( iboundaryFace ); // fill the flag-to-faceIdList map

        for (UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace)
        {

            boundaryDofGlobalIdVector.clearVector();
            boundaryDofGlobalIdVector.resize( M_numTotalDofPerFaceVector[iFESpace] );

            // updating finite element information
            M_currentBdFEPtrVector[iFESpace]->updateMeas( M_meshPtr->bElement( iboundaryFace ) );

            // ===================================================
            // Vertex based Dof
            // ===================================================
            if ( M_numDofPerVertexVector[iFESpace] )
            {

                // loop on face vertices
                for ( ID iVertexPerFace = 1; iVertexPerFace <= M_numVerticesPerFace; ++iVertexPerFace )
                {
                	// local vertex number (in element)
                    iVertexLocalId = elementGeometricShape_Type::fToP( iFaceLocalId, iVertexPerFace );

                    // Loop number of Dof per vertex
                    for ( ID localDof = 1; localDof <= M_numDofPerVertexVector[iFESpace]; ++localDof )
                    {
                        dofLocalId = ( iVertexPerFace - 1 ) * M_numDofPerVertexVector[iFESpace] + localDof ; // local Dof
                        dofGlobalId = M_dofPtrVector[iFESpace]->localToGlobal(
                        		iFirstAdjacentElement, ( iVertexLocalId - 1 ) * M_numDofPerVertexVector[iFESpace] + localDof ); // global Dof
                        dofGlobalIdVectorIterator = find(
                        		M_dofGlobalIdVector[iFESpace].begin(), M_dofGlobalIdVector[iFESpace].end(), dofGlobalId );
                        if ( dofGlobalIdVectorIterator == M_dofGlobalIdVector[iFESpace].end() )
                        { // the dofGlobalId has been encountered for the first time
                            boundaryDofGlobalIdVector( dofLocalId ) = numBoundaryDofVector[iFESpace];
                            M_dofGlobalIdVector[iFESpace].push_back( dofGlobalId ); // local to boundary global on this face
                            numBoundaryDofVector[iFESpace]++;
                        }
                        else
                        { // the dofGlobalId has been already inserted in the M_dofGlobalIdVector vector
                            dofAuxiliaryId = ( ID ) ( ( dofGlobalIdVectorIterator - M_dofGlobalIdVector[iFESpace].begin() ) ) + 1;
                            boundaryDofGlobalIdVector( dofLocalId ) = dofAuxiliaryId; // local to boundary global on this face
                        }
                    }
                }
            }
            // ===================================================
            // Edge based Dof
            // ===================================================
            if ( M_numDofPerEdgeVector[iFESpace] )
            {
                // loop on face edges
                for ( ID iEdgePerFace = 1; iEdgePerFace <= M_numVerticesPerFace; ++iEdgePerFace )
                {
                	// local edge number (in element)
                    iEdgeLocalId = elementGeometricShape_Type::fToE( iFaceLocalId, iEdgePerFace ).first;
                    // Loop number of Dof per edge
                    for ( ID localDof = 1; localDof <= M_numDofPerEdgeVector[iFESpace]; ++localDof )
                    {
                        dofLocalId = M_numVertexDofPerFaceVector[iFESpace] +
                        		( iEdgePerFace - 1 ) * M_numDofPerEdgeVector[iFESpace] + localDof ; // local Dof
                        dofGlobalId = M_dofPtrVector[iFESpace]->localToGlobal(
                        		iFirstAdjacentElement, M_numVertexDofPerElement[iFESpace] + ( iEdgeLocalId - 1 ) *
                        		                       M_numDofPerEdgeVector[iFESpace] + localDof ); // global Dof
                        dofGlobalIdVectorIterator = find(
                        		M_dofGlobalIdVector[iFESpace].begin(), M_dofGlobalIdVector[iFESpace].end(), dofGlobalId );
                        if ( dofGlobalIdVectorIterator == M_dofGlobalIdVector[iFESpace].end() )
                        { // the dofGlobalId has been encountered for the first time
                            boundaryDofGlobalIdVector( dofLocalId ) = numBoundaryDofVector[iFESpace];
                            M_dofGlobalIdVector[iFESpace].push_back( dofGlobalId ); // local to boundary global on this face
                            numBoundaryDofVector[iFESpace]++;
                        }
                        else
                        { // the dofGlobalId has been already inserted in the M_dofGlobalIdVector vector
                            dofAuxiliaryId = ( ID ) ( dofGlobalIdVectorIterator - M_dofGlobalIdVector[iFESpace].begin() ) + 1;
                            boundaryDofGlobalIdVector( dofLocalId ) = dofAuxiliaryId; // local to boundary global on this face
                        }
                    }
                }
            }
            // ===================================================
            // Face based Dof
            // ===================================================
            // Loop on number of Dof per face
            for ( ID localDof = 1; localDof <= M_numDofPerFaceVector[iFESpace]; ++localDof )
            {
            	// local Dof
                dofLocalId = M_numEdgeDofPerFaceVector[iFESpace] + M_numVertexDofPerFaceVector[iFESpace] + localDof;
                dofGlobalId = M_dofPtrVector[iFESpace]->localToGlobal(
                		iFirstAdjacentElement, M_numEdgeDofPerElementVector[iFESpace] + M_numVertexDofPerElement[iFESpace] +
                		( iFaceLocalId - 1 ) * M_numDofPerFaceVector[iFESpace] + localDof ); // global Dof
                dofGlobalIdVectorIterator = find(
                		M_dofGlobalIdVector[iFESpace].begin(), M_dofGlobalIdVector[iFESpace].end(), dofGlobalId );
                if ( dofGlobalIdVectorIterator == M_dofGlobalIdVector[iFESpace].end() )
                { // the dofGlobalId has been encountered for the first time
                    boundaryDofGlobalIdVector( dofLocalId ) = numBoundaryDofVector[iFESpace];
                    M_dofGlobalIdVector[iFESpace].push_back( dofGlobalId ); // local to boundary global on this face
                    numBoundaryDofVector[iFESpace]++;
                }
                else
                { // the dofGlobalId has been already inserted in the M_dofGlobalIdVector vector
                    dofAuxiliaryId = ( ID ) ( dofGlobalIdVectorIterator - M_dofGlobalIdVector[iFESpace].begin() ) + 1;
                    boundaryDofGlobalIdVector( dofLocalId ) = dofAuxiliaryId; // local to boundary global on this face
                }
            }

            M_vectorNumberingPerFacetVector[iFESpace].push_back( boundaryDofGlobalIdVector );
        }
    }

    for (UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace)
    {
        // each processor holds information on HIS OWN patches
        M_numBoundaryDofVector[iFESpace] = M_dofGlobalIdVector[iFESpace].size();

        M_patchMeasureVector[iFESpace].resize( M_numBoundaryDofVector[iFESpace] );
        for ( Vector::iterator it = M_patchMeasureVector[iFESpace].begin(); it<M_patchMeasureVector[iFESpace].end(); it++ )
            *it = 0.0;

        M_patchNormalVector[iFESpace].resize( M_numBoundaryDofVector[iFESpace] * NDIM );
        for ( Vector::iterator it = M_patchNormalVector[iFESpace].begin(); it<M_patchNormalVector[iFESpace].end(); it++ )
            *it = 0.0;

        M_patchIntegratedPhiVector[iFESpace].resize( M_numBoundaryDofVector[iFESpace] );
        for ( Vector::iterator it = M_patchIntegratedPhiVector[iFESpace].begin();
        		it<M_patchIntegratedPhiVector[iFESpace].end(); it++ )
            *it = 0.0;
    }
}

///////////////////////////////////////////////


// Area of faces with a certain marker
template<typename MeshType>
Real PostProc<MeshType>::area( const entityFlag_Type& flag )
{
    // Each processor computes the area across his own flagged faces --> areaScatter
    // At the end I'll reduce the process areas --> area
    Real areaScatter(0.0), area(0.);

    std::list<ID> faceList( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    //
    // Loop on flagged processor faces
    //
    for (Iterator j=faceList.begin(); j != faceList.end(); ++j)
    {

        M_currentBdFEPtrVector[0]->updateMeas( M_meshPtr->bElement( *j ) );  // updating finite element information

        areaScatter += M_currentBdFEPtrVector[0]->measure();

    }

    // reducing per-processor information
    M_epetraMapPtr->Comm().SumAll( &areaScatter, &area, 1 );

    return area;
}


// flux of vector field "field" through faces with a certain marker
template<typename MeshType>
template<typename VectorType>
Real PostProc<MeshType>::flux( const VectorType& field, const entityFlag_Type& flag, UInt feSpace,
                           UInt nDim )
{
    // Each processor computes the flux across his own flagged faces --> fluxScatter
    // At the end I'll reduce the process fluxes --> flux
    Real fluxScatter(0.0), flux(0.);

    // I need the global Dof ID to query the vector
    // dofVectorIndex is the index of the dof in the data structure of PostProc class
    // dofGlobalId is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt dofVectorIndex, dofGlobalId;

    // list of flagged faces on current processor
    std::list<ID> faceList( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of field in the current face
    Vector localFieldVector(nDim * M_numTotalDofPerFaceVector[feSpace]);

    // Loop on faceList
    for (Iterator j=faceList.begin(); j != faceList.end(); ++j)
    {

        // Updating quadrature data on the current face
        M_currentBdFEPtrVector[feSpace]->updateMeasNormalQuadPt(M_meshPtr->bElement(*j));

        // Quadrature formula
        // Loop on quadrature points
        for (Int iq=0; iq< M_currentBdFEPtrVector[feSpace]->nbQuadPt; ++iq)
        {

            // Dot product
            // Loop on components
            for (UInt iComponent =0; iComponent<nDim; ++iComponent)
            {

                // Interpolation
                // Loop on local dof
                for (ID iDof=1; iDof<=M_numTotalDofPerFaceVector[feSpace]; ++iDof)
                {

                    // Extracting nodal values of field in the current face
                    dofVectorIndex = M_vectorNumberingPerFacetVector[feSpace][ ( UInt ) *j - 1 ][ iDof - 1 ];
                    dofGlobalId = M_dofGlobalIdVector[feSpace][dofVectorIndex-1]; // this is in the GLOBAL mesh

                    localFieldVector[iComponent*M_numTotalDofPerFaceVector[feSpace]+iDof-1] =
                    		field[iComponent*M_numTotalDofVector[feSpace]+dofGlobalId];

                    fluxScatter += M_currentBdFEPtrVector[feSpace]->weightMeas(iq)
                                    * localFieldVector[iComponent*M_numTotalDofPerFaceVector[feSpace]+iDof-1]
                                    * M_currentBdFEPtrVector[feSpace]->phi(Int(iDof-1),iq)
                                    * M_currentBdFEPtrVector[feSpace]->normal(Int(iComponent),iq);
                }
            }
        }
    }
    // Reducing per-processor values
    M_epetraMapPtr->Comm().SumAll( &fluxScatter, &flux, 1 );

    return flux;
}


// Average value of field on faces with a certain marker
template<typename MeshType>
template<typename VectorType>
Vector PostProc<MeshType>::average( const VectorType& field, const entityFlag_Type& flag,
                                UInt feSpace, UInt nDim )
{
    // Each processor computes the average value on his own flagged faces --> fieldAverageScatter
    // At the end I'll reduce the process values --> fieldAverage
    Vector fieldAverageScatter(nDim), fieldAverage(nDim), localField(nDim);
    // basic policy for type VectorType: operator[] available
    for ( UInt iComponent=0; iComponent < nDim; ++iComponent )
    {
        fieldAverageScatter[iComponent] = 0.;
        fieldAverage[iComponent] = 0.;
        localField[iComponent] = 0.;
    }

    // The total area of the considered faces
    Real areaScatter(0.), area;

    // I need the global Dof ID to query the Oseen solution vector
    // dofVectorIndex is the id of the dof in the data structure of PostProc class
    // dofGlobalId is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt dofVectorIndex, dofGlobalId;

    // list of flagged faces on current processor
    std::list<ID> faceList( M_boundaryMarkerToFacetIdMap[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of field in the current face
    Vector localFieldVector(M_numTotalDofPerFaceVector[feSpace]);

    // Loop on faces
    for (Iterator j=faceList.begin(); j != faceList.end(); ++j)
    {

        // basic policy for type VectorType: operator[] available
        for ( UInt iComponent=0; iComponent < nDim; ++iComponent ) localField[iComponent] = 0.;

        // Updating quadrature data on the current face
        M_currentBdFEPtrVector[feSpace]->updateMeasNormalQuadPt(M_meshPtr->bElement(*j));

        // Loop on components
        for (UInt iComponent =0; iComponent<nDim; ++iComponent)
        {

            // Quadrature formula
            // Loop on quadrature points
            for (Int iq=0; iq< M_currentBdFEPtrVector[feSpace]->nbQuadPt; ++iq)
            {

                // Interpolation
                // Loop on local dof
                for (ID iDof=1; iDof<=M_numTotalDofPerFaceVector[feSpace]; ++iDof)
                {

                    // Extracting nodal values of field in the current face
                    dofVectorIndex = M_vectorNumberingPerFacetVector[feSpace][ ( UInt ) *j - 1 ][ iDof - 1 ];
                    dofGlobalId = M_dofGlobalIdVector[feSpace][dofVectorIndex-1]; // this is in the GLOBAL mesh

                    // basic policy for type VectorType: operator[] available
                    localFieldVector[iDof-1] = field[iComponent*M_numTotalDofVector[feSpace]+dofGlobalId];

                    localField[iComponent] += M_currentBdFEPtrVector[feSpace]->weightMeas(iq)
                                    * localFieldVector[iDof-1] * M_currentBdFEPtrVector[feSpace]->phi(Int(iDof-1),iq);
                }
            }
            // Computing the field integral over the boundary faces
            fieldAverageScatter[iComponent] += localField[iComponent];
        }

        // Computing the area
        areaScatter += M_currentBdFEPtrVector[feSpace]->measure();
    }

    M_epetraMapPtr->Comm().SumAll( &areaScatter, &area, 1 );

    // Reducing per-processor values
    for ( UInt iComponent=0; iComponent < nDim; ++iComponent )
    {
//      fieldAverageScatter[iComponent] /= area;
        M_epetraMapPtr->Comm().SumAll( &fieldAverageScatter[iComponent], &fieldAverage[iComponent], 1 );
    }

    return fieldAverage / area;
}


// Area of patches on the boundary
template<typename MeshType>
void PostProc<MeshType>::set_area()
{
    for ( UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace )
    {

        // area of the mesh face
        Real localArea;

        // index of the considered dof in this class' vectors
        ID dofVectorIndex;

        // ===================================================
        // Loop on boundary faces
        // ===================================================
        for ( ID iboundaryFace = 1 ; iboundaryFace <= M_numBoundaryFaces; ++iboundaryFace )
        {
        	// updating finite element information
            M_currentBdFEPtrVector[iFESpace]->updateMeas( M_meshPtr->bElement( iboundaryFace ) );

            localArea = M_currentBdFEPtrVector[iFESpace]->measure();
            // Loop on the total Dof per Face
            for ( ID iDof = 1; iDof <= M_numTotalDofPerFaceVector[iFESpace]; ++iDof )
            {
                // Extracting local ID of iDof
                dofVectorIndex = M_vectorNumberingPerFacetVector[iFESpace][ iboundaryFace - 1 ][ iDof - 1 ];
                M_patchMeasureVector[iFESpace][dofVectorIndex-1] += localArea;
            }
        }
    }
}

template <typename MeshType>
void PostProc<MeshType>::show_area()
{
    for ( UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace )
    {

        std::cout << "\n***** Post Proc: Area of the patches *****"
                  << "\n\tfor variable " << iFESpace << std::endl;
        ID counter = 1;

        for ( Vector::iterator it = M_patchMeasureVector[iFESpace].begin(); it<M_patchMeasureVector[iFESpace].end(); it++ )
        {
            std::cout << "Boundary Dof: " << counter
                      << ", corresponding to Global Dof: " << M_dofGlobalIdVector[iFESpace][ counter - 1 ]
                      << " has patch area: " << *it << std::endl;
            counter++;
        }
    }
}

template<typename MeshType>
void PostProc<MeshType>::show_bdLtoG()
{
    for ( UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace )
    {

        Int counter = 0;
        std::cout << "\n***** Post Proc: Bd Local To Global *****"
                  << "\n\tfor variable " << iFESpace << std::endl;
        std::cout << M_vectorNumberingPerFacetVector[iFESpace].size() << std::endl;
        for ( std::vector<SimpleVect<ID> >::iterator it1 = M_vectorNumberingPerFacetVector[iFESpace].begin();
                it1<M_vectorNumberingPerFacetVector[iFESpace].end(); it1++ )
        {
            counter++;
            std::cout << "Bd Face " << counter << std::endl;
            for ( SimpleVect<ID>::iterator it2 = it1->begin(); it2<it1->end(); it2++ )
            {
                std::cout << *it2 << ",";
            }
            std::cout << std::endl;
        }

        std::cout << "***** Post Proc: From Boundary Faces to Global Dof *****" << std::endl;
        std::cout << M_dofGlobalIdVector[iFESpace].size() << std::endl;

        for ( std::vector<ID>::iterator it3 = M_dofGlobalIdVector[iFESpace].begin();
                it3<M_dofGlobalIdVector[iFESpace].end(); it3++ )
        {
            std::cout << "Index :" << it3 - M_dofGlobalIdVector[iFESpace].begin()
                      << ", Global Dof: " << *it3 << std::endl;
        }
    }
}

/////////////////////////////////////////////////

///////////////////////////////////////////////


// Normal vectors of patches on the boundary
template<typename MeshType>
void PostProc<MeshType>::set_normal()
{
    for ( UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace )
    {

        // for each patch, the average of each component of the normal vector
        Real sum;

        // index of the considered dof in this class' vectors
        ID dofVectorIndex;

        // ===================================================
        // Loop on boundary faces
        // ===================================================
        for ( ID iboundaryFace = 1 ; iboundaryFace <= M_numBoundaryFaces; ++iboundaryFace )
        {
            // updating finite element information
            M_currentBdFEPtrVector[iFESpace]->updateMeasNormal( M_meshPtr->bElement( iboundaryFace ) );

            // Loop on the components
            for ( Int iComponent = 0; iComponent < NDIM; iComponent++ )
            {
                sum = 0.;
                // Loop on the quadrature points
                for ( Int iQuadraturePoint = 0; iQuadraturePoint < M_currentBdFEPtrVector[iFESpace]->nbQuadPt;
                		++iQuadraturePoint )
                {
                    sum += M_currentBdFEPtrVector[iFESpace]->normal( iComponent, iQuadraturePoint ) *
                    		M_currentBdFEPtrVector[iFESpace]->weightMeas( iQuadraturePoint );
                }
                for ( ID iDof = 1; iDof <= M_numTotalDofPerFaceVector[iFESpace]; ++iDof )
                {
                    // Extracting local ID of iDof
                    dofVectorIndex = M_vectorNumberingPerFacetVector[iFESpace][ iboundaryFace - 1 ][ iDof - 1 ];
                    M_patchNormalVector[iFESpace][ iComponent * M_numBoundaryDofVector[iFESpace] + dofVectorIndex - 1 ] += sum;
                }
            }
        }
        // Normalization of the averaged normals with the patch area
        for ( UInt iBoundaryDof = 0; iBoundaryDof < M_numBoundaryDofVector[iFESpace]; ++iBoundaryDof )
        {
            Real localArea = M_patchMeasureVector[iFESpace][ iBoundaryDof ];
            for ( Int icc = 0; icc < NDIM; icc++ )
                M_patchNormalVector[iFESpace][ icc * M_numBoundaryDofVector[iFESpace] + iBoundaryDof ] *= 1. / localArea;
        }
    }
}


template <typename MeshType>
void PostProc<MeshType>::show_normal()
{
    for ( UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace )
    {

        std::cout << "\n***** Post Proc: Normal vector on the patches *****"
                  << "\n\tfor variable " << iFESpace << std::endl;

        ID counter = 1;

        for ( Vector::iterator it = M_patchMeasureVector[iFESpace].begin(); it<M_patchMeasureVector[iFESpace].end(); it++ )
        {
            std::cout << "Boundary Dof: " << counter
                      << ", corresponding to Global Dof: " << M_dofGlobalIdVector[iFESpace][ counter - 1 ]
                      << " has patch area: " << *it << std::endl;
            std::cout << "and normal components " ;
            for ( Int iComponent = 0; iComponent<NDIM; iComponent++ )
                std::cout <<
                M_patchNormalVector[iFESpace][ iComponent * M_numBoundaryDofVector[iFESpace] + counter - 1 ] << " ";

            std::cout << std::endl;
            counter++;
        }
    }

    std::cout << "End SHOW NORMAL" << std::endl;
}

//////////////////////////////////////////////////
//
///////////////////////////////////////////////////

// Vector with the integral of the shape functions on the patches on the boundary
template<typename MeshType>
void PostProc<MeshType>::set_phi()
{
    for ( UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace )
    {

        // sum contributions from each face of the patch
        Real sum;

        // ID of the considered dof in this class' vectors
        ID dofVectorIndex;

        // ===================================================
        // Loop on boundary faces
        // ===================================================
        for ( ID iboundaryFace = 1 ; iboundaryFace <= M_numBoundaryFaces; ++iboundaryFace )
        {
        	// updating finite element information
            M_currentBdFEPtrVector[iFESpace]->updateMeas( M_meshPtr->bElement( iboundaryFace ) );

            for ( ID iDof = 1; iDof <= M_numTotalDofPerFaceVector[iFESpace]; ++iDof )
            {
                sum = 0.0;
                // global dof
                dofVectorIndex = M_vectorNumberingPerFacetVector[iFESpace][ ( UInt ) iboundaryFace - 1 ][ iDof - 1 ];

                // Loop on the quadrature points
                for ( Int iQuadraturePoint = 0; iQuadraturePoint < M_currentBdFEPtrVector[iFESpace]->nbQuadPt;
                		++iQuadraturePoint )
                {
                    sum += M_currentBdFEPtrVector[iFESpace]->phi( ( Int ) ( iDof - 1 ), iQuadraturePoint )
                    		* M_currentBdFEPtrVector[iFESpace]->weightMeas( iQuadraturePoint );
                }
                M_patchIntegratedPhiVector[iFESpace][ dofVectorIndex - 1 ] += sum;
            }
        }
    }
}


template <typename MeshType>
void PostProc<MeshType>::show_phi()
{
    for ( UInt iFESpace=0; iFESpace<M_numFESpaces; ++iFESpace )
    {

        std::cout << "\n***** Post Proc: Average phi on the patches *****"
                  << "\n\tfor variable " << iFESpace << std::endl;

        ID counter = 0;

        for ( Vector::iterator it = M_patchMeasureVector[iFESpace].begin(); it<M_patchMeasureVector[iFESpace].end(); it++ )
        {
            std::cout << "Boundary Dof: " << counter + 1
                      << ", corresponding to Global Dof: " << M_dofGlobalIdVector[iFESpace][ counter - 1 ]
                      << " has patch area: " << *it << std::endl;
            std::cout << "and average phi  " << M_patchIntegratedPhiVector[iFESpace][ counter ] << std::endl ;
            counter++;
        }
    }
    std::cout << "End SHOW PHI" << std::endl;
}

// the following part is definitely not ready yet TP 09/2008
#if 0
////////////////////////////////
////////////////////////////////
////////////////////////////////
template<typename MeshType>
template<typename VectorType>
VectorType PostProc<MeshType>::compute_sstress( const VectorType& r, UInt ncomp, bool residual = true ) const
{
    ASSERT( ncomp = NDIM, "Error: Shear stress computation possible only for vector unknowns" );

    // prepare the vectors
    VectorType stress( M_numBoundaryDofVector[0] * NDIM );
    stress.clear();
    VectorType nstress( M_numBoundaryDofVector[0] * NDIM );
    nstress.clear();
    VectorType sstress( M_numBoundaryDofVector[0] * NDIM );
    sstress.clear();
    ID counter, dofGlobalId;

    // number of DOFs for each component
    UInt dim = r.size() / ncomp;

    // helper structures to avoid "if" statement
    // if residual==true, vector r has ncomp*M_numBoundaryDofVector components
    // if residual==false, vector r has ncomp*M_numTotalDofVector components
    Vector coef(2);
    std::vector<ID> index(2);

    // loop on locally stored DOFs
    for ( counter = 0; counter < M_numBoundaryDofVector[0]; ++counter )
    {
        // find global ID for boundary dof
        dofGlobalId = M_dofGlobalIdVector[0][counter]; // this is in the GLOBAL mesh

        coef[0] = 1.;
        index[0] = counter;
        coef[1] = 1./ M_patchIntegratedPhiVector[0][counter];
        index[1] = dofGlobalId - 1;

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
VectorType PostProc<MeshType>::compute_stress( const VectorType& grad, const UInt& dim,
                                  const Real& visc ) const
{
    VectorType stress( M_numBoundaryDofVector[0] * NDIM );
    stress.clear();

    ID dofGlobalId;

    // cycle over boundary dof
    for ( ID bound_dof = 0; bound_dof < M_numBoundaryDofVector[0]; ++bound_dof )
    {
        // find global ID for boundary dof
        dofGlobalId = M_dofGlobalIdVector[0][bound_dof]; // this is in the GLOBAL mesh

        // cycle over stress components
        for ( UInt scomp=0; scomp<NDIM; ++scomp )

            // cycle over normal components
            for ( UInt ncomp = 0; ncomp < NDIM; ++ncomp )
            {
                stress[ bound_dof + scomp*M_numBoundaryDofVector[0] ] += visc *
                                                          // grad!
                                                          ( grad[ dofGlobalId-1 + ncomp*dim + scomp*nDimensions*dim ] +
                                                            // transpose grad!
                                                            grad[ dofGlobalId-1 + scomp*dim + ncomp*nDimensions*dim ] )
                                                          * M_patchNormalVector[ bound_dof + ncomp*M_numBoundaryDofVector[0] ];
            }
    }


    return stress;
} // compute_stress
#endif
} // namespace LifeV

#endif /* POSTPROC_H */

