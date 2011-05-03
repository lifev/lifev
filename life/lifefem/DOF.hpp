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
    @brief Degrees of freedom, the class that provides the localtoglobal table

    @author M.A. Fernandez & Luca Formaggia
    @date 00-07-2002

    @contributor Vincent Martin
                 Mohamed Belhadj
                 Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    The numbering of the degrees of freedom follow the order
    Vertices - Edges - Faces - Volumes
    If more than one degree of freedon is associated to a geometrical entity, the
    global numbering associated to two "consecutive" degree of freedom  on the given
    entity is consecutive, and the order is according to the entity orientation.

 */

#ifndef _DOF_HH
#define _DOF_HH

#include <life/lifearray/ArraySimple.hpp>

#include <life/lifecore/LifeV.hpp>

#include <life/lifefem/DOFLocalPattern.hpp>

#include <life/lifemesh/ElementShapes.hpp>

#include <algorithm>
#include <map>

namespace LifeV
{
/*! Local-to-global table

This class provides the localtoglobal table that relates the local DOF of
a finite element to its global numbering. It needs a DOFLocalPattern in
order to obtain all the necessary information about the local pattern. In
fact it stores a copy of it so to make the local pattern available, if
needed.

It is useless until is has not been set up on a specific RegionMesh. This is accomplished either by
passing the mesh to the constructor, or calling the method DOF::update().

\note The methods bulds the table for ALL degrees of freedom, i.e. it does not handle any essential
boundary condition.

Now the class include also a local-to-global table with DOF grouped by (internal) face that was implemented in the old versions into the dofByFace.hpp and dofByFace.cpp files created by D. A. Di Pietro in 2004
*/
class DOF
{
public:

    //! @name Public Types
    //@{

    //@}


    //! @name Constructor & Destructor
    //@{

    /*! The minimal constructor
      \param fePattern is the DOFLocalPattern on which the ref FE is built
      \param Offset: the smallest DOF numbering. It might be used if we want the
      degrees of freedom numbering start from a specific value.
    */
    DOF( const DOFLocalPattern& fePattern);

    //! Copy constructor
    DOF( const DOF & dof2 );

    //! Constructor accepting a mesh as parameter
    /*!
      \param mesh a RegionMesh3D
      \param _fe is the DOFLocalPattern on which the ref FE is built
      \param Offset: the smallest DOF numbering. It might be used if we want the
      degrees of freedom numbering start from a specific value.
    */
    template <typename MeshType>
    DOF( MeshType& mesh, const DOFLocalPattern& fePattern);

    //@}


    //! @name Methods
    //@{

    //! Build the localToGlobal table
    /*!
      \param mesh A RegionMesh3D
      Updates the LocaltoGlobal array
    */
    template <typename MeshType>
    void update( MeshType & );

    /*!
      Returns the global numbering of a DOF, given an boundary facet and the
      local numbering
      \param faceId the boundary facet ID
      \param localDOF the local DOF numbering
      \return The global numbering of the DOF
    */
     ID localToGlobalMapByBdFacet(const ID& facetId, const ID& localDof ) const;

    inline VectorSimple<ID> localToGlobalMapOnBdFacet(const ID& facetId) const{
    	return M_localToGlobalByBdFacet[facetId];
    }

    //! Ouput
    void showMe( std::ostream & out = std::cout, bool verbose = false ) const;
    void showMeByFace(std::ostream& out = std::cout, bool verbose = false) const;

    //@}


    //! @name Set Methods
    //@{

    void setTotalDof(const UInt& totalDof)
    {
        M_totalDof = totalDof;
    }

    //@}


    //! @name Get Methods
    //@{

    //! The total number of Dof
    const UInt& numTotalDof() const
    {
        return M_totalDof;
    }

    //! The number of local DOF (nodes) in the finite element
    const UInt& numLocalDof() const
    {
        return M_elementDofPattern.nbLocalDof();
    }

    //! Return the specified entries of the localToGlobal table
    /*!
      Returns the global numbering of a DOF, given an element and the local numbering
      \param ELId the element ID
      \param localNode the local DOF numbering
      \return The numbering of the DOF
    */
    const ID& localToGlobalMap( const ID ElId, const ID localNode) const
	{
		return M_localToGlobal( localNode, ElId );
	}

    //! Number of elements in mesh
    const UInt& numElements() const
    {
        return M_numElement;
    }

    //! Number of faces in the mesh
    const UInt& numFaces() const
    {
        return M_nbFacets;
    }

    //! Number of local vertices (in a elment)
    const UInt& numLocalVertices() const
    {
        return M_nbLocalPeaks;
    }

    //! Number of local edges (in a elment)
    const UInt& numLocalEdges() const
    {
        return M_nbLocalRidges;
    }

    //! Number of local faces (in a elment)
    const UInt& numLocalFaces() const
    {
        return M_nbLocalFacets;
    }

    //! Number of Local DofByFace
    const UInt& numLocalDofByFace() const
    {
        ASSERT_PRE( (M_numLocalDofByFacet>0) , "This data are not available for this reference element");
        return M_numLocalDofByFacet;
    }

    //! Getter for the localDofPattern
    const DOFLocalPattern& localDofPattern() const
    {
        return M_elementDofPattern;
    }


    //@}



private:

    typedef ArraySimple<UInt> Container_Type;

    typedef ID ( *facetToPointPtr_Type )( ID const& localFace, ID const& point );

    //! The pattern of the local degrees of freedom.
    const DOFLocalPattern& M_elementDofPattern;

    // Total number of degrees of freedom
    UInt M_totalDof;

    // Number of elements in the considered mesh
    UInt M_numElement;


    UInt M_nbLocalPeaks;
    UInt M_nbLocalRidges;
    UInt M_nbLocalFacets;

    // The local to global table
    Container_Type M_localToGlobal;

    // number of faces in the mesh
    UInt M_nbFacets;

    // The local to global table based on the boundary facets
    std::vector<VectorSimple<ID> > M_localToGlobalByBdFacet;

    // local array that maps the local dof of the
    facetToPointPtr_Type M_facetToPoint;

    // face to the local dof the the element
    UInt M_numLocalDofByFacet;

    // Just 5 counters
    UInt M_dofPositionByEntity[ 5 ];
};


// ===================================================
// Constructors & Destructor
// ===================================================

//! Constructor that builds the localToglobal table
template <typename MeshType>
DOF::DOF( MeshType& mesh, const DOFLocalPattern& fePattern):
        M_elementDofPattern       ( fePattern ),
        M_totalDof( 0 ),
        M_numElement     ( 0 ),
        M_nbLocalPeaks      ( 0 ),
        M_nbLocalRidges      ( 0 ),
        M_nbLocalFacets      ( 0 ),
        M_localToGlobal     (),
        M_nbFacets( 0 ),
        M_localToGlobalByFacet(),
        M_globalToLocalByFacet()
{
    //Getting the face
    switch ( fePattern.nbLocalDof() )
    {
    case 1: //P0 Q0
    	M_numLocalDofByFace = 0;
    	break;
    case 2:
        // No M_faceToPoint (it is 1D)
        M_numLocalDofByFace = 1;
        break;
    case 4:
        M_faceToPoint = LinearTetra::faceToPoint;
        M_numLocalDofByFace = 3;
        break;
    case 5:
        M_faceToPoint = LinearTetraBubble::faceToPoint;
        M_numLocalDofByFace = 3;
        break;
    case 10:
        M_faceToPoint = QuadraticTetra::faceToPoint;
        M_numLocalDofByFace = 6;
        break;
    case 8:
        M_faceToPoint = LinearHexa::faceToPoint;
        M_numLocalDofByFace = 4;
        break;
    case 27:
        M_faceToPoint = QuadraticHexa::faceToPoint;
        M_numLocalDofByFace = 27;
        break;
    default:
        std::cout << "Warning: This refFE is not available for the dof by face." << std::endl;
        M_numLocalDofByFace = 0;
        break;
    }

    for ( UInt i = 0; i < 5; ++i )
        M_dofPositionByEntity[ i ] = 0;
    update( mesh );
}

// ===================================================
// Methods
// ===================================================

//! Build the localToGlobal table
template <typename MeshType>
void DOF::update( MeshType& mesh )
{

    typedef typename MeshType::elementShape_Type geoShape_Type;
    typedef typename geoShape_Type::GeoBShape geoBShape_Type;

    // Some useful local variables, to save some typing
    UInt nbLocalDofPerPeak = M_elementDofPattern.nbDofPerPeak();
    UInt nbLocalDofPerRidge = M_elementDofPattern.nbDofPerRidge();
    UInt nbLocalDofPerFacet = M_elementDofPattern.nbDofPerFacet();
    UInt nbLocalDofPerElement = M_elementDofPattern.nbDofPerElement();

    M_nbLocalPeaks = geoShape_Type::S_numPeaks;
    M_nbLocalRidges = geoShape_Type::S_numRidges;
    M_nbLocalFacets = geoShape_Type::S_numFacets;

    M_numElement = mesh.numElements();

    UInt nbGlobalElements = mesh.numGlobalElements();
    UInt nbGlobalRidges = mesh.numGlobalRidges();
    UInt nbGlobalPeaks = mesh.numGlobalPeaks();
    UInt nbGlobalFacets = mesh.numGlobalFacets();

    UInt nbDofElemPeaks = nbLocalDofPerPeak * M_nbLocalPeaks; // number of vertex's DOF on a Element
    UInt nbDofElemRidges = nbLocalDofPerRidge * M_nbLocalRidges; // number of edge's DOF on a Element


    UInt i, l, ie;

    // Total number of degree of freedom for each element

    UInt nldof = nbLocalDofPerElement
        		+ nbLocalDofPerRidge * M_nbLocalRidges
    			+ nbLocalDofPerPeak * M_nbLocalPeaks
    			+ nbLocalDofPerFacet * M_nbLocalFacets;

    // Consistency check
    ASSERT_PRE( nldof == UInt( M_elementDofPattern.nbLocalDof() ), "Something wrong in FE specification" ) ;

    // Global total of degrees of freedom
	M_totalDof = nbGlobalElements * nbLocalDofPerElement
		+ nbGlobalRidges * nbLocalDofPerRidge
		+ nbGlobalPeaks * nbLocalDofPerPeak
		+ nbGlobalFacets * nbLocalDofPerFacet;


    // Reshape the container to fit the needs
    M_localToGlobal.reshape( nldof, M_numElement );

    // Make sure the mesh has everything needed


    UInt gcount( 0 );
    UInt lcount;
    UInt lc;

    bool update_ridges( nbLocalDofPerRidge != 0 && ! mesh.hasLocalRidges() );
    if ( update_ridges )
        mesh.updateElementRidges();

    bool update_facets( nbLocalDofPerFacet != 0 && ! mesh.hasLocalFacets() );
    if ( update_facets )
        mesh.updateElementFacets();

    // Peak Based DOFs
    if ( nbLocalDofPerPeak > 0 )
        for ( ie = 0; ie < M_numElement; ++ie )//for each element
        {
            lc = 0;
            for ( i = 0; i < M_nbLocalPeaks; ++i )//for each vertex in the element
                for ( l = 0; l < nbLocalDofPerPeak; ++l )//for each degree of freedom per vertex
                {
                    // label of the ith point of the mesh element
                    M_localToGlobal( lc++, ie ) = gcount + mesh.element( ie ).point( i ).id() * nbLocalDofPerPeak + l;
               }
        }
    gcount += nbLocalDofPerPeak * nbGlobalPeaks;//dof per vertex * total # vertices
    lcount = nbLocalDofPerPeak * M_nbLocalPeaks;

    // Ridge Based DOFs
	if ( nbLocalDofPerRidge > 0 )
		for ( ie = 0; ie < M_numElement; ++ie )
		{
			lc = lcount;
			for ( i = 0; i < M_nbLocalRidges; ++i )
			{
				UInt eID = mesh.ridge(mesh.localRidgeId(ie, i)).id();
				for ( l = 0; l < nbLocalDofPerRidge; ++l )
					M_localToGlobal( lc++, ie ) = gcount + eID * nbLocalDofPerRidge + l;
			}
		}
	gcount += nbGlobalRidges * nbLocalDofPerRidge;
	lcount += nbLocalDofPerRidge * M_nbLocalRidges;

    //Facet based DOFs
    if ( nbLocalDofPerFacet > 0 )
        for ( ie = 0; ie < M_numElement; ++ie )
        {
            lc = lcount;

            for ( i = 0; i < M_nbLocalFacets; ++i )
            {
            	UInt fID = mesh.facet( mesh.localFacetId( ie, i ) ).id();
                for ( l = 0; l < nbLocalDofPerFacet; ++l )
                    M_localToGlobal( lc++, ie ) = gcount + fID * nbLocalDofPerFacet + l;
            }
        }
    gcount += nbGlobalFacets * nbLocalDofPerFacet;
    lcount += nbLocalDofPerFacet * M_nbLocalFacets;
    // Element  Based DOFs
    if ( nbLocalDofPerElement > 0 )
        for ( ie = 0; ie < M_numElement; ++ie )
        {
            lc = lcount;
            for ( l = 0; l < nbLocalDofPerElement; ++l )
                M_localToGlobal( lc++, ie ) = gcount + mesh.element( ie ).id() * nbLocalDofPerElement + l;
       }
    gcount += nbGlobalElements * nbLocalDofPerElement;

    UInt nBElemRidges = geoBShape_Type::S_numRidges; // Number of boundary facet's vertices

    UInt nBElemFacets = geoBShape_Type::S_numFacets;    // Number of boundary facet's edges

    ASSERT_POS( gcount == M_totalDof , "Something wrong in Dof Setup " << gcount  << " " << M_totalDof ) ;


    VectorSimple<ID> globalDOFOnBdFacet(nbLocalDofPerPeak*nBElemRidges + nBElemFacets*nbLocalDofPerRidge + nbLocalDofPerFacet);
    M_localToGlobalByBdFacet.resize(mesh.numBFacets());

    for ( ID iBoundaryFacet = 0 ; iBoundaryFacet < mesh.numBFacets(); ++iBoundaryFacet )
	{
		ID iAdjacentElem = mesh.bFacet( iBoundaryFacet ).firstAdjacentElementIdentity();  // id of the element adjacent to the face
		ID iElemBFacet = mesh.bFacet( iBoundaryFacet ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

		UInt lDof(0), dofOffset; //local DOF on boundary element

		//loop on Dofs associated with peaks
		if(nbLocalDofPerPeak)
			for ( ID iBElemRidge = 0; iBElemRidge < nBElemRidges; ++iBElemRidge )
			{
				ID iElemPeak = geoShape_Type::facetToPeak( iElemBFacet, iBElemRidge ); // local vertex number (in element)
				dofOffset = iElemPeak * nbLocalDofPerPeak;
				for ( ID l = 0; l < nbLocalDofPerPeak; ++l )
					globalDOFOnBdFacet[ lDof++ ] = M_localToGlobal( dofOffset + l, iAdjacentElem );
			}

		//loop on Dofs associated with Ridges
		if(nbLocalDofPerRidge)
			for ( ID iBElemFacets = 0; iBElemFacets < nBElemFacets; ++iBElemFacets )
			{
				ID iElemRidge = geoShape_Type::facetToRidge( iElemBFacet, iBElemFacets ); // local edge number (in element)
				dofOffset = nbDofElemPeaks + iElemRidge * nbLocalDofPerRidge;
				for ( ID l = 0; l < nbLocalDofPerRidge; ++l )
					globalDOFOnBdFacet[ lDof++ ] = M_localToGlobal( dofOffset + l, iAdjacentElem ); // global Dof
			}

		//loop on Dofs associated with facets
		dofOffset = nbDofElemPeaks + nbDofElemRidges + iElemBFacet * nbLocalDofPerFacet;
		for ( ID l = 0; l < nbLocalDofPerFacet; ++l )
			globalDOFOnBdFacet[ lDof++ ] = M_localToGlobal( dofOffset + l, iAdjacentElem ); // global Dof


		M_localToGlobalByBdFacet[iBoundaryFacet] = globalDOFOnBdFacet;
	}


    if ( update_ridges )
             mesh.cleanElementRidges();
	 if ( update_facets )
		 mesh.cleanElementFacets();
}

}
#endif
