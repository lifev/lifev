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

#include <algorithm>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/ArraySimple.hpp>

#include <lifev/core/fem/DOFLocalPattern.hpp>

#include <lifev/core/mesh/ElementShapes.hpp>

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

Now the class include also a local-to-global table with DOF grouped by (internal) face
that was implemented in the old versions into the dofByFace.hpp and dofByFace.cpp files created by D. A. Di Pietro in 2004

@note The dof numbering refers to the global numbering, not the numbering local to a partition of a partitioned mesh
@todo This class must be bettered. The logic by which the dof table is built may fail for certain type of elements.
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
      \param mesh a RegionMesh
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
      \param mesh A RegionMesh
      Updates the LocaltoGlobal array
    */
    template <typename MeshType>
    void update( MeshType & );

    //! Build the globalElements list
    /*!
      @param mesh A RegionMesh
    */
    template <typename MeshType>
    std::vector<Int> globalElements( MeshType & mesh );

    /*!
      Returns the global numbering of a DOF, given an boundary facet and the
      local numbering
      \param faceId the boundary facet ID
      \param localDOF the local DOF numbering
      \return The global numbering of the DOF
    */
     ID localToGlobalMapByBdFacet(const ID& facetId, const ID& localDof ) const;

    inline std::vector<ID> localToGlobalMapOnBdFacet(const ID& facetId) const{
    	return M_localToGlobalByBdFacet[facetId];
    }

    //! Ouput
    void showMe( std::ostream & out = std::cout, bool verbose = false ) const;
    void showMeByBdFacet(std::ostream& out = std::cout, bool verbose = false) const;

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

    // The local to global table based on the boundary facets
    std::vector<std::vector<ID> > M_localToGlobalByBdFacet;

    // local array that maps the local dof of the
    facetToPointPtr_Type M_facetToPoint;

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
        M_localToGlobal     ()
{
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
    bool update_ridges( nbLocalDofPerRidge != 0 && ! mesh.hasLocalRidges() && (MeshType::S_geoDimensions == 3));
    if ( update_ridges )
        mesh.updateElementRidges();

    bool update_facets( nbLocalDofPerFacet != 0 && ! mesh.hasLocalFacets() );
    if ( update_facets )
        mesh.updateElementFacets();

    UInt gcount( 0 );
    UInt lcount;
    UInt lc;

    // Peak Based DOFs
    M_dofPositionByEntity[ 0 ] = gcount;
    if ( nbLocalDofPerPeak > 0 )
        for ( ie = 0; ie < M_numElement; ++ie )//for each element
        {
            lc = 0;
            for ( i = 0; i < M_nbLocalPeaks; ++i )//for each vertex in the element
            {
                ID pID = mesh.element( ie ).point( i ).id();
                for ( l = 0; l < nbLocalDofPerPeak; ++l )//for each degree of freedom per vertex
                {
                    // label of the ith point of the mesh element
                    UInt dof = gcount +  pID * nbLocalDofPerPeak + l;
                    ASSERT( dof != NotAnId, "the dof is not properly set" );
                    M_localToGlobal( lc++, ie ) = dof;
                }
            }
        }

    // Ridge Based DOFs
    gcount += nbLocalDofPerPeak * nbGlobalPeaks;//dof per vertex * total # vertices
    lcount = nbLocalDofPerPeak * M_nbLocalPeaks;
    M_dofPositionByEntity[ 1 ] = gcount;
    
    if ( nbLocalDofPerRidge > 0 )
        for ( ie = 0; ie < M_numElement; ++ie )
        {
            lc = lcount;
            for ( i = 0; i < M_nbLocalRidges; ++i )
            {
                UInt rID = mesh.ridge(mesh.localRidgeId(ie, i)).id();
                for ( l = 0; l < nbLocalDofPerRidge; ++l )
                {
                    UInt dof = gcount +  rID * nbLocalDofPerRidge + l;
                    ASSERT( dof != NotAnId, "the dof is not properly set" );
                    M_localToGlobal( lc++, ie ) = dof;
                }
            }
        }

    //Facet based DOFs
    gcount += nbGlobalRidges * nbLocalDofPerRidge;
    lcount += nbLocalDofPerRidge * M_nbLocalRidges;
    M_dofPositionByEntity[ 2 ] = gcount;

    if ( nbLocalDofPerFacet > 0 )
        for ( ie = 0; ie < M_numElement; ++ie )
        {
            lc = lcount;

            for ( i = 0; i < M_nbLocalFacets; ++i )
            {
                UInt fID = mesh.facet( mesh.localFacetId( ie, i ) ).id();
                for ( l = 0; l < nbLocalDofPerFacet; ++l )
                {
                    UInt dof = gcount +  fID * nbLocalDofPerFacet + l;
                    ASSERT( dof != NotAnId, "the dof is not properly set" );
                    M_localToGlobal( lc++, ie ) = dof;
                }
            }
        }

    // Element  Based DOFs
    gcount += nbGlobalFacets * nbLocalDofPerFacet;
    lcount += nbLocalDofPerFacet * M_nbLocalFacets;

    M_dofPositionByEntity[ 3 ] = gcount;
    if ( nbLocalDofPerElement > 0 )
        for ( ie = 0; ie < M_numElement; ++ie )
        {
            lc = lcount;
            ID eID = mesh.element( ie ).id();
            for ( l = 0; l < nbLocalDofPerElement; ++l )
            {
                UInt dof = gcount +  eID * nbLocalDofPerElement + l;
                ASSERT( dof != NotAnId, "the dof is not properly set" );
                M_localToGlobal( lc++, ie ) = dof;
            }
        }
    gcount += nbGlobalElements * nbLocalDofPerElement;
    M_dofPositionByEntity[ 4 ] = gcount;

    ASSERT_POS( gcount == M_totalDof , "Something wrong in Dof Setup " << gcount  << " " << M_totalDof ) ;

	//Building map of global DOF on boundary facets
    UInt nBElemRidges = geoBShape_Type::S_numRidges; // Number of boundary facet's vertices
    UInt nBElemFacets = geoBShape_Type::S_numFacets;    // Number of boundary facet's edges
    
    std::vector<ID> globalDOFOnBdFacet(nbLocalDofPerPeak*nBElemRidges + nBElemFacets*nbLocalDofPerRidge + nbLocalDofPerFacet);
    M_localToGlobalByBdFacet.resize(mesh.numBoundaryFacets());

    for ( ID iBoundaryFacet = 0 ; iBoundaryFacet < mesh.numBoundaryFacets(); ++iBoundaryFacet )
	{
		ID iAdjacentElem = mesh.boundaryFacet( iBoundaryFacet ).firstAdjacentElementIdentity();  // id of the element adjacent to the face
		ID iElemBFacet = mesh.boundaryFacet( iBoundaryFacet ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

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

template <typename MeshType>
std::vector<Int> DOF::globalElements( MeshType& mesh )
{
    std::set<Int> dofNumberSet;
    // Gather all dofs local to the given mesh (dofs use global numbering)
    // The set ensures no repetition
    for (UInt elementId=0; elementId < mesh.numElements(); ++elementId )
        for (UInt localDof=0; localDof < this->numLocalDof();++localDof )
            dofNumberSet.insert( static_cast<Int>( this->localToGlobalMap(elementId,localDof ) ) );
    // dump the set into a vector for adjacency
    // to save memory I use copy() and not the vector constructor directly
    std::vector<Int> myGlobalElements(dofNumberSet.size());
    std::copy(dofNumberSet.begin(),dofNumberSet.end(),myGlobalElements.begin());
    // Save memory
    dofNumberSet.clear();

    return myGlobalElements;
}

}
#endif
