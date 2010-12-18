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
    @brief   Special routines to read meshes and special structures for
  sides and faces handling

  Classes BareFace and BareEdge have been created to give an UNIQUE
  Representation for mesh faces and edges and thus allow the construction
  of global tables or fields.

    @author Luca Formaggia <luca.formaggia@polimi.it>
    @contributor Luca Bertagna <lbertag@emory.edu>
    @contributor Lucia Mirabella <lucia.mirabell@gmail.com>
    @maintainer Lucia Mirabella <lucia.mirabell@gmail.com>

    @date 19-08-1999

    Classes BareFace and BareEdge have been created to give an UNIQUE internal
	representation for mesh faces and edges, allowing thus the construction of
	Dof objects (which are naturally linked to mesh entities).

	\par Introduction

	One of the paradigms chosen for the development of this	library is the fact
	that degrees of freedom (Dof) are linked to geometrical	entities.  Now if
	we have degrees of freedom associated, for instance, to	an Edge (like in
	a P2 Tetra) in order to build the global numbering of the Dof and the
	association between local (element-wise) and global numbering, we need to
	identify edges and give them a unique ID. Yet, we may not want to
	build a full Edge (GeoElement2D) object: we only need
	the ID of the edge and a way of computing the ID's of the degrees of
	freedom on the edge, all the remaining data of a full Edge object is not
	necessarily needed.

	Another related problem is how to uniquely identify a face or an edge in the mesh.

	The dilemma has been resolved by creating the concept of a BareEdge and
	BareFace (bare geometry items).  A bare geometry item is formed by the
	minimal information required to uniquely identify it, namely 2
	<tt>Point</tt>'s ID 's for an edge and 3 <tt>Point</tt>'s ID 's for the
	Faces (it is enough also for Quad faces!). We build the bare items by
	looping through the elements and obviously we make sure that the BareItem
	ID is consistent with that of the corresponding ``full item'' if the latter
	has been instantiated.

	Another <em>very important</em> issue is that of orientation.  There are
	different ways of considering orientation of a Face or an Edge. The first is
	the <em>local</em> orientation of a Face or Edge of the reference finite
	element. This is conventionally chosen when designing the finite
	elements. For the faces, we have adopted the convention that the local
	positive orientation is such that the face normal calculated with the right
	hand rule is <em>outwardly</em> oriented. As for the edges, the local
	orientation for a 3D element is more arbitrary.

	However, for a 2D element, the positive orientation of an Edge is the one
	which is in accordance with the right hand rule applied to that element.

	When a Face or an edge is <em>active</em>, i.e. is effectively stored in
	the mesh, then there is another obvious orientation, of global rather than
	local nature: that induced by the way the edge or face is stored. For
	boundary elements (faced in 3D or edges in 2D) it is compulsory that the
	orientation of the stored item be consistent with the convention chosen for
	the orientation of the domain boundary. More precisely, boundary elements
	are stored so that the normal (always calculated following the right hand
	rule) is outward with respect to the domain.

	However, there is the need of defining a <em>global</em> orientation also
	for <em>non active</em> entities. This because we do not always store all
	faces and all edges. We need then to choose a unique way to identify the
	orientation of an Edge or of a Face <em>independently </em> from the fact
	that they are active or not. We will call this orientation the <em>natural </em>
	orientation. We have chosen the following convention for natural orientation
	of faces and edges

	<ul>
	<li>The positive natural orientation of an  <em>Edge</em> is given by \f$V_{min} \rightarrow V_{max} \f$,
        \f$V_{min}\f$ being the Edge Vertex with smallest ID</li>

	<li>The positive natural orientation of a  <em>Face</em> is given by the cicle
        \f$V_{min} \rightarrow V_2\rightarrow V_3 \f$,
        \f$V_{min}\f$ being the Face Vertex with smallest ID, \f$V_2\f$ the second smallest
        and \f$V_2\f$ the thirsd smallest.</li>
	</ul>

	Note that the latter definition applies both to triangular and to quad faces.

	\warning If I want to associate boundary conditions I need the active
	entity, since BareEdges do not store Marker data. (This is why the
	RegionMesh classes treat boundary items in a rather special way).


 */

#ifndef BAREITEMS_H
#define BAREITEMS_H 1

#include<utility>
#include<vector>
#include<map>
#include<algorithm>
#include<iostream>

#include <life/lifecore/life.hpp>

namespace LifeV
{
#ifndef _LIFEV_HH_
/*! \typedef typedef unsigned int ID
\brief type used for Identifiers (Integral type in the range [1, MAX_INT]).
All principal items handled by the library have an identified, which is an
unsigned integer <em>greater or equal to one</em>
 */
// more correct version
typedef size_t ID;
// original version
//typedef unsigned int ID;
//! I define UInt=unsigned int. This allow to use this stuff also outside lifeV
// more correct version
//typedef size_t UInt;
// original version
typedef unsigned int UInt;
#endif

//! The Edge basis class
/*! It contains the attributes common to all Edges. In particular, it
contains the two ID's (first and second) of the points at the two ends
of the edge.
\invariant first < second
 */
struct BareEdge
{
    //! @name Constructor & Destructor
    //@{
    //! Empty Constructor
    BareEdge() : first( 0 ), second( 0 )
    {}
    ;
    //! Constructor that takes the ID's as parameter
    /*!
        @param i ID of the first node of the edge
        @param j ID of the second node of the edge
     */
    BareEdge( ID i, ID j ) : first( i ), second( j )
    {}
    ;
    //@}
    ID first; //!< First ID which defines the Edge
    ID second; //!< Second ID which defines the Edge
};

//! The base Face class
/*! It contains the attributes common to all Faces. In particular, a it contains three
  ID's (first, second and third) of points at face vertex, chosen so
  to uniquely identify the face.
\invariant first<second<third.
 */
struct BareFace
{
    //! @name Constructor & Destructor
    //@{
    //! Empty Constructor
    BareFace() : first( 0 ), second( 0 ), third( 0 )
    {}
    ;
    //! Constructor that takes the ID's as parameter
    /*!
        @param i ID of the first node of the face
        @param j ID of the second node of the face
        @param k ID of the third node of the face
     */
    //!
    BareFace( ID i, ID j, ID k ) : first( i ), second( j ), third( k )
    {}
    ;
    //! Constructor that takes a BareEdge object and an ID. The face is then identified by the ID of the Points on the BareEdge + the point with identified by id
    /*!
        @param id of the third point defining the face
        @param id of the edgedefining the face

     */
    BareFace( ID id, const BareEdge & edge ) : first( id ), second( edge.first ), third( edge.second )
    {}
    ;
    ID first;  //!< First ID which defines the BareFace
    ID second; //!< Second ID which defines the BareFace
    ID third;  //!< Third ID which defines the BareFace
};


//! \defgroup BareItemsBuilder Global functions to build Bare Items.

/*! \ingroup BareItemsBuilder
  \brief It creates a BareEdge end returns the orientation of the created edge with respect
  to the given data.
  @param i is a Point ID
  @param j is a Point ID
  @return a pair composed of the created BareEdge and a boolean value which is False if orientation has been changed, True otherwise

  The BareEdge that will be built is the one passing by <tt>i</tt> and
  <tt>j</tt>. The orientation is i->j if the returned parameter is a
  true.

  \pre i and j >0, i!=j */
inline
std::pair<BareEdge, bool>
makeBareEdge( ID const i, ID const j )
{
    if ( i < j )
    {
        return std::make_pair( BareEdge( i, j ), true );
    }
    else
    {
        return std::make_pair( BareEdge( j, i ), false );
        ;
    }
}

/*! \ingroup BareItemsBuilder
  \brief It creates a BareEdge, ignoring orientation.
  @param i is a Point ID
  @param j is a Point ID
  @return the BareEdge created

  The BareEdge that will be built is the one passing by <tt>i</tt> and <tt>j</tt>
  A lighter version of MakeBareEdge, to be used
  if orientation flag is not needed;
  \pre i and j >0, i!=j
 */
inline
BareEdge
setBareEdge( ID const i, ID const j )
{
    BareEdge bareEdge;
    bareEdge.first = i < j ? i : j;
    bareEdge.second = ( i - bareEdge.first ) + j;
    return bareEdge;
}

/*! \ingroup BareItemsBuilder
  \brief It creates a non-standard BareEdge.
  @param i is a Point ID
  @param j is a Point ID
  @return the BareEdge created

  \pre i and j >0

  Yet another  lighter version of MakeBareEdge, without orientation, To be used for
  non-oriented graphs.

  \warning It produces a BareEdge which does not comply with the invariant of the
  class (first < second). It must be used only if the BareEdge class is NOT used to uniquely identify edges.
 */
inline
BareEdge
setBareEdgeNo( ID const i, ID const j )
{
    return BareEdge( i, j );
}


/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from three Point ID's
  \param i is a Point ID
  \param j is a Point ID
  \param k is a Point ID
  @return a pair composed of the created BareFace and a boolean value which is False if orientation has been changed, True otherwise

  To be used for triangular faces.
  \pre i, j and k >0. i!=j!=k

 */
std::pair<BareFace, bool> makeBareFace( ID const i, ID const j, ID const k );

/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from four Point ID's. To be used with Quad faces.
  \param i is a Point ID
  \param j is a Point ID
  \param k is a Point ID
  \param l is a Point ID
  @return a pair composed of the created BareFace and a boolean value which is False if orientation has been changed, True otherwise

  \remarks For quad faces the construction process is more complex. We start from
  the smallest vertex and we take the first three vertices in the
  sequence. We then proceede as for the triangles.

 */
std::pair<BareFace, bool> makeBareFace( ID const i, ID const j, ID const k, ID const l );

/*! \defgroup comparison Comparison Operators
  Operators for comparing BareItems
 */

/*! \ingroup comparison
    \brief inequality
 */
inline
bool
operator!=( const BareEdge & edge1 , const BareEdge & edge2 )
{
    return edge1.first != edge2.first || edge1.second != edge2.second;
}

/*! \ingroup comparison
     \brief equality
 */
inline
bool
operator==( const BareEdge & edge1 , const BareEdge & edge2 )
{
    return edge1.first == edge2.first && edge1.second == edge2.second;
}

/*! \ingroup comparison
    \brief greater than
 */
inline
bool
operator>( const BareEdge & edge1 , const BareEdge & edge2 )
{
    return edge2.first > edge1.first || ( edge2.first == edge1.first && edge2.second > edge1.second );
}

/*! \ingroup comparison
    \brief greater-equal than
 */
inline
bool
operator>=( const BareEdge & edge1 , const BareEdge & edge2 )
{
    return edge1 == edge2 || edge1 > edge2;
}

/*! \ingroup comparison
    \brief less than
 */
inline
bool
operator<( const BareEdge & edge1 , const BareEdge & edge2 )
{
    return edge2.first < edge1.first || ( edge2.first == edge1.first && edge2.second < edge1.second );
}

/*! \ingroup comparison
    \brief less-equal than
 */
inline
bool
operator<=( const BareEdge & edge1 , const BareEdge & edge2 )
{
    return edge1 == edge2 || edge1 < edge2;
    ;
}

/*! \ingroup comparison
  \brief inequality
 */
inline
bool
operator!=( const BareFace & face1 , const BareFace & face2 )
{
    return face1.first != face2.first || face1.second != face2.second || face1.third != face2.third;
}


/*! \ingroup comparison
  \brief equality
 */
inline
bool
operator==( const BareFace & face1 , const BareFace & face2 )
{
    return face1.first == face2.first && face1.second == face2.second && face1.third == face2.third;
}

/*! \ingroup comparison
  \brief General functor for lexicographic comparison
 */
template <typename T>
struct cmpBareItem;

/*! \ingroup comparison
   \brief Specialized functor for Edges
 */
template <>
struct cmpBareItem<BareEdge> //!< The actual comparison operator
{
    bool operator() ( const BareEdge & edge1, const BareEdge & edge2 ) const
    {
        return edge2.first > edge1.first || ( edge2.first == edge1.first && edge2.second > edge1.second );
    }
};

/*! \ingroup comparison
  \brief Specialized functor for Faces
 */
template <>
struct cmpBareItem<BareFace>
{
    bool operator() ( const BareFace & edge1, const BareFace & edge2 ) const
    {
        if ( edge2.first > edge1.first )
            return true;
        if ( edge2.first == edge1.first )
        {
            if ( edge2.second > edge1.second )
                return true;
            if ( edge2.second == edge1.second )
                return edge2.third > edge1.third;
        }
        return false;
    }
};

//! BareItemsHandler class - Class to handle bare edges and faces construction
/*!
    @author Luca Formaggia
    @see

	This class handles mesh bare edges and faces construction. Used only in mesh builders
	A BareItemsHandler is a specialisation of a STL map which holds the pair
	formed by a bareitem and its ID.
	The ID  is automatically generated if one uses the method addIfNotThere
 */
template <typename BareItemType>
class BareItemsHandler: public std::map<BareItemType, ID, cmpBareItem<BareItemType> >
{
public:
    //! @name Public Types
    //@{
    typedef BareItemType 												bareItem_Type;
    typedef std::map<bareItem_Type, UInt, cmpBareItem<bareItem_Type> > 	container_Type;
    typedef typename container_Type::iterator 							containerIterator_Type;
    typedef typename container_Type::const_iterator 					containerConstIterator_Type;
    typedef std::pair<const bareItem_Type, UInt> 						value_Type;
    //@}

    //! @name Constructors & Destructor
    //@{
    //! Empty Constructor
    BareItemsHandler();
    //@}

    //! @name Methods
    //@{

    //! Method to ask if an item already exists
    /*!
        @param item Item we are looking for
        @return True if the item has been found, False otherwise
     */
    bool isThere( bareItem_Type const & item) const;

    //! Method that adds a BareItem if it is not already there and automatically generates the ID
    /*!
        @param item Item to be added
    	@return a pair composed of the ID of the added item and a boolean value which is True if the item has been successfully added and False otherwise
     */
    std::pair<ID, bool> addIfNotThere( bareItem_Type const & item );

    //! Method that adds a bareItem_Type if it is not already there and assigns the ID
    /*!
        @param item Item to be added
        @param id ID to be assigned to the item
    	@return a pair composed of the ID of the added item and a boolean value which is True if the item has been successfully added and False otherwise
     */
    std::pair<ID, bool> addIfNotThere( bareItem_Type const & item, const ID id );

    //! Method that removes a bareItem_Type if it is there (the ID is then lost)
    /*!
        @param item Item to be removed
    	@return True if the item has been erased and False otherwise
     */
    bool deleteIfThere( bareItem_Type const & item);

    //! Method that removes a bareItem_Type if it is there (the ID is then lost)
    /*!
        @deprecated
        @param item Item to be removed
    	@return True if the item has been erased and False otherwise
     */
    bool isThereDel( bareItem_Type const & item);

    //! Method that counts how many items are stored
    /*!
    	@return the number of entities actually stored
     */
    UInt howMany() const;

    //! Method that returns the maximum id currently in use
    /*!
    	@return the maximum id currently in use
     */
    UInt maxId() const
    {
        return M_idCount;
    }

    //! Method that writes info in output
    /*!
     */
    void showMe() const;

    //! Method that returns the ID of a BareItem. It returns 0 if the item doesn't exist
    /*!
        @param item Item we are looking for
        @return ID of the item. 0 if the item doesn't exist
     */
    ID id( bareItem_Type const & item ) const;

    //@}

    //! @name Set Methods
    //@{
    //! Method that returns the ID of a BareItem. It returns 0 if the item doesn't exist
    /*!
        @param item Item to modify
        @param id new ID to assign to item
    	@return True if the item has been found and modified, False otherwise
     */
    bool setId( bareItem_Type const & item, const ID& id );

    //! Method that returns the ID of a BareItem. It returns 0 if the item doesn't exist
    /*!
     	@deprecated
        @param item Item to modify
        @param id new ID to assign to item
    	@return True if the item has been found and modified, False otherwise
     */
    bool setId( bareItem_Type const & item, ID const id );
    //@}

    //! @name Get Methods
    //@{


    //@}
private:
    UInt M_idCount;
};

/*********************************************************************************
               IMPLEMENTATIONS
 *********************************************************************************/

// ===================================================
// Constructors & Destructor
// ===================================================
template <class BareItemType>
BareItemsHandler<BareItemType>::BareItemsHandler() :
        M_idCount( 0 )
{ }

// ===================================================
// Methods
// ===================================================
template <class BareItemType>
inline
bool
BareItemsHandler<BareItemType>::isThere( const bareItem_Type & item ) const
{
    return find( item ) != container_Type::end();
}

template <class BareItemType>
inline
std::pair<ID, bool>
BareItemsHandler<BareItemType>::addIfNotThere( const bareItem_Type & item )
{
    std::pair<typename BareItemsHandler<BareItemType>::containerIterator_Type, bool> i( insert( std::make_pair( item, M_idCount + 1 ) ) );
    if ( i.second )
        ++M_idCount;
    return std::make_pair( ( i.first )->second, i.second );
}

template <class BareItemType>
inline
std::pair<ID, bool>
BareItemsHandler<BareItemType>::addIfNotThere( const bareItem_Type & item, const ID id )
{
    std::pair<typename BareItemsHandler<BareItemType>::containerIterator_Type, bool> i( insert( std::make_pair( item, id ) ) );
    ( i.first ) ->second = id; // Set new id in any case.
    return std::make_pair( id, i.second ); // for consistency with other version.
}
template <class BareItemType>
bool
BareItemsHandler<BareItemType>::deleteIfThere( bareItem_Type const & item )
{
    return erase( item ) != 0;
}

template <class BareItemType>
bool  __attribute__ ((__deprecated__))
BareItemsHandler<BareItemType>::isThereDel( bareItem_Type const & item )
{
    return erase( item ) != 0;
}

template <class BareItemType>
inline
UInt
BareItemsHandler<BareItemType>::howMany() const
{
    return container_Type::size();
}

template <typename BareItemType>
inline
void BareItemsHandler<BareItemType>::showMe() const
{
    std::cout << "BareItemsHandler: " << std::endl;
    std::cout << "Number of Items stored: " << this->size() << std::endl;
    std::cout << "Max Id stored         : " << this->maxId() << std::endl;
    std::cout << "End of Information";
}

template <class BareItemType>
inline
ID
BareItemsHandler<BareItemType>::id( const bareItem_Type & item ) const
{
    containerConstIterator_Type i = this->find( item );
    if ( i != container_Type::end() )
        return i->second;
    else
        return 0;
}

// ===================================================
// Set Methods
// ===================================================
template <class BareItemType>
inline
bool  __attribute__ ((__deprecated__))
BareItemsHandler<BareItemType>::setId( const bareItem_Type & item, const ID& id )
{
    containerConstIterator_Type i = find( item );
    if ( i != container_Type::end() )
    {
        i->second = id;
        return true;
    }
    else
    {
        return false;
    }

}


template <class BareItemType>
inline
bool
BareItemsHandler<BareItemType>::setId( const bareItem_Type & item, ID const id )
{
    containerConstIterator_Type i = find( item );
    if ( i != container_Type::end() )
    {
        i->second = id;
        return true;
    }
    else
    {
        return false;
    }

}

// ===================================================
// Get Methods
// ===================================================


}
#endif
