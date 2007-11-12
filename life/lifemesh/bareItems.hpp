/*
This file is part of the LifeV library
Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*! \file bareItems.h
\brief Special structures for handling mesh faces and sides
\version 0.0 Experimental   19/8/99. Luca Formaggia

Classes BareFace and BareEdge have been created to give an UNIQUE internal
representation for mesh faces and edges, allowing thus the construction of
Dof objects (which are naturally linked to mesh entities).

\par Introduction One of the paradigms chosen for the development of this
library is the fact that degrees of freedom (Dof) are linked to geometrical
entities.  Now if we have degrees of freedom associated, for instance, to
an Edge (like in a P2 Tetra) in order to build the global numbering of the
Dof and the association between local (element-wise) and global numbering,
I need to identify edges and give them a unique id_type. Yet, I may not want
want to build a full Edge (GeoElement2D) object: after all that I need is
the id_type of the edge and a way of computing the id_type's of the degrees of
freedom on the edge, all the remaining data of a full Edge object is not
necessarily needed.

Another related problem is how to uniquely identify a face or an edge in the mesh.


The dilemma has been resolved by creating the concept of a BareEdge and
BareFace (bare geometry items).  A bare geometry item is formed by the
minimal information required to uniquely identify it, namely 2
<tt>Point</tt>'s id_type 's for an edge and 3 <tt>Point</tt>'s id_type 's for the
Faces (it is enough also for Quad faces!). We build the bare items by
looping through the elements and obviously we make sure that the BareItem
id_type is consistent with that of the corresponding ``full item'' if the latter
has been instantiated.

Another <em>very important</em> issue is that of orientation.  There are
different ways of considerin orientation of a Face or an Edge. The first is
the <em>local</em> orientation of a Face or Egde of the reference finite
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
are stored so that the normal (always calculated following the righ hand
rule) is outward with respect to the domain.

However, there is the need of defining a <em>global</em> orientation also
for <em>non active</em> entities. This because we do not always store all
faces and all edges. We need then to choose a unique way to identify the
orientation of an Edge or of a Face <em>independently </em> from the fact
that they are active or not. We will call this orientation the <em>natural </em>
orientation. We have chosen the following
convention for natural orientation of faces and edges

<ul>
<li>The positive natural orientation of an  <em>Edge</em> is given by \f$V_{min} \rightarrow V_{max} \f$,
    \f$V_{min}\f$ being the Edge Vertex with smallest id_type</li>

<li>The positive natural orientation of a  <em>Face</em> is given by the cicle
    \f$V_{min} \rightarrow V_2\rightarrow V_3 \f$,
    \f$V_{min}\f$ being the Face Vertex with smallest id_type, \f$V_2\f$ the second smallest
    and \f$V_2\f$ the thirsd smallest.</li>
</ul>

Note that the latter definition applies both to triangular and to quad faces.

\warning If I want to associate boundary conditions I need the active
entity, since BareEdges do not store Marker data. (This is why the
RegionMesh classes treat boundary items in a rather special way).



\par Usage

Since BareItems are used only internally, we are (temporarily) omitting a
detailed documentation.
*/

#ifndef _MESHBAREITEMS_HH_
#define _MESHBAREITEMS_HH_
/*!
  Special routines to read meshes and special structures for
  sides and faces handling

  Classes BareFace and BareEdge have been created to give an UNIQUE
  Representation for mesh faces and edges and thus allow the construction
  of global tables or Fields.

*/
#include<utility>
#include<vector>
#include<map>
#include<algorithm>
#include<iostream>

#include <life/lifecore/life.hpp>

namespace LifeV
{
#ifndef _LIFEV_HH_
/*! \typedef typedef unsigned int id_type
\brief type used for Identifiers (Integral type in the range [1, MAX_INT]).
All principal items handled by the library have an identified, which is an
unsigned integer <em>greater or equal to one</em>
*/
// more correct version
typedef size_t id_type;
// original version
//typedef unsigned int id_type;
//! I define UInt=unsigned int. This allow to use this stuff also outside lifeV
// more correct version
//typedef size_t UInt;
// original version
typedef unsigned int UInt;
#endif

//! The Edge basis class
/*! It contains the attributes common to all Edges In particular, it
contains the two id_type's (first and second) of the points at the two ends
of the edge.
\invariant first < second
*/
struct BareEdge
{
    BareEdge() : first( 0 ), second( 0 )
    {}
    ; /*!< Standard Constructor*/
    BareEdge( id_type i, id_type j ) : first( i ), second( j )
    {}
    ; /*!< constructor taking the id_type's */
    id_type first; //!< First id_type which defines the Edge
    id_type second; //!< Second id_type which defines the Edge
};

//! The base Face class
/*! Contains the common parameters of all Edges A Face contains three
  id_type's (first, second and third) of points at face vertex, chosen so
  to uniquely identify the face.
\invariant first<second<third.

*/
struct BareFace
{
    BareFace() : first( 0 ), second( 0 ), third( 0 )
    {}
    ; //!< Default constructor
    //! Constructor that takes the id_type's as parameter
    BareFace( id_type i, id_type j, id_type k ) : first( i ), second( j ), third( k )
    {}
    ;
    /*! Constructor that takes a BareEdge Object and an id_type. The face is
      then identified by the id_type of the Points on the BareEdge +
      <tt>i</tt>*/
    BareFace( id_type i, const BareEdge & e ) : first( i ), second( e.first ), third( e.second )
    {}
    ;
    id_type first; //!< First id_type which defines the BareFace
    id_type second; //!< Second id_type which defines the BareFace
    id_type third; //!< Third id_type which defines the BareFace
};


//! \defgroup BareItemsBuilder Global functions to build Bare Items.

/*! \ingroup BareItemsBuilder
  \brief It creates a BareEdge end returns the orientation of the created edge with respect
  to the given data.
  \param  bool is false if orientation  has been changed.
  \param i is a Point id_type
  \param j is a Point id_type

  The BareEdge that will be built is the one passing by <tt>i</tt> and
  <tt>j</tt>. The orientation is i->j if the returned parameter is a
  true.

  \pre i and j >0, i!=j */
inline
std::pair<BareEdge, bool>
makeBareEdge( id_type const i, id_type const j )
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
  \param i is a Point id_type
  \param j is a Point id_type
  \brief It creates a BareEdge, ignoring orientation.

  The BareEdge that will be built is the one passing by <tt>i</tt> and <tt>j</tt>
  A lighter version of MakeBareEdge, to be used
  if orientation flag is not needed;
  \pre i and j >0, i!=j
*/
inline
BareEdge
setBareEdge( id_type const i, id_type const j )
{
    BareEdge be;
    be.first = i < j ? i : j;
    be.second = ( i - be.first ) + j;
    return be;
}

/*! \ingroup BareItemsBuilder
  \param  bool is false if orientation  has been changed.
  \param i is a Point id_type
  \param j is a Point id_type
  \pre i and j >0
  \brief It creates a non-standard BareEdge.

  Yet another  lighter version of MakeBareEdge, without orientation, To be used for
  non-oriented graphs.

  \warning It produces a BareEdge which does not comply with the invariant of the
  class (first < second). It must be used only if the BareEdge class is NOT used to uniquely identify edges.
 */
inline
BareEdge
setBareEdgeNo( id_type const i, id_type const j )
{
    return BareEdge( i, j );
}


/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from three Point id_type's
  \param  bool is false if orientation  has been changed.
  \param i is a Point id_type
  \param j is a Point id_type
  \param k is a Point id_type

  To be used for triangular faces.
  \pre i, j and k >0. i!=j!=k

*/
std::pair<BareFace, bool> makeBareFace( id_type const i, id_type const j, id_type const k );

/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from four Point id_type's
  \param  bool is false if orientation  has been changed.
  \param i is a Point id_type
  \param j is a Point id_type
  \param k is a Point id_type
  \param l is a Point id_type

  To be used for triangular faces.
  \pre i, j and k >0. i!=j!=k

*/
std::pair<BareFace, bool> makeBareFace( id_type const i, id_type const j, id_type const k, id_type const l );

/*! \defgroup comparison Comparison Operators
  Operators for comparing BareItems
*/

/*! \ingroup comparison
    inequality
*/
inline
bool
operator!=( const BareEdge & p1 , const BareEdge & p2 )
{
    return p1.first != p2.first || p1.second != p2.second;
}

/*! \ingroup comparison
     equality
*/
inline
bool
operator==( const BareEdge & p1 , const BareEdge & p2 )
{
    return p1.first == p2.first && p1.second == p2.second;
}

/*! \ingroup comparison
    greater than
*/
inline
bool
operator>( const BareEdge & e1 , const BareEdge & e2 )
{
    return e2.first > e1.first || ( e2.first == e1.first && e2.second > e1.second );
}

/*! \ingroup comparison
     greater-equal than
*/
inline
bool
operator>=( const BareEdge & e1 , const BareEdge & e2 )
{
    return e1 == e2 || e1 > e2;
}

/*! \ingroup comparison
    less than
*/
inline
bool
operator<( const BareEdge & e1 , const BareEdge & e2 )
{
    return e2.first < e1.first || ( e2.first == e1.first && e2.second < e1.second );
}

/*! \ingroup comparison
     less-equal than
*/
inline
bool
operator<=( const BareEdge & e1 , const BareEdge & e2 )
{
    return e1 == e2 || e1 < e2;
    ;
}

/*! \ingroup comparison*/
inline
bool
operator!=( const BareFace & p1 , const BareFace & p2 )
{
    return p1.first != p2.first || p1.second != p2.second || p1.third != p2.third;
}


/*! \ingroup comparison*/
inline
bool
operator==( const BareFace & p1 , const BareFace & p2 )
{
    return p1.first == p2.first && p1.second == p2.second && p1.third == p2.third;
}

/*! \ingroup comparison
  General functor for lexicographic comparison
*/
template <typename T>
struct cmpBareItem;

/*! \ingroup comparison
   Specialised functor for Edges
*/
// Specialisations
template <>
struct cmpBareItem<BareEdge> //!< The actual comparison operator
{
    bool operator() ( const BareEdge & e1, const BareEdge & e2 ) const
    {
        return e2.first > e1.first || ( e2.first == e1.first && e2.second > e1.second );
    }
};

/*! \ingroup comparison
  Specialised functor for Faces
*/
template <>
struct cmpBareItem<BareFace>
{
    bool operator() ( const BareFace & e1, const BareFace & e2 ) const
    {
        if ( e2.first > e1.first )
            return true;
        if ( e2.first == e1.first )
        {
            if ( e2.second > e1.second )
                return true;
            if ( e2.second == e1.second )
                return e2.third > e1.third;
        }
        return false;
    }
};

//! Handles Bare Items
/*! This class handles mesh bare edges and faces construction. Used only in mesh builders
 A BareItemsHandler is a specialisation of a STL map which holds the pair
 formed by a bareitem and its id_type.
 The id_type  is automatically generated if one uses the method  addIfNotThere
*/
template <typename BareItem>
class BareItemsHandler: public std::map<BareItem, id_type, cmpBareItem<BareItem> >
{
public:
    typedef std::map<BareItem, UInt, cmpBareItem<BareItem> > container;
    typedef typename container::iterator iterator;
    typedef typename container::const_iterator const_iterator;
    typedef std::pair<const BareItem, UInt> value_type;

    BareItemsHandler();
    bool isThere( BareItem const & ) const; //!< is the item there? I just ask
    id_type id( BareItem const & ) const; //!< Returns id_type of a BareItem. 0 if not there
    bool setId( BareItem const & item, id_type const i ); //!< To modify id_type of bareitem item in the list
    std::pair<id_type, bool> addIfNotThere( BareItem const & ); //!< if not there adds it, the item id_type is autogenerated
    std::pair<id_type, bool> addIfNotThere( BareItem const &, const id_type id ); //!<if not there adds it, and sets id_type id
    bool isThereDel( BareItem const & ); //!< if it is there take it out (Id is lost)

    UInt howMany() const; //!< The # of entities ones actually stored.
    UInt maxId() const
    {
        return _idCount;
    } //!< Max id_type currently in use
    void showMe() const; //!< Writes info in output
private:
    UInt _idCount;
};

/*********************************************************************************
               IMPLEMENTATIONS
*********************************************************************************/
//
/*! \defgroup Helper Some helper functions
*/

//!\ingroup Helper
template <typename BareItem>
inline
id_type getId( std::pair<BareItem, id_type> const & i )
{
    return i.second;
}

//!\ingroup Helper
template <typename BareItem>
inline
BareItem getItem( std::pair<BareItem, id_type> const & i )
{
    return i.first;
}


//                   BareItemsHandler
template <class BareItem>
BareItemsHandler<BareItem>::BareItemsHandler() :
        _idCount( 0 )
{ }


template <class BareItem>
inline
bool
BareItemsHandler<BareItem>::isThere( const BareItem & s ) const
{
    return find( s ) != container::end();
}

template <class BareItem>
inline
bool
BareItemsHandler<BareItem>::setId( const BareItem & s, id_type const id )
{
    const_iterator i = find( s );
    if ( i != container::end() )
    {
        i->second = id;
        return true;
    }
    else
    {
        return false;
    }

}

template <class BareItem>
inline
id_type
BareItemsHandler<BareItem>::id( const BareItem & s ) const
{
    const_iterator i = this->find( s );
    if ( i != container::end() )
        return i->second;
    else
        return 0;
}

template <class BareItem>
inline
std::pair<id_type, bool>
BareItemsHandler<BareItem>::addIfNotThere( const BareItem & s )
{
    std::pair<typename BareItemsHandler<BareItem>::iterator, bool> i( insert( std::make_pair( s, _idCount + 1 ) ) );
    if ( i.second )
        ++_idCount;
    return std::make_pair( ( i.first )->second, i.second );
}

template <class BareItem>
inline
std::pair<id_type, bool>
BareItemsHandler<BareItem>::addIfNotThere( const BareItem & s, const id_type id )
{
    std::pair<typename BareItemsHandler<BareItem>::iterator, bool> i( insert( std::make_pair( s, id ) ) );
    ( i.first ) ->second = id; // Set new id in any case.
    return std::make_pair( id, i.second ); // for consistency with other version.
}

template <class BareItem>
inline
UInt
BareItemsHandler<BareItem>::howMany() const
{
    return container::size();
}

template <class BareItem>
bool
BareItemsHandler<BareItem>::isThereDel( BareItem const & s )
{
    return erase( s ) != 0;
}

template <typename BareItem>
inline
void BareItemsHandler<BareItem>::showMe() const
{
    std::cout << "BareItemsHandler: " << std::endl;
    std::cout << "Number of Items stored: " << this->size() << std::endl;
    std::cout << "Max Id stored         : " << this->maxId() << std::endl;
    std::cout << "End of Information";
}
}
#endif

