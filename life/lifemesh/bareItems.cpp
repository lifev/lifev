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

#include <life/lifecore/life.hpp>
#include <life/lifemesh/bareItems.hpp>

namespace LifeV
{
std::pair<BareFace, bool>
makeBareFace( ID const i, ID const j, ID const k )
{
    if ( i < j && i < k )
    {
        if ( j < k )
        {
            return std::make_pair( BareFace( i, j, k ), true );
        }
        else
        {
            return std::make_pair( BareFace( i, k, j ), false );
        }
    }
    else if ( j < k && j < i )
    {
        if ( k < i )
        {
            return std::make_pair( BareFace( j, k, i ), true );
        }
        else
        {
            return std::make_pair( BareFace( j, i, k ), false );
        }
    }
    else
    {
        if ( i < j )
        {
            return std::make_pair( BareFace( k, i, j ), true );
        }
        else
        {
            return std::make_pair( BareFace( k, j, i ), false );
        }
    }
}


std::pair<BareFace, bool>
makeBareFace( ID const i, ID const j, ID const k, ID const l )
{
    std::vector<ID> helper( 4 );
    helper[ 0 ] = i;
    helper[ 1 ] = j;
    helper[ 2 ] = k;
    helper[ 3 ] = l;
    std::vector<ID>::iterator vi = std::max_element( helper.begin(), helper.end() );
    std::rotate( helper.begin(), vi, helper.end() );
    return makeBareFace( helper[ 1 ], helper[ 2 ], helper[ 3 ] );
}
}
