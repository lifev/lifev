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
     @brief

     @date 10/2011
     @author A. Cervone <ant.cervone@gmail.com>
 */

#ifndef _NEIGHBORMARKER_H_
#define _NEIGHBORMARKER_H_ 1

#include <boost/unordered_set.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/Marker.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

// LifeV namespace
namespace LifeV
{

typedef boost::unordered_set<ID> neighbors_Type;
typedef std::vector<neighbors_Type> neighborList_Type;


//! @class NeighborMarker
/*! This class extends the default marker adding containers to store all adjacency informations.
 */

template <typename FlagPolicy = MarkerIDStandardPolicy>
class NeighborMarker: public Marker<FlagPolicy>
{
public:

    typedef std::set<ID>                   neighbors_Type;
    typedef neighbors_Type::iterator       neighborIterator_Type;
    typedef neighbors_Type::const_iterator neighborConstIterator_Type;

    NeighborMarker () : Marker<FlagPolicy>() {}

    explicit NeighborMarker ( markerID_Type& p ) : Marker<FlagPolicy> ( p ) {}

    neighbors_Type& pointNeighbors ()
    {
        return M_pointNeighbors;
    }
    //    neighbors_Type & edgeNeighbors () { return M_edgeNeighbors; }
    //    neighbors_Type & faceNeighbors () { return M_faceNeighbors; }
    //    neighbors_Type & elementNeighbors () { return M_elementNeighbors; }

    void setPointNeighbors ( neighbors_Type const& pointNeighbors )
    {
        M_pointNeighbors = pointNeighbors;
    }
    //    void setEdgeNeighbors ( neighbors_Type const & edgeNeighbors ) { M_edgeNeighbors = edgeNeighbors; }
    //    void setFaceNeighbors ( neighbors_Type const & faceNeighbors ) { M_faceNeighbors = faceNeighbors; }
    //    void setElementNeighbors ( neighbors_Type const & elementNeighbors ) { M_elementNeighbors = elementNeighbors; }

protected:
    neighbors_Type M_pointNeighbors;
    //    neighbors_Type M_edgeNeighbors;
    //    neighbors_Type M_faceNeighbors;
    //    neighbors_Type M_elementNeighbors;
};

//! @class NeighborMarkerCommon
/*! Trait class to collect all mesh entity markers.
 *  The only difference with the default one is a NeighborMarker for points.
 */

template <class MT>
class NeighborMarkerCommon
{
public:

    //! @name Public Types
    //@{
    //! The marker used for the Points
    typedef NeighborMarker<MT> pointMarker_Type;
    //! The marker used for the Edges
    typedef Marker<MT> edgeMarker_Type;
    //! The marker used for the Faces
    typedef Marker<MT> faceMarker_Type;
    //! The marker used for the Volumes
    typedef Marker<MT> volumeMarker_Type;
    //! The marker used for the Regions
    typedef Marker<MT> regionMarker_Type;
    //@}
};

//! The NeighborMarkerCommon: uses all defaults except for Points
typedef NeighborMarkerCommon<MarkerIDStandardPolicy> neighborMarkerCommon_Type;

//! this routine generates point neighbors for the given mesh
/*! the routine assumes that the mesh is not yet partitioned or reordered
 *  (i.e. the local id and the global id are the same).
 *  if this is not true the method should be changed to use a more
 *  expensive STL find on the mesh points to get the correct point that has
 *  the given global id or construct a globalToLocal map beforehand.
 */
template <typename MeshType>
void createPointNeighbors ( MeshType& mesh )
{
    // TODO: ASSERT_COMPILE_TIME that MeshType::pointMarker == NeighborMarker
    // this guarantees that the pointNeighbors structure is available.

    // generate point neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < mesh.numEdges(); ie++ )
    {
        ID id0 = mesh.edge ( ie ).point ( 0 ).id();
        ID id1 = mesh.edge ( ie ).point ( 1 ).id();

        ASSERT ( mesh.point ( id0 ).id() == id0 , "the mesh has been reordered, the point must be found" );
        ASSERT ( mesh.point ( id1 ).id() == id1 , "the mesh has been reordered, the point must be found" );

        mesh.point ( id0 ).pointNeighbors().insert ( id1 );
        mesh.point ( id1 ).pointNeighbors().insert ( id0 );
    }
}

template <typename MeshType>
void createPointNeighbors ( MeshType const& mesh, neighborList_Type& neighborList )
{
    neighborList.resize ( mesh.numGlobalPoints() );
    // generate point neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < mesh.numEdges(); ie++ )
    {
        ID id0 = mesh.edge ( ie ).point ( 0 ).id();
        ID id1 = mesh.edge ( ie ).point ( 1 ).id();

        ASSERT ( mesh.point ( id0 ).id() == id0 && mesh.point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        neighborList[ id0 ].insert ( id1 );
        neighborList[ id1 ].insert ( id0 );
    }
}

} // namespace LifeV

#endif // _NEIGHBORMARKER_H_

