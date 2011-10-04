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

#include <life/lifecore/LifeV.hpp>

#include <life/lifemesh/Marker.hpp>

// LifeV namespace
namespace LifeV
{

//! @class NeighborMarker
/*! This class extends the default marker adding containers to store all adjacency informations.
 */

template <typename FlagPolicy = EntityFlagStandardPolicy>
class NeighborMarker: public Marker<FlagPolicy>
{
public:

    typedef std::set<ID>                   neighbors_Type;
    typedef neighbors_Type::iterator       neighborIterator_Type;
    typedef neighbors_Type::const_iterator neighborConstIterator_Type;

    NeighborMarker (): Marker<FlagPolicy>() {}

    explicit NeighborMarker ( markerID_Type & p ): Marker<FlagPolicy>( p ) {}

    neighbors_Type & nodeNeighbors () { return M_nodeNeighbors; }
//    neighbors_Type & edgeNeighbors () { return M_edgeNeighbors; }
//    neighbors_Type & faceNeighbors () { return M_faceNeighbors; }
//    neighbors_Type & elementNeighbors () { return M_elementNeighbors; }

    void setNodeNeighbors ( neighbors_Type const & nodeNeighbors ) { M_nodeNeighbors = nodeNeighbors; }
//    void setEdgeNeighbors ( neighbors_Type const & edgeNeighbors ) { M_edgeNeighbors = edgeNeighbors; }
//    void setFaceNeighbors ( neighbors_Type const & faceNeighbors ) { M_faceNeighbors = faceNeighbors; }
//    void setElementNeighbors ( neighbors_Type const & elementNeighbors ) { M_elementNeighbors = elementNeighbors; }

protected:
    neighbors_Type M_nodeNeighbors;
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
typedef NeighborMarkerCommon<EntityFlagStandardPolicy> neighborMarkerCommon_Type;

}

#endif // _DATATRACK_H_

