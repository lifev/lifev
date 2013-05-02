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
    @brief A simple implementations of Markers

    @date 00-00-0000
    @author Luca Formaggia
    @contributor Luca Bertagna <lbertag@emody.edu>

    This is the simplest implementation of the markers, which just adopts the
    basis marker classes defined in marker_base.h.

    Specialised markers can be implemented using this file as a
    reference.
*/

#ifndef MARKERDEFINITIONS_H
#define MARKERDEFINITIONS_H 1

#include <lifev/core/mesh/Marker.hpp>

namespace LifeV
{
//! MarkerCommon - A trait class that defines the markers used in RegionMesh
/*!
 * It takes as template parameter the policy for the entityFlag treatment.
 * It is a concrete class that may be redefined by the user if different markers are
 * needed. It is passed as template argument to the RegionMesh classes
 */

template
<class MT>
class MarkerCommon
{
public:

    //! @name Public Types
    //@{
    //! The marker used for the Points
    typedef Marker<MT> pointMarker_Type;
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

//! The simplest MarkerCommon: uses all defaults
typedef MarkerCommon<MarkerIDStandardPolicy> defaultMarkerCommon_Type;

} // Namespace LifeV

#endif /* MARKERDEFINITIONS_H */
