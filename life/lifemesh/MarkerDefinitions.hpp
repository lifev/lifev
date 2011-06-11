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

#include <life/lifemesh/Marker.hpp>

namespace LifeV
{
//! MarkerCommon - A trait class that defines the markers used in RegionMesh
/*!
 * It takes as template parameter the policy for the entityFlag treatment. In fact
 * it is a concrete class that may be redefined by the used if different markers are
 * needed
 */

template
<class MT>
class MarkerCommon
{
public:

    //! @name Public Types
    //@{
    //! The policy used to treat entityFlags
    typedef MT entityFlagPolicy_Type;
    //! The marker used for the Points
    typedef Marker<MT> pointMarker_Type;
    //! The marker used for the Edges
    typedef Marker<MT> edgeMarker_Type;
    //! The marker used for the faces
    typedef Marker<MT> faceMarker_Type;
    //! The marker used for the volumes
    typedef Marker<MT> volumeMarker_Type;
    //! The marker used for the regions
    typedef Marker<MT> regionMarker_Type;

    // Old typedefs to delete

 //   typedef MT MarkerTraits;
    typedef Marker<MT> PointMarker;
    typedef Marker<MT> EdgeMarker;
    typedef Marker<MT> FaceMarker;
    typedef Marker<MT> VolumeMarker;
    typedef Marker<MT> RegionMarker;

    //@}
};

//! The simples MarkerCommon: uses all defaults
typedef MarkerCommon<EntityFlagStandardPolicy> defaultMarkerCommon_Type;

//!Expose NULLFLAG
static const entityFlag_Type S_NULLFLAG = EntityFlagStandardPolicy::S_NULLFLAG;


} // Namespace LifeV

#endif /* MARKERDEFINITIONS_H */
