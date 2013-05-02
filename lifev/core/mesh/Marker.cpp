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
    @brief Implementations for Marker.hpp

    @contributor Luca Bertagna <lbertag@emory.edu>
    @date 00-00-0000

 */

#include <limits>
#include <lifev/core/mesh/Marker.hpp>

namespace LifeV
{

//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************

///////////////////////
// MarkerIDStandardPolicy //
///////////////////////

const markerID_Type MarkerIDStandardPolicy::S_NULLMARKERID =
    std::numeric_limits<Int>::max();

//MM: if you modify these changes here recheck function readNetgenMesh
//        because it uses this changes

markerID_Type MarkerIDStandardPolicy::strongerMarkerID ( markerID_Type const& markerID1, markerID_Type const& markerID2 )
{
    if ( markerID1 == S_NULLMARKERID )
    {
        return markerID2;
    }
    if ( markerID2 == S_NULLMARKERID )
    {
        return markerID1;
    }
    return markerID1 > markerID2 ? markerID1 : markerID2 ;
}

markerID_Type MarkerIDStandardPolicy::weakerMarkerID ( markerID_Type const& markerID1, markerID_Type const& markerID2 )
{
    if ( markerID1 == S_NULLMARKERID )
    {
        return markerID2;
    }
    if ( markerID2 == S_NULLMARKERID )
    {
        return markerID1;
    }
    return markerID1 < markerID2 ? markerID1 : markerID2 ;
}

bool MarkerIDStandardPolicy::equalMarkerID (const markerID_Type& markerID1, const markerID_Type& markerID2)
{
    return markerID1 == markerID2;
}

} // Namespace LifeV
