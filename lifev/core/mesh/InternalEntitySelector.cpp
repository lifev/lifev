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
    @brief This file implements the standard selector for internal entities

    @author
    @contributor Nur Aiman Fadel <nur.fadel@mail.polimi.it>
    @maintainer Nur Aiman Fadel <nur.fadel@mail.polimi.it>

    @date

    A more detailed description of the file (if necessary)
 */

#include <lifev/core/mesh/InternalEntitySelector.hpp>

namespace LifeV
{

/*! @todo This class was meant to separate internal from boundary flags. With the
 *  new way of selecting boundary entities this will be useless!
 */
const markerID_Type InternalEntitySelector::defMarkFlag ( markerID_Type ( 100000 ) );

// ===================================================
// Constructors & Destructor
// ===================================================
InternalEntitySelector::InternalEntitySelector() :
    M_watermarkFlag ( defMarkFlag )
{}

InternalEntitySelector::InternalEntitySelector ( const markerID_Type& w ) :
    M_watermarkFlag ( w )
{}

SetFlagAccordingToWatermarks::SetFlagAccordingToWatermarks ( const flag_Type& flagToSet,
                                                             const std::vector<markerID_Type>& watermarks,
                                                             const flagPolicy_ptr& flagPolicy) :
    M_flagToSet (flagToSet),
    M_watermarks (watermarks),
    M_flagPolicy (flagPolicy)
{
    std::sort (M_watermarks.begin(), M_watermarks.end() );
}

// ===================================================
// Operators
// ===================================================
bool InternalEntitySelector::operator() ( markerID_Type const& test ) const
{
    return ( test == markerID_Type ( 0 ) || test > M_watermarkFlag );
}

//========================================================
// Methods
//========================================================
void SetFlagAccordingToMarkerRanges::insert ( rangeID_Type const& key, flag_Type flag )
{
    M_map[ key ] = flag;
}

std::pair<flag_Type, bool> SetFlagAccordingToMarkerRanges::findFlag ( markerID_Type const& m ) const
{
    if ( !M_map.empty() )
    {
        const_iterator_Type it = M_map.upper_bound ( std::make_pair ( m, markerID_Type ( 0 ) ) );
        if ( it-- != M_map.begin() ) // go back one
        {
            markerID_Type first  = it->first.first;
            markerID_Type second = it->first.second;
            if ( m >= first && m <= second )
            {
                return std::make_pair ( it->second, true );
            }
        }
    }
    return std::make_pair ( flag_Type ( 0 ), false );
}

} // Namespace LifeV
