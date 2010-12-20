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

#include <climits>
#include <life/lifemesh/Marker.hpp>

namespace LifeV
{

//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************

///////////////////////
// MarkerTraits //
///////////////////////

const MarkerTraits::entityFlag_Type MarkerTraits::S_NULLFLAG = LONG_MIN;

//MM: if you modify these changes here recheck function readNetgenMesh
//        because it uses this changes

MarkerTraits::entityFlag_Type MarkerTraits::strongerFlag( entityFlag_Type const & flag1, entityFlag_Type const & flag2 )
{
    return flag1 > flag2 ? flag1 : flag2 ;
}

MarkerTraits::entityFlag_Type MarkerTraits::weakerFlag( entityFlag_Type const & flag1, entityFlag_Type const & flag2 )
{
    if ( flag1 == S_NULLFLAG )
        return flag2;
    if ( flag2 == S_NULLFLAG )
        return flag1;
    return flag1 < flag2 ? flag1 : flag2 ;
}

bool MarkerTraits::EqualFlags(const entityFlag_Type& flag1, const entityFlag_Type& flag2)
{
    return flag1 == flag2;
}

} // Namespace LifeV
