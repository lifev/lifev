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
#include "mesh_util_base.hpp"

namespace LifeV
{
GetCoordComponent::GetCoordComponent() : comp( -1 )
{}

GetCoordComponent::GetCoordComponent( int i ) : comp( i )
{}

void GetCoordComponent::operator() ( Real const x, Real const y, Real const z, Real ret[ 3 ] ) const
{
    switch ( comp )
    {
    case( 0 ) :
                    ret[ 0 ] = x;
        ret[ 1 ] = 0.0;
        ret[ 2 ] = 0.0;
        break;
    case( 1 ) :
                    ret[ 0 ] = 0.0;
        ret[ 1 ] = y;
        ret[ 2 ] = 0.0;
        break;
    case( 2 ) :
                    ret[ 0 ] = 0.0;
        ret[ 1 ] = 0.0;
        ret[ 2 ] = z;
        break;
    default:
        ret[ 0 ] = x;
        ret[ 1 ] = y;
        ret[ 2 ] = z;
    }
}

void GetOnes::operator() ( Real const /*x*/, Real const /*y*/, Real const /*z*/, Real ret[ 3 ] )
const
{
    ret[ 0 ] = 1.0;
    ret[ 1 ] = 1.0;
    ret[ 2 ] = 1.0;
}
}
