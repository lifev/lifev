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
    @brief Contains methods which generate structured meshes.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @contributor -
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 16-04-2010

    Such methods will be usefull in order to test problems at different
    scales.
 */

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>

namespace LifeV
{

markerID_Type regularMeshPointPosition ( const UInt& i_x,
                                         const UInt& i_y,
                                         const UInt& i_z,
                                         const UInt& n_x,
                                         const UInt& n_y,
                                         const UInt& n_z )
{
    // We will use a binary representation to
    // find the position of the point
    // 000000 = in the cube        ( 0)
    // 000001 = on the x-z plane   ( 1)
    // 000010 = on the plane x=n_x ( 2)
    // 000100 = on the plane y=n_y ( 4)
    // 001000 = on the plane y-z   ( 8)
    // 010000 = on the plane x-y   (16)
    // 100000 = on the plane z=n_z (32)
    UInt pointPosition (0);
    if ( i_y == 0 )
    {
        pointPosition = pointPosition |  1;
    }
    if ( i_x == n_x - 1)
    {
        pointPosition = pointPosition |  2;
    }
    if ( i_y == n_y - 1)
    {
        pointPosition = pointPosition |  4;
    }
    if ( i_x == 0 )
    {
        pointPosition = pointPosition |  8;
    }
    if ( i_z == 0 )
    {
        pointPosition = pointPosition | 16;
    }
    if ( i_z == n_z - 1 )
    {
        pointPosition = pointPosition | 32;
    }

    switch ( pointPosition )
    {
            // We are inside
        case 0:
            return 0;
            // We are on 1 face
        case 1:
            return 1;
        case 2:
            return 2;
        case 4:
            return 3;
        case 8:
            return 4;
        case 16:
            return 5;
        case 32:
            return 6;
            // We are on an edge
        case 17:
            return 7;
        case 18:
            return 8;
        case 20:
            return 9;
        case 24:
            return 10;
        case 9:
            return 11;
        case 3:
            return 12;
        case 6:
            return 13;
        case 12:
            return 14;
        case 33:
            return 15;
        case 34:
            return 16;
        case 36:
            return 17;
        case 40:
            return 18;
            // We are on a corner
        case 25:
            return 19;
        case 19:
            return 20;
        case 22:
            return 21;
        case 28:
            return 22;
        case 41:
            return 23;
        case 35:
            return 24;
        case 38:
            return 25;
        case 44:
            return 26;
        default:
            return 0;
    }
}


} // Namespace LifeV
