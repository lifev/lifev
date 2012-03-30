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
    @brief Contains methods which generate 2D structured meshes.

    @author Iori Guido <guido.iori@mail.polimi.it>
    @contributor -

    @date 23-05-2011

 */

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>

namespace LifeV
{

markerID_Type regularMeshPointPosition2D( const UInt& i_x,
                                     const UInt& i_y,
                                     const UInt& n_x,
                                     const UInt& n_y )
{
    // We will use a binary representation to
    // find the position of the point
    // 0000 = in the rectangle        ( 0)
    // 0001 = on the y=0 axis         ( 1)
    // 0010 = on the y=l_y axis       ( 2)
    // 0100 = on the x=l_x axis       ( 4)
    // 1000 = on the x=0 axis         ( 8)
    
    UInt pointPosition(0);
    if ( i_y == 0 )       { pointPosition = pointPosition |  1; }
    if ( i_y == n_y - 1)  { pointPosition = pointPosition |  2; }
    if ( i_x == n_x - 1)  { pointPosition = pointPosition |  4; }
    if ( i_x == 0 )       { pointPosition = pointPosition |  8; }

    switch ( pointPosition )
    {
        // We are inside
    case 0:
        return 0;
        // We are on 1 edge
    case 1:
        return 1;
    case 2:
        return 2;
    case 4:
        return 3;
    case 8:
        return 4;
        // We are on a corner
    case 5:
        return 5;
    case 6:
        return 6;
    case 9:
        return 7;
    case 10:
        return 8;
    
    default:
        return 0;
    }
    
    /*
    // The boundary of the structured mesh
    // labels:

    INTERNAL = 0;

    // Edges
    BOTTOMEDGE  =  1;
    LEFTEDGE    =  2;
    RIGHTEDGE   =  4;
    TOPEDGE     =  8;

    // Corners
    BOTTOMLEFTCORNER   = 9;
    BOTTOMRIGHTCORNER  = 5;
    TOPLEFTCORNER      = 10;
    TOPRIGHTCORNER     = 6;
    */
    
}


} // Namespace LifeV
