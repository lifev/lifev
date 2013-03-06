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

#include <lifev/core/mesh/RegionMesh2DStructured.hpp>

namespace LifeV
{

markerID_Type regularMeshPointPosition2D ( const UInt& i_x,
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

    UInt pointPosition (0);

    if ( i_y == 0 )
    {
        pointPosition = pointPosition |  1;
    }
    if ( i_y == n_y - 1)
    {
        pointPosition = pointPosition |  2;
    }
    if ( i_x == n_x - 1)
    {
        pointPosition = pointPosition |  4;
    }
    if ( i_x == 0 )
    {
        pointPosition = pointPosition |  8;
    }

    /*
    // The boundary of the structured mesh
    // labels:

    INTERNAL = 0;

    // Edges
    LEFTEDGE     =  8;
    BOTTOMEDGE   =  1;
    TOPEDGE      =  2;
    RIGHTEDGE    =  4;

    // Corners
    BOTTOMLEFTCORNER   = 9;
    BOTTOMRIGHTCORNER  = 5;
    TOPLEFTCORNER      = 10;
    TOPRIGHTCORNER     = 6;
    */
    switch ( pointPosition )
    {
            // We are inside
        case 0:
            return Structured2DLabel::INTERNAL;
            // We are on 1 edge
        case 1:
            return Structured2DLabel::BOTTOM;
        case 2:
            return Structured2DLabel::TOP;
        case 4:
            return Structured2DLabel::RIGHT;
        case 8:
            return Structured2DLabel::LEFT;
            // We are on a corner
        case 5:
            return Structured2DLabel::BOTTOM_RIGHT;
        case 6:
            return Structured2DLabel::TOP_RIGHT;
        case 9:
            return Structured2DLabel::BOTTOM_LEFT;
        case 10:
            return Structured2DLabel::TOP_LEFT;
        default:
            return Structured2DLabel::INTERNAL;
    }

}

} // Namespace LifeV
