/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Daniele Antonio Di Pietro <dipietro@unibg.it>
Date: 26-1-2005

Copyright (C) 2005 Università degli Studi di Bergamo

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
/**
   \file LevelSetHandler_utils.hpp
   \author Daniele Antonio Di Pietro <dipietro@unibg.it>
   \date 1-26-2005
*/
#ifndef _LEVELSETSOLVERUTILS_H_
#define _LEVELSETSOLVERUTILS_H_

namespace LifeV {
    /**
       \Compute point-to-point distance
    */

    template<class Point>
    inline Real point_to_point_distance(Point& P1, Point& P2) {
        return boost::numeric::ublas::norm_2(P2 - P1);
    }

    /**
       \Compute point-to-plane distance.
       \P    : the point
       \P0, n: point and normal vector defining the plane
  
    */

    template<class Point, class Vector>
    inline Real point_to_plane_distance(Point& P, Point& P0, Vector& n) {
        Real d = - inner_prod(n, P0);
     
        return fabs( ( inner_prod(n, P) + d ) / boost::numeric::ublas::norm_2(n) );
    }
}

#endif
