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

#include <functional>

namespace LifeV {
    
    /**
       \Compute cross product of two vectors
    */

    template<class Vector>
    Vector cross_prod(Vector& v1, Vector& v2) {
        Vector v;
        v[0] = v1[1] * v2[2] - v2[1] * v1[2];
        v[1] = - v1[0] * v2[2] + v2[0] * v1[2];
        v[2] = v1[0] * v2[1] - v2[0] * v1[1];
        
        return v;
    }

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

    /**
       \Compute a point's projection on a plane. It takes the following 
       \arguments:
       \P  : point to project on the plane
       \P0 : point the plane passes through
       \n  : normal to the plane
       \Q  : point projection on the plane
    */

    template<class Point, class Vector>
    void point_projection_on_plane(Point& P, Point& P0, Vector& n, Point& Q) {
        Real d = - inner_prod(n, P0);

        Real D = pow(boost::numeric::ublas::norm_2(n), 2);

        Real a = n[0];
        Real b = n[1];
        Real c = n[2];

        Real B2 = b * P[0] - a * P[1];
        Real B3 = c * P[0] - a * P[2];

        Q[0] = (- a * d + b * B2 + c * B3) / D;
        Q[1] = (- a * b * d - a * a * B2 + b * c * B3 - c * c * B2) / (a * D);
        Q[2] = (- a * c * d - a * a * B3 - b * b * B3 + b * c * B2) / (a * D);
    }

    /**
       \Determine whether a point belongs to a (triangular) face.
    */

    template<class Point, class FACE>
    bool point_is_on_face(Point& P, FACE& f) {

        Real a11 = f.point(2)[0] - f.point(1)[0];
        Real a12 = f.point(3)[0] - f.point(1)[0];
        Real a21 = f.point(2)[1] - f.point(1)[1];
        Real a22 = f.point(3)[1] - f.point(1)[1];

        Real b1 = P[0] - f.point(1)[0];
        Real b2 = P[1] - f.point(1)[1];

        Real D = a11 * a22 - a12 * a21;

        Real xi = ( b2 * a12 - b1 * a22 ) / D;
        Real eta = ( b2 * a11 - b1 * a12 ) / D;

        return (xi >= 0 && xi <= 1 && eta >= 0 && eta <= 1 - xi);
    }

    /**
       \Find the point where u is zero on an edge
       \P1  : first segment node
       \P2  : second segment node
       \P   : the intersection point with the interface
    */

    template<class Point>
    inline void find_zero_on_edge(Point& P1, Real u1, Point& P2, Real u2, Point& P) {
        Real s = fabs(u1) / (fabs(u1) + fabs(u2));

        P = s * P2 + (1 - s) * P1;
    }

    /**
       \A functor to determine whether two points coincide up to a certain tolerance
    */

    template<typename Point> class points_coincide : public std::binary_function<Point, Point, bool> {
    public:
        points_coincide() { _M_toll = 1e-6; }
        points_coincide(Real toll) : _M_toll(toll) {}
        inline bool operator()(const Point& P1, const Point& P2) const {
            return point_to_point_distance(P1, P2) < _M_toll;
        }
    private:
        Real _M_toll;
    };

    /**
       \Point type conversion from scalar x, y, z-like to node_type-like format
    */

    template<class PointType1, class PointType2>
    void convert_point_type(PointType1& Ppt1, PointType2& Ppt2){
        Ppt1[0] = Ppt2.x();
        Ppt1[1] = Ppt2.y();
        Ppt1[2] = Ppt2.z();
    }

}
#endif
