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
