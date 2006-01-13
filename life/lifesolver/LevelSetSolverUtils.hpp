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
   \file LevelSetSolverUtils.hpp
   \author Daniele Antonio Di Pietro <dipietro@unibg.it>
   \date 1-26-2005
*/
#ifndef _LEVELSETSOLVERUTILS_HPP_
#define _LEVELSETSOLVERUTILS_HPP_

#include <functional>

namespace LifeV {

/*! \return the signum of a number */
inline Real signum(Real x) {
    return x >= 0 ? 1. : -1.;
}

//! Compute cross product of two vectors
template<class Vector>
Vector crossProd(const Vector& v1, const Vector& v2) {
    Vector v;
    v[0] =   v1[1] * v2[2] - v2[1] * v1[2];
    v[1] = - v1[0] * v2[2] + v2[0] * v1[2];
    v[2] =   v1[0] * v2[1] - v2[0] * v1[1];

    return v;
}

//! Compute point-to-point distance
template<class Point>
inline Real pointToPointDistance(const Point& P1, const Point& P2) {
    return boost::numeric::ublas::norm_2(P2 - P1);
}

/*!
  Compute point-to-plane distance.
  \param P the point
  \param P0 point defining the plane
  \param n normal vector defining the plane
*/
template<class Point, class Vector>
inline Real pointToPlaneDistance(Point& P, Point& P0, Vector& n) {
    return fabs( ( inner_prod(n, P - P0) ) / boost::numeric::ublas::norm_2(n) );
}

/*!
  Compute a point's projection on a plane. It takes the following arguments:
  @param P point to project on the plane
  @param P0 point the plane passes through
  @param n normal to the plane
  @param Q point projection on the plane
  @param d the distance between point P and the plane (P0, n)
*/
template<class Point, class Vector>
void pointProjectionOnPlane(const Point& P, const Point& P0,
                            const Vector& n, Point& Q, Real& d) {
    Vector versor_n = n / boost::numeric::ublas::norm_2(n);
    d = pointToPlaneDistance(P, P0, versor_n);
    Q = P - signum(inner_prod(versor_n, P - P0)) * d * versor_n;
}

/*!
  Find the point where u is zero on an edge
  @param P1 first segment node
  @param u1 value of u in P1
  @param P2 second segment node
  @param u2 value of u in P2
  @param P the intersection point with the interface
*/
    template<class Point>
    inline void findZeroOnEdge(const Point& P1, const Real u1,
                               const Point& P2, const Real u2, Point& P) {
        Real s = fabs(u1) / (fabs(u1) + fabs(u2));

        P = s * P2 + (1 - s) * P1;
    }

/*!
  A functor to determine whether two points coincide up to a certain tolerance
*/
template<typename Point> class pointsCoincide
    : public std::binary_function<Point, Point, bool> {
public:
    pointsCoincide() { _M_toll = 1e-6; }
    pointsCoincide(Real toll) : _M_toll(toll) {}
    inline bool operator()(const Point& P1, const Point& P2) const {
        return pointToPointDistance(P1, P2) < _M_toll;
    }
private:
    Real _M_toll;
};

//! Point type conversion from scalar x, y, z-like to node_type-like format
template<class PointType1, class PointType2>
void convertPointType(PointType1& Ppt1, const PointType2& Ppt2){
    Ppt1[0] = Ppt2.x();
    Ppt1[1] = Ppt2.y();
    Ppt1[2] = Ppt2.z();
}

} // namespace LifeV

#endif /* _LEVELSETSOLVERUTILS_HPP_ */
