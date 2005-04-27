/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2004-11-12

 Copyright (C) 2004 EPFL

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
   \file ud_functions.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-11-12
*/

#include <ud_functions.hpp>

namespace LifeV {

Real EthierSteinmanSteady::f(const Real& t, const Real& x, const Real& y,
                             const Real& z, const ID& i) {
    switch(i) {
        case 1:
            return
                -2.0*nu*exp(a*x-a*z+b*y-b*z)*b*b*b
                +2.0*nu*exp(a*z-a*y+b*x-b*y)*a*b*b
                +2.0*nu*exp(a*z-a*y+b*x-b*y)*a*a*b
                -2.0*nu*exp(a*x-a*z+b*y-b*z)*b*a*a
                +2.0*nu*exp(a*z-a*y+b*x-b*y)*a*a*a
                -2.0*nu*exp(a*x-a*z+b*y-b*z)*b*b*a
                + sigma*uexact(t,x,y,z,i);
            break;
        case 2:
            return
                -2.0*nu*exp(a*y-a*x+b*z-b*x)*b*b*b
                +2.0*nu*exp(a*x-a*z+b*y-b*z)*b*a*a
                +2.0*nu*exp(a*x-a*z+b*y-b*z)*a*a*a
                +2.0*nu*exp(a*x-a*z+b*y-b*z)*b*b*a
                -2.0*nu*exp(a*y-a*x+b*z-b*x)*b*b*a
                -2.0*nu*exp(a*y-a*x+b*z-b*x)*b*a*a
                + sigma*uexact(t,x,y,z,i);
            break;
        case 3:
            return
                -2.0*nu*exp(a*z-a*y+b*x-b*y)*b*b*b
                +2.0*nu*exp(a*y-a*x+b*z-b*x)*a*a*a
                -2.0*nu*exp(a*z-a*y+b*x-b*y)*a*b*b
                -2.0*nu*exp(a*z-a*y+b*x-b*y)*a*a*b
                +2.0*nu*exp(a*y-a*x+b*z-b*x)*b*b*a
                +2.0*nu*exp(a*y-a*x+b*z-b*x)*b*a*a
                + sigma*uexact(t,x,y,z,i);
            break;
    }
    exit(1);
}

Real EthierSteinmanSteady::uexact(const Real& /* t */,
                                  const Real& x, const Real& y,
                                  const Real& z, const ID& i) {
    switch(i) {
        case 1:
            return
                b*exp(a*(x-z)+b*(y-z))-
                a*exp(a*(z-y)+b*(x-y));
            break;
        case 2:
            return
                b*exp(a*(y-x)+b*(z-x))-
                a*exp(a*(x-z)+b*(y-z));
            break;
        case 3:
            return
                b*exp(a*(z-y)+b*(x-y))-
                a*exp(a*(y-x)+b*(z-x));
            break;
    }

    exit(1);
}

Real EthierSteinmanSteady::ux( const Real& /* t */,
                               const Real& x, const Real& y,
                               const Real& z, const ID& i) {
    switch(i) {
        case 1:
            return
                a * b*exp(a*(x-z)+b*(y-z))-
                b * a*exp(a*(z-y)+b*(x-y));
            break;
        case 2:
            return
                (-a-b) * b*exp(a*(y-x)+b*(z-x))-
                 a * a*exp(a*(x-z)+b*(y-z));
            break;
        case 3:
            return
                b * b*exp(a*(z-y)+b*(x-y))-
                (-a-b) * a*exp(a*(y-x)+b*(z-x));
            break;
    }

    exit(1);
}

Real EthierSteinmanSteady::uy( const Real& /* t */,
                               const Real& x, const Real& y,
                               const Real& z, const ID& i) {
    switch(i) {
        case 1:
            return
                b * b*exp(a*(x-z)+b*(y-z))-
                (-a-b) * a*exp(a*(z-y)+b*(x-y));
            break;
        case 2:
            return
                a * b*exp(a*(y-x)+b*(z-x))-
                b * a*exp(a*(x-z)+b*(y-z));
            break;
        case 3:
            return
                (-a-b) * b*exp(a*(z-y)+b*(x-y))-
                a * a*exp(a*(y-x)+b*(z-x));
            break;
    }

    exit(1);
}

Real EthierSteinmanSteady::uz( const Real& /* t */,
                               const Real& x, const Real& y,
                               const Real& z, const ID& i) {
    switch(i) {
        case 1:
            return
                (-a-b) * b*exp(a*(x-z)+b*(y-z))-
                a * a*exp(a*(z-y)+b*(x-y));
            break;
        case 2:
            return
                b * b*exp(a*(y-x)+b*(z-x))-
                (-a-b) * a*exp(a*(x-z)+b*(y-z));
            break;
        case 3:
            return
                a * b*exp(a*(z-y)+b*(x-y))-
                b * a*exp(a*(y-x)+b*(z-x));
            break;
    }

    exit(1);
}


Real EthierSteinmanSteady::pexact(const Real& /* t */,
                                  const Real& x, const Real& y,
                                  const Real& z, const ID& /* i */) {
    return (a*a+b*b+a*b)*(exp(a*(x-y)+b*(x-z))+
                          exp(a*(y-z)+b*(y-x))+
                          exp(a*(z-x)+b*(z-y)));

}

Real EthierSteinmanSteady::xexact(const Real& t, const Real& x, const Real& y,
                                  const Real& z, const ID& i) {

    switch(i) {
        case 1:
        case 2:
        case 3:
            return uexact(t, x, y, z, i);
            break;
        case 4:
            return pexact(t, x, y, z, 1);
            break;
        default:
            exit(1);
    }
}


// Initial velocity
Real EthierSteinmanSteady::x0(const Real& t, const Real& x, const Real& y,
                              const Real& z, const ID& i)
{
    // return 0.0;
    return xexact(t,x,y,z,i);
}

Real EthierSteinmanSteady::fNeumann(const Real& t, const Real& x,
                                    const Real& y,
                                    const Real& z, const ID& i)
{
    Real nx=0.;
    Real ny=0.;
    Real nz=-1.;
    switch(i) {
        case 1:
            return pexact(t, x, y, z, 1) * nx
                - mu * ( ux(t, x, y, z, 1) * nx * 2 +
                         ux(t, x, y, z, 2) * ny +
                         ux(t, x, y, z, 3) * nz +
                         uy(t, x, y, z, 1) * ny +
                         uz(t, x, y, z, 1) * nz );
        case 2:
            return pexact(t, x, y, z, 1) * ny
                - mu * ( uy(t, x, y, z, 1) * nx +
                         uy(t, x, y, z, 2) * ny * 2 +
                         uy(t, x, y, z, 3) * nz +
                         ux(t, x, y, z, 2) * nx +
                         uz(t, x, y, z, 2) * nz );
        case 3:
            return pexact(t, x, y, z, 1) * nz
                - mu * ( uz(t, x, y, z, 1) * nx +
                         uz(t, x, y, z, 2) * ny +
                         uz(t, x, y, z, 3) * nz * 2 +
                         ux(t, x, y, z, 3) * nx +
                         uy(t, x, y, z, 3) * ny);
        default:
            exit(1);
    }
}

void EthierSteinmanSteady::setParamsFromGetPot( const GetPot& dataFile )
{
    a = dataFile( "fluid/problem/a", 0.75 );
    b = dataFile( "fluid/problem/b", 0.75 );
    mu = dataFile( "fluid/physics/viscosity", 0.0001 );
    nu = mu / dataFile( "fluid/physics/density", 1. );
}

Real EthierSteinmanSteady::nu;
Real EthierSteinmanSteady::mu;
Real EthierSteinmanSteady::sigma;
Real EthierSteinmanSteady::a;
Real EthierSteinmanSteady::b;

} // namespace LifeV
