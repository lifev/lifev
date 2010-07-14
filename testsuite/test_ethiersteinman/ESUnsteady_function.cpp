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
   \file ethierSteinman.cpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-11-12
*/

#include <ESUnsteady_function.hpp>

namespace LifeV
{

Real EthierSteinmanUnsteady::f( const Real& /* t */,
                                const Real& /* x */,
                                const Real& /* y */,
                                const Real& /* z */,
                                const ID& /* i */ ) { return 0; }

Real EthierSteinmanUnsteady::xexact( const Real& t,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const ID& i )
{
    Real tau = nu * t;

    switch(i) {
        case 1:
            return -a * exp(-d*d*tau) *
                ( exp(a*x) * sin(a*y+d*z) +
                  exp(a*z) * cos(a*x+d*y) );
            break;
        case 2:
            return -a * exp(-d*d*tau) *
                ( exp(a*y) * sin(a*z+d*x) +
                  exp(a*x) * cos(a*y+d*z) );
            break;
        case 3:
            return -a * exp(-d*d*tau) *
                ( exp(a*z) * sin(a*x+d*y) +
                  exp(a*y) * cos(a*z+d*x) );
            break;
        case 4:
            return -a*a / 2 * exp(-2*d*d*tau) *
                ( exp(2*a*x) + exp(2*a*y) + exp(2*a*z) +
                  2 * sin(a*x+d*y) * cos(a*z+d*x) * exp(a*(y+z)) +
                  2 * sin(a*y+d*z) * cos(a*x+d*y) * exp(a*(z+x)) +
                  2 * sin(a*z+d*x) * cos(a*y+d*z) * exp(a*(x+y)) );
            break;
        default:
            exit(1);
    }
}

Real EthierSteinmanUnsteady::uexact( const Real& t,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const ID& i)
{
    if (i < 4)
        return xexact(t, x, y, z, i);
    else
        return 0.;
}

Real EthierSteinmanUnsteady::pexact( const Real& t,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const ID& /* i */ )
{
    return xexact(t, x, y, z, 4);
}

// Initial velocity
Real EthierSteinmanUnsteady::x0( const Real& /* t */, const Real& x, const Real& y,
                                 const Real& z, const ID& i )
{
    return xexact(0.,x, y, z, i);
}

// Derivative of u with respect to the time
Real EthierSteinmanUnsteady::uderexact( const Real& t,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const ID& i)
{

    if (i < 4)
        return - nu*d*d*xexact(t, x, y, z, i);
    else
        return 0.;
}






// derivatives for neumann
Real EthierSteinmanUnsteady::ux( const Real& t, const Real& x, const Real& y,
                                 const Real& z, const ID& i )
{
    Real tau = nu * t;

    switch(i) {
        case 1:
            return -a * exp(-d*d*tau) *
                ( a * exp(a*x) * sin(a*y+d*z) -
                  a * exp(a*z) * sin(a*x+d*y) );
            break;
        case 2:
            return -a * exp(-d*d*tau) *
                ( d * exp(a*y) * cos(a*z+d*x) +
                  a * exp(a*x) * cos(a*y+d*z) );
            break;
        case 3:
            return -a * exp(-d*d*tau) *
                ( a * exp(a*z) * cos(a*x+d*y) -
                  d * exp(a*y) * sin(a*z+d*x) );
            break;
        default:
            exit(1);
    }
}

Real EthierSteinmanUnsteady::uy( const Real& t, const Real& x, const Real& y,
                                 const Real& z, const ID& i )
{
    Real tau = nu * t;

    switch(i) {
        case 1:
            return -a * exp(-d*d*tau) *
                ( a * exp(a*x) * cos(a*y+d*z) -
                  d * exp(a*z) * sin(a*x+d*y) );
            break;
        case 2:
            return -a * exp(-d*d*tau) *
                ( a * exp(a*y) * sin(a*z+d*x) -
                  a * exp(a*x) * sin(a*y+d*z) );
            break;
        case 3:
            return -a * exp(-d*d*tau) *
                ( d * exp(a*z) * cos(a*x+d*y) +
                  a * exp(a*y) * cos(a*z+d*x) );
            break;
        default:
            exit(1);
    }
}

Real EthierSteinmanUnsteady::uz( const Real& t, const Real& x, const Real& y,
                                 const Real& z, const ID& i )
{
    Real tau = nu * t;

    switch(i) {
        case 1:
            return -a * exp(-d*d*tau) *
                ( d * exp(a*x) * cos(a*y+d*z) +
                  a * exp(a*z) * cos(a*x+d*y) );
            break;
        case 2:
            return -a * exp(-d*d*tau) *
                ( a * exp(a*y) * cos(a*z+d*x) -
                  d * exp(a*x) * sin(a*y+d*z) );
            break;
        case 3:
            return -a * exp(-d*d*tau) *
                ( a * exp(a*z) * sin(a*x+d*y) -
                  a * exp(a*y) * sin(a*z+d*x) );
            break;
        default:
            exit(1);
    }
}

Real EthierSteinmanUnsteady::fNeumann( const Real& t,
                                       const Real& x,
                                       const Real& y,
                                       const Real& z,
                                       const ID& i )
{
    Real nx = 0.;
    Real ny = 0.;
    Real nz = 0.;
    if        ( x == -1. ) {
        nx = -1.;
    } else if ( x ==  1. ) {
        nx =  1.;
    } else if ( y == -1. ) {
        ny = -1.;
    } else if ( y ==  1. ) {
        ny =  1.;
    } else if ( z == -1. ) {
        nz = -1.;
    } else if ( z ==  1. ) {
        nz =  1.;
    } else {
        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                  << std::endl;
    }

    switch(i) {
        case 1:
            return - pexact(t, x, y, z, 1) * nx
                + mu * ( ux(t, x, y, z, 1) * nx * 2 +
                         ux(t, x, y, z, 2) * ny +
                         ux(t, x, y, z, 3) * nz +
                         uy(t, x, y, z, 1) * ny +
                         uz(t, x, y, z, 1) * nz );
        case 2:
            return - pexact(t, x, y, z, 1) * ny
                + mu * ( uy(t, x, y, z, 1) * nx +
                         uy(t, x, y, z, 2) * ny * 2 +
                         uy(t, x, y, z, 3) * nz +
                         ux(t, x, y, z, 2) * nx +
                         uz(t, x, y, z, 2) * nz );
        case 3:
            return - pexact(t, x, y, z, 1) * nz
                + mu * ( uz(t, x, y, z, 1) * nx +
                         uz(t, x, y, z, 2) * ny +
                         uz(t, x, y, z, 3) * nz * 2 +
                         ux(t, x, y, z, 3) * nx +
                         uy(t, x, y, z, 3) * ny);
        default:
            exit(1);
    }

    /*
    switch(i) {
        case 1:
            return - mu/nu*pexact(t, x, y, z, 1) * nx
                + mu * ( ux(t, x, y, z, 1) * nx +
                         uy(t, x, y, z, 1) * ny +
                         uz(t, x, y, z, 1) * nz );
        case 2:
            return - mu/nu*pexact(t, x, y, z, 1) * ny
                + mu * ( ux(t, x, y, z, 2) * nx +
                         uy(t, x, y, z, 2) * ny +
                         uz(t, x, y, z, 2) * nz );
        case 3:
            return - mu/nu*pexact(t, x, y, z, 1) * nz
                + mu * ( ux(t, x, y, z, 3) * nx +
                         uy(t, x, y, z, 3) * ny +
                         uz(t, x, y, z, 3) * nz);
        default:
            exit(1);
    }
    */
}

void EthierSteinmanUnsteady::setParamsFromGetPot( const GetPot& dataFile )
{
    a = dataFile( "fluid/problem/a", 1. );
    d = dataFile( "fluid/problem/d", 1. );
    mu = dataFile( "fluid/physics/viscosity", 1. );
    nu = mu / dataFile( "fluid/physics/density", 1. );
}

Real EthierSteinmanUnsteady::a;
Real EthierSteinmanUnsteady::d;
Real EthierSteinmanUnsteady::nu;
Real EthierSteinmanUnsteady::mu;

} // namespace LifeV
