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

Real EthierSteinmanSteady::uexact(const Real& t, const Real& x, const Real& y,
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


Real EthierSteinmanSteady::pexact(const Real& t, const Real& x, const Real& y,
                                  const Real& z, const ID& i) {
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
Real EthierSteinmanSteady::u0(const Real& t, const Real& x, const Real& y,
                              const Real& z, const ID& i)
{
    // return 0.0;
    return uexact(t,x,y,z,i);
}

Real EthierSteinmanSteady::fNeumann(const Real& t, const Real& x,
                                    const Real& y,
                                    const Real& z, const ID& i)
{
    return 0.0;
}

void EthierSteinmanSteady::setParamsFromGetPot( const GetPot& dataFile )
{
    a = dataFile( "fluid/problem/a", 0.75 );
    b = dataFile( "fluid/problem/b", 0.75 );
    nu =
        dataFile( "fluid/physics/viscosity", 0.0001 ) /
        dataFile( "fluid/physics/density", 1. );
}

Real EthierSteinmanSteady::nu;
Real EthierSteinmanSteady::sigma;
Real EthierSteinmanSteady::a;
Real EthierSteinmanSteady::b;

} // namespace LifeV
