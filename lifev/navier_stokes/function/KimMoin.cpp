/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Mauro Perego <mauro@mathcs.emory.edu>
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
   \file KimMoin.cpp
   \author Mauro Perego <mauro@mathcs.emory.edu>
   \date 2009-10-02
*/
#include <lifev/core/LifeV.hpp>

#include <lifev/navier_stokes/function/KimMoin.hpp>

namespace LifeV
{
const Real Pi = 3.14159265358979323846264338328;

Real KimMoin::uexact( const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
    Real e = std::exp(-8*Pi*Pi*a*a*nu*t);
    switch(i) {
        case 0:  //u_1
            return -std::cos(2*Pi*a*x)*std::sin(2*Pi*a*y)*e;
        case 1:  //u_2
            return std::sin(2*Pi*a*x)*std::cos(2*Pi*a*y)*e;
        default:
            exit(1);
    }
   //dummy return as required by IBM xlC compiler
   return 1.;
}

Real KimMoin::pexact( const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& /* i */ )
{
    Real e2 = rho * std::exp(-16*a*a*Pi*Pi*nu*t);
    return -(std::cos(4*Pi*a*x)+std::cos(4*Pi*a*y))*e2/4;
}

Real KimMoin::uderexact( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return -8*Pi*Pi*a*a*nu*uexact(t, x, y, z, i);
}

Real KimMoin::grad_u( const UInt& icoor, const Real& t, const Real& x, const Real& y, const Real& /*z*/, const ID& i )
{
    Real e = std::exp(-8*Pi*Pi*a*a*nu*t);
    switch(icoor) {
        case 0:    // u_x
            switch(i) {
                    case 0:
                        return 2*Pi*e*a*std::sin(2*Pi*a*x)*std::sin(2*Pi*a*y);
                    case 1:
                        return 2*Pi*e*a*std::cos(2*Pi*a*x)*std::cos(2*Pi*a*y);
                    default:
                        exit(1);
             }
        case 1:   // u_y
            switch(i) {
                    case 0:
                        return -2*Pi*e*a*std::cos(2*Pi*a*x)*std::cos(2*Pi*a*y);
                    case 1:
                        return -2*Pi*e*a*std::sin(2*Pi*a*x)*std::sin(2*Pi*a*y);
                    default:
                        exit(1);
            }
        default:
            exit(1);
    }
   //dummy return as required by IBM xlC compiler
   return 1.;
}

Real KimMoin::f( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& /* i */ ) { return 0; }

Real KimMoin::xexact( const Real& t,
                    const Real& x,
                    const Real& y,
                    const Real& z,
                    const ID& i )
{
    switch(i) {
        case 0:  //u_1
        case 1:  //u_2
            return uexact(t, x, y, z, i);
        case 2:  //pressure
            return pexact(t, x, y, z, 1);
            break;
        default:
            exit(1);
    }
   //dummy return as required by IBM xlC compiler
   return 1.;
}




Real KimMoin::x0( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    return xexact(t,x,y,z,i);
}

// Initial velocity
Real KimMoin::u0( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    return uexact(t,x,y,z,i);
}

// Initial pressure
Real KimMoin::p0( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/ )
{
    return pexact(t,x,y,z,1);
}




//we suppose that the problem geometry is the square [0,1]x[0,1]
Real KimMoin::fNeumann( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    Real n[2] = {0, 0}; Real out=0;
    if        ( x == 0. ) {
        n[0] = -1.;
    } else if ( x ==  1. ) {
        n[0] =  1.;
    } else if ( y == 0. ) {
        n[1] = -1.;
    } else if ( y ==  1. ) {
        n[1] =  1.;
    } else {
        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                  << std::endl;
    }

    for (UInt k =0; k< 2; k++)  //mu gradu n
        out += mu* grad_u(k, t, x, y, z, i)*n[k];

    if(flag_strain)
        for (UInt k =0; k< 2; k++)  //mu gradu^T n
            out += mu* grad_u(i, t, x, y, z, k)*n[k];

    out -= pexact(t, x, y, z, i) * n[i];


    return out;
}

//we suppose that the problem geometry is the square [0,1]x[0,1]
Real KimMoin::normalVector( const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i )
{
    Real n[2] = {0, 0};
    if        ( x == 0. ) {
        n[0] = -1.;
    } else if ( x ==  1. ) {
        n[0] =  1.;
    } else if ( y == 0. ) {
        n[1] = -1.;
    } else if ( y ==  1. ) {
        n[1] =  1.;
    } else {
//        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
//                  << std::endl;
    }

    return n[i];
}

Real KimMoin::fShearStress( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    Real out=0;

    for (UInt k =0; k< nDimensions; k++)  //mu gradu n
        out += mu* grad_u(k, t, x, y, z, i)*normalVector( t, x, y, z, k );

    if(flag_strain)
        for (UInt k =0; k< nDimensions; k++)  //mu gradu^T n
            out += mu* grad_u(i, t, x, y, z, k)*normalVector( t, x, y, z, k );

    return  out;
}

Real KimMoin::fWallShearStress( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    // wss = s_n - ( s_n \cdot n ) n
    Real wss=0;
    Real s_n_n(0.);
    // this is the i-component of the normal stress
    wss = fShearStress(t, x, y, z, i);

    for (UInt k =0; k< nDimensions; k++) // ( s_n \cdot n )
        s_n_n += fShearStress(t, x, y, z, k) * normalVector( t, x, y, z, k );

    wss -= s_n_n * normalVector( t, x, y, z, i );

    return  wss;
}

void KimMoin::setParamsFromGetPot( const GetPot& dataFile )
{
    mu = dataFile( "fluid/physics/viscosity", 1. );
    rho = dataFile( "fluid/physics/density", 1. );
    flag_strain = dataFile( "fluid/space_discretization/stiff_strain", false );
    nu = mu / rho;
    a = dataFile( "fluid/problem/a", 1. );
}

Real KimMoin::nu;
Real KimMoin::rho;
Real KimMoin::mu;
Real KimMoin::a;
bool  KimMoin::flag_strain;

} // namespace LifeV
