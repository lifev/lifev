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
    @brief Fully 3D solution for the unsteady Navier-Stokes equations

    @author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 12-11-2004

    see [Ethier and Steinman] for more details
 */

#include "EthierSteinmanUnsteady.hpp"

namespace LifeV
{

Real EthierSteinmanUnsteady::f( const Real& /* t */,
                                const Real& /* x */,
                                const Real& /* y */,
                                const Real& /* z */,
                                const ID& /* i */ ) const { return 0; }

Real EthierSteinmanUnsteady::xexact( const Real& t,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const ID& i ) const
{
    Real tau = nu * t;

    switch (i)
    {
    case 0:
        return -a * exp(-d*d*tau) *
               ( exp(a*x) * sin(a*y+d*z) +
                 exp(a*z) * cos(a*x+d*y) );
        break;
    case 1:
        return -a * exp(-d*d*tau) *
               ( exp(a*y) * sin(a*z+d*x) +
                 exp(a*x) * cos(a*y+d*z) );
        break;
    case 2:
        return -a * exp(-d*d*tau) *
               ( exp(a*z) * sin(a*x+d*y) +
                 exp(a*y) * cos(a*z+d*x) );
        break;
    case 3:
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
                                     const ID& i) const
{
    if (i < 3)
        return xexact(t, x, y, z, i);
    else
        return 0.;
}

Real EthierSteinmanUnsteady::pexact( const Real& t,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const ID& /* i */ ) const
{
    return xexact(t, x, y, z, 3);
}

// Initial velocity
Real EthierSteinmanUnsteady::x0( const Real& /* t */, const Real& x, const Real& y,
                                 const Real& z, const ID& i ) const
{
    return xexact(0.,x, y, z, i);
}

// Derivative of u with respect to the time
Real EthierSteinmanUnsteady::uderexact( const Real& t,
                                        const Real& x,
                                        const Real& y,
                                        const Real& z,
                                        const ID& i) const
{

    if (i < 3)
        return - nu*d*d*xexact(t, x, y, z, i);
    else
        return 0.;
}


EthierSteinmanUnsteady::EthierSteinmanUnsteady()
:a(1.), d(1.), mu(1.), density(1.), nu(mu/density)
{

}

EthierSteinmanUnsteady::EthierSteinmanUnsteady(const EthierSteinmanUnsteady& ESU)
:a(ESU.a), d(ESU.d), mu(ESU.mu), density(ESU.density), nu(mu/density)
{

}

EthierSteinmanUnsteady::~EthierSteinmanUnsteady()
{

}



// derivatives for neumann
Real EthierSteinmanUnsteady::ux( const Real& t, const Real& x, const Real& y,
                                 const Real& z, const ID& i ) const
{
    Real tau = nu * t;

    switch (i)
    {
    case 0:
        return -a * exp(-d*d*tau) *
               ( a * exp(a*x) * sin(a*y+d*z) -
                 a * exp(a*z) * sin(a*x+d*y) );
        break;
    case 1:
        return -a * exp(-d*d*tau) *
               ( d * exp(a*y) * cos(a*z+d*x) +
                 a * exp(a*x) * cos(a*y+d*z) );
        break;
    case 2:
        return -a * exp(-d*d*tau) *
               ( a * exp(a*z) * cos(a*x+d*y) -
                 d * exp(a*y) * sin(a*z+d*x) );
        break;
    default:
        exit(1);
    }
}

Real EthierSteinmanUnsteady::uy( const Real& t, const Real& x, const Real& y,
                                 const Real& z, const ID& i ) const
{
    Real tau = nu * t;

    switch (i)
    {
    case 0:
        return -a * exp(-d*d*tau) *
               ( a * exp(a*x) * cos(a*y+d*z) -
                 d * exp(a*z) * sin(a*x+d*y) );
        break;
    case 1:
        return -a * exp(-d*d*tau) *
               ( a * exp(a*y) * sin(a*z+d*x) -
                 a * exp(a*x) * sin(a*y+d*z) );
        break;
    case 2:
        return -a * exp(-d*d*tau) *
               ( d * exp(a*z) * cos(a*x+d*y) +
                 a * exp(a*y) * cos(a*z+d*x) );
        break;
    default:
        exit(1);
    }
}

Real EthierSteinmanUnsteady::uz( const Real& t, const Real& x, const Real& y,
                                 const Real& z, const ID& i ) const
{
    Real tau = nu * t;

    switch (i)
    {
    case 0:
        return -a * exp(-d*d*tau) *
               ( d * exp(a*x) * cos(a*y+d*z) +
                 a * exp(a*z) * cos(a*x+d*y) );
        break;
    case 1:
        return -a * exp(-d*d*tau) *
               ( a * exp(a*y) * cos(a*z+d*x) -
                 d * exp(a*x) * sin(a*y+d*z) );
        break;
    case 2:
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
                                       const ID& i ) const
{
    Real nx = 0.;
    Real ny = 0.;
    Real nz = 0.;
    if        ( x == -1. )
    {
        nx = -1.;
    }
    else if ( x ==  1. )
    {
        nx =  1.;
    }
    else if ( y == -1. )
    {
        ny = -1.;
    }
    else if ( y ==  1. )
    {
        ny =  1.;
    }
    else if ( z == -1. )
    {
        nz = -1.;
    }
    else if ( z ==  1. )
    {
        nz =  1.;
    }
    else
    {
        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                  << std::endl;
    }

    switch (i)
    {
    case 0:
        return - pexact(t, x, y, z, 0) * nx
               + mu * ( ux(t, x, y, z, 0) * nx * 2 +
                        ux(t, x, y, z, 1) * ny +
                        ux(t, x, y, z, 2) * nz +
                        uy(t, x, y, z, 0) * ny +
                        uz(t, x, y, z, 0) * nz );
    case 1:
        return - pexact(t, x, y, z, 0) * ny
               + mu * ( uy(t, x, y, z, 0) * nx +
                        uy(t, x, y, z, 1) * ny * 2 +
                        uy(t, x, y, z, 2) * nz +
                        ux(t, x, y, z, 1) * nx +
                        uz(t, x, y, z, 1) * nz );
    case 2:
        return - pexact(t, x, y, z, 0) * nz
               + mu * ( uz(t, x, y, z, 0) * nx +
                        uz(t, x, y, z, 1) * ny +
                        uz(t, x, y, z, 2) * nz * 2 +
                        ux(t, x, y, z, 2) * nx +
                        uy(t, x, y, z, 2) * ny);
    default:
        exit(1);
    }

    /*
    switch(i) {
        case 1:
            return - mu/nu*pexact(t, x, y, z, 0) * nx
                + mu * ( ux(t, x, y, z, 0) * nx +
                         uy(t, x, y, z, 0) * ny +
                         uz(t, x, y, z, 0) * nz );
        case 2:
            return - mu/nu*pexact(t, x, y, z, 0) * ny
                + mu * ( ux(t, x, y, z, 1) * nx +
                         uy(t, x, y, z, 1) * ny +
                         uz(t, x, y, z, 1) * nz );
        case 3:
            return - mu/nu*pexact(t, x, y, z, 0) * nz
                + mu * ( ux(t, x, y, z, 2) * nx +
                         uy(t, x, y, z, 2) * ny +
                         uz(t, x, y, z, 2) * nz);
        default:
            exit(1);
    }
    */
}

void EthierSteinmanUnsteady::setA(const Real& a)
{
    this->a = a;
}
void EthierSteinmanUnsteady::setD(const Real& d)
{
    this->d = d;
}
void EthierSteinmanUnsteady::setViscosity(const Real& mu)
{
    this->mu = mu;
    this->nu = mu/density;
}
void EthierSteinmanUnsteady::setDensity(const Real& density)
{
    this->density = density;
    this->nu = mu/density;
}

void setParamsFromGetPot( EthierSteinmanUnsteady& ESU, const GetPot& dataFile )
{
    ESU.setA(dataFile( "fluid/problem/a", 1. ));
    ESU.setD(dataFile( "fluid/problem/d", 1. ));
    ESU.setViscosity(dataFile( "fluid/physics/viscosity", 1. ));
    ESU.setDensity(dataFile( "fluid/physics/density", 1. ));
}

} // namespace LifeV
