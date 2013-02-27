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
    @file Ross-Ethier Steinmann Analytical Solution
    @brief A short description of the test content

    @author      Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @contributor Mauro Perego         <mauro@mathcs.emory.edu>
    @contributor Umberto Villa        <uvilla@emory.edu>
    @maintainer  Umberto Villa        <uvilla@emory.edu>

    @date 2004-11-12

    Implementation.
*/

#include <lifev/core/LifeV.hpp>

#include <lifev/navier_stokes/function/RossEthierSteinmanDec.hpp>

namespace LifeV
{


Real RossEthierSteinmanUnsteadyDec::uexact ( const Real& t,
                                             const Real& x,
                                             const Real& y,
                                             const Real& z,
                                             const ID& i)
{
    Real e = std::exp (-S_d * S_d * S_nu * t);

    switch (i)
    {
        case 0:
            return -S_a * e * ( std::exp (S_a * x) * std::sin (S_a * y + S_d * z) + std::exp (S_a * z) * std::cos (S_a * x + S_d * y) );
        case 1:
            return -S_a * e * ( std::exp (S_a * y) * std::sin (S_a * z + S_d * x) + std::exp (S_a * x) * std::cos (S_a * y + S_d * z) );
        case 2:
            return -S_a * e * ( std::exp (S_a * z) * std::sin (S_a * x + S_d * y) + std::exp (S_a * y) * std::cos (S_a * z + S_d * x) );
        default:
            exit (1);
    }
    return 1.;
}

Real RossEthierSteinmanUnsteadyDec::uderexact ( const Real& t,
                                                const Real& x,
                                                const Real& y,
                                                const Real& z,
                                                const ID& i)
{

    if (i < 3)
    {
        return -S_d * S_d * S_nu * xexact (t, x, y, z, i);
    }
    else
    {
        return 0.;
    }
}

Real RossEthierSteinmanUnsteadyDec::pexact ( const Real& t,
                                             const Real& x,
                                             const Real& y,
                                             const Real& z,
                                             const ID& /* i */ )
{
    return - S_rho * S_a * S_a / 2. * std::exp (-2.*S_d * S_d * S_nu * t) *
           ( std::exp (2.*S_a * x) + std::exp (2.*S_a * y) + std::exp (2.*S_a * z) +
             2. * std::sin (S_a * x + S_d * y) * std::cos (S_a * z + S_d * x) * std::exp (S_a * (y + z) ) +
             2. * std::sin (S_a * y + S_d * z) * std::cos (S_a * x + S_d * y) * std::exp (S_a * (z + x) ) +
             2. * std::sin (S_a * z + S_d * x) * std::cos (S_a * y + S_d * z) * std::exp (S_a * (x + y) ) );
}

Real RossEthierSteinmanUnsteadyDec::grad_u ( const UInt& icoor, const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    Real e = std::exp (-S_d * S_d * S_nu * t);
    switch (icoor)
    {
        case 0:    // u_x
            switch (i)
            {
                case 0:
                    return -S_a * e * ( S_a * std::exp (S_a * x) * std::sin (S_a * y + S_d * z) - S_a * std::exp (S_a * z) * std::sin (S_a * x + S_d * y) );
                case 1:
                    return -S_a * e * ( S_d * std::exp (S_a * y) * std::cos (S_a * z + S_d * x) + S_a * std::exp (S_a * x) * std::cos (S_a * y + S_d * z) );
                case 2:
                    return -S_a * e * ( S_a * std::exp (S_a * z) * std::cos (S_a * x + S_d * y) - S_d * std::exp (S_a * y) * std::sin (S_a * z + S_d * x) );
                default:
                    exit (1);
            }
        case 1:   // u_y
            switch (i)
            {
                case 0:
                    return -S_a * e * ( S_a * std::exp (S_a * x) * std::cos (S_a * y + S_d * z) - S_d * std::exp (S_a * z) * std::sin (S_a * x + S_d * y) );
                case 1:
                    return -S_a * e * ( S_a * std::exp (S_a * y) * std::sin (S_a * z + S_d * x) - S_a * std::exp (S_a * x) * std::sin (S_a * y + S_d * z) );
                case 2:
                    return -S_a * e * ( S_d * std::exp (S_a * z) * std::cos (S_a * x + S_d * y) + S_a * std::exp (S_a * y) * std::cos (S_a * z + S_d * x) );
                default:
                    exit (1);
            }
        case 2:
            switch (i)
            {
                case 0:
                    return -S_a * e * ( S_d * std::exp (S_a * x) * std::cos (S_a * y + S_d * z) +  S_a * std::exp (S_a * z) * std::cos (S_a * x + S_d * y) );
                case 1:
                    return -S_a * e * ( S_a * std::exp (S_a * y) * std::cos (S_a * z + S_d * x) - S_d * std::exp (S_a * x) * std::sin (S_a * y + S_d * z) );
                case 2:
                    return -S_a * e * ( S_a * std::exp (S_a * z) * std::sin (S_a * x + S_d * y) - S_a * std::exp (S_a * y) * std::sin (S_a * z + S_d * x) );
                default:
                    exit (1);
            }
        default:
            exit (1);
    }
    return 1.;
}

Real RossEthierSteinmanUnsteadyDec::f ( const Real& /* t */,
                                        const Real& /* x */,
                                        const Real& /* y */,
                                        const Real& /* z */,
                                        const ID& /* i */ )
{
    return 0.;
}

Real RossEthierSteinmanUnsteadyDec::xexact ( const Real& t,
                                             const Real& x,
                                             const Real& y,
                                             const Real& z,
                                             const ID& i )
{
    switch (i)
    {
        case 0:
        case 1:
        case 2:
            return uexact (t, x, y, z, i);
        case 3:
            return pexact (t, x, y, z, 1);
        default:
            exit (1);
    }
    return 1.;
}



// Initial velocity
Real RossEthierSteinmanUnsteadyDec::x0 ( const Real& t, const Real& x, const Real& y,
                                         const Real& z, const ID& i )
{
    return xexact (t, x, y, z, i);
}



//we suppose that the problem geometry is the cube [-1,1]x[-1,1]x[-1,1].
Real RossEthierSteinmanUnsteadyDec::fNeumann ( const Real& t,
                                               const Real& x,
                                               const Real& y,
                                               const Real& z,
                                               const ID& i )
{
    Real n[3] = {0., 0., 0.};
    Real out = 0.;
    if        ( x == -1. )
    {
        n[0] = -1.;
    }
    else if ( x ==  1. )
    {
        n[0] =  1.;
    }
    else if ( y == -1. )
    {
        n[1] = -1.;
    }
    else if ( y ==  1. )
    {
        n[1] =  1.;
    }
    else if ( z == -1. )
    {
        n[2] = -1.;
    }
    else if ( z ==  1. )
    {
        n[2] =  1.;
    }
    else
    {
        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                  << std::endl;
    }

    for (UInt k = 0; k < 3; k++) //mu grad_u n
    {
        out += S_mu * grad_u (k, t, x, y, z, i) * n[k];
    }

    if (S_flagStrain)
        for (UInt k = 0; k < 3; k++) //mu grad_u^T n
        {
            out += S_mu * grad_u (i, t, x, y, z, k) * n[k];
        }

    out -= pexact (t, x, y, z, i) * n[i]; //grad_p n
    return out;
}


//we suppose that the problem geometry is the cube [-1,1]x[-1,1]x[-1,1].
Real RossEthierSteinmanUnsteadyDec::normalVector ( const Real& /*t*/,
                                                   const Real& x,
                                                   const Real& y,
                                                   const Real& z,
                                                   const ID& i )
{
    Real n[3] = {0., 0., 0.};
    if        ( x == -1. )
    {
        n[0] = -1.;
    }
    else if ( x ==  1. )
    {
        n[0] =  1.;
    }
    else if ( y == -1. )
    {
        n[1] = -1.;
    }
    else if ( y ==  1. )
    {
        n[1] =  1.;
    }
    else if ( z == -1. )
    {
        n[2] = -1.;
    }
    else if ( z ==  1. )
    {
        n[2] =  1.;
    }
    else
    {
        std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                  << std::endl;
    }
    return n[i];

}


//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
Real RossEthierSteinmanUnsteadyDec::fShearStress ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    Real out = 0.;

    for (UInt k = 0; k < nDimensions; k++) //mu gradu n
    {
        out += S_mu * grad_u (k, t, x, y, z, i) * normalVector ( t, x, y, z, k );
    }

    if (S_flagStrain)
        for (UInt k = 0; k < nDimensions; k++) //mu gradu^T n
        {
            out += S_mu * grad_u (i, t, x, y, z, k) * normalVector ( t, x, y, z, k );
        }

    return  out;
}

//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
Real RossEthierSteinmanUnsteadyDec::fWallShearStress ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    // wss = s_n - ( s_n \cdot n ) n
    Real wss = 0.;
    Real s_n_n (0.);
    // this is the i-component of the normal stress
    wss = fShearStress (t, x, y, z, i);

    for (UInt k = 0; k < nDimensions; k++) // ( s_n \cdot n )
    {
        s_n_n += fShearStress (t, x, y, z, k) * normalVector ( t, x, y, z, k );
    }

    wss -= s_n_n * normalVector ( t, x, y, z, i );

    return  wss;
}

void RossEthierSteinmanUnsteadyDec::setParamsFromGetPot ( const GetPot& dataFile )
{
    S_a = dataFile ( "fluid/problem/a", 1. ) ;
    S_d = dataFile ( "fluid/problem/d", 1. ) ;
    S_mu = dataFile ( "fluid/physics/viscosity", 1. );
    S_rho = dataFile ( "fluid/physics/density", 1. );
    S_nu = S_mu / S_rho;
    S_flagStrain = dataFile ( "fluid/physics/flag_strain", 0 );
}
void RossEthierSteinmanUnsteadyDec::setA (const Real& aValue)
{
    S_a = aValue;
}
void RossEthierSteinmanUnsteadyDec::setD (const Real& dValue)
{
    S_d = dValue;
}
void RossEthierSteinmanUnsteadyDec::setViscosity (const Real& mu)
{
    S_mu = mu;
    S_nu = S_mu / S_rho;
}
void RossEthierSteinmanUnsteadyDec::setDensity (const Real& rho)
{
    S_rho = rho;
    S_nu = S_mu / S_rho;
}
void RossEthierSteinmanUnsteadyDec::setFlagStrain (const Int& flagValue)
{
    S_flagStrain = flagValue;
}

Real RossEthierSteinmanUnsteadyDec::S_a;
Real RossEthierSteinmanUnsteadyDec::S_d;
Real RossEthierSteinmanUnsteadyDec::S_mu;
Real RossEthierSteinmanUnsteadyDec::S_rho;
Real RossEthierSteinmanUnsteadyDec::S_nu;
int  RossEthierSteinmanUnsteadyDec::S_flagStrain;

} // namespace LifeV
