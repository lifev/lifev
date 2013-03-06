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

#include <lifev/navier_stokes/function/RossEthierSteinmanInc.hpp>

namespace LifeV
{


Real RossEthierSteinmanUnsteadyInc::uexact ( const Real& t,
                                             const Real& x,
                                             const Real& y,
                                             const Real& z,
                                             const ID& i)
{
    Real e = std::exp (2.* (S_a * S_a + S_d * S_d + S_a * S_d) * S_nu * t);
    switch (i)
    {
        case 0:
            return e * ( S_d * std::exp (S_a * (x - z) + S_d * (y - z) ) - S_a * std::exp (S_a * (z - y) + S_d * (x - y) ) );
        case 1:
            return e * ( S_d * std::exp (S_a * (y - x) + S_d * (z - x) ) - S_a * std::exp (S_a * (x - z) + S_d * (y - z) ) );
        case 2:
            return e * ( S_d * std::exp (S_a * (z - y) + S_d * (x - y) ) - S_a * std::exp (S_a * (y - x) + S_d * (z - x) ) );
        default:
            exit (1);
    }
    return 1.;
}

Real RossEthierSteinmanUnsteadyInc::uderexact ( const Real& t,
                                                const Real& x,
                                                const Real& y,
                                                const Real& z,
                                                const ID& i)
{

    if (i < 3)
    {
        return 2.* (S_a * S_a + S_d * S_d + S_a * S_d) * S_nu * xexact (t, x, y, z, i);
    }
    else
    {
        return 0.;
    }
}

Real RossEthierSteinmanUnsteadyInc::pexact ( const Real& t,
                                             const Real& x,
                                             const Real& y,
                                             const Real& z,
                                             const ID& /* i */ )
{
    return S_rho * (S_a * S_a + S_d * S_d + S_a * S_d) * std::exp (4.* (S_a * S_a + S_d * S_d + S_a * S_d) * S_nu * t) *
           (std::exp (S_a * (x - y) + S_d * (x - z) ) + std::exp (S_a * (y - z) + S_d * (y - x) ) + std::exp (S_a * (z - x) + S_d * (z - y) )  );
}

Real RossEthierSteinmanUnsteadyInc::grad_u ( const UInt& icoor, const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    Real e = std::exp (2.* (S_a * S_a + S_d * S_d + S_a * S_d) * S_nu * t);
    switch (icoor)
    {
        case 0:    // u_x
            switch (i)
            {
                case 0:
                    return S_a * S_d * e * ( std::exp (S_a * (x - z) + S_d * (y - z) ) - std::exp (S_a * (z - y) + S_d * (x - y) ) );
                case 1:
                    return e * ( - (S_a + S_d) * S_d * std::exp (S_a * (y - x) + S_d * (z - x) ) - S_a * S_a * std::exp (S_a * (x - z) + S_d * (y - z) ) );
                case 2:
                    return e * ( S_d * S_d * std::exp (S_a * (z - y) + S_d * (x - y) ) + S_a * (S_a + S_d) * std::exp (S_a * (y - x) + S_d * (z - x) ) );
                default:
                    exit (1);
            }
        case 1:   // u_y
            switch (i)
            {
                case 0:
                    return e * ( S_d * S_d * std::exp (S_a * (x - z) + S_d * (y - z) ) + S_a * (S_a + S_d) * std::exp (S_a * (z - y) + S_d * (x - y) ) );
                case 1:
                    return e * ( S_a * S_d * std::exp (S_a * (y - x) + S_d * (z - x) ) - S_a * S_d * std::exp (S_a * (x - z) + S_d * (y - z) ) );
                case 2:
                    return e * ( -S_d * (S_a + S_d) * std::exp (S_a * (z - y) + S_d * (x - y) ) - S_a * S_a * std::exp (S_a * (y - x) + S_d * (z - x) ) );
                default:
                    exit (1);
            }
        case 2:
            switch (i)
            {
                case 0:
                    return e * ( -S_d * (S_a + S_d) * std::exp (S_a * (x - z) + S_d * (y - z) ) - S_a * S_a * std::exp (S_a * (z - y) + S_d * (x - y) ) );
                case 1:
                    return e * ( S_d * S_d * std::exp (S_a * (y - x) + S_d * (z - x) ) + S_a * (S_a + S_d) * std::exp (S_a * (x - z) + S_d * (y - z) ) );
                case 2:
                    return e * ( S_a * S_d * std::exp (S_a * (z - y) + S_d * (x - y) ) - S_a * S_d * std::exp (S_a * (y - x) + S_d * (z - x) ) );
                default:
                    exit (1);
            }
        default:
            exit (1);
    }
    return 1.;
}

Real RossEthierSteinmanUnsteadyInc::f ( const Real& /* t */,
                                        const Real& /* x */,
                                        const Real& /* y */,
                                        const Real& /* z */,
                                        const ID& /* i */ )
{
    return 0.;
}

Real RossEthierSteinmanUnsteadyInc::xexact ( const Real& t,
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
Real RossEthierSteinmanUnsteadyInc::x0 ( const Real& t, const Real& x, const Real& y,
                                         const Real& z, const ID& i )
{
    return xexact (t, x, y, z, i);
}



//we suppose that the problem geometry is the cube [-1,1]x[-1,1]x[-1,1].
Real RossEthierSteinmanUnsteadyInc::fNeumann ( const Real& t,
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


void RossEthierSteinmanUnsteadyInc::setParamsFromGetPot ( const GetPot& dataFile )
{
    S_a = dataFile ( "fluid/problem/a", 1. ) ;
    S_d = dataFile ( "fluid/problem/d", 1. ) ;
    S_mu = dataFile ( "fluid/physics/viscosity", 1. );
    S_rho =  dataFile ( "fluid/physics/density", 1. );
    S_nu = S_mu / S_rho;
    S_flagStrain = dataFile ( "fluid/physics/flag_strain", 0 );
}
void RossEthierSteinmanUnsteadyInc::setA (const Real& aValue)
{
    S_a = aValue;
}
void RossEthierSteinmanUnsteadyInc::setD (const Real& dValue)
{
    S_d = dValue;
}
void RossEthierSteinmanUnsteadyInc::setViscosity (const Real& mu)
{
    S_mu = mu;
    S_nu = S_mu / S_rho;
}
void RossEthierSteinmanUnsteadyInc::setDensity (const Real& rho)
{
    S_rho = rho;
    S_nu = S_mu / S_rho;
}
void RossEthierSteinmanUnsteadyInc::setFlagStrain (const Int& flagValue)
{
    S_flagStrain = flagValue;
}

Real RossEthierSteinmanUnsteadyInc::S_a;
Real RossEthierSteinmanUnsteadyInc::S_d;
Real RossEthierSteinmanUnsteadyInc::S_mu;
Real RossEthierSteinmanUnsteadyInc::S_rho;
Real RossEthierSteinmanUnsteadyInc::S_nu;
Int  RossEthierSteinmanUnsteadyInc::S_flagStrain;

} // namespace LifeV
