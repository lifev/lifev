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
    @brief Womersley Analytical Solution

    @author      Mauro Perego  <mauro@mathcs.emory.edu>
    @contributor Umberto Villa <uvilla@emory.edu>
    @maintainer  Umberto Villa <uvilla@emory.edu>

    @date 11-12-2009

    Implementation
 */

#include <lifeconfig.h>
#include <life/lifefunctions/Womersley.hpp>
#include <life/lifefunctions/bessel/bessel.hpp>

namespace LifeV
{

Real Womersley::uexact( const Real& t, const Real& /*x*/, const Real& y, const Real& z, const ID& i)
{
	Real r=sqrt(z*z+y*y);
	std::complex<Real> z2, b2;
	z2 = 2.*r/S_D*S_z1;
	cbessjy01(z2, b2, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);
	Real u = real(S_A/S_L/S_rho/S_wi*(1.-b2/S_b1)*exp(S_wi*t));
    switch(i) {
        case 1:  //u_1
            return u;//-4*x*y*y; //u-4*x*y*y;
        case 2:  //u_2
            return 0;//y*y*y;
        case 3:
        	return 0;
        default:
            exit(1);
    }
}

Real Womersley::pexact( const Real& t, const Real& x, const Real& /*y*/, const Real& /*z*/, const ID& /* i */ )
{
    return S_A/S_L*(S_L-x)*cos(S_w*t);
}

Real Womersley::grad_u( const UInt& icoor, const Real& t, const Real& /*x*/, const Real& y, const Real& z, const ID& i )
{
	Real r=sqrt(y*y+z*z);
	std::complex<Real> z2, b2;
	z2 = 2.*r/S_D*S_z1;
	cbessjy01(z2, b2, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);
	b2 = -2./S_D*S_z1*S_cj0p;
	Real u_r = real(S_A/S_L/S_rho/S_wi*+b2/S_b1*exp(S_wi*t));

	switch(icoor) {
		case 1:
			switch(i) {
			case 1:
				return 0;
			case 2:
				return 0;
			case 3:
				return 0;
			default:
				exit(1);
			 }
		case 2:   // u_y
			switch(i) {
	        case 1:
	            return u_r/r*y;
	        case 2:
	            return 0;
	        case 3:
	        	return 0;
	        default:
	            exit(1);
			}
		case 3:
			switch(i) {
	        case 1:
	            return u_r/r*z;
	        case 2:
	            return 0;
	        case 3:
	        	return 0;
	        default:
	            exit(1);
			}
		default:
			exit(1);
	}
}

Real Womersley::f( const Real& /* t */, const Real&  /*x*/ , const Real&  /*y*/ , const Real& /* z */, const ID& i ) {
	switch(i) {
	        case 1:
	            return 0;
	        case 2:
	            return 0;
	        case 3:
	        	return 0;
	        default:
	            exit(1);
	    }
	return 0; }

Real Womersley::xexact( const Real& t,
					const Real& x,
                    const Real& y,
                    const Real& z,
                    const ID& i )
{
    switch(i) {
        case 1:  //u_1
        case 2:  //u_2
        	return uexact(t, x, y, z, i);
        case 3:  //pressure
            return pexact(t, x, y, z, 1);
            break;
        default:
            exit(1);
    }
}




Real Womersley::x0( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    return xexact(t,x,y,z,i);
}

// Initial velocity
Real Womersley::u0( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    return uexact(t,x,y,z,i);
}

// Initial pressure
Real Womersley::p0( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/ )
{
    return pexact(t,x,y,z,1);
}

//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
Real Womersley::fNeumann( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
	Real r=sqrt(y*y+z*z);
    Real n[3] = {0, 0, 0}; Real out=0;
    if        ( x < 1e-6/S_L ) {
        n[0] = -1.;
    } else if ( x >  S_L*(1-1e-6) ) {
        n[0] =  1.;
    } else if ( r > S_D/2.*0.9 ) {
        n[1] = y/r;
        n[2] = z/r;
    } else {
        // std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }

    for (UInt k =0; k< nDimensions; k++)  //mu gradu n
        out += S_mu* grad_u(k+1, t, x, y, z, i)*n[k];

    if(S_flagStrain)
    	for (UInt k =0; k< nDimensions; k++)  //mu gradu^T n
    		out += S_mu* grad_u(i, t, x, y, z, k+1)*n[k];

    out -= pexact(t, x, y, z, i) * n[i-1];


    return  out;
}

Real Womersley::normalVector( const Real& /*t*/, const Real& x, const Real& y,
        const Real& z, const ID& i )
{
    Real r=sqrt(y*y+z*z);
    Real n[3] = {0, 0, 0};
    if        ( x < 1e-6/S_L ) {
        n[0] = -1.;
    } else if ( x >  S_L*(1-1e-6) ) {
        n[0] =  1.;
    } else if ( r > S_D/2.*0.9 ) {
        n[1] = y/r;
        n[2] = z/r;
    } else {
        // std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }
    return n[i-1];
};

//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
Real Womersley::fShearStress( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    Real r=sqrt(y*y+z*z);
    Real n[3] = {0, 0, 0}; Real out=0;
    if        ( x < 1e-6/S_L ) {
        n[0] = -1.;
    } else if ( x >  S_L*(1-1e-6) ) {
        n[0] =  1.;
    } else if ( r > S_D/2.*0.9 ) {
        n[1] = y/r;
        n[2] = z/r;
    } else {
        // std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }

    for (UInt k =0; k< nDimensions; k++)  //mu gradu n
        out += S_mu* grad_u(k+1, t, x, y, z, i)*n[k];

    if(S_flagStrain)
        for (UInt k =0; k< nDimensions; k++)  //mu gradu^T n
            out += S_mu* grad_u(i, t, x, y, z, k+1)*n[k];

    return  out;
}

//we suppose that the problem geometry is the cylinder having axis x, origin (0,0,0), diameter D and height L
Real Womersley::fWallShearStress( const Real& t, const Real& x, const Real& y, const Real& z, const ID& i )
{
    Real r=sqrt(y*y+z*z);
    Real n[3] = {0, 0, 0};
    if        ( x < 1e-6/S_L ) {
        n[0] = -1.;
    } else if ( x >  S_L*(1-1e-6) ) {
        n[0] =  1.;
    } else if ( r > S_D/2.*0.9 ) {
        n[1] = y/r;
        n[2] = z/r;
    } else {
        //std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
        //          << std::endl;
    }

    // wss = s_n - ( s_n \cdot n ) n
    Real wss=0;
    Real s_n_n(0.);
    // this is the i-component of the normal stress
    wss = fShearStress(t, x, y, z, i);

    for (UInt k =0; k< nDimensions; k++) // ( s_n \cdot n )
        s_n_n += fShearStress(t, x, y, z, k+1) * n[k];

    wss -= s_n_n * n[i-1];

    return  wss;
}

void Womersley::setParamsFromGetPot( const GetPot& dataFile )
{
    S_mu = dataFile( "fluid/physics/viscosity", 1. );
    S_flagStrain = dataFile( "fluid/physics/flag_strain", 0 );
    S_nu = S_mu / dataFile( "fluid/physics/density", 1. );
    S_D = dataFile( "fluid/problem/D", 1. );
    S_T = dataFile( "fluid/problem/T", 1. );
    S_rho = dataFile( "fluid/physics/density", 1. );
    S_W0 =S_D/2*sqrt(2*Pi/S_nu/S_T);
    S_L = dataFile( "fluid/problem/L", 1. );
    S_A = dataFile( "fluid/problem/A", 16*S_mu/S_D/S_D*S_L*S_L);
    S_ii = std::complex<Real>(0.,1.);
    S_w = 2.*Pi/S_T;
    S_wi = S_w*S_ii;
    S_z1 = S_W0*pow(S_ii,1.5);
    cbessjy01(S_z1, S_b1, S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p);

}

void Womersley::showMe()
{
	std::cout << "Kynetic viscosity " << S_nu << std::endl;
	std::cout << "Pipe radius " << S_D/2. << " Pipe lenght " << S_L << std::endl;
	std::cout << "Oscillation period " << S_T << std::endl;
	std::cout << "Pressure Drop " << S_A << std::endl;
	std::cout << "Womersley Number " << S_D/2.*sqrt(S_w/S_nu) << std::endl;
}

Real Womersley::S_nu;
Real Womersley::S_mu;
Int  Womersley::S_flagStrain;
Real Womersley::S_D;
Real Womersley::S_rho;
Real Womersley::S_T;
Real Womersley::S_W0;
Real Womersley::S_L;
Real Womersley::S_A;
Real Womersley::S_w;
std::complex<Real> Womersley::S_cj1, Womersley::S_cy0, Womersley::S_cy1,
	Womersley::S_cj0p, Womersley::S_cj1p, Womersley::S_cy0p, Womersley::S_cy1p ,
	Womersley::S_z1, Womersley::S_b1 , Womersley::S_ii, Womersley::S_wi;
} // namespace LifeV
