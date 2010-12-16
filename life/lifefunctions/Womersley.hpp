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

    Analytic solution of Womersley for unsteady Navier-Stokes 3D on the cylinder having axis x, origin (0,0,0), diameter D and height L.
	Solution of incompressible NS equation in a cylindrical vessel with a sinusoidal pressure drop (deltaP = A cos(wt) )
	between the inflow and the outflow and no-slip conditions on the vessel wall.
	this solution works also in the 2D-axisymmetric formulation (geometry [0, L]x[0, D/2]).
	<br>
	The Womersley number \alpha arises in the solution of the linearized Navier Stokes equations
	for oscillatory flow (presumed to be laminar and incompressible) in a tube.
	When \alpha is small (1 or less), it means the frequency of pulsations is sufficiently low
	that a parabolic velocity profile has time to develop during each cycle,
	and the flow will be very nearly in phase with the pressure gradient,
	and will be given to a good approximation by Poiseuille's law,
	using the instantaneous pressure gradient.
	When \alpha is large (10 or more), it means the frequency of pulsations is sufficiently large
	that the velocity profile is relatively flat or plug-like,
	and the mean flow lags the pressure gradient by about 90 degrees.

 */

#ifndef __WOMERSLEY_HPP
#define __WOMERSLEY_HPP 1

#include <life/lifearray/tab.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>



namespace LifeV
{

class Womersley
{
public:
    static Real f( const Real& t, const Real& x, const Real& y,
                   const Real& z, const ID& i );

    static Real xexact( const Real& t, const Real& x, const Real& y,
                        const Real& z, const ID& i );
    static Real uexact( const Real& t, const Real& x, const Real& y,
                        const Real& z, const ID& i );
    static Real pexact( const Real& t, const Real& x, const Real& y,
                        const Real& z, const ID& i );

    // Initial velocity
    static Real x0( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );

    static Real u0( const Real& t, const Real& x, const Real& y,
                        const Real& z, const ID& i );

    static Real p0( const Real& t, const Real& x, const Real& y,
                        const Real& z, const ID& i );

    static Real grad_u( const UInt& icoor, const Real& t, const Real& x, const Real& y,
                                     const Real& z, const ID& i );

    static Real fNeumann( const Real& t, const Real& x, const Real& y,
                          const Real& z, const ID& i );

    static Real normalVector( const Real& t, const Real& x, const Real& y,
                          const Real& z, const ID& i );

    static Real fShearStress( const Real& t, const Real& x, const Real& y,
                              const Real& z, const ID& i );

    static Real fWallShearStress( const Real& t, const Real& x, const Real& y,
                              const Real& z, const ID& i );
    static void setParamsFromGetPot( const GetPot& dataFile );
    static void showMe();


private:

    static Int  S_flagStrain;
    static Real S_mu;
    static Real S_nu;
    static Real S_D;
    static Real S_T;
    static Real S_rho;
    static Real S_W0;
    static Real S_L;
    static Real S_A;
    static Real S_w;
    static std::complex<Real> S_cj1, S_cy0, S_cy1, S_cj0p, S_cj1p, S_cy0p, S_cy1p , S_z1, S_b1 , S_ii, S_wi;
}; // class Womersley

} // namespace LifeV

#endif /* __WOMERSLEY_HPP */
