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

	esponentially increasing analytic solution of Ross-Ethier Steinmann for unsteady Navier-Stokes 3D on the cube [0,1]x[0,1]x[0,1].
	@see Exact fully 3D Navier-Stokes solutions for benchmarking - C. Ross Ethier, D. A. Steinman
 */

/*
 *
AU: C. Ross Ethier
AU: D. A. Steinman
TI: Exact fully 3D Navier-Stokes solutions for benchmarking
SO: International Journal for Numerical Methods in Fluids
VL: 19
NO: 5
PG: 369-375
YR: 1994
CP: Copyright © 1994 John Wiley & Sons, Ltd
ON: 1097-0363
PN: 0271-2091
AD: Department of Mechanical Engineering, 5 King's College Road, University of Toronto, Toronto, Ontario M5S IA4, Canada
DOI: 10.1002/fld.1650190502
US: http://dx.doi.org/10.1002/fld.1650190502
 *
 */

#ifndef __ROSS_ETHIER_STEINMAN_DEC_HPP
#define __ROSS_ETHIER_STEINMAN_DEC_HPP 1

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>

namespace LifeV
{

class RossEthierSteinmanUnsteadyDec
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
    static Real grad_u( const UInt& icoor, const Real& t, const Real& x, const Real& y,
                            const Real& z, const ID& i );

    // Initial velocity
    static Real x0( const Real& t, const Real& x, const Real& y,
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



private:

    static Real S_a;
    static Real S_d;
    static Real S_mu;
    static Real S_rho;
    static Real S_nu;
    static Int  S_flagStrain;
}; // class RossEthierSteinmanUnsteadyDec

} // namespace LifeV

#endif /* __ROSS_ETHIER_STEINMAN_DEC_HPP */
