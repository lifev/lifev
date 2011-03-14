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

#ifndef ETHIERSTEINMANUNSTEADY_HPP
#define ETHIERSTEINMANUNSTEADY_HPP 1

#include <life/lifecore/LifeV.hpp>
#include <life/lifefilters/GetPot.hpp>

namespace LifeV
{

class EthierSteinmanUnsteady
{
private:
    // derivatives for neumann
    Real ux( const Real& t, const Real& x, const Real& y,
             const Real& z, const ID& i ) const;
    Real uy( const Real& t, const Real& x, const Real& y,
             const Real& z, const ID& i ) const;
    Real uz( const Real& t, const Real& x, const Real& y,
             const Real& z, const ID& i ) const;

    Real a;
    Real d;
    Real mu;
    Real density;
    Real nu;

public:
    EthierSteinmanUnsteady();
    EthierSteinmanUnsteady(const EthierSteinmanUnsteady& ESU);
    ~EthierSteinmanUnsteady();

    Real f        ( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i ) const;
    Real xexact   ( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i ) const;
    Real uexact   ( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i ) const;
    Real uderexact( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i ) const;
    Real pexact   ( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i ) const;

    // Initial velocity
    Real x0( const Real& t, const Real& x, const Real& y,
             const Real& z, const ID& i ) const;

    Real fNeumann( const Real& t, const Real& x, const Real& y,
                   const Real& z, const ID& i ) const;

    void setA(const Real& a);
    void setD(const Real& d);
    void setViscosity(const Real& mu);
    void setDensity(const Real& density);

}; // class EthierSteinmanUnsteady

void setParamsFromGetPot( EthierSteinmanUnsteady& ESU, const GetPot& dataFile );

} // namespace LifeV

#endif /* ETHIERSTEINMANUNSTEADY_HPP */
