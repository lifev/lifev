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
   \file ethierSteinman.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-11-12
*/
#ifndef __ESUNSTEADY_HPP
#define __ESUNSTEADY_HPP 1

#include <life/lifecore/life.hpp>
#include <life/lifefilters/GetPot.hpp>

namespace LifeV
{

class EthierSteinmanUnsteady
{
public:
    static Real f        ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i );
    static Real xexact   ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i );
    static Real uexact   ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i );
    static Real uderexact( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i );
    static Real pexact   ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i );

    // Initial velocity
    static Real x0( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );

    static Real fNeumann( const Real& t, const Real& x, const Real& y,
                          const Real& z, const ID& i );
    static void setParamsFromGetPot( const GetPot& dataFile );

private:
    // derivatives for neumann
    static Real ux( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );
    static Real uy( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );
    static Real uz( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );

    static Real a;
    static Real d;
    static Real mu;
    static Real nu;
}; // class EthierSteinmanUnsteady

} // namespace LifeV

#endif /* __ETHIER_STEINMAN_HPP */
