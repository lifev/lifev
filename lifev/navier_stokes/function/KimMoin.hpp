/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Mauro Perego <mauro@mathcs.emory.edu>
      Date: 2009-10-02

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

//! analytic solution of Kim & Moin for unsteady Navier-Stokes 2D on the square [0,1]x[0,1].

#ifndef __KIM_MOIN_HPP
#define __KIM_MOIN_HPP 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>

namespace LifeV
{

class KimMoin
{
public:
    static Real f ( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );

    static Real xexact ( const Real& t, const Real& x, const Real& y,
                         const Real& z, const ID& i );
    static Real uexact ( const Real& t, const Real& x, const Real& y,
                         const Real& z, const ID& i );
    static Real uderexact ( const Real& t, const Real& x, const Real& y,
                            const Real& z, const ID& i );
    static Real pexact ( const Real& t, const Real& x, const Real& y,
                         const Real& z, const ID& i );

    // Initial velocity
    static Real x0 ( const Real& t, const Real& x, const Real& y,
                     const Real& z, const ID& i );

    static Real u0 ( const Real& t, const Real& x, const Real& y,
                     const Real& z, const ID& i );

    static Real p0 ( const Real& t, const Real& x, const Real& y,
                     const Real& z, const ID& i );

    static Real grad_u ( const UInt& icoor, const Real& t, const Real& x, const Real& y,
                         const Real& z, const ID& i );

    static Real fNeumann ( const Real& t, const Real& x, const Real& y,
                           const Real& z, const ID& i );

    static Real normalVector ( const Real& t, const Real& x, const Real& y,
                               const Real& z, const ID& i );

    static Real fShearStress ( const Real& t, const Real& x, const Real& y,
                               const Real& z, const ID& i );

    static Real fWallShearStress ( const Real& t, const Real& x, const Real& y,
                                   const Real& z, const ID& i );
    static void setParamsFromGetPot ( const GetPot& dataFile );

    static bool flag_strain;


private:

    static Real mu;
    static Real rho;
    static Real nu;
    static Real a;
}; // class KimMoin

} // namespace LifeV

#endif /* __KIM_MOIN_HPP */
