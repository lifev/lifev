/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

namespace LifeV
{
    Real g(const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& /* i */) {
        return 0.;
    }

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_type;

    struct srcfct : public source_type {
        Real operator()(Real const& /* x */, Real const& /* y */, Real const& /* z */, ID const& /* comp_id */) {
            return 0.;
        }
    };

    Real sphere(const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */) {
        return sqrt( x * x + y * y + z * z ) - .4;
    }

    Real zero(const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& /* i */) {
        return 0.;
    }

    Real gravity(const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& i) {
        const Real g = - 9.81;
        if(i == 3)
            return g;
        else
            return 0.;
    }
}
