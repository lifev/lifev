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
    @brief

    @author
    @date 00-00-0000

 */

#ifndef UDF_HPP
#define UDF_HPP

#include <life/lifecore/life.hpp>

namespace LifeV
{
Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real u1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real fZero(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

// Initial velocity
Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
Real p0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real u2(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);


// Initial displacement and velocity
Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

Real w0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
Real PhysFlux(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i);
Real aortaPhysPress(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i);

}

#endif
