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
    @file Laplacian Analytical Solution
    @brief This file contains the exact solution and some additive function for testing codes

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 2011-08-03

    Exact solution for the problem \f$-\Delta\mathbf{u}=\mathbf{f} on the cube [0,1]x[0,1]x[0,1].
 */

#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP 1

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

class Laplacian
{
public:
    static Real f ( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );
    static Real uexact ( const Real& t, const Real& x, const Real& y,
                         const Real& z, const ID& i );

    static Real duexactdx ( const Real& t, const Real& x, const Real& y,
                            const Real& z, const ID& i );

    static Real duexactdy ( const Real& t, const Real& x, const Real& y,
                            const Real& z, const ID& i );

    static Real duexactdz ( const Real& t, const Real& x, const Real& y,
                            const Real& z, const ID& i );

    static void setModes ( const Int& xMode, const Int& yMode, const Int& zMode );

private:
    static Int  M_xMode;
    static Int  M_yMode;
    static Int  M_zMode;

}; // class Laplacian

} // namespace LifeV

#endif /* LAPLACIAN_HPP */
