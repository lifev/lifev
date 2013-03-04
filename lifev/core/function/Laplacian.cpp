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

    Implementation.
*/

#include <lifev/core/function/Laplacian.hpp>

namespace LifeV
{

Real Laplacian::f ( const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */ )
{
    return  (M_xMode * M_xMode + M_yMode * M_yMode + M_zMode * M_zMode) * 4 * M_PI * M_PI * std::sin (M_xMode * 2 * M_PI * x) * std::sin (M_yMode * 2 * M_PI * y) * std::sin (M_zMode * 2 * M_PI * z);
}

Real Laplacian::uexact ( const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */ )
{
    return  std::sin (M_xMode * 2 * M_PI * x) * std::sin (M_yMode * 2 * M_PI * y) * std::sin (M_zMode * 2 * M_PI * z);
}
Real Laplacian::duexactdx ( const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */ )
{
    return  M_xMode * 2 * M_PI * std::cos (M_xMode * 2 * M_PI * x) * std::sin (M_yMode * 2 * M_PI * y) * std::sin (M_zMode * 2 * M_PI * z);
}

Real Laplacian::duexactdy ( const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */ )
{
    return  M_yMode * 2 * M_PI * std::sin (M_xMode * 2 * M_PI * x) * std::cos (M_yMode * 2 * M_PI * y) * std::sin (M_zMode * 2 * M_PI * z);
}

Real Laplacian::duexactdz ( const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */ )
{
    return  M_zMode * 2 * M_PI * std::sin (M_xMode * 2 * M_PI * x) * std::sin (M_yMode * 2 * M_PI * y) * std::cos (M_zMode * 2 * M_PI * z);
}

void Laplacian::setModes ( const Int& xMode, const Int& yMode, const Int& zMode )
{
    M_xMode = xMode;
    M_yMode = yMode;
    M_zMode = zMode;
}

Int  Laplacian::M_xMode = 1;
Int  Laplacian::M_yMode = 1;
Int  Laplacian::M_zMode = 1;

} // namespace LifeV

