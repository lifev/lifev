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
 *  @file
 *  @brief File containing the boundary conditions for the Monolithic Test
 *
 *  @date 2009-04-09
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @maintainer Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#ifndef UDFNS_HPP
#define UDFNS_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>

namespace LifeV
{

Real fZero (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
	return 0.0;
}

Real inflow (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real Q_hat = 5.0;
    Real Tr    = 0.01;
    Real Q     = 0;
    
    if (t<=Tr)
    {
        Q = Q_hat/2.0*(1.0 - std::cos(t*M_PI/Tr));
    }
    else
    {
        Q = Q_hat;
    }
    
    Real fluidRadiusSquared = 0.5*0.5;
    Real A = M_PI * fluidRadiusSquared;
    
	switch (i)
	{
	case 0:
		return 0.0;
		break;
	case 1:
		return 0.0;
		break;
	case 2:
		return 2.0*Q/A*(fluidRadiusSquared-(x*x+y*y))/(fluidRadiusSquared);
		break;
	}
	return 0;
}

}



#endif
