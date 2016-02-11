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

Real fPressure (const Real& time, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
	Real i_HeartBeat = 0.0;
	Real T_heartbeat = 0.8;
	Real Q;

	if ( time < T_heartbeat )
	{
		i_HeartBeat = 0.0;
	}
	else if ( time >= T_heartbeat && time < 2*T_heartbeat )
	{
		i_HeartBeat = 1.0;
	}
	else if ( time >= 2*T_heartbeat && time < 3*T_heartbeat )
	{
		i_HeartBeat = 2.0;
	}
	else if ( time >= 3*T_heartbeat && time < 4*T_heartbeat )
	{
		i_HeartBeat = 3.0;
	}

	if ( (time >= 0.05 && time <= 0.42) || (time >= (0.05+T_heartbeat) && time <= (0.42+T_heartbeat) ) || (time >= (0.05+2*T_heartbeat) && time <= (0.42+2*T_heartbeat) ) || (time >= (0.05+3*T_heartbeat) && time <= (0.42+3*T_heartbeat) ) )
	{
		// Q_in = 2.422818092859456e+8*std::pow(time,8)-4.764207344433996e+8*std::pow(time,7) + 3.993883831476327e+8*std::pow(time,6) -1.867066900011057e+8*std::pow(time,5) +0.533079809563519e+8*std::pow(time,4) -0.094581323616832e+8*std::pow(time,3) +0.009804512311267e+8*std::pow(time,2) -0.000482942399225e+8*time+0.000008651437192e+8;
		Q = 2.117637666632775e+04*std::pow(time-i_HeartBeat*T_heartbeat,6)-3.370930726888496e+04*std::pow(time-i_HeartBeat*T_heartbeat,5)+2.133377678002176e+04*std::pow(time-i_HeartBeat*T_heartbeat,4)-6.666366536069445e+03*std::pow(time-i_HeartBeat*T_heartbeat,3)+1.011772959679957e+03*std::pow(time-i_HeartBeat*T_heartbeat,2)-6.023975547926423e+01*(time-i_HeartBeat*T_heartbeat)+1.192718364532979e+00;
	}
	else
	{
		Q = 0.0;
	}

	Real pressure = 1500.0/2.51*(299.5*Q);

	return pressure;
}

Real inflow (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real Q_hat = 1;
    Real Tr    = 0.1;
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
