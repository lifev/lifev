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
		// old
		//Q = 2.117637666632775e+04*std::pow(time-i_HeartBeat*T_heartbeat,6)-3.370930726888496e+04*std::pow(time-i_HeartBeat*T_heartbeat,5)+2.133377678002176e+04*std::pow(time-i_HeartBeat*T_heartbeat,4)-6.666366536069445e+03*std::pow(time-i_HeartBeat*T_heartbeat,3)+1.011772959679957e+03*std::pow(time-i_HeartBeat*T_heartbeat,2)-6.023975547926423e+01*(time-i_HeartBeat*T_heartbeat)+1.192718364532979e+00;

		Q = -2.314569820334801e+09*std::pow(time-i_HeartBeat*T_heartbeat,9) +
				4.952537061974133e+09*std::pow(time-i_HeartBeat*T_heartbeat,8) -
				4.532060231242586e+09*std::pow(time-i_HeartBeat*T_heartbeat,7) +
				2.325743716202249e+09*std::pow(time-i_HeartBeat*T_heartbeat,6) -
				7.387577876374097e+08*std::pow(time-i_HeartBeat*T_heartbeat,5) +
				1.514516710083440e+08*std::pow(time-i_HeartBeat*T_heartbeat,4) -
				2.018053394181958e+07*std::pow(time-i_HeartBeat*T_heartbeat,3) +
				1.667954643625200e+06*std::pow(time-i_HeartBeat*T_heartbeat,2) -
				7.160662399848596e+04*(time-i_HeartBeat*T_heartbeat) +
				1.184312187078482e+03;
		Q = Q/394;
	}
	else
	{
		Q = 0.0;
	}

	Real Q_inflow = 394*Q;
	Real Q_flag4  = 20.25*Q; // left_common_carotid
	Real Q_flag5  = 21.82*Q; // right_common_carotid
	Real Q_flag6  = 1.43*Q;  // right_vertebral
	Real Q_flag7  = 24.77*Q; // right_subclavian
	Real Q_flag8  = 4.69*Q;  // left_vertebral
	Real Q_flag9  = 21.54*Q; // left_subclavian

	Real pressureValue = 1500.0/2.51*(Q_inflow - Q_flag4 - Q_flag5 - Q_flag6 - Q_flag7 - Q_flag8 - Q_flag9);

	return pressureValue;
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
