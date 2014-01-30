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
 */

#include "flowConditions.hpp"

#define PI 3.141592653589793

namespace LifeV
{
FlowConditions::FlowConditions() :
    pi (3.141592635),
    M_outflux (0),
    M_resistance (0),
    M_hydrostaticP (0),
    M_outP (0),
    M_flag (0),
    M_name ( ),
    conditionNumber (0)
{
    outputVector.push_back (0);
    conditionNumber = FlowConditions::outputVector.size() - 1;
}


void FlowConditions::initParameters ( const int flag,
                                      const Real resistance,
                                      const Real hydrostatic,
                                      const std::string name)
{
    // In this function we set up the pressure and the resistance value
    // The resistance BCs are applied as Neumann BCs using the flux at
    // the previous time

    M_resistance   = resistance;
    M_hydrostaticP = hydrostatic;
    M_flag         = flag;
    M_name         = name;

}

void FlowConditions::renewParameters ( OseenSolver<RegionMesh<LinearTetra> > &  solver,
                                       const VectorEpetra& solution)
{

    // Compute the flux using the solution on the desired flag
    M_outflux = solver.flux( M_flag, solution);

    M_outP = 1.0 * ( M_resistance * M_outflux + M_hydrostaticP );

    FlowConditions::outputVector[conditionNumber] = M_outP;
}




Real FlowConditions::fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real FlowConditions::outPressure0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[0];
}

Real FlowConditions::outPressure1 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[1];
}

Real FlowConditions::outPressure2 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[2];
}
Real FlowConditions::outPressure3 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[3];
}

Real FlowConditions::outPressure4 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[4];
}


Real FlowConditions::outPressure5 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[5];
}

Real FlowConditions::outPressure6 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[6];
}

std::vector<Real> FlowConditions::outputVector;
}
