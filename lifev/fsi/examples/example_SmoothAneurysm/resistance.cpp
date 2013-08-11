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

#include "resistance.hpp"

#define PI 3.141592653589793

namespace LifeV
{

ResistanceBCs::ResistanceBCs() :
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
    conditionNumber = ResistanceBCs::outputVector.size() - 1;
}


void ResistanceBCs::initParameters ( const int flag,
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

void ResistanceBCs::renewParameters ( OseenSolverShapeDerivative<RegionMesh<LinearTetra> > &  solver,
                                      const VectorEpetra& solution,
                                      const Real time)
{

    M_resistance = 0.0;
    M_resistance = computeResistance( time );

    // Compute the flux using the solution on the desired flag
    M_outflux = solver.flux( M_flag, solution);

    M_outP = 1.0 * ( M_resistance * M_outflux + M_hydrostaticP );

    solver.getDisplayer().leaderPrint ( " ****************** Resistance BCs infos ***************************x\n" );
    solver.getDisplayer().leaderPrint ( " Flow rate = " , M_outflux );
    solver.getDisplayer().leaderPrint ( " \n" );
    solver.getDisplayer().leaderPrint ( " Area Inlet = " , solver.area(2)  );
    solver.getDisplayer().leaderPrint ( " \n" );
    solver.getDisplayer().leaderPrint ( " Area Outlet = ", solver.area(3) );
    solver.getDisplayer().leaderPrint ( " \n" );
    solver.getDisplayer().leaderPrint ( " Hydrostatic pressure   = " , M_hydrostaticP );
    solver.getDisplayer().leaderPrint ( " \n" );
    solver.getDisplayer().leaderPrint ( " Resistance   = " , M_resistance );
    solver.getDisplayer().leaderPrint ( " \n" );
    solver.getDisplayer().leaderPrint ( " Outflow pressure   = " , M_outP );
    solver.getDisplayer().leaderPrint ( " \n" );
    solver.getDisplayer().leaderPrint ( " ****************** Resistance BCs infos ***************************" );

    ResistanceBCs::outputVector[conditionNumber] = M_outP;
}

Real ResistanceBCs::computeResistance (const Real t )
{
    Real resistance(0);

    Real highestResistance( 501950 );
    Real totalTime = 1.012;
    Real halfTime = totalTime / 3.0;

    Real m = ( 0.8 * highestResistance ) / ( totalTime - halfTime );

    if ( t <= halfTime )
      resistance =   ( ( highestResistance / 5 ) / ( halfTime ) ) * t ;

    if ( t > halfTime && t <= totalTime)
      resistance = m * (t - halfTime ) + highestResistance / 5 ;

    if ( t > totalTime )
      resistance = highestResistance;

    // Real a = ( highestResistance / 2 ) * ( 1/ ( halfTime * halfTime ) );

    // if ( t <= halfTime )
    //     resistance = a * t*t;

    // if ( t > halfTime )
    //     resistance = - a * (t - totalTime)*(t - totalTime) + highestResistance;


    return resistance;
}




Real ResistanceBCs::fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real ResistanceBCs::outPressure0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -ResistanceBCs::outputVector[0];
}

Real ResistanceBCs::outPressure1 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -ResistanceBCs::outputVector[1];
}

Real ResistanceBCs::outPressure2 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -ResistanceBCs::outputVector[2];
}
Real ResistanceBCs::outPressure3 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -ResistanceBCs::outputVector[3];
}

Real ResistanceBCs::outPressure4 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -ResistanceBCs::outputVector[4];
}


Real ResistanceBCs::outPressure5 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -ResistanceBCs::outputVector[5];
}

Real ResistanceBCs::outPressure6 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -ResistanceBCs::outputVector[6];
}

std::vector<Real> ResistanceBCs::outputVector;
}
