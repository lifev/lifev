/* -*- mode: c++ -*-

   This file is part of the LifeV library

   Author(s): Paolo Crosetto <paolo.crosetto@epfl.ch>
              Simone Deparis <simone.deparis@epfl.ch>
   Date: 2009-06-03

   Copyright (C) 2009 EPFL

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
   \file flowConditions.cpp
   \author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   \date 2009-06-03
*/


#include "flowConditions.hpp"
#include "Epetra_SerialDenseVector.h"

namespace LifeV
{
FlowConditions::FlowConditions():
    pi(3.141592635),
    bcOnFluid(true),
    M_outflux(0),
    M_influx(0),
    M_outP(0),
    M_area0(0),
    M_inRadius0(0),
    M_outRadius0(0),
    M_inDeltaRadius(0),
    M_outDeltaRadius(0),
    M_beta(0),
    M_rhos(0),
    conditionNumber(0)
{
    outputVector.push_back(0);
    conditionNumber= FlowConditions::outputVector.size()-1;
}


void FlowConditions::initParameters( FSIOperator&  oper,
                                     const int&    outflowFlag)
{

    Epetra_SerialDenseVector fluidQuantities(1); // M_area0
    Epetra_SerialDenseVector solidQuantities(2); // M_beta and M_rhos

    if (oper.isFluid())
    {
        fluidQuantities(0) = oper.fluid().area(outflowFlag);
    }

    M_area0      = fluidQuantities(0);
    M_outRadius0 = std::sqrt(M_area0/pi);
    M_inRadius0 = M_outRadius0;

    oper.displayer().leaderPrint( "  Outflow BC : area0     = ", M_area0 );
    oper.displayer().leaderPrint( "  Outflow BC : radius    = ", M_outRadius0 );


    if (oper.isSolid())
    {
        solidQuantities(0) =  ( ( oper.solid().thickness()*oper.solid().young()     ) /
                                ( 1 - oper.solid().poisson()*oper.solid().poisson() )*
                                pi/M_area0 );

        solidQuantities(1) = oper.solid().rho();

        oper.displayer().leaderPrint( "  Outflow BC : thickness = " , oper.solid().thickness() );
        oper.displayer().leaderPrint( "  Outflow BC : young     = " , oper.solid().young() );
        oper.displayer().leaderPrint( "  Outflow BC : poisson   = " , oper.solid().poisson() );

    }

    //oper.worldComm().Broadcast( solidQuantities.Values(), solidQuantities.Length(),
    //oper.getSolidLeaderId() );


    M_beta  = solidQuantities(0);
    M_rhos  = solidQuantities(1);
    oper.displayer().leaderPrint( "  Outflow BC : beta      = " , M_beta );
    oper.displayer().leaderPrint( "  Outflow BC : rho       = " , M_rhos );


}

void FlowConditions::renewParameters ( FSISolver&  oper_,
                                       const int&    outflowFlag)
{

    Epetra_SerialDenseVector fluidQuantities(2); // Flux and Area
    //Epetra_SerialDenseVector solidQuantities(0); // M_beta and M_rhos
    FSIOperator* oper(oper_.FSIOper().get());

    if (oper->isFluid())
    {
        fluidQuantities(0) = oper->fluid().flux(outflowFlag, oper_.displacement());
        fluidQuantities(1) = oper->fluid().area(outflowFlag);
    }

    oper->worldComm()->Broadcast( fluidQuantities.Values(), fluidQuantities.Length(),
                                 oper->getFluidLeaderId() );


    Real qn;
    Real area;

    qn   = fluidQuantities(0);
    area = fluidQuantities(1);

    // Setting parameters for our simulation:
    // if imposing the absorbing boundary condition through the pressure:
    if (bcOnFluid)
    {
        M_outP =  pow((M_rhos/(2.*sqrt(2))*qn/area + sqrt(M_beta*sqrt(M_area0))),2)
            - M_beta*sqrt(M_area0);
        FlowConditions::outputVector[conditionNumber]=M_outP;

        oper->displayer().leaderPrint( " Flow rate = " , qn );
        oper->displayer().leaderPrint( " outflow pressure   = " , M_outP );

        M_outDeltaRadius = 0;


    } else {
        // if imposing the absorbing boundary condition through change in radius: --> Not ready
#ifdef  TESTING
        M_outP = Pout;

        area = qn * std::sqrt(M_rhos) / ( (2.*std::sqrt(2)) *
                                          std::sqrt( M_outP + M_beta*sqrt(M_area0) ) - std::sqrt( M_beta*sqrt(M_area0) ) );

        assert(area >= 0 );
        if (area < 1e-8*M_area0) area = M_area0;

        M_outDeltaRadius = std::sqrt( area / pi  ) - M_outRadius0;

        oper->displayer().leaderPrint( " outflow A = " , area );
        oper->displayer().leaderPrint( " outflow dr = " , M_outDeltaRadius );
        oper->displayer().leaderPrint( " Flow rate = " , qn );
        oper->displayer().leaderPrint( " outflow pressure   = " , M_outP );
#endif

    }

    // for now applying absBC only at outflow
    M_inDeltaRadius = 0;
    //    M_inP = Pin;

}




Real FlowConditions::fZero(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real FlowConditions::outPressure0(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -FlowConditions::outputVector[0];
}

Real FlowConditions::outPressure1(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -FlowConditions::outputVector[1];
}

Real FlowConditions::outPressure2(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -FlowConditions::outputVector[2];
}
Real FlowConditions::outPressure3(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -FlowConditions::outputVector[4];
}

Real FlowConditions::outPressure5(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -FlowConditions::outputVector[5];
}

Real FlowConditions::outPressure6(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -FlowConditions::outputVector[6];
}


Real FlowConditions::inDeltaRadius (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    if (i == 3) return 0;

    Real r ( sqrt(x*x + y*y) );

    switch(i) {
    case 1:
        return M_inDeltaRadius * x/r;
    case 2:
        return M_inDeltaRadius * y/r;
    default:
        ERROR_MSG("This entry is not allowed: flowConditions.hpp");
    };
    return 0.;
}

Real FlowConditions::outDeltaRadius(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    if (i == 3) return 0;

    Real r ( sqrt(x*x + y*y) );

    switch(i) {
    case 1:
        return M_outDeltaRadius * x/r;
    case 2:
        return M_outDeltaRadius * y/r;
    default:
        ERROR_MSG("This entry is not allowed: flowConditions.hpp");
    };
    return 0.;
}

std::vector<Real> FlowConditions::outputVector;
}
