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


void FlowConditions::initParameters( FSIOperator&  Oper,
                                     const int&    outflowFlag)
{

  UInt flag = 1;
    Epetra_SerialDenseVector fluidQuantities(1); // M_area0
    Epetra_SerialDenseVector solidQuantities(2); // M_beta and M_rhos

    /*
        std::cout << "Here I am " << std::endl;
        int n1;
        std::cin >> n1;
    */


    if (Oper.isFluid())
    {

        fluidQuantities(0) = Oper.fluid().area(outflowFlag);

    }

    M_area0      = fluidQuantities(0);
    M_outRadius0 = std::sqrt(M_area0/pi);
    M_inRadius0 = M_outRadius0;

    Oper.displayer().leaderPrint( "  Outflow BC : area0     = ", M_area0 );
    Oper.displayer().leaderPrint( "  Outflow BC : radius    = ", M_outRadius0 );
    if (Oper.isSolid())
    {
      solidQuantities(0) =  ( ( Oper.solid().thickness()*Oper.solid().young( flag )     )/( 1 - Oper.solid().poisson( flag )*Oper.solid().poisson( flag ) )*pi/M_area0 );

        solidQuantities(1) = Oper.solid().rho( );

        Oper.displayer().leaderPrint( "  Outflow BC : thickness = " , Oper.solid().thickness() );
        Oper.displayer().leaderPrint( "  Outflow BC : young     = " , Oper.solid().young( flag ) );
        Oper.displayer().leaderPrint( "  Outflow BC : poisson   = " , Oper.solid().poisson( flag ) );

    }

    //Oper.worldComm().Broadcast( solidQuantities.Values(), solidQuantities.Length(),
    //Oper.getSolidLeaderId() );


    M_beta  = solidQuantities(0);
    M_rhos  = solidQuantities(1);
    Oper.displayer().leaderPrint( "  Outflow BC : beta      = " , M_beta );
    Oper.displayer().leaderPrint( "  Outflow BC : rho       = " , M_rhos );


}

void FlowConditions::renewParameters ( FSISolver&  oper_,
                                       const int&    outflowFlag)
{

    Epetra_SerialDenseVector fluidQuantities(2); // Flux and Area
    //Epetra_SerialDenseVector solidQuantities(0); // M_beta and M_rhos
    FSIOperator* Oper(oper_.FSIOper().get());

    if (Oper->isFluid())
    {
        fluidQuantities(0) = Oper->fluid().flux(outflowFlag, *Oper->fluid().solution());
        fluidQuantities(1) = Oper->fluid().area(outflowFlag);
    }

    Oper->worldComm()->Broadcast( fluidQuantities.Values(), fluidQuantities.Length(),
                                  Oper->getFluidLeaderId() );


    Real qn;
    Real area;
    Real area0;

    qn   = fluidQuantities(0);
    area = fluidQuantities(1);
    area0 = 0.193529;
    //Fluid density
    Real density = 1.0;
    UInt flag   = 1;

    // Setting parameters for our simulation:
    // if imposing the absorbing boundary condition through the pressure:
    if (bcOnFluid)
    {

        // Moura et al.
          //Alexandra's Abc
    Real exp  = 5/4;
    Real beta = ( std::sqrt(PI) * Oper->solid().thickness() * Oper->solid().young( flag ) ) / (1 - Oper->solid().poisson( flag ) * Oper->solid().poisson( flag ) );
    Real R    = ( std::sqrt(Oper->solid().rho( ) * beta ) ) / ( std::sqrt(2.0) * std::pow(area0,exp) );

    M_outP       = R * qn;

        // Nobile & Vergara
        // M_outP =  pow((sqrt(density)/(2.*sqrt(2))*qn/area + sqrt(M_beta*sqrt(M_area0))),2)
        //           - M_beta*sqrt(M_area0);
        FlowConditions::outputVector[conditionNumber]=M_outP;

        Oper->displayer().leaderPrint( " Flow rate = " , qn );
        Oper->displayer().leaderPrint( " outflow pressure   = " , M_outP );

        M_outDeltaRadius = 0;


    }
    else
    {
        // if imposing the absorbing boundary condition through change in radius: --> Not ready
#ifdef  TESTING
        M_outP = Pout;

        area = qn * std::sqrt(M_rhos) / ( (2.*std::sqrt(2)) *
                                          std::sqrt( M_outP + M_beta*sqrt(M_area0) ) - std::sqrt( M_beta*sqrt(M_area0) ) );

        assert(area >= 0 );
        if (area < 1e-8*M_area0) area = M_area0;

        M_outDeltaRadius = std::sqrt( area / pi  ) - M_outRadius0;

        Oper->displayer().leaderPrint( " outflow A = " , area );
        Oper->displayer().leaderPrint( " outflow dr = " , M_outDeltaRadius );
        Oper->displayer().leaderPrint( " Flow rate = " , qn );
        Oper->displayer().leaderPrint( " outflow pressure   = " , M_outP );
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

Real FlowConditions::outPressure0(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[0];
}

Real FlowConditions::outPressure1(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[1];
}

Real FlowConditions::outPressure2(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[2];
}
Real FlowConditions::outPressure3(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[3];
}

Real FlowConditions::outPressure4(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[4];
}


Real FlowConditions::outPressure5(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[5];
}

Real FlowConditions::outPressure6(const Real&/*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -FlowConditions::outputVector[6];
}


Real FlowConditions::inDeltaRadius (const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
    if (i == 2) return 0;

    Real r ( sqrt(x*x + y*y) );

    switch (i)
    {
    case 0:
        return M_inDeltaRadius * x/r;
    case 1:
        return M_inDeltaRadius * y/r;
    default:
        ERROR_MSG("This entry is not allowed: flowConditions.hpp");
    };
    return 0.;
}

Real FlowConditions::outDeltaRadius(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const ID& i)
{
    if (i == 2) return 0;

    Real r ( sqrt(x*x + y*y) );

    switch (i)
    {
    case 0:
        return M_outDeltaRadius * x/r;
    case 1:
        return M_outDeltaRadius * y/r;
    default:
        ERROR_MSG("This entry is not allowed: flowConditions.hpp");
    };
    return 0.;
}

std::vector<Real> FlowConditions::outputVector;
}
