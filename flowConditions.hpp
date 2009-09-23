/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Paolo Crosetto <crosetto@iacspc70.epfl.ch>
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
   \file flowConditions.hpp
   \author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   \date 2009-06-03
 */
#ifndef __FLOWCONDITIONS_HPP
#define __FLOWCONDITIONS_HPP

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSIOperator.hpp>

namespace LifeV
{
class FlowConditions
{
public:

    FlowConditions();

    void initParameters      ( FSIOperator&  oper,
                                    const int&    outflowFlag);

    void renewParameters     ( FSISolver&  oper,
                                    const int&    outflowFlag);

    Real fZero               (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    static Real outPressure0         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure1         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure2         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure3         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure4         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure5         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    static Real outPressure6         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    Real inDeltaRadius          (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);
    Real outDeltaRadius          (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);


private:
    Real pi;

    bool bcOnFluid;

    Real M_outflux;
    Real M_influx;
    Real M_outP;

    Real M_area0;
    Real M_inRadius0;
    Real M_outRadius0;
    Real M_inDeltaRadius;
    Real M_outDeltaRadius;

    Real M_beta;
    Real M_rhos;
    static    std::vector<Real> outputVector;
    UInt conditionNumber;
};
}

#endif /* __FLOWCONDITIONS_HPP */

