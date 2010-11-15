//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/**
   \file lumpedHeart.cpp
   \author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   \date 2009-06-03
*/


#include <math.h>
#include "lumpedHeart.hpp"
#define TESTING
#define PI 3.14159265

namespace LifeV
{

void LumpedHeart::initParameters( FSIOperator&  oper,
                                  const std::string&    FileName )
{
    M_ODEscheme.initialize_unk(0.);

    GetPot       dataFile(FileName);
    M_dt         =dataFile("problem/timestep", 0.001);
    M_T_max      =dataFile("problem/T_max", 0.8);
    M_E_max      =dataFile("problem/E_max", 1333.22);
    M_V_0        =dataFile("problem/V_0", 2);
    M_Vt_ao      =dataFile("problem/Vt_ao", 5);
    M_RV_art     =dataFile("problem/RV_art", 10);
    M_RA_V       =dataFile("problem/RA_V", 10);
    M_LV_art     =dataFile("problem/LV_art", 0.6879);
    M_LA_V       =dataFile("problem/LA_V", 0.6670);
    M_pressure   =dataFile("problem/p0", 1000);
    M_Tpb        =dataFile("problem/Tpb", 0.1);
    M_Tpw        =dataFile("problem/Tpw", 0.55);
}

Real& LumpedHeart::outPressure         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return LumpedHeart::M_pressure;
}

void LumpedHeart::renewParameters ( FSIOperator&  oper, const int& flag, const Real& time , const Real& flux)
{
    M_intFlux += flux*M_dt;
    //should have a different sign, but it is assigned as a Normal bc so we take the opposite
    LumpedHeart::M_pressure = (-flux*(M_RV_art)-M_elastance(time)*(M_Vt_ao-M_intFlux-M_V_0))+(M_LV_art)*(-flux/M_dt+M_ODEscheme.time_der(M_dt));
    //    flux     = (M_ODEScheme.time_der(M_dt)-M_Pressure+M_elastance*(M_Vt_ao-M_intFlux-M_V_0))/(M_LV_art+M_elastance*M_dt+M_RV_art);
    M_ODEscheme.shift_right(flux);
}

Real LumpedHeart::M_elastance(const Real& t)
{
    Real ret = 1 - cos((t-M_Tpb)/M_Tpw*2*PI);
    if(t<M_Tpb)
        return (0.0);
    else
        if(t<M_Tpb+M_Tpw)
            return ret*1333.22;
        else
            return (0.0);
}


Real LumpedHeart::fZero(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real LumpedHeart::M_pressure(0);
}
