/* -*- mode: c++ -*-
   This program is part of the LifeV library

   Author(s): Alexandra Moura <moura@mate.polimi.it>
       Date: 2004-11-05


   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/**
   \file meanPressure.hpp
   \author Alexandra Moura <moura@mate.polimi.it>
   \date 2004-11-05
 */

#ifndef _MEANPRESSURE
#define _MEANPRESSURE

#include "C3d0d.hpp"

/**

    Mean Pressure strategy for coupling 3D (NS - non compliant) and 0D (lumped parameter model) models

        - 3D passes the flux to 0D, whereas 0D passes the mean pressure to 3D

*/

namespace LifeV
{

  class meanPressure : public C3d0d {

  public:

    //Constructors
    meanPressure(GetPot& data_file, BCHandler& BChx_u, int interfaceFlag);

    //Destructor
    ~meanPressure();

    //Member Functions   

    void initialize(const Function& u0, const Function& p0, Real t0, Real dt);

    void timeAdvance(source_type const& source, Real const& time);

    void iterate(Real const& time);

    //Mutators 

    //Setters
    
    //virtual NavierStokesWithFlux< NS > & nsFlux() = 0;


  protected: 

    int       _M_interfaceFlag;

    UInt      _M_dimFluid;
    Vector    _M_vectorPressure;

    Real      _M_Q;
    Real      _M_deltaP;

  private: 
  
  };  
    
//************************************************************************************
//                                   IMPLEMENTATION                                   
//************************************************************************************


//Constructors
meanPressure::meanPressure(GetPot& data_file, BCHandler& BCh_u, int interfaceFlag)
  :
  C3d0d(data_file, BCh_u),
  _M_interfaceFlag(interfaceFlag),
  _M_dimFluid(_M_ns->uDof().numTotalDof()),
  _M_vectorPressure(_M_dimFluid)
{  
  BCVector _M_bcvecPressure(_M_vectorPressure, _M_dimFluid, 1);
  _M_ns->bcHandler().modifyBC(_M_interfaceFlag, _M_bcvecPressure);
}

//
void 
meanPressure::initialize(const Function& u0, const Function& p0, Real t0, Real dt)
{
  _M_ns->initialize(u0, p0, t0, dt);
}

//
void 
meanPressure::timeAdvance(source_type const& source, Real const& time)
{

  _M_Q      = _M_ns->flux(_M_interfaceFlag);
  _M_deltaP = _M_network.getPressureFromQ(time,_M_Q);
  
  _M_vectorPressure = ScalarVector( _M_dimFluid, -_M_deltaP );

  _M_outfile.open("res_Q.m", std::ios::app);
  _M_outfile << "     " << -_M_Q << ";" << std::endl;
  _M_outfile.close();

  _M_outfile.open("res_DP.m", std::ios::app);
  _M_outfile << "      " <<  _M_deltaP << ";" << std::endl;
  _M_outfile.close();

  _M_ns->timeAdvance(source, time);
}


//
void 
meanPressure::iterate(Real const& time)
{
  _M_ns->iterate(time);

}

// Destructor
meanPressure::~meanPressure()
{
}


} // namespace LifeV

#endif
