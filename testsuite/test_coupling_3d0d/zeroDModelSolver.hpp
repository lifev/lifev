/* -*- mode: c++ -*-
   This program is part of the LifeV library 
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
/*!
  \file zeroDModelSolver.hpp
  \author Alexandra Moura
  \date 09/2004
  \version 1.0

  \brief This file contains a solver for a very specific case of the 0D model 
         (small electrical network whitch forcing term is a coseno) to be coupled
         with a 3D model (NS), i.e., a part of the network is substituted by the 3D model 
*/

#ifndef _ZERODMODELSOLVER_H_
#define _ZERODMODELSOLVER_H_
#include <string>
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include "RNM.hpp"

#include "chrono.hpp"
#include "GetPot.hpp"

#define pi 3.14159265

namespace LifeV
{
template<typename SOLVER>
class zeroDModelSolver // getPressureFromQ
{

public:
  
  //! Constructor
  /*!
    \param solver (containing a member flux(int)) 
    \param time step
  */
  zeroDModelSolver(SOLVER  & solver, Real & dt);

  //! Computes the pressure of the network, to give to the 3D NS
  //! \param current time stepx
  Real pression(Real & time);

private:

  SOLVER  * _M_solver;
  Real _M_dt;
 
  Real _M_C1, _M_C2, _M_C3, _M_C4; //Capacities in the network
  Real _M_L5, _M_L6, _M_L7, _M_L8; //Induntances in the network
  Real _M_R5, _M_R6, _M_R7, _M_R8; //Resitences in the network
  Real _M_U;                              //Forcing term
  Real _M_Q;                              //Flux from the 3D NS

  Real x[8], y[8];  

  std::ofstream _M_outfile;

};

//***************************************************************************************//
//                                     IMPLEMENTATION                                    //
//***************************************************************************************//

template<typename SOLVER>
zeroDModelSolver<SOLVER>::zeroDModelSolver(SOLVER  & solver, Real & dt):
  //getPressureFromQ<SOLVER>::getPressureFromQ(SOLVER  & solver, Real & dt):
  _M_solver(&solver),
  _M_dt(dt),

  _M_C1(0.05), 
  _M_C2(0.5), 
  _M_C3(0.02), 
  _M_C4(0.001),
  _M_L5(0.5), 
  _M_L6(0.1),  
  _M_L8(0.1),
  _M_R5(5), 
  _M_R6(6), 
  _M_R8(5.001)

{
  x[0]=0.0178;
  x[1]=0.3740;
  x[2]=0.3346;
  x[3]=-0.0309;
  x[4]=-0.0908;
  x[5]=-0.0011;
  x[6]=0.0211;
  x[7]=0.0257;

  _M_outfile.open("res_Q.m", std::ios::app);
  _M_outfile << "Q = [ " << std::endl;
  _M_outfile.close();

  _M_outfile.open("res_DP.m", std::ios::app);
  _M_outfile << "DP = [ " << std::endl;
  _M_outfile.close();

}

template<typename SOLVER>
Real
zeroDModelSolver<SOLVER>::pression(Real & time)
  //getPressureFromQ<SOLVER>::pression(Real & time)
{ 
  _M_Q=_M_solver->flux(1); //flux coming from the 3D Model
  _M_U=0.1+cos(2*pi*time);

  y[0]=x[0]+(1/_M_C1*x[4]-1/_M_C1*x[7])*_M_dt;
  y[1]=x[1]+(-1/_M_C2*x[4]+1/_M_C2*x[5])*_M_dt;
  y[2]=x[2]+(-1/_M_C3*x[5]+1/_M_C3*x[6])*_M_dt;
  y[3]=x[3]+(1/_M_C4*x[7]-1/_M_C4*x[6])*_M_dt;
  y[4]=x[4]+(-1/_M_L5*x[0]+1/_M_L5*x[1]-_M_R5/_M_L5*x[4]+_M_U/_M_L5)*_M_dt;
  y[5]=x[5]+(-1/_M_L6*x[1]+1/_M_L6*x[2]-_M_R6/_M_L6*x[5]-_M_U/_M_L6)*_M_dt;
  y[6]=-_M_Q;  //The flux is now given from the 3D Model -- y[6]=x[6] +(x[3]/L-delta_P/L-x[2]/L)*_dt
  y[7]=x[7]+(1/_M_L8*x[0]-1/_M_L8*x[3]-_M_R8/_M_L8*x[7])*_M_dt;
    
  for(UInt i=0; i<8; ++i){
    x[i]=y[i];
  };

  _M_outfile.open("res_Q.m", std::ios::app);
  _M_outfile << "     " << y[6] << ";" << std::endl;
  _M_outfile.close();

  _M_outfile.open("res_DP.m", std::ios::app);
  _M_outfile << "      " <<  y[3]-y[2] << ";" << std::endl;
  _M_outfile.close();

  return y[3]-y[2];
}
}

#endif
