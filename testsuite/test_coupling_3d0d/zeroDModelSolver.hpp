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
  \file oneDModelSolver.hpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains a solver class for the 1D model. 
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


//#include <iostream>

namespace LifeV
{
template<typename SOLVER>
class getPressureFromQ
{

public:
  getPressureFromQ(SOLVER  & solver, GetPot & data, Real Length, 
		   Real Radius, Real Thick);

  Real pression(Real & time);

private:
  Real C, R, L; //Variables related to the 3D Part
  Real C1, C2, C3, C4; //Capacities in the network
  Real L5, L6, L7, L8; //Induntances in the network
  Real R5, R6, R7, R8; //Resitences in the network
  Real U; //Forcing term
  Real dt;
  Real x[8], y[8];
  SOLVER  * _solver;
  Real pi;

  std::ofstream outfile;

};

//***************************************************************//
//                         IMPLEMENTATION                        //
//***************************************************************//

template<typename SOLVER>
getPressureFromQ<SOLVER>::getPressureFromQ(SOLVER  & solver, GetPot & dfile, Real Length, 
					   Real Radius, Real Thick/*, Real Pinit*/):
  _solver(&solver)
{
  Real rho       = dfile("fluid/physics/density",1.);
  Real mu        = dfile("fluid/physics/viscosity",1.);
  Real E         = dfile("solid/physics/young",1.);
  dt             = dfile("fluid/discretization/timestep",0.); 
  if (Thick!=0)
    C=3*3.1415*pow(Radius,3)*Length/(2*E*Thick);  
  R= rho*8*3.1415*mu*Length/(pow(3.1415*Radius*Radius,2));
  L= (rho*Length)/(3.1415*pow(Radius,2));

  C1=0.05; 
  C2=0.5; 
  C3=0.02; 
  C4=0.001;
  L5=0.5; 
  L6=0.1;  
  L8=0.1;
  R5=5; R6=6; 
  R8=5.001;
 
  pi=acos(-1.);

  x[0]=0.0178;
  x[1]=0.3740;
  x[2]=0.3346;
  x[3]=-0.0309;
  x[4]=-0.0908;
  x[5]=-0.0011;
  x[6]=0.0211;
  x[7]=0.0257;

  outfile.open("res_Q.m", std::ios::app);
  outfile << "Q = [ " << std::endl;
  outfile.close();

  outfile.open("res_DP.m", std::ios::app);
  outfile << "DP = [ " << std::endl;
  outfile.close();

}

template<typename SOLVER>
Real
getPressureFromQ<SOLVER>::pression(Real & time)
{ 
  Real Q=_solver->flux(1);
  U=0.1+cos(2*pi*time);

  y[0]=x[0]+(1/C1*x[4]-1/C1*x[7])*dt;
  y[1]=x[1]+(-1/C2*x[4]+1/C2*x[5])*dt;
  y[2]=x[2]+(-1/C3*x[5]+1/C3*x[6])*dt;
  y[3]=x[3]+(1/C4*x[7]-1/C4*x[6])*dt;
  y[4]=x[4]+(-1/L5*x[0]+1/L5*x[1]-R5/L5*x[4]+U/L5)*dt;
  y[5]=x[5]+(-1/L6*x[1]+1/L6*x[2]-R6/L6*x[5]-U/L6)*dt;
  //y[6]=x[6] +(x[3]/L-delta_P/L-x[2]/L)*dt;
  y[6]=Q;
  y[7]=x[7]+(1/L8*x[0]-1/L8*x[3]-R8/L8*x[7])*dt;
    
  for(UInt i=0; i<8; ++i){
    x[i]=y[i];
  };

  outfile.open("res_Q.m", std::ios::app);
  outfile << "     " << y[7] << ";" << std::endl;
  outfile.close();

  outfile.open("res_DP.m", std::ios::app);
  outfile << "      " <<  y[3]-y[2] << ";" << std::endl;
  outfile.close();

  return y[3]-y[2];
}
}

#endif
