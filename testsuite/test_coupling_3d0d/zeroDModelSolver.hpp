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
class zeroDModelSolver
{

public:
  
  //! Constructor
  /*!
    \param solver
    \param data file
    \param Lenght of the 3D tube
    \param Radius of the 3D tube
    \param Thicknes of the wall (only in the case where 3D is compliant - FSI)
  */
  zeroDModelSolver(GetPot & data_file, Real Length, Real Radius, Real Thick);

  //! Computes the pressure of the network, given the flux from 3D 
  //! \param current time step
  //! \param flux from 3d
  Real getPressureFromQ(Real & time, Real & Q);

  //! Computes the flux of the network, given the pression from 3D 
  //! \param current time step
  //! \param pressure from 3D
  Real getQFromPressure(Real & time, Real & deltaP);

  //! Indicates the coupling strategy (acording to data file)
  bool isMeanPressProb();

private:

  Real _M_dt;
  Real _M_rho;
  Real _M_nu;
  Real _M_E;
  UInt _M_isMeanPressProb;
 
  Real _M_C1, _M_C2, _M_C3, _M_C4; // Capacities in the network
  Real _M_L5, _M_L6, _M_L7, _M_L8; // Induntances in the network
  Real _M_R5, _M_R6, _M_R7, _M_R8; // Resitences in the network
  Real _M_U;                       // Forcing term
  
  Real _M_L, _M_R, _M_C;           // Parameters for the electrical net analog of the 3D 
                                   // Model (here substituted for the 3D Model itself)

  Real x[8], y[8];  

  std::ofstream _M_outfile;

};

//***************************************************************************************//
//                                     IMPLEMENTATION                                    //
//***************************************************************************************//

zeroDModelSolver::zeroDModelSolver(GetPot & data_file, Real Length, Real Radius, Real Thick):

  _M_C1(0.05), 
  _M_C2(0.5), 
  _M_C3(0.02), 
  _M_C4(0.001),
  _M_L5(0.5), 
  _M_L6(0.1),  
  _M_L8(0.1),
  _M_R5(5.), 
  _M_R6(6.), 
  _M_R8(5.001)

{
  _M_dt  =data_file("fluid/discretization/timestep", 0. );
  _M_rho =data_file("fluid/physics/density", 1. );
  _M_nu  =data_file("fluid/physics/viscosity", 1. );
  _M_E   =data_file("solid/physics/young", 1. );
  _M_isMeanPressProb = data_file("fluid/miscellaneous/coupStrategy", 1 );

  if (Thick!=0)
    _M_C=3*pi*pow(Radius,3)*Length/(2*_M_E*Thick);  
  _M_R= _M_rho*8*pi*_M_nu*Length/(pow(pi*Radius*Radius,2));
  _M_L= (_M_rho*Length)/(pi*pow(Radius,2));

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

Real zeroDModelSolver::getPressureFromQ(Real & time, Real & Q)
{ 
  Real _M_Q = Q;                // Flux coming from the 3D
  _M_U=0.1+cos(2*pi*time);      // Forcing term

  y[0]=x[0]+(1/_M_C1*x[4]-1/_M_C1*x[7])*_M_dt;
  y[1]=x[1]+(-1/_M_C2*x[4]+1/_M_C2*x[5])*_M_dt;
  y[2]=x[2]+(-1/_M_C3*x[5]+1/_M_C3*x[6])*_M_dt;
  y[3]=x[3]+(1/_M_C4*x[7]-1/_M_C4*x[6])*_M_dt;
  y[4]=x[4]+(-1/_M_L5*x[0]+1/_M_L5*x[1]-_M_R5/_M_L5*x[4]+_M_U/_M_L5)*_M_dt;
  y[5]=x[5]+(-1/_M_L6*x[1]+1/_M_L6*x[2]-_M_R6/_M_L6*x[5]-_M_U/_M_L6)*_M_dt;
  y[6]=-_M_Q;  //The flux is now given from the 3D Model 
  y[7]=x[7]+(1/_M_L8*x[0]-1/_M_L8*x[3]-_M_R8/_M_L8*x[7])*_M_dt;
    
  for(UInt i=0; i<8; ++i){
    x[i]=y[i];
  };

  _M_outfile.open("res_Q.m", std::ios::app);
  _M_outfile << "     " << -_M_Q << ";" << std::endl;
  _M_outfile.close();

  _M_outfile.open("res_DP.m", std::ios::app);
  _M_outfile << "      " <<  y[3]-y[2] << ";" << std::endl;
  _M_outfile.close();

  return y[3]-y[2];
}


Real zeroDModelSolver::getQFromPressure(Real & time, Real & deltaP)
{ 
   Real _M_L_int  = 2*_M_L; // Induntance necessary on the interface for the flow rate problem
   Real _M_deltaP = deltaP; // Pressure coming from the 3D Model
   _M_U=0.1+cos(2*pi*time); // Forcing term

  y[0]=x[0]+(1/_M_C1*x[4]-1/_M_C1*x[7])*_M_dt;
  y[1]=x[1]+(-1/_M_C2*x[4]+1/_M_C2*x[5])*_M_dt;
  y[2]=x[2]+(-1/_M_C3*x[5]+1/_M_C3*x[6])*_M_dt;
  y[3]=x[3]+(1/_M_C4*x[7]-1/_M_C4*x[6])*_M_dt;
  y[4]=x[4]+(-1/_M_L5*x[0]+1/_M_L5*x[1]-_M_R5/_M_L5*x[4]+_M_U/_M_L5)*_M_dt;
  y[5]=x[5]+(-1/_M_L6*x[1]+1/_M_L6*x[2]-_M_R6/_M_L6*x[5]-_M_U/_M_L6)*_M_dt;
  y[6]=x[6] +(x[3]/_M_L_int-_M_deltaP/_M_L_int-x[2]/_M_L_int)*_M_dt; // P is given from the 3D
  y[7]=x[7]+(1/_M_L8*x[0]-1/_M_L8*x[3]-_M_R8/_M_L8*x[7])*_M_dt;
    
  for(UInt i=0; i<8; ++i){
    x[i]=y[i];
  };

  _M_outfile.open("res_Q.m", std::ios::app);
  _M_outfile << "     " << -y[6] << ";" << std::endl;
  _M_outfile.close();

  _M_outfile.open("res_DP.m", std::ios::app);
  _M_outfile << "      " <<  _M_deltaP << ";" << std::endl;
  _M_outfile.close();

  return -y[6];
}

bool zeroDModelSolver::isMeanPressProb()
{ 

  if (_M_isMeanPressProb == 1)
    return true;
  else
    return false;

}
 
}

#endif
