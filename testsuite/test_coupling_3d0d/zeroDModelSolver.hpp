/* -*- mode: c++ -*-
   This program is part of the LifeV library

   Author(s): Alexandra Moura <moura@mate.polimi.it>
   Date: 09/2004

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
    zeroDModelSolver(GetPot const& data_file);

    /*! Constructor
      \param zeroDModelSolver object
     */
    zeroDModelSolver(zeroDModelSolver const& network);

    //! Computes the pressure of the network, given the flux from 3D
    //! \param current time step
    //! \param flux from 3d
    Real getPressureFromQ(Real const& time, Real const& Q);

    //! Computes the pressure of the network, given the flux from 3D
    //! \param current time step
    //! \param flux from 3d
    Real getPressureFromQFSI(Real const& time, Real const& Qin, Real const& Qout);

    //! Computes the flux of the network, given the pression from 3D
    //! \param current time step
    //! \param pressure from 3D
    Real getQFromPressure(Real const& time, Real const& deltaP);

    // Ratio betwen the reference radius and the new radius due to the pressure 
    // (algebraic relationship)
    Real radius();

private:

    Real _M_Length;
    Real _M_Radius0;
    Real _M_Thickness;
    Real _M_dt;
    Real _M_rho;
    Real _M_nu;
    Real _M_E;
    Real _M_Poisson;

    Real _M_C1, _M_C2, _M_C3, _M_C4; // Capacities in the network
    Real _M_L5, _M_L6, _M_L7, _M_L8; // Induntances in the network
    Real _M_R5, _M_R6, _M_R7, _M_R8; // Resitences in the network
    Real _M_U;                       // Forcing term

    Real _M_L, _M_R, _M_C;           // Parameters for the electrical net analog of the 3D
                                     // Model (here substituted for the 3D Model itself)

    Real _M_L_int;                   // For the flux imposed problem
                                     // Induntance necessary on the interface for the flow rate problem
    Real x[9], y[9];

    Real _M_Q;
    Real _M_deltaP;
    Real _M_Qin, _M_Qout;

    Real _M_Radius;

    std::ofstream _M_outfile;

};

//***************************************************************************************//
//                                     IMPLEMENTATION                                    //
//***************************************************************************************//

zeroDModelSolver::zeroDModelSolver(GetPot const& data_file):

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
    _M_Length     = data_file( "fluid/geometry/length", 0. );
    _M_Radius0    = data_file( "fluid/geometry/radius", 1. );
    _M_Thickness  = data_file( "fluid/geometry/thickness", 0. );

    _M_dt      = data_file("fluid/discretization/timestep", 0. );
    _M_rho     = data_file("fluid/physics/density", 1. );
    _M_nu      = data_file("fluid/physics/viscosity", 1. );
    _M_E       = data_file("solid/physics/young", 1. );
    _M_Poisson = data_file("solid/physics/poisson", 1. );

    if (_M_Thickness!=0)
        _M_C=3*pi*pow(_M_Radius0,3)*_M_Length/(2*_M_E*_M_Thickness);
    _M_R= _M_rho*8*pi*_M_nu*_M_Length/(pow(pi*_M_Radius0*_M_Radius0,2));
    _M_L= (_M_rho*_M_Length)/(pi*pow(_M_Radius0,2));

    x[0]=0.0178;
    x[1]=0.3740;
    x[2]=0.3346;
    x[3]=-0.0309;
    x[4]=-0.0908;
    x[5]=-0.0011;
    x[6]=0.0211;
    x[7]=0.0257;
    x[8]=0.0211;

}

zeroDModelSolver::zeroDModelSolver(zeroDModelSolver const& __network):

    _M_Length(__network._M_Length),
    _M_Radius0(__network._M_Radius),
    _M_Thickness(__network._M_Thickness),

    _M_dt(__network._M_dt),
    _M_rho(__network._M_rho),
    _M_nu(__network._M_nu),
    _M_E(__network._M_E),
    _M_Poisson(__network._M_Poisson),

    _M_C1(0.05),
    _M_C2(0.5),
    _M_C3(0.02),
    _M_C4(0.001),
    _M_L5(0.5),
    _M_L6(0.1),
    _M_L8(0.1),
    _M_R5(5.),
    _M_R6(6.),
    _M_R8(5.001),

    _M_L(__network._M_L),
    _M_R(__network._M_R),
    _M_C(__network._M_C)

{
    x[0]=0.0178;
    x[1]=0.3740;
    x[2]=0.3346;
    x[3]=-0.0309;
    x[4]=-0.0908;
    x[5]=-0.0011;
    x[6]=0.0211;
    x[7]=0.0257;
    x[8]=0.0211;

}


Real zeroDModelSolver::getPressureFromQ(Real const& time, Real const& Q)
{
    _M_Q = Q;                     // Flux coming from the 3D
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

    return y[3]-y[2];
}


Real zeroDModelSolver::getQFromPressure(Real const& time, Real const& deltaP)
{
    _M_L_int  = 2*_M_L;       // Induntance necessary on the interface for the flow rate problem
    _M_deltaP = deltaP;       // Pressure coming from the 3D Model
    _M_U=0.1+cos(2*pi*time);  // Forcing term

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

    return -y[6];
}

Real zeroDModelSolver::getPressureFromQFSI(Real const& time, Real const& Qin, Real const& Qout)
{
    _M_Qin = Qin;              // Flux coming from the 3D at inflow (flag=2)
    _M_Qout = Qout;            // Flux coming from the 3D at outflow (flag=3)
    _M_U=0.1+cos(2*pi*time);   // Forcing term

    y[0]=x[0]+(1/_M_C1*x[4]-1/_M_C1*x[7])*_M_dt;
    y[1]=x[1]+(-1/_M_C2*x[4]+1/_M_C2*x[5])*_M_dt;
    y[2]=x[2]+(-1/_M_C3*x[5]+1/_M_C3*x[6])*_M_dt;
    y[3]=x[3]+(1/_M_C4*x[7]-1/_M_C4*x[6])*_M_dt;
    y[4]=x[4]+(-1/_M_L5*x[0]+1/_M_L5*x[1]-_M_R5/_M_L5*x[4]+_M_U/_M_L5)*_M_dt;
    y[5]=x[5]+(-1/_M_L6*x[1]+1/_M_L6*x[2]-_M_R6/_M_L6*x[5]-_M_U/_M_L6)*_M_dt;
    y[6]=-_M_Qin;  //The inflow flux is now given from the 3D Model
    y[7]=x[7]+(1/_M_L8*x[0]-1/_M_L8*x[3]-_M_R8/_M_L8*x[7])*_M_dt;
    y[8]=_M_Qout;  //The outflow flux is now given from the 3D Model

    for(UInt i=0; i<9; ++i){
        x[i]=y[i];
    };

    _M_deltaP = y[3]-y[2];

    _M_Radius = _M_Radius0 + _M_deltaP*(pow(_M_Radius0,2)*(1-pow(_M_Poisson,2))/(_M_Thickness*_M_E));

    return _M_deltaP;
}

Real
zeroDModelSolver::radius()
{
    return _M_Radius/_M_Radius0;
}

}

#endif
