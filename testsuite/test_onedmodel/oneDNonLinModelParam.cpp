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
  \file oneDNonLinModelParam.cpp
  \author Vincent Martin
  \date 09/2004
  \version 1.0

  \brief File containing a class for the parameter

*/
#include "oneDNonLinModelParam.hpp"

#include <cmath>

namespace LifeV
{
 
OneDNonLinModelParam::OneDNonLinModelParam(const GetPot& dfile) :
  //! we suppose so far that all params are constant along the vessel
  //! otherwise replace _M_paramSize(1) by _M_paramSize(_M_dimDof)
  _M_paramSize(1),
  _M_Area0(_M_paramSize),
  _M_AlphaCoriolis(_M_paramSize),
  _M_PressBeta0(_M_paramSize),
  _M_PressBeta1(_M_paramSize),
  _M_FrictionKr(_M_paramSize)
{

  //-------------------------------------------
  //! Initialisation of the parameter variables
  //-------------------------------------------
  for (UInt iz = 0; iz < _M_paramSize ; iz ++ ) {
    _M_Area0[iz] = dfile("parameters/Area0",1.);
    _M_AlphaCoriolis[iz] =  dfile("parameters/alphaCor",1.);
    _M_PressBeta0[iz] = dfile("parameters/beta0",1.);
    _M_PressBeta1[iz] = dfile("parameters/beta1",0.5);
    _M_FrictionKr = dfile("parameters/Kr",1.);
  }
  _M_DensityRho = dfile("parameters/rho",1.);

}

//! so far we assume that all params are constant along the vessel
//! otherwise, replace [0] by [ii]
Real OneDNonLinModelParam::Area0(const UInt& ii) const {
  return _M_Area0[0]; 
}

Real OneDNonLinModelParam::AlphaCor(const UInt& ii) const {
  return _M_AlphaCoriolis[0];
}

Real OneDNonLinModelParam::Beta0(const UInt& ii) const {
  return _M_PressBeta0[0];
}

Real OneDNonLinModelParam::Beta1(const UInt& ii) const {
  return _M_PressBeta1[0];
}


Real OneDNonLinModelParam::FrictionKr(const UInt& ii) const {
  return _M_FrictionKr[0];
}

Real OneDNonLinModelParam::DensityRho() const {
  return _M_DensityRho;
}

//! initialisation from physical values
void OneDNonLinModelParam::initParam(const Real& Young_modulus) 
{  
  //  Real Young_modulus = 4.e6;
  Real thickness  = 0.065;
  Real reference_radius = 0.75;
  Real density   = 1.;  //???
  Real viscosity = 0.;  //???
  Real ksi       = 0.;  //???

  Real Coriolis_coeff = 1;

  //-------------------------------------------
  //! Initialisation of the parameter variables
  //-------------------------------------------
  for (UInt iz = 0; iz < _M_paramSize ; iz ++ ) {
    //! A0
    _M_Area0[iz] = M_PI*reference_radius * reference_radius;
    //! alpha
    _M_AlphaCoriolis[iz] = Coriolis_coeff;
    //! beta0
    _M_PressBeta0[iz] = thickness * Young_modulus * sqrt(M_PI) / 
      ( sqrt(_M_Area0[iz]) * (1 - ksi * ksi) );
  //! beta1 
    _M_PressBeta1[iz] = 0.5;
  //! Kr  
  _M_FrictionKr = 8 * M_PI * viscosity;
  }

  //! rho
  _M_DensityRho = density;

}

void OneDNonLinModelParam::showMeData(std::ostream& c) const
{
  //! parameters
  c << "\t[parameters]\n";
  c << "alphaCor = " << _M_AlphaCoriolis << "\n";
  c << "beta0 = " << _M_PressBeta0 << "\n";
  c << "beta1 = " << _M_PressBeta1 << "\n";
  c << "Kr       = " << _M_FrictionKr << "\n";
  c << "Area0    = " << _M_Area0 << "\n";
  c << "rho      = " << _M_DensityRho << "\n" << std::endl;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++


//+++++++++++++++++++++++++++++++++++++++++++++++++++

LinearSimpleParam::LinearSimpleParam(const GetPot& dfile) :
  //! we suppose so far that all params are constant along the vessel
  //! otherwise replace _M_paramSize(1) by _M_paramSize(_M_dimDof)
  _M_paramSize(1),
  _M_Flux11(_M_paramSize),
  _M_Flux12(_M_paramSize),
  _M_Flux21(_M_paramSize),
  _M_Flux22(_M_paramSize),
  _M_celerity1(_M_paramSize),
  _M_celerity2(_M_paramSize),
  _M_celer1_left_eigVector1(_M_paramSize),
  _M_celer1_left_eigVector2(_M_paramSize),
  _M_celer2_left_eigVector1(_M_paramSize),
  _M_celer2_left_eigVector2(_M_paramSize),
  _M_Source10(_M_paramSize),
  _M_Source20(_M_paramSize),
  _M_Source11(_M_paramSize),
  _M_Source12(_M_paramSize),
  _M_Source21(_M_paramSize),
  _M_Source22(_M_paramSize)
{

  //-------------------------------------------
  //! Initialisation of the parameter variables
  //-------------------------------------------
  for (UInt iz = 0; iz < _M_paramSize ; iz ++ ) {
    _M_Flux11[iz] = dfile("parameters/flux11",1.);
    _M_Flux12[iz] = dfile("parameters/flux12",0.);
    _M_Flux21[iz] = dfile("parameters/flux21",0.);
    _M_Flux22[iz] = dfile("parameters/flux22",1.);
    _M_celerity1[iz] = dfile("parameters/celer1",1.);
    _M_celerity2[iz] = dfile("parameters/celer2",1.);
    _M_celer1_left_eigVector1[iz] = dfile("parameters/left_eigvec11",1.);
    _M_celer1_left_eigVector2[iz] = dfile("parameters/left_eigvec12",0.);
    _M_celer2_left_eigVector1[iz] = dfile("parameters/left_eigvec21",0.);
    _M_celer2_left_eigVector2[iz] = dfile("parameters/left_eigvec22",1.);
    _M_Source10[iz] = dfile("parameters/source10",0.);
    _M_Source20[iz] = dfile("parameters/source20",0.);
    _M_Source11[iz] = dfile("parameters/source11",0.);
    _M_Source12[iz] = dfile("parameters/source12",0.);
    _M_Source21[iz] = dfile("parameters/source21",0.);
    _M_Source22[iz] = dfile("parameters/source22",0.);
  }

}

Real LinearSimpleParam::Flux11(const UInt& ii) const {
  return _M_Flux11[0];
}
Real LinearSimpleParam::Flux12(const UInt& ii) const {
  return _M_Flux12[0];
}
Real LinearSimpleParam::Flux21(const UInt& ii) const {
  return _M_Flux21[0];
}
Real LinearSimpleParam::Flux22(const UInt& ii) const {
  return _M_Flux22[0];
}

Real LinearSimpleParam::Celerity1(const UInt& ii) const {
  return _M_celerity1[0];
}

Real LinearSimpleParam::Celerity2(const UInt& ii) const {
  return _M_celerity2[0];
}

Real LinearSimpleParam::LeftEigenVector11(const UInt& ii) const {
  return _M_celer1_left_eigVector1[0];
}

Real LinearSimpleParam::LeftEigenVector12(const UInt& ii) const {
  return _M_celer1_left_eigVector2[0];
}

Real LinearSimpleParam::LeftEigenVector21(const UInt& ii) const {
  return _M_celer2_left_eigVector1[0];
}

Real LinearSimpleParam::LeftEigenVector22(const UInt& ii) const {
  return _M_celer2_left_eigVector2[0];
}

Real LinearSimpleParam::Source10(const UInt& ii) const {
  return _M_Source10[0];
}
Real LinearSimpleParam::Source20(const UInt& ii) const {
  return _M_Source20[0];
}
Real LinearSimpleParam::Source11(const UInt& ii) const {
  return _M_Source11[0];
}
Real LinearSimpleParam::Source12(const UInt& ii) const {
  return _M_Source12[0];
}
Real LinearSimpleParam::Source21(const UInt& ii) const {
  return _M_Source21[0];
}
Real LinearSimpleParam::Source22(const UInt& ii) const {
  return _M_Source22[0];
}

void LinearSimpleParam::showMeData(std::ostream& c) const
{
  //! parameters
  c << "\t[parameters]\n";
  c << "flux11      = " << _M_Flux11 << "\n";
  c << "flux12      = " << _M_Flux12 << "\n";
  c << "flux21      = " << _M_Flux21 << "\n";
  c << "flux22      = " << _M_Flux22 << "\n";
  c << "celer1      = " << _M_celerity1 << "\n";
  c << "celer2      = " << _M_celerity2 << "\n";
  c << "eigenvector11  = " << _M_celer1_left_eigVector1 << "\n";
  c << "eigenvector12  = " << _M_celer1_left_eigVector2 << "\n";
  c << "eigenvector21  = " << _M_celer2_left_eigVector1 << "\n";
  c << "eigenvector22  = " << _M_celer2_left_eigVector2 << "\n";
  c << "source10      = " << _M_Source10 << "\n";
  c << "source20      = " << _M_Source20 << "\n";
  c << "source11      = " << _M_Source11 << "\n";
  c << "source12      = " << _M_Source12 << "\n";
  c << "source21      = " << _M_Source21 << "\n";
  c << "source22      = " << _M_Source22 << "\n";
  c << std::endl;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++


}
