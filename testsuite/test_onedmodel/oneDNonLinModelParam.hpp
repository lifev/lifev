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
  \file oneDNonLinModelParam.hpp
  \author Vincent Martin
  \date 09/2004
  \version 1.0

  \brief File containing a class for the parameter

*/
#ifndef _PARAMETER_H_
#define _PARAMETER_H_
#include <string>
#include <iostream>
#include "lifeV.hpp"
#include "GetPot.hpp"

#include "RNM.hpp"

namespace LifeV
{
/*!
  \class OneDNonLinModelParam

  Class which holds all the physical parameters necessary for the
  1D Non Linear blood flow model.


  parameters: 
     Area0, alpha, beta0, beta1, Kr, rho.

  dA/dt + dQ/dz = 0
  dQ/dt + d/dz(alpha * Q^2/A) + A/rho * dP/dz + Kr * Q/A = 0

  with 
  P - P_ext = beta0 [ ( A / Area0 )^{beta1} - 1 ]

*/
class OneDNonLinModelParam
{
public : 
  
  //! constructor
  OneDNonLinModelParam(const GetPot& dfile);

  //! return the values 
  Real Area0(const UInt& ii) const;
  Real AlphaCor(const UInt& ii) const;
  Real Beta0(const UInt& ii) const;
  Real Beta1(const UInt& ii) const;
  Real FrictionKr(const UInt& ii) const;
  Real DensityRho() const;

  //! initialisation from physical values
  void initParam( const Real& Young_modulus );
  
  //! output
  void showMeData(std::ostream& c) const;

private :

  //! size of the parameter vectors (=1 if they are constant along the vessel)
  UInt _M_paramSize;

  //! reference area (often called A0) 
  KN<Real> _M_Area0;
  
  //! coriolis coefficient (often called alpha) 
  KN<Real> _M_AlphaCoriolis;

  /*!-----------------------------------------------------
    BETA0 AND BETA1
   -------------------------------------------------------
    pressure function factor and power parameters
    (often called beta0 (or beta) and beta1)

   BEWARE: there are at least 2 or 3 different ways of defining it!!!
   
   -------------------------------------------------------
   CONVENTIONS used here: 
   -------------------------------------------------------
   1/ Parameter homogeneous to a pressure:
   P - P_ext = PressBeta0 [ ( A / Area0 )^{PressBeta1} - 1 ]

   This PressBeta0 is homogeneous to a pressure. 
   In most cases PressBeta1 is taken equal to 1/2. 
   
   PressBeta0 = ( \sqrt{\pi} h_0 E ) / ( ( 1 - \ksi^2 ) * \sqrt{Area0} )
   -------------------------------------------------------
   (
   other conventions: (not used!!)
   1/ from Formaggia and Veneziani (p. 1.10, MOX report No 21 - june 2003)
   P - P_ext = \tilde{\beta_0} ( \sqrt{A} - \sqrt{A_0} ) / A_0 
   with 
   \beta0 = ( \sqrt{\pi} h_0 E ) / ( 1 - \ksi^2 )
   
   link with PressBeta0:

      \tilde{\beta_0} = PressBeta0 * \sqrt{A_0}

   2/ Auxiliary Parameter often used in the model1D code (by J-FG or D Lamponi)
      (ONLY when beta1=1/2 !!) 
   P - P_ext = 2 * rho * AuxBeta ( \sqrt{A} - \sqrt{A_0} )
   
   link with PressBeta0:

      AuxBeta = PressBeta0 * PressBeta1 / ( rho * Area0^(PressBeta1) ) 

        or whenever PressBeta1 = 1/2 :
      AuxBeta = PressBeta0 / ( 2 * rho * \sqrt{A_0} ) 
   )
   -------------------------------------------------------
  */

  //! P - P_ext = PressBeta0 [ ( A / Area0 )^{PressBeta1} - 1 ]
  KN<Real> _M_PressBeta0; //! homogeneous to a pressure
  KN<Real> _M_PressBeta1; //! power coeff (>0, often=1/2)


  // friction parameter Kr 
  KN<Real> _M_FrictionKr;

  //! density rho (always taken constant along the vessel)
  Real _M_DensityRho;

};

//+++++++++++++++++++++++++++++++++++++++++++++++++++
/*!
  \class LinearSimpleParam

  Class which holds all the physical parameters necessary for the
  1D  Linear hyperbolic equation.


  parameters: 
     F11, F12, F21, F22,
     celerity1, celerity2
  
     
  dU1/dt + F11 dU1/dz + F12 dU2/dz = 0
  dU2/dt + F21 dU1/dz + F22 dU2/dz = 0

  The flux matrix F = [F11, F12 ; F21 F22] has the eigenvalues

  celerity1, celerity2.


*/
class LinearSimpleParam {
public:
  //! constructor
  LinearSimpleParam(const GetPot& dfile);

  //! return the values 
  Real Flux11(const UInt& ii) const; 
  Real Flux12(const UInt& ii) const; 
  Real Flux21(const UInt& ii) const; 
  Real Flux22(const UInt& ii) const; 

  Real Celerity1(const UInt& ii) const;
  Real Celerity2(const UInt& ii) const;

  Real LeftEigenVector11(const UInt& ii) const;
  Real LeftEigenVector12(const UInt& ii) const;
  Real LeftEigenVector21(const UInt& ii) const;
  Real LeftEigenVector22(const UInt& ii) const;

  Real Source10(const UInt& ii) const; 
  Real Source20(const UInt& ii) const; 
  Real Source11(const UInt& ii) const; 
  Real Source12(const UInt& ii) const; 
  Real Source21(const UInt& ii) const; 
  Real Source22(const UInt& ii) const; 

  //! output
  void showMeData(std::ostream& c) const;
 
private :

  //! size of the parameter vectors (=1 if they are constant along the vessel)
  UInt _M_paramSize;

  //! flux matrix
  KN<Real> _M_Flux11;
  KN<Real> _M_Flux12;
  KN<Real> _M_Flux21;
  KN<Real> _M_Flux22;

  //! celerities of the linear problem (eigenvalues of the flux matrix)
  KN<Real> _M_celerity1;
  KN<Real> _M_celerity2;

  //! eigenvector for first eigenvalue
  KN<Real> _M_celer1_left_eigVector1;
  KN<Real> _M_celer1_left_eigVector2;
  //! eigenvector for second eigenvalue
  KN<Real> _M_celer2_left_eigVector1;
  KN<Real> _M_celer2_left_eigVector2;

  //! source matrix
  KN<Real> _M_Source10;
  KN<Real> _M_Source20;
  KN<Real> _M_Source11;
  KN<Real> _M_Source12;
  KN<Real> _M_Source21;
  KN<Real> _M_Source22;

};


}
#endif
