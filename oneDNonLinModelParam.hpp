/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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
#include <fstream>

#include <life/lifemesh/basicOneDMesh.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/tab.hpp>

#define M_ROBERTSON_CORRECTION 1 // 0.75

namespace LifeV
{

/*!
  \class BloodFlowParam

  Class which holds all the physical parameters necessary for the
  1D blood flow model.


  parameters:
  Area0, alpha, beta0, beta1, Kr, rho.

  Euler equations
  dA/dt + dQ/dz = 0
  dQ/dt + d/dz(alpha * Q^2/A) + A/rho * dP/dz + Kr * Q/A = 0

  with
  P - P_ext = beta0 [ ( A / Area0 )^{beta1} - 1 ]

*/
class BloodFlowParam
{
public :

    //! constructor
    BloodFlowParam(const GetPot& dfile, const std::string& section = "" );

    //! return the values
    Real Area0(const UInt& ii) const;
    Real Beta0(const UInt& ii) const;
    Real Beta1(const UInt& ii) const;
    Real dArea0dz(const UInt& ii) const;
    Real dBeta0dz(const UInt& ii) const;
    Real dBeta1dz(const UInt& ii) const;
    Real Celerity0(const UInt& ii) const;
    Real FrictionKr(const UInt& ii) const;
    Real DensityRho() const;
    Real DensityWall() const;
    Real Gamma() const;
    Real CoeffA() const;
    Real Thickness() const;
    Real XLeft() const;
    Real XRight() const;
    Real Length() const;
    UInt ParamSize() const;

    //! initialisation from physical values
    void initParam( const GetPot& dfile );

    //! compute the total pressure : Pt = P + rho/2 * (Q/A)^2
    Real totalPressure(const Real& _A, const Real& _Q,
                       const UInt& indz = 0) const;

    // compute the pressure (with viscoelastic term):
    // P = beta0 * ( ( _A / Area0 )^beta1 - 1 ) +
    //     + 1/(2*sqrt(pi*A)) * gamma * dA / dt
    Vector pressure(const Real& _A, const Real& _A_n, const Real& _A_nm1, const Real& dt,
                    const UInt& indz, const UInt& steps = 1, const bool& visco = 1, const bool& linearized = 1) const;

    //! compute the pressure : beta0 * ( ( _A / Area0 )^beta1 - 1 )
    Real pressure(const Real& _A, const UInt& indz = 0) const;

    //! compute the derivative of the pressure with respect to A
    Real pressureDiff(const Real& _A, const UInt& indz = 0) const;

    //! Compute A from P
    Real A_from_P(const Real& P, const UInt& indz=0) const;

    //! compute the derivative of total pressure with respect to A and Q
    Real totalPressureDiff( const Real& _A, const Real& _Q,
                            const ID& ii,
                            const UInt& indz = 0) const;

    //! These routines change the elastic modulus along the vessel
    /*!
       \brief Make the vessel stiffer on the left side of interval [xl, xr]

       When x < alpha - delta/2, the Young modulus is E * factor
       When x > alpha + delta/2, the Young modulus is E
       When alpha - delta/2 < x < alpha + delta/2, the Young modulus changes
        smoothly from the larger to the smaller value, according to a
        polynomial law of order n

       The grid size can be adapted (yesadaptive=1) in the nieghborhood of alpha,
       where the spatial derivative of the parameter will be maximum.
       However, the grid size is not allowed to be smaller than min_deltax

     */
    void stiffenVesselLeft( const Real& xl, const Real& xr,
                            const Real& factor, const Real& alpha,
                            const Real& delta, const Real& n,
                            const Real& min_deltax=1, const UInt& yesAdaptive=0 );

    /*!
       \brief Make the vessel stiffer on the right side of interval [xl, xr]

       \sa stiffenVesselLeft
    */
    void stiffenVesselRight( const Real& xl, const Real& xr,
                             const Real& factor, const Real& alpha,
                             const Real& delta, const Real& n,
                             const Real& min_deltax=1, const UInt& yesAdaptive=0  );

    //! output
    virtual void showMe(std::ostream& c = std::cout) const;

    //! destructor
    virtual ~BloodFlowParam(){}

protected :

    //! size of the parameter vectors (=1 if they are constant along the vessel)
    UInt _M_paramSize;

    //! length of the tube TP 1/05
    Real _M_length, _M_xleft, _M_xright;

    //! Resistence for the ResistiveLoad bc
    //      Real _M_terminalResistence;

    //! reference area (often called A0)
    Vector _M_Area0;
    Vector _M_dArea0dz;

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
    Vector _M_PressBeta0;    //! homogeneous to a pressure
    Vector _M_dPressBeta0dz; //! homogeneous to a pressure
    Vector _M_PressBeta1;    //! power coeff (>0, often=1/2)
    Vector _M_dPressBeta1dz; //! power coeff (>0, often=1/2)


    // friction parameter Kr
    Vector _M_FrictionKr;

    //! density rho (always taken constant along the vessel)
    Real _M_DensityRho;
    Real _M_DensityWall;
    Real _M_Thickness;
    Real _M_Gamma;
    Real _M_CoeffA;
    Real _M_PerfusionCoeff;
    Real _M_muscle;
    Real _M_muscle_theta;
    Real _M_muscle_period;
    Real _M_muscle_x1;
    Real _M_muscle_x2;
    Real _M_muscle_x3;
    Real _M_muscle_x4;
    Real _M_muscle_tau1;
    Real _M_muscle_tau2;
    Real _M_muscle_startT;
};



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
class OneDNonLinModelParam :
        public BloodFlowParam
{
public :

    //! constructor
    OneDNonLinModelParam(const GetPot& dfile, const std::string& section = "" );

    //! return the values
    Real AlphaCor(const UInt& ii) const;
    Real dAlphaCordz(const UInt& ii) const;
    //    Real Celerity0(const UInt& ii) const;

    //! Compute W from U
    void W_from_U( Real& _W1, Real& _W2,
                   const Real& _U1, const Real& _U2,
                   const UInt& indz ) const;

    //! Compute U from W
    void U_from_W( Real& _U1, Real& _U2,
                   const Real& _W1, const Real& _W2,
                   const UInt& indz ) const;

    //! compute the pressure as a function of (W1, W2)
    Real pressure_W(const Real& _W1, const Real& _W2, const UInt& indz = 0) const;

    /*! compute W1 or W2 given the pressure:
      W1 - W2 = 4 * sqrt( beta0 / (beta1 * rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
    */
    Real W_from_P(const Real& _P, const Real& _W, const ID& ii, const UInt& indz) const;

    Real W_from_Q(const Real& _Q, const Real& _W_n, const Real& _W, const ID& ii, const UInt& indz) const;

    //! compute the derivative of pressure with respect to W1 and W2
    Real pressure_WDiff( const Real& _W1, const Real& _W2,
                         const ID& ii,
                         const UInt& indz = 0) const;

    //! output
    void showMe(std::ostream& c = std::cout) const;

private :

    //! coriolis coefficient (often called alpha)
    Vector _M_AlphaCoriolis;
    Vector _M_dAlphaCoriolisdz;
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
    LinearSimpleParam(const GetPot& dfile );

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
    virtual void showMe(std::ostream& c = std::cout) const;

    //! destructor
    virtual ~LinearSimpleParam(){}

protected :

    //! size of the parameter vectors (=1 if they are constant along the vessel)
    UInt _M_paramSize;

    //! flux matrix
    Vector _M_Flux11;
    Vector _M_Flux12;
    Vector _M_Flux21;
    Vector _M_Flux22;

    //! celerities of the linear problem (eigenvalues of the flux matrix)
    Vector _M_celerity1;
    Vector _M_celerity2;

    //! eigenvector for first eigenvalue
    Vector _M_celer1_left_eigVector1;
    Vector _M_celer1_left_eigVector2;
    //! eigenvector for second eigenvalue
    Vector _M_celer2_left_eigVector1;
    Vector _M_celer2_left_eigVector2;

    //! source matrix
    Vector _M_Source10;
    Vector _M_Source20;
    Vector _M_Source11;
    Vector _M_Source12;
    Vector _M_Source21;
    Vector _M_Source22;

};



//+++++++++++++++++++++++++++++++++++++++++++++++++++
/*!
  \class LinearizedParam

  This class contains the parameters to solve a linearized
  blood flow problem.


  parameters:
  Area0, alpha, beta0, beta1, Kr, rho.

  dA/dt + dQ/dz = 0
  dQ/dt + A_0/rho * dP/dz + Kr/A_0 * Q = 0

  with
  P - P_ext = beta0 [ ( A / Area0 )^{beta1} - 1 ]

  which means
  dP / dz = beta0 beta1 ( A / Area0 )^{beta1 - 1} dA / dz
  \simeq beta0 beta1 dA / dz

*/
class LinearizedParam :
        public BloodFlowParam, public LinearSimpleParam
{
public:
    //! constructor
    LinearizedParam(const GetPot& dfile );

    //! return the values
    //    Real Celerity0(const UInt& /*ii*/) const;

    //! Compute W from U
    void W_from_U( Real& _W1, Real& _W2,
                   const Real& _U1, const Real& _U2,
                   const UInt& indz ) const;

    //! Compute U from W
    void U_from_W( Real& _U1, Real& _U2,
                   const Real& _W1, const Real& _W2,
                   const UInt& indz ) const;

    //! compute the pressure as a function of (W1, W2)
    Real pressure_W(const Real& _W1, const Real& _W2, const UInt& indz = 0) const;

    //! compute the derivative of pressure with respect to W1 and W2
    Real pressure_WDiff( const Real& _W1, const Real& _W2,
                         const ID& ii,
                         const UInt& indz = 0) const;

    /*! compute W1 or W2 given the pressure:
      W1 - W2 = 4 * sqrt( beta0 / (beta1 * rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
    */
    Real W_from_P(const Real& _P, const Real& _W, const ID& ii, const UInt& indz) const;

    Real W_from_Q(const Real& _Q, const Real& _W_n, const Real& _W, const ID& ii, const UInt& indz) const;

    //! output
    void showMe(std::ostream& c = std::cout) const;

};


}
#endif
