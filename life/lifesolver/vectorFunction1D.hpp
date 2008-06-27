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
#ifndef _VECTORFUNCTION1D_H_
#define _VECTORFUNCTION1D_H_

#include <life/lifecore/life.hpp>
#include <life/lifesolver/oneDNonLinModelParam.hpp>

namespace LifeV
{

class VectorFunction1D
{
public:

    //! do nothing constructor
    VectorFunction1D() {}

    virtual ~VectorFunction1D() = 0 ;

    //! i=1 -> F1, i=2 -> F2
    //! indz is the index position (not really beautiful...) VM 09/04
    virtual Real operator()(const Real& x1, const Real& x2,
                            const ID& i,
                            const UInt& indz = 0) const = 0;

    /*! Jacobian matrix dFi/dxj :
      diff(1,1) = dF1/dx1    diff(1,2) = dF1/dx2
      diff(2,1) = dF2/dx1    diff(2,2) = dF2/dx2
    */
    virtual Real diff(const Real& x1, const Real& x2,
                      const ID& i, const ID& j,
                      const UInt& indz = 0) const = 0;

    /*! Eigenvalues of the Jacobian matrix dFi/dxj :
      jacEigenValue(1) = lambda1
      jacEigenValue(2) = lambda2
    */
    virtual Real jacEigenValue(const Real& x1, const Real& x2,
                               const ID& i,
                               const UInt& indz = 0) const = 0;

    /*! Second derivative tensor d2Fi/(dxj dxk)
      diff2(1,1,1) = d2F1/dx1dx1    diff2(1,1,2) = d2F1/dx1dx2
      diff2(1,2,1) = d2F1/dx2dx1    diff2(1,2,2) = d2F1/dx2dx2

      diff2(2,1,1) = d2F2/dx1dx1    diff2(2,1,2) = d2F2/dx1dx2
      diff2(2,2,1) = d2F2/dx2dx1    diff2(2,2,2) = d2F2/dx2dx2

      with d2Fi/dx1dx2 = d2Fi/dx2dx1 .
    */
    virtual Real diff2(const Real& x1, const Real& x2,
                       const ID& i, const ID& j, const ID& k,
                       const UInt& indz = 0) const = 0;

};

//+++++++++++++++++++++++++++++++++++++++++++++++++++
/*!
  Class containing the non-linear flux function F
  of the 1D hyperbolic problem

  dU/dt + dF(U)/dz + B(U) = 0

  with U=[A,Q]^T
*/
class NonLinearFluxFun1D
{

private:
    //! the parameters that are used in the function
    const OneDNonLinModelParam& _M_oneDParam;

public:

  typedef std::pair< Real, Real > Vec2D;

    //! constructor
    NonLinearFluxFun1D(const OneDNonLinModelParam& onedparam);

    //! do nothing destructor
    // ~NonLinearFluxFun1D() {}

    /*! F = [Q,
      alpha*Q^2/A + beta0*beta1/(rho*(beta1+1)*A0^beta1) * A^(beta1+1) ]

      \param indz : is the index position for the parameters
      when they are space dependent.
      This is NOT pretty. I should try to
      remove this dependency. VM 09/04
    */
    Real operator()(const Real& _A, const Real& _Q,
                    const ID& ii,
                    const UInt& indz = 0) const ;

    //! Jacobian matrix Hij = dFi/dxj
    Real diff(const Real& _A, const Real& _Q,
              const ID& ii, const ID& jj,
              const UInt& indz = 0) const;

    /*! Eigenvalues and eigenvectors of the Jacobian matrix dFi/dxj
      \param eigi is the ith eigen value of the matrix dF/dx (i=1,2).
      \param lefteigvecij is the jth component of the left eigen vector
      associated to eigi. (i,j=1,2)
    */
    void jacobian_EigenValues_Vectors(const Real& _A, const Real& _Q,
                                      Real& eig1, Real& eig2,
                                      Real& lefteigvec11, Real& lefteigvec12,
                                      Real& lefteigvec21, Real& lefteigvec22,
                                      const UInt& indz = 0 ) const;

 //! Second derivative tensor d2Fi/(dxj dxk)
    Real diff2(const Real& _A, const Real& _Q,
               const ID& ii, const ID& jj, const ID& kk,
               const UInt& indz = 0) const;

};

/*!
  Class containing the non-linear source function B
  of the 1D hyperbolic problem

  dU/dt + dF(U)/dz + B(U) = 0

  with U=[A,Q]^T
*/
class NonLinearSourceFun1D
{

private:
    //! the parameters that are used in the function
    const OneDNonLinModelParam& _M_oneDParam;

public:

    //! constructor
    NonLinearSourceFun1D(const OneDNonLinModelParam& onedparam);

    //! do nothing destructor
    // ~NonLinearSourceFun1D() {}

    /*! B = [0, B2]^T

    with B2 such that
    B2 =   Kr*Q/A
    - beta1 * beta0/( rho*(beta1+1) ) * (A/A0)^(beta1+1)  * dA0/dz
    + A0/rho * [ 1/(beta1+1) * (A/A0)^(beta1+1) - A/A0 ]  * dbeta0/dz
    + A0    * beta0/( rho*(beta1+1) ) * (A/A0)^(beta1+1)
    * [ log(A/A0) - 1/(beta1+1) ]                       * dbeta1/dz

    \param indz : is the index position for the parameter
    */
    Real operator()(const Real& _A, const Real& _Q,
                    const ID& ii,
                    const UInt& indz = 0) const ;

    //! Jacobian matrix dBi/dxj
    Real diff(const Real& _A, const Real& _Q,
              const ID& ii, const ID& jj,
              const UInt& indz = 0) const;

    //! Second derivative tensor d2Bi/(dxj dxk)
    Real diff2(const Real& _A, const Real& _Q,
               const ID& ii, const ID& jj, const ID& kk,
               const UInt& indz = 0) const;

    /*! Sql = [Sql1, Sql2]^T

    Sql source term of the equation under its quasi-linear
    formulation :

    dU/dt + H(U) dU/dz + Sql(U) = 0

    \param indz : is the index position for the parameter
    */
    Real QuasiLinearSource(const Real& _U1, const Real& _U2,
                           const ID& ii,
                           const UInt& indz = 0) const ;

};
//+++++++++++++++++++++++++++++++++++++++++++++++++++


//+++++++++++++++++++++++++++++++++++++++++++++++++++
/*!
  Class containing the linear flux function F
  of the 1D hyperbolic problem

  dU/dt + dF(U)/dz + S(U) = 0


  with U=[U1, U2]^T and F(U) = [F11 U1 + F12 U2 ; F21 U1 + F22 U2]

  Fij are constant.

*/
class LinearSimpleFluxFun1D
{

private:
    //! the parameters that are used in the function
    const LinearSimpleParam& _M_oneDParam;

public:

    //! constructor
    LinearSimpleFluxFun1D(const LinearSimpleParam& onedparam);

    //! do nothing destructor
    // ~LinearSimpleFluxFun1D() {}

    /*! F = F(U) = [F11 U1 + F12 U2 , F21 U1 + F22 U2]^T

    \param indz : is the index position for the parameters
    when they are space dependent.
    This is NOT pretty. I should try to
    remove this dependency. VM 09/04
    */
    Real operator()(const Real& _U1, const Real& _U2,
                    const ID& ii,
                    const UInt& indz = 0) const ;

    //! Jacobian matrix dFi/dxj
    Real diff(const Real& _U1, const Real& _U2,
              const ID& ii, const ID& jj,
              const UInt& indz = 0) const;

    /*! Eigenvalues and eigenvectors of the Jacobian matrix dFi/dxj
      \param eigi is the ith eigen value of the matrix dF/dx (i=1,2).
      \param lefteigvecij is the jth component of the left eigen vector
      associated to eigi. (i,j=1,2)
    */
    void jacobian_EigenValues_Vectors(const Real& _U1, const Real& _U2,
                                      Real& eig1, Real& eig2,
                                      Real& lefteigvec11, Real& lefteigvec12,
                                      Real& lefteigvec21, Real& lefteigvec22,
                                      const UInt& indz = 0 ) const;

    //! Second derivative tensor d2Fi/(dxj dxk)
    Real diff2(const Real& _U1, const Real& _U2,
               const ID& ii, const ID& jj, const ID& kk,
               const UInt& indz = 0) const;
};

/*!
  Class containing the linear source function S
  of the 1D hyperbolic problem

  dU/dt + dF(U)/dz + S(U) = 0

  with U=[U1,U2]^T
*/
class LinearSimpleSourceFun1D
{

private:
    //! the parameters that are used in the function
    const LinearSimpleParam& _M_oneDParam;

public:

    //! constructor
    LinearSimpleSourceFun1D(const LinearSimpleParam& onedparam);

    //! do nothing destructor
    // ~LinearSimpleSourceFun1D() {}

    /*! S = [S1, S2]^T

    S1 = S10 + S11 U1 + S12 U2
    S2 = S20 + S21 U1 + S22 U2

    \param indz : is the index position for the parameter
    */
    Real operator()(const Real& _U1, const Real& _U2,
                    const ID& ii,
                    const UInt& indz = 0) const ;

    //! Jacobian matrix dBi/dxj
    Real diff(const Real& _U1, const Real& _U2,
              const ID& ii, const ID& jj,
              const UInt& indz = 0) const;

    //! Second derivative tensor d2Bi/(dxj dxk)
    Real diff2(const Real& _U1, const Real& _U2,
               const ID& ii, const ID& jj, const ID& kk,
               const UInt& indz = 0) const;

    /*! Sql = [Sql1, Sql2]^T

    Sql source term of the equation under its quasi-linear
    formulation :

    dU/dt + H(U) dU/dz + Sql(U) = 0

    Here H is constant w.r. to U.
    And Sql = S(U), because there is no variation of
    the coefficients.

    \param indz : is the index position for the parameter
    */
    Real QuasiLinearSource(const Real& _U1, const Real& _U2,
                           const ID& ii,
                           const UInt& indz = 0) const ;

};
//+++++++++++++++++++++++++++++++++++++++++++++++++++

}

#endif
