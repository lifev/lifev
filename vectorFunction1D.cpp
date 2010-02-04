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
#include "vectorFunction1D.hpp"

#include <cmath>

namespace LifeV
{

//---------------------------------------------
//! FLUX FUNCTION AND DERIVATIVES
//---------------------------------------------
NonLinearFluxFun1D::NonLinearFluxFun1D(const OneDNonLinModelParam& onedparam) :
    _M_oneDParam(onedparam)
{}

//! ii=1 -> F1, ii=2 -> F2
Real NonLinearFluxFun1D::operator()(const Real& _A, const Real& _Q,
                                    const ID& ii,
                                    const UInt& indz) const
{
    Real F2;
    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1;

    if( ii == 1 ) { //! F1
        return _Q;
    }
    if( ii == 2 ) { //! F2

        Area0 = _M_oneDParam.Area0(indz);
        alphaCor = M_ROBERTSON_CORRECTION * _M_oneDParam.AlphaCor(indz);
        beta0 = _M_oneDParam.Beta0(indz);
        beta1 = _M_oneDParam.Beta1(indz);
        rho   = _M_oneDParam.DensityRho();

        AoverA0 = _A / Area0;
        AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

        F2 = alphaCor * _Q * _Q / _A ;
        F2 += M_ROBERTSON_CORRECTION * beta0 * beta1 / ( rho * ( beta1 + 1) ) * _A * AoverA0POWbeta1;

        // added Robertson, Zakaria correction factor 3/4 (TP 05/07)
        //        return _M_Robertson_correction * F2;
        return F2;
    }
    ERROR_MSG("The flux function has only 2 components.");
    return -1.;
}

/*! Jacobian matrix dFi/dxj :
  diff(1,1) = dF1/dx1    diff(1,2) = dF1/dx2
  diff(2,1) = dF2/dx1    diff(2,2) = dF2/dx2
*/
Real NonLinearFluxFun1D::diff(const Real& _A, const Real& _Q,
                              const ID& ii, const ID& jj,
                              const UInt& indz) const
{
    Real dF2dA;
    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1;

    if( ii == 1 && jj == 1 ) { //! dF1/dA
        return 0.;
    }
    if( ii == 1 && jj == 2 ) { //! dF1/dQ
        return 1.;
    }
    if( ii == 2 && jj == 1 ) { //! dF2/dA

        Area0 = _M_oneDParam.Area0(indz);
        alphaCor = M_ROBERTSON_CORRECTION * _M_oneDParam.AlphaCor(indz);
        beta0 = _M_oneDParam.Beta0(indz);
        beta1 = _M_oneDParam.Beta1(indz);
        rho   = _M_oneDParam.DensityRho();

        AoverA0 = _A / Area0;
        AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

        dF2dA = - alphaCor * _Q * _Q / ( _A * _A );
        dF2dA += M_ROBERTSON_CORRECTION * beta0 * beta1 / rho * AoverA0POWbeta1;

        // added Robertson, Zakaria correction factor 3/4 (TP 05/07)
        return dF2dA;
    }
    if( ii == 2 && jj == 2 ) { //! dF2/dQ
        alphaCor = M_ROBERTSON_CORRECTION * _M_oneDParam.AlphaCor(indz);
        // added Robertson, Zakaria correction factor 3/4 (TP 05/07)
        return 2 * alphaCor * _Q / _A;
    }

    ERROR_MSG("Flux's differential function has only 4 components.");
    return -1.;
}

//! Eigenvalues and eigenvectors of the Jacobian matrix dFi/dxj :
void NonLinearFluxFun1D::
jacobian_EigenValues_Vectors(const Real& _A,
                             const Real& _Q,
                             Real&       eig1,
                             Real&       eig2,
                             Real&       lefteigvec11,
                             Real&       lefteigvec12,
                             Real&       lefteigvec21,
                             Real&       lefteigvec22,
                             const UInt& indz ) const
{
    Debug(6312) << "[NonLinearFluxFun1D]::jabocian_EigenValues_Vectors\n";

    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1, QoverA;

    Real celeralpha;

    Area0           = _M_oneDParam.Area0(indz - 1);
    alphaCor        = M_ROBERTSON_CORRECTION * _M_oneDParam.AlphaCor(indz - 1);
    beta0           = _M_oneDParam.Beta0(indz - 1);
    beta1           = _M_oneDParam.Beta1(indz - 1);
    rho             = _M_oneDParam.DensityRho();

    AoverA0         = _A / Area0;
    AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

    QoverA          = _Q / _A;

    celeralpha      = alphaCor * ( alphaCor - 1) * QoverA * QoverA;
    celeralpha     += M_ROBERTSON_CORRECTION * beta0 * beta1 / rho * AoverA0POWbeta1;
    celeralpha      = std::sqrt( celeralpha );

    /*
      std::cout << "\nArea in compute eigenvalues/vectors " << _A
      << "\nFlux in compute eigenvalues/vectors " << _Q
      << "\nArea0 in compute eigenvalues/vectors " << Area0
      << "\nalphaCor in compute eigenvalues/vectors " << alphaCor
      << "\nbeta0 in compute eigenvalues/vectors " << beta0
      << "\nbeta1 in compute eigenvalues/vectors " << beta1
      << "\nrho in compute eigenvalues/vectors " << rho
      << "\nAoverA0 in compute eigenvalues/vectors " << AoverA0
      << "\nAoverA0POWbeta1 in compute eigenvalues/vectors " << AoverA0POWbeta1
      << "\nQoverA in compute eigenvalues/vectors " << QoverA
      << "\nceleralpha in compute eigenvalues/vectors " << celeralpha
      << "\nRobertson_correction in compute eigenvalues/vectors " << M_ROBERTSON_CORRECTION
      << std::endl;
    */

    //! eigen values
    eig1 =   celeralpha + alphaCor * QoverA;
    eig2 = - celeralpha + alphaCor * QoverA;
    //! eigen vectors
    lefteigvec11 = - eig2 / _A;
    lefteigvec12 = 1. / _A;
    lefteigvec21 = - eig1 / _A;
    lefteigvec22 = 1. / _A;
}


/*! Second derivative tensor d2Fi/(dxj dxk)
  diff2(1,1,1) = d2F1/dx1dx1    diff2(1,1,2) = d2F1/dx1dx2
  diff2(1,2,1) = d2F1/dx2dx1    diff2(1,2,2) = d2F1/dx2dx2

  diff2(2,1,1) = d2F2/dx1dx1    diff2(2,1,2) = d2F2/dx1dx2
  diff2(2,2,1) = d2F2/dx2dx1    diff2(2,2,2) = d2F2/dx2dx2

  with d2Fi/dx1dx2 = d2Fi/dx2dx1 .
*/
Real NonLinearFluxFun1D::diff2(const Real& _A, const Real& _Q,
                               const ID& ii, const ID& jj, const ID& kk,
                               const UInt& indz) const
{
    Real d2F2dA2;
    Real Area0, alphaCor, beta0, beta1, rho;
    Real AoverA0, AoverA0POWbeta1divA;

    //! diff second of F1 is always 0.
    if( ii == 1 ) { //! d2F1/dUjdUk = 0.
        if( ( jj == 1 || jj == 2 ) && ( kk == 1 || kk == 2 ) ) {
            return 0.;
        }
    }
    if( ii == 2 ) {
        if( jj == 1 && kk == 1 ) { //! d2F2/dA2

            Area0 = _M_oneDParam.Area0(indz);
            alphaCor = M_ROBERTSON_CORRECTION * _M_oneDParam.AlphaCor(indz);
            beta0 = _M_oneDParam.Beta0(indz);
            beta1 = _M_oneDParam.Beta1(indz);
            rho   = _M_oneDParam.DensityRho();

            AoverA0 = _A / Area0;
            AoverA0POWbeta1divA = std::pow( AoverA0, beta1 ) / _A;

            d2F2dA2 = 2 * alphaCor * _Q * _Q / ( _A * _A * _A );
            d2F2dA2 += M_ROBERTSON_CORRECTION * beta0 * beta1 * beta1 / rho * AoverA0POWbeta1divA;

            return d2F2dA2;
        }
        //! cross terms (equal)
        if( (jj == 1 && kk == 2) || (jj == 2 && kk == 1) ) { //! d2F2/dAdQ=d2F2/dQdA
            alphaCor = M_ROBERTSON_CORRECTION * _M_oneDParam.AlphaCor(indz);
            return - 2 * alphaCor * _Q / ( _A * _A );
        }
        if( jj == 2 && kk == 2 ) { //! d2F2/dQ2
            alphaCor = M_ROBERTSON_CORRECTION * _M_oneDParam.AlphaCor(indz);
            return 2 * alphaCor / _A;
        }
    }
    ERROR_MSG("Flux's second differential function has only 8 components.");
    return -1.;
}

//---------------------------------------------
//! SOURCE FUNCTION AND DERIVATIVES
//---------------------------------------------
NonLinearSourceFun1D::NonLinearSourceFun1D(const OneDNonLinModelParam& onedparam) :
    _M_oneDParam(onedparam)
{}

//! i=1 -> F1, i=2 -> F2
Real NonLinearSourceFun1D::operator()(const Real& _A, const Real& _Q,
                                      const ID& ii,
                                      const UInt& indz) const
{
    Real B2(0.);
    Real Area0, beta0, beta1, rho, Kr;
    Real beta1plus1, AoverA0, AoverA0POWbeta1plus1;
    Real tmp;


    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;

    if( ii == 1 ) { //! B1
        return 0.;
    }
    if( ii == 2 ) { //! B2

        Area0 = _M_oneDParam.Area0(indz);
        beta0 = _M_oneDParam.Beta0(indz);
        beta1 = _M_oneDParam.Beta1(indz);
        dArea0dz = _M_oneDParam.dArea0dz(indz);
        dbeta0dz = _M_oneDParam.dBeta0dz(indz);
        dbeta1dz = _M_oneDParam.dBeta1dz(indz);
        Kr    = _M_oneDParam.FrictionKr(indz);
        rho   = _M_oneDParam.DensityRho();

        beta1plus1 = beta1 + 1;
        AoverA0 = _A / Area0;
        AoverA0POWbeta1plus1 = std::pow( AoverA0, beta1plus1 );

        tmp = beta0 / ( rho * beta1plus1 ) * AoverA0POWbeta1plus1;
        //! friction term
        B2 = Kr * _Q /_A;

        //! term with the derivative of A0 with respect to z
        B2 += - tmp * beta1 * dArea0dz;

        //! term with the derivative of beta0 with respect to z
        B2 += Area0 / rho * ( AoverA0POWbeta1plus1 / beta1plus1  -  AoverA0 )
            * dbeta0dz;

        //! term with the derivative of beta1 with respect to z
        B2 += tmp * Area0 * ( std::log( AoverA0 ) - 1. / beta1plus1 )
            * dbeta1dz;

        return M_ROBERTSON_CORRECTION * B2;
    }
    ERROR_MSG("The flux function has only 2 components.");
    return -1.;
}

//! Jacobian matrix dBi/dxj
Real NonLinearSourceFun1D::diff(const Real& _A, const Real& _Q,
                                const ID& ii, const ID& jj,
                                const UInt& indz) const
{
    Real dB2dA;
    Real Area0, beta0, beta1, rho, Kr;
    Real AoverA0, AoverA0POWbeta1;
    Real tmp;

    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;

    //! B1 = 0 so...
    if( ii == 1 ) {
        if( jj == 1 || jj == 2 ) { //! dB2/dUj = 0
            return 0.;
        }
    }
    if( ii == 2 && jj == 1 ) { //! dB2/dA

        Area0 = _M_oneDParam.Area0(indz);
        beta0 = _M_oneDParam.Beta0(indz);
        beta1 = _M_oneDParam.Beta1(indz);
        dArea0dz = _M_oneDParam.dArea0dz(indz);
        dbeta0dz = _M_oneDParam.dBeta0dz(indz);
        dbeta1dz = _M_oneDParam.dBeta1dz(indz);
        Kr    = _M_oneDParam.FrictionKr(indz);
        rho   = _M_oneDParam.DensityRho();

        AoverA0 = _A / Area0;
        AoverA0POWbeta1 = std::pow( AoverA0, beta1 );

        tmp = beta0 / rho * AoverA0POWbeta1;

        //! friction term
        dB2dA = - Kr * _Q / ( _A * _A );

        //! term with the derivative of A0 with respect to z
        dB2dA += - tmp * beta1 / Area0 * dArea0dz;

        //! term with the derivative of beta0 with respect to z
        dB2dA += ( AoverA0POWbeta1 - 1. ) / rho * dbeta0dz;

        //! term with the derivative of beta1 with respect to z
        dB2dA += tmp * std::log( AoverA0 ) * dbeta1dz;

        return M_ROBERTSON_CORRECTION * dB2dA;
    }
    if( ii == 2 && jj == 2 ) { //! dB2/dQ
        Kr    = _M_oneDParam.FrictionKr(indz);
        return M_ROBERTSON_CORRECTION * Kr / _A;
    }

    ERROR_MSG("Source's differential function has only 4 components.");
    return -1.;
}

//! Second derivative tensor d2Bi/(dxj dxk)
Real NonLinearSourceFun1D::diff2(const Real& _A, const Real& _Q,
                                 const ID& ii, const ID& jj, const ID& kk,
                                 const UInt& indz) const
{
    Real d2B2dA2;
    Real Area0, beta0, beta1, rho, Kr;
    Real AoverA0, AoverA0POWbeta1divA;
    Real tmp;

    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;

    //! B1 = 0 so ...
    if( ii == 1 ) { //! d2B1/dUjdUk = 0.
        if( ( jj == 1 || jj == 2 ) && ( kk == 1 || kk == 2 ) ) {
            return 0.;
        }
    }
    if( ii == 2 ) {
        if( jj == 1 && kk == 1 ) { //! d2B2/dA2
            //! this term is not strictly necessary as it is always multiplied by 0.
            //! but for the sake of generality...

            Area0 = _M_oneDParam.Area0(indz);
            beta0 = _M_oneDParam.Beta0(indz);
            beta1 = _M_oneDParam.Beta1(indz);
            dArea0dz = _M_oneDParam.dArea0dz(indz);
            dbeta0dz = _M_oneDParam.dBeta0dz(indz);
            dbeta1dz = _M_oneDParam.dBeta1dz(indz);
            Kr    = _M_oneDParam.FrictionKr(indz);
            rho   = _M_oneDParam.DensityRho();

            AoverA0 = _A / Area0;
            AoverA0POWbeta1divA = std::pow( AoverA0, beta1 ) / _A;

            tmp = beta0 / rho * AoverA0POWbeta1divA;

            //! friction term
            d2B2dA2 = 2 * Kr * _Q / ( _A * _A * _A);

            //! term with the derivative of A0 with respect to z
            d2B2dA2 += - tmp * beta1 * beta1 / Area0 * dArea0dz;

            //! term with the derivative of beta0 with respect to z
            d2B2dA2 += beta1 / rho * AoverA0POWbeta1divA * dbeta0dz;

            //! term with the derivative of beta1 with respect to z
            d2B2dA2 += tmp * ( beta1 * std::log( AoverA0 ) + 1. ) * dbeta1dz;

            return M_ROBERTSON_CORRECTION * d2B2dA2;
        }
        //! cross terms (equal)
        if( (jj == 1 && kk == 2) || (jj == 2 && kk == 1) ) {//! d2B2/dAdQ=d2B2/dQdA
            Kr    = _M_oneDParam.FrictionKr(indz);
            return - M_ROBERTSON_CORRECTION * Kr / ( _A * _A );
        }
        if( jj == 2 && kk == 2 ) { //! d2B2/dQ2
            return 0.;
        }
    }

    ERROR_MSG("Source's second differential function has only 8 components.");
    return -1.;
}

Real NonLinearSourceFun1D::
QuasiLinearSource(const Real& _A, const Real& _Q,
                  const ID& ii,
                  const UInt& indz) const
{
    Real Sql2;
    Real Area0, beta0, beta1, rho, Kr;
    Real AoverA0, AoverA0POWbeta1timesA;
    Real tmp;

    Real dArea0dz = 0.;
    Real dbeta0dz = 0.;
    Real dbeta1dz = 0.;
    Real dalphadz = 0.;

    if( ii == 1 ) { //! Sql1
        return 0.;
    }
    if( ii == 2 ) { //! Sql2

        Area0 = _M_oneDParam.Area0(indz);
        beta0 = _M_oneDParam.Beta0(indz);
        beta1 = _M_oneDParam.Beta1(indz);
        dArea0dz = _M_oneDParam.dArea0dz(indz);
        dbeta0dz = _M_oneDParam.dBeta0dz(indz);
        dbeta1dz = _M_oneDParam.dBeta1dz(indz);
        dalphadz = _M_oneDParam.dAlphaCordz(indz);
        Kr    = _M_oneDParam.FrictionKr(indz);
        rho   = _M_oneDParam.DensityRho();

        AoverA0               = _A/Area0;
        AoverA0POWbeta1timesA = std::pow( AoverA0, beta1 ) * _A;

        tmp  = beta0 / rho * AoverA0POWbeta1timesA;
        //! friction term
        Sql2 = Kr * _Q/_A;
        //! term with the derivative of A0 with respect to z
        Sql2 += - tmp * beta1 / Area0 * dArea0dz;
        //! term with the derivative of beta0 with respect to z
        Sql2 += 1 / rho * ( AoverA0POWbeta1timesA  -  _A ) * dbeta0dz;
        //! term with the derivative of beta1 with respect to z
        Sql2 += tmp * std::log( AoverA0 ) * dbeta1dz;
        //! term with the derivative of alpha with respect to z
        Sql2 += _Q * _Q / _A * dalphadz;

        return M_ROBERTSON_CORRECTION * Sql2;
    }
    ERROR_MSG("The QL source function has only 2 components.");
    return -1.;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++
}
