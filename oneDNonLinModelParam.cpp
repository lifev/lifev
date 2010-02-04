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
  \file oneDNonLinModelParam.cpp
  \author Vincent Martin
  \date 09/2004
  \version 1.0

  \brief File containing a class for the parameter

*/
#include <lifemc/lifesolver/oneDNonLinModelParam.hpp>

#include <cmath>

namespace LifeV
{

BloodFlowParam::BloodFlowParam(const GetPot& dfile, const std::string& section) :
        _M_paramSize      (dfile((section + "/discretization/nb_elem").data(), 0 ) + 1),
        _M_length         (0.),
        _M_xleft          (dfile((section + "/discretization/x_left").data(), 0. )),
        _M_xright         (dfile((section + "/discretization/x_right").data(), 0. )),
        _M_Area0          (_M_paramSize),
        _M_dArea0dz       (_M_paramSize),
        _M_PressBeta0     (_M_paramSize),
        _M_dPressBeta0dz  (_M_paramSize),
        _M_PressBeta1     (_M_paramSize),
        _M_dPressBeta1dz  (_M_paramSize),
        _M_FrictionKr     (_M_paramSize)
{
    _M_length = _M_xright - _M_xleft;


    Real _A0       = dfile((section + "/parameters/Area0").data(), M_PI);
    Real _beta0    = dfile((section + "/parameters/beta0").data(), 1.e6);
    Real _beta1    = dfile((section + "/parameters/beta1").data(), 0.5);
    Real _kr       = dfile((section + "/parameters/Kr").data(),    1.);
    Real _dA0dz    = 0.;
    Real _dbeta0dz = 0.;
    Real _dbeta1dz = 0.;


    std::cout << "NB ELEM = " << _M_paramSize << std::endl;
    //-------------------------------------------
    //! Initialisation of the parameter variables
    //-------------------------------------------
//     _M_Area0         = Vector( _M_paramSize, (double) _A0 );
//     _M_dArea0dz      = Vector( _M_paramSize, _dA0dz );
//     _M_PressBeta0    = Vector( _M_paramSize, _beta0 );
//     _M_dPressBeta0dz = Vector( _M_paramSize, _dbeta0dz );
//     _M_PressBeta1    = Vector( _M_paramSize, _beta1 );
//     _M_dPressBeta1dz = Vector( _M_paramSize, _dbeta1dz );
//     _M_FrictionKr    = Vector( _M_paramSize, _kr );

    for (int ii = 0; ii < _M_paramSize; ++ii)
      {
          _M_Area0[ii]         = _A0;
          _M_dArea0dz[ii]      = _dA0dz;
          _M_PressBeta0[ii]    = _beta0;
          _M_dPressBeta0dz[ii] = _dbeta0dz;
          _M_PressBeta1[ii]    = _beta1;
          _M_dPressBeta1dz[ii] = _dbeta1dz;
          _M_FrictionKr[ii]    = _kr;
      }

    _M_DensityRho  = dfile("parameters/rho",1.);

    _M_DensityWall = dfile("parameters/rho_w",1.);

    _M_Thickness   = dfile("parameters/thickness",0.);

    _M_Gamma       = dfile("parameters/gamma",0.);

    _M_CoeffA      = dfile("parameters/coeffA",0.);

    if( dfile("parameters/use_physical_values",false) )
        {
            Debug( 6320 ) << "[BloodFlowParam] initializing from physical values\n";
            initParam(dfile);
        }

    //          Debug( 6320 ) << "[BloodFlowParam] Robertson correction coefficient = "
    //                    << M_ROBERTSON_CORRECTION << "\n";
}

Real BloodFlowParam::Area0(const UInt& ii) const {
    return _M_Area0[ii];
}

Real BloodFlowParam::Beta0(const UInt& ii) const {
    return _M_PressBeta0[ii];
}

Real BloodFlowParam::Beta1(const UInt& ii) const {
    return _M_PressBeta1[ii];
}

Real BloodFlowParam::dArea0dz(const UInt& ii) const {
    return _M_dArea0dz[ii];
}

Real BloodFlowParam::dBeta0dz(const UInt& ii) const {
    return _M_dPressBeta0dz[ii];
}

Real BloodFlowParam::dBeta1dz(const UInt& ii) const {
    return _M_dPressBeta1dz[ii];
}

Real BloodFlowParam::Celerity0(const UInt& ii) const {
    Real celerity0( _M_PressBeta0[ii] * _M_PressBeta1[ii] / _M_DensityRho );
    return std::sqrt( celerity0 );
}

Real BloodFlowParam::FrictionKr(const UInt& ii) const {
    return _M_FrictionKr[ii];
}

Real BloodFlowParam::DensityRho() const {
    return _M_DensityRho;
}

Real BloodFlowParam::DensityWall() const {
    return _M_DensityWall;
}

Real BloodFlowParam::Thickness() const {
    return _M_Thickness;
}

Real BloodFlowParam::Gamma() const {
    return _M_Gamma;
}

Real BloodFlowParam::CoeffA() const {
    return _M_CoeffA;
}

Real BloodFlowParam::XLeft() const {
    return _M_xleft;
}

Real BloodFlowParam::XRight() const {
    return _M_xright;
}

Real BloodFlowParam::Length() const {
    return _M_length;
}

UInt BloodFlowParam::ParamSize() const {
    return _M_paramSize;
}

//! initialisation from physical values
void BloodFlowParam::initParam(const GetPot& dfile)
{


    Real Young_modulus    = dfile("1d_physics/young"          , 5.e6);
    Real thickness        = dfile("1d_physics/thickness"      , 0.05);
    Real reference_radius = dfile("1d_physics/radius"         , 1.);
    Real viscosity        = dfile("1d_physics/viscosity"      , 0.035);  //???
    Real ksi              = dfile("1d_physics/ksi"            , 0.5);  //???
    Real friction_factor  = dfile("1d_physics/friction_factor", 8.);  //???
    bool thick_vessel     = dfile("1d_physics/thick_vessel"   , 0);

    Real _A0( M_PI*reference_radius*reference_radius );

    Real _beta0, _beta1;
    if( thick_vessel ){ // see Amadori, Ferrari, Formaggia (MOX report 86)
        //! beta0
        _beta0 = - thickness*Young_modulus*sqrt(M_PI) /
            ( sqrt(_A0) *
              ( (1 - ksi * ksi)
                + ksi * (1 + ksi) * (thickness * sqrt(M_PI)
                                     / sqrt(_A0) )
                )
              );
        //! beta1
        _beta1 = - 0.5;
    }
    else{
        //! beta0
        _beta0 = thickness * Young_modulus * sqrt(M_PI) /
            ( sqrt(_A0) * (1 - ksi * ksi) );
        //! beta1
        _beta1 = 0.5;
    }

    Real _kr( friction_factor * M_PI * viscosity );

    //-------------------------------------------
    //! Initialisation of the parameter variables
    //-------------------------------------------
    //! A0
    _M_Area0      = Vector( _M_paramSize );
    //! beta0
    _M_PressBeta0 = Vector( _M_paramSize );
    //! beta1
    _M_PressBeta1 = Vector( _M_paramSize );
    //! Kr
    _M_FrictionKr = Vector( _M_paramSize );

    for (int ii = 0; ii < _M_paramSize; ++ii)
      {
          _M_Area0[ii]         = _A0;
          _M_PressBeta0[ii]    = _beta0;
          _M_PressBeta1[ii]    = _beta1;
          _M_FrictionKr[ii]    = _kr;
      }

    _M_Thickness  = thickness;

}

// compute the pressure (with viscoelastic term):
// P = beta0 * ( ( _A / Area0 )^beta1 - 1 ) +
//     + 1/(2*sqrt(pi*A)) * gamma * dA / dt
Vector BloodFlowParam::
pressure(const Real& _A,
         const Real& _A_n,
         const Real& _A_nm1,
         const Real& dt,
         const UInt& indz,
         const UInt& steps,
         const bool& visco,
         const bool& linearized ) const
{
    Vector _a(2), _b(2), _c(2), area(2), result(4);

    _a(0) =  1.;
    _b(0) = -1.;
    _c(0) =  0.;

    _a(1) =  3./2.;
    _b(1) = -2.;
    _c(1) =  1./2.;

    area(0) = _A;
    area(1) = _M_Area0[indz];

    Real _pi( 4*std::atan(1) );

    result(3) = ( _a(steps) * _A + _b(steps) * _A_n + _c(steps) * _A_nm1 ) / dt; //> dA/dt

    result(2) = _M_Gamma * result(3) / ( 2*sqrt(_pi*area(linearized)) );               //> visc_component

    result(1) = pressure( _A, indz );                                                  //> elast_component

    result(0) = result(1) + visco * result(2);

    return result;
}

// compute the pressure : beta0 * ( ( _A / Area0 )^beta1 - 1 )
Real BloodFlowParam::
pressure(const Real& _A, const UInt& indz) const
{
    return ( _M_PressBeta0[indz]
             * ( std::pow( _A/_M_Area0[indz], _M_PressBeta1[indz] ) - 1 ) );
}

/*! Derivative of pressure as a function of A
  dP(A)/dA = beta1 * beta0 * ( _A / Area0 )^beta1 / A
*/
Real BloodFlowParam::
pressureDiff(const Real& _A, const UInt& indz) const
{
    Real AoverA0POWbeta1( std::pow( _A / _M_Area0[indz], _M_PressBeta1[indz] ) );

    //std::cout << indz << " -> " <<  _M_PressBeta0[indz] << " " <<  _M_PressBeta1[indz] << " " <<  AoverA0POWbeta1 << " " << _A << std::endl;
    return _M_PressBeta0[indz] * _M_PressBeta1[indz] * AoverA0POWbeta1 / _A;
}

//! compute area given the pressure: A = A0 * ( P / beta0 + 1 )^(1/beta1)
//! To be used in initialization! when time derivative of A is supposed null
Real BloodFlowParam::
A_from_P(const Real& _P, const UInt& indz) const
{
    return ( _M_Area0[indz]
             * std::pow( _P / _M_PressBeta0[indz] + 1, 1/_M_PressBeta1[indz] )  );
}


/*! Total pressure (used for interface conditions)
  Pt = P + rho/2 * (Q/A)^2
*/
Real BloodFlowParam::
totalPressure(const Real& _A, const Real& _Q,
              const UInt& indz) const
{
    Real vel( _Q / _A );

    return ( pressure( _A, indz ) + (_M_DensityRho/2) * vel * vel );
}

/*! Derivative of Total pressure (used for interface conditions)
  dPt/dU_ii = dP/dU_ii + rho/2 * d(Q/A)^2/dU_ii
*/
Real BloodFlowParam::
totalPressureDiff( const Real& _A, const Real& _Q,
                   const ID& ii,
                   const UInt& indz) const
{
    Real vel( _Q / _A );

    if( ii == 1 ) { //! dPt/dA

        Real dPtdA( pressureDiff( _A, indz )
                    - _M_DensityRho * vel * vel / _A );

        return dPtdA;
    }
    if( ii == 2 ) { //! dPt/dQ
        return ( _M_DensityRho * vel / _A );
    }
    ERROR_MSG("Total pressure's differential function has only 2 components.");
    return -1.;
}

void BloodFlowParam::stiffenVesselLeft( const Real& xl, const Real& xr,
                                        const Real& factor, const Real& alpha,
                                        const Real& delta, const Real& n,
                                        const Real& min_deltax, const UInt& yesAdaptive )
{
    /* Stiffen Left boundary with a fifth order polynomial law

       if (alpha-delta/2) <= x < alpha
       coeff = ( (alpha + delta/2) - x )^5 * 2^4 / delta^5;

       if alpha <= x <= alpha + delta/2
       coeff = 1 - ( x - (alpha - delta/2))^5 * 2^4 / delta^5;

    */
    if (yesAdaptive)
        {
            Real ratio, n_elem_delta,n_elem_r,n_elem_l;

            UInt iz=0, alpha_iz;

            //      alpha_iz = static_cast<UInt>( alpha / (xr-xl) * static_cast<Real>( _M_paramSize-1 ) );
            alpha_iz = static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) +
                ( (_M_paramSize - 1) -
                  static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) ) -
                  static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) ) / 2;

            //      n_elem_r = static_cast<Real>( (_M_paramSize-1) - alpha_iz );
            n_elem_r = ( (_M_paramSize-1) - alpha_iz ) -
                static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) );

            //      n_elem_l = static_cast<Real>( alpha_iz );
            n_elem_l = alpha_iz -
                static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) );

            n_elem_delta = static_cast<Real>(_M_paramSize - 1) / (xr - xl) * delta;

            //      n_elem_delta = n_elem_r + n_elem_l;


            Real x_current,deltax,deltax_adaptive,deltax_uniform;

            x_current = alpha;

            do
                {
                    //! beta0
                    // fifth order
                    ratio=(( (alpha + delta/2) - x_current ) / delta);

                    _M_dPressBeta0dz[alpha_iz+iz] = _M_PressBeta0[alpha_iz+iz] *
                        ( factor * (- n / delta) *
                          ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
                    _M_dPressBeta0dz[alpha_iz-iz] = _M_dPressBeta0dz[alpha_iz+iz];

                    _M_PressBeta0[alpha_iz+iz] = _M_PressBeta0[alpha_iz+iz] *
                        ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n)));
                    _M_PressBeta0[alpha_iz-iz] = _M_PressBeta0[alpha_iz+iz] /
                        ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n)))
                        * ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) );

                    // first order
                    //        _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] * ( -factor * (n / delta) *
                    //                       ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
                    //_M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor * ratio );


                    deltax_adaptive = ( -1/n_elem_delta ) *
                        ( 1 / ( (-n/delta) * pow(2,(n-1)) *
                                pow( ratio , (n-1) )
                                )
                          );

                    deltax_uniform = ( (alpha+delta/2) - x_current) / ( n_elem_r - iz );

                    iz++;

                    deltax = ( ( deltax_adaptive < deltax_uniform ) && (iz < n_elem_r) )
                        ? deltax_adaptive : deltax_uniform;

                    //( xr - xl ) / _M_edgeList.size();
                    ASSERT_PRE( deltax > 0 ,
                                "The left point is on the right..." );

                    x_current += deltax;

                }
            while ( ( x_current < ( alpha + delta/2 ) ) && ( (alpha_iz - (iz - 1)) > 0) );

            if ( ( alpha_iz - (iz - 1)) > 0)
                {
                    do
                        {
                            _M_PressBeta0[alpha_iz-iz] = _M_PressBeta0[alpha_iz-iz] *
                                ( 1 + factor );
                            iz++;
                        }
                    while ( (alpha_iz - (iz - 1)) > 0 );

                    //      _M_PressBeta0[0] = _M_PressBeta0[0] *
                    //  ( 1 + factor );
                }
            else
                std::cout << "[stiffenVesselRight] error! out of left boundary" << std::endl;
        }

    else
        {
            UInt iz=0;

            Real ratio, x_current=xl, deltax;

            deltax=(xr-xl)/static_cast<Real>(_M_paramSize-1);

            while ( (x_current < (alpha - delta/2)) && (iz < _M_paramSize) )
                {
                    _M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor );
                    iz++;
                    x_current+=deltax;
                }

            while ( (x_current < alpha) && (iz < _M_paramSize) )
                {
                    ratio=(( x_current - (alpha-delta/2) ) / delta);

                    _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] *
                        ( factor * (- n / delta) *
                          ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );

                    _M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor *
                                                              ( 1 - pow(2,(n-1)) * pow(ratio,n)) );
                    iz++;
                    x_current+=deltax;
                }

            while ( (x_current < (alpha+delta/2)) && (iz < _M_paramSize) )
                {
                    ratio=(( (alpha+delta/2) - x_current ) / delta);

                    _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] *
                        ( factor * ( -n / delta) *
                          ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );

                    _M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor *
                                                              ( pow(2,(n-1)) * pow(ratio,n)) );
                    iz++;
                    x_current+=deltax;
                }
        }

}

void BloodFlowParam::stiffenVesselRight( const Real& xl, const Real& xr,
                                         const Real& factor, const Real& alpha,
                                         const Real& delta, const Real& n,
                                         const Real& min_deltax, const UInt& yesAdaptive )
{
  Debug( 6320 ) << "stiffenVesselright ...\n";
    /* Stiffen Left boundary with a fifth order polynomial law

       if (alpha-delta/2) <= x < alpha
       coeff = ( x - (alpha - delta/2) )^5 * 2^4 / delta^5;

       if alpha <= x <= alpha + delta/2
       coeff = 1 + ( x - (alpha + delta/2))^5 * 2^4 / delta^5;

    */
    if (yesAdaptive)
        {
            Real ratio, n_elem_delta,n_elem_r,n_elem_l;

            UInt iz=0, alpha_iz;

            //      alpha_iz = static_cast<UInt>( alpha / (xr-xl) * ( static_cast<Real>( _M_paramSize-1 ) ) );
            alpha_iz = static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) +
                ( (_M_paramSize - 1) -
                  static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) ) -
                  static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) ) ) / 2;

            n_elem_delta = static_cast<Real>(_M_paramSize - 1) / (xr - xl) * delta;

            //      n_elem_r = static_cast<Real>( (_M_paramSize-1) - alpha_iz );
            n_elem_r = ( (_M_paramSize-1) - alpha_iz ) -
                static_cast<int>( std::floor( (xr - (alpha+delta/2)) / min_deltax + 0.5 ) );

            //      n_elem_l = static_cast<Real>( alpha_iz );
            n_elem_l = alpha_iz -
                static_cast<int>( std::floor( (alpha-delta/2) / min_deltax + 0.5 ) );

            Real x_current,deltax,deltax_adaptive,deltax_uniform;

            x_current = alpha;

            do
                {
                    //! beta0
                    // fifth order
                    ratio=(( (alpha + delta/2) - x_current) / delta);

                    _M_dPressBeta0dz[alpha_iz+iz] = _M_PressBeta0[alpha_iz+iz] *
                        ( factor * ( n / delta) *
                          ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );

                    _M_dPressBeta0dz[alpha_iz-iz] = _M_dPressBeta0dz[alpha_iz+iz];

                    _M_PressBeta0[alpha_iz+iz] = _M_PressBeta0[alpha_iz+iz] *
                        ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) );

                    _M_PressBeta0[alpha_iz-iz] = _M_PressBeta0[alpha_iz+iz] /
                        ( 1 + factor * ( 1 - ( pow(2,(n-1)) * pow(ratio,n) ) ) )
                        * ( 1 + factor * ( pow(2,(n-1)) * pow(ratio,n) ) );

                    // first order
                    //        _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] * ( -factor * (n / delta) *
                    //                       ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
                    //_M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor * ratio );

                    deltax_adaptive = ( -1/n_elem_delta ) *
                        ( 1 / ( (-n/delta) * pow(2,(n-1)) *
                                pow( ratio , (n-1) )
                                )
                          );

                    deltax_uniform = ( (alpha+delta/2) - x_current) / ( n_elem_r - iz );

                    iz++;

                    deltax = ( ( deltax_adaptive < deltax_uniform ) && (iz < n_elem_r) )
                        ? deltax_adaptive : deltax_uniform;

                    //( xr - xl ) / _M_edgeList.size();
                    ASSERT_PRE( deltax > 0 ,
                                "The left point is on the right..." );

                    x_current += deltax;

                }
            while ( x_current < ( alpha + delta/2 ) && ( (alpha_iz - (iz - 1)) > 0) );

            if ( ( alpha_iz + iz ) <= (_M_paramSize -1) )
                {
                    do
                        {
                            _M_PressBeta0[alpha_iz+iz] = _M_PressBeta0[alpha_iz+iz] *
                                ( 1 + factor );
                            iz++;
                        }
                    while ( (alpha_iz + iz - 1) < (_M_paramSize -1) );

                    //      _M_PressBeta0[0] = _M_PressBeta0[0] *
                    //  ( 1 + factor );
                }
            else
                std::cout << "\n[stiffenVesselRight] error! out of right boundary" << std::endl;
        }
    else
        {
            UInt iz=_M_paramSize-1;

            Real ratio, x_current=xr, deltax;

            deltax=(xr-xl)/static_cast<Real>(_M_paramSize-1);

            while ( (x_current > (alpha+delta/2)) && ((iz+1) > 0) )
                {
                    _M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor );
                    iz--;
                    x_current-=deltax;
                }

            while ( (x_current > alpha) && ((iz+1) > 0 ) )
                {
                    ratio=(( (alpha+delta/2) - x_current ) / delta);

                    _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] *
                        ( factor * ( n / delta) *
                          ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );

                    _M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor *
                                                              ( 1 - pow(2,(n-1)) * pow(ratio,n)) );
                    iz--;
                    x_current-=deltax;
                }

            while ( (x_current > (alpha-delta/2)) && ((iz+1) > 0) )
                {
                    ratio=(( x_current - (alpha-delta/2) ) / delta);

                    _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] *
                        ( factor * ( n / delta) *
                          ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );

                    _M_PressBeta0[iz] = _M_PressBeta0[iz] * ( 1 + factor *
                                                              ( pow(2,(n-1)) * pow(ratio,n)) );
                    iz--;
                    x_current-=deltax;
                }

        }

}


void BloodFlowParam::showMe(std::ostream& c) const
{
    //! parameters
    c << "\t[parameters]\n";
    c << "Area0\t= " << _M_Area0 << "\n";
    c << "beta0\t= " << _M_PressBeta0 << "\n";
    c << "beta1\t= " << _M_PressBeta1 << "\n";
    c << "dArea0dz    = " << _M_dArea0dz << "\n";
    c << "dbeta0dz = " << _M_dPressBeta0dz << "\n";
    c << "dbeta1dz = " << _M_dPressBeta1dz << "\n";
    c << "fluid mass density\t= " << _M_DensityRho << "\n";
    c << "celerity\t= " << Celerity0(0) << "\n" << std::endl;
    c << "friction parameter (Kr)\t= " << _M_FrictionKr << "\n";
    c << "wall mass density\t= " << _M_DensityWall << "\n";
    c << "viscoelastic modulus\t= " << _M_Gamma << "\n";
    c << "inertial modulus\t= " << _M_CoeffA << "\n";
    c << "thickness\t= " << _M_Thickness << "\n";
    c << "left boundary abscissa\t= " << _M_xleft << "\n";
    c << "right boundary abscissa\t= " << _M_xright << "\n";
    c << "length\t= " << _M_length << "\n";
    c << "size of the parameter vectors\t= " << _M_paramSize << "\n";
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++


OneDNonLinModelParam::OneDNonLinModelParam(const GetPot& dfile, const std::string& section ) :
        BloodFlowParam     ( dfile, section ),
        _M_AlphaCoriolis   (this->_M_paramSize),
        _M_dAlphaCoriolisdz(this->_M_paramSize)
{
    Real _alpha = dfile("parameters/alphaCor",1./M_ROBERTSON_CORRECTION);
    Real _dalphadz = 0.;

    //-------------------------------------------
    //! Initialisation of the parameter variables
    //-------------------------------------------
    _M_AlphaCoriolis = ScalarVector( _M_paramSize, _alpha );
    _M_dAlphaCoriolisdz = ScalarVector( _M_paramSize, _dalphadz );

    Debug( 6320 ) << "[OneDNonLinModelParam] Robertson correction coefficient = "
                  << M_ROBERTSON_CORRECTION << "\n";
}

Real OneDNonLinModelParam::AlphaCor(const UInt& ii) const {
    return _M_AlphaCoriolis[ii];
}

Real OneDNonLinModelParam::dAlphaCordz(const UInt& ii) const {
    return _M_dAlphaCoriolisdz[ii];
}

// Reimann Invariants corresponding to data (Q, A) at node indz
/*
  W1,2 = (Q / A) +- (2 / beta1) * sqrt(chi) * (celerity - celerity0)

  being chi the correction coefficient proposed by A. Robertson and H. Zakaria
*/
void
OneDNonLinModelParam::W_from_U( Real& _W1, Real& _W2,
                                const Real& _U1, const Real& _U2,
                                const UInt& indz ) const
{
    Real QoverA  = _U2 / _U1;

    Real AoverA0  = _U1 / _M_Area0[indz];

    Real celerity( Celerity0(indz) * std::sqrt( std::pow( AoverA0, _M_PressBeta1[indz] ) ) );

    Real add( std::sqrt( M_ROBERTSON_CORRECTION ) * ( celerity - Celerity0(indz) )
              * 2 / _M_PressBeta1[indz] );

    _W1 = QoverA + add;

    _W2 = QoverA - add;

}

// Physical variables corresponding to (W1, W2) at node indz
/*
  A = A0 * ( rho / (beta0 * beta1) )^(1/beta1)
  * ( beta1 / (4 * sqrt(chi) ) * (W1 - W2) + celerity0 )^(2/beta1)

  Q = A (W1 + W2) / 2
*/
void
OneDNonLinModelParam::U_from_W( Real& _U1, Real& _U2,
                                const Real& _W1, const Real& _W2,
                                const UInt& indz ) const
{
    Real rhooverbeta0beta1 ( _M_DensityRho
                             / ( _M_PressBeta0[indz] * _M_PressBeta1[indz] ) );

    Real beta1over4SQRTchi( _M_PressBeta1[indz] / ( std::sqrt(M_ROBERTSON_CORRECTION) * 4 ) );

    _U1 = _M_Area0[indz]
        * std::pow( rhooverbeta0beta1,
                    (1/_M_PressBeta1[indz]) )
        * std::pow( beta1over4SQRTchi * (_W1 - _W2)
                    + Celerity0(indz),
                    (2/_M_PressBeta1[indz]) );

    _U2 = _U1 * ( _W1 + _W2 ) / 2;

}

//! compute the pressure as a fuction of W1, W2:
/*
  P = beta0 * ( rho / (beta0 * beta1) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )^2 - 1 )
*/
Real OneDNonLinModelParam::
pressure_W(const Real& _W1, const Real& _W2,
           const UInt& indz) const
{
    Real rhooverbeta0beta1 ( _M_DensityRho
                             / ( _M_PressBeta0[indz] * _M_PressBeta1[indz] ) );

    Real beta1over4SQRTchi( _M_PressBeta1[indz] / ( std::sqrt(M_ROBERTSON_CORRECTION) * 4 ) );

    return ( _M_PressBeta0[indz]
    		 * ( rhooverbeta0beta1
                 * std::pow( beta1over4SQRTchi * (_W1 - _W2)
                             + Celerity0(indz), 2 )
                 - 1 )
             );
}

/*! Derivative of pressure as a function of (W1, W2)
  dP(W1,W2)/dW_1 = rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
  dP(W1,W2)/dW_2 = - rho / (2 * sqrt(chi)) * ( beta1 / (4 * sqrt(chi)) * (W1 - W2) + celerity0 )
*/
Real OneDNonLinModelParam::
pressure_WDiff( const Real& _W1, const Real& _W2,
                const ID& ii,
                const UInt& indz) const
{
    Real rhoover2SQRTchi ( _M_DensityRho / ( std::sqrt(M_ROBERTSON_CORRECTION) * 2 ) );

    Real beta1over4SQRTchi( _M_PressBeta1[indz] / ( std::sqrt(M_ROBERTSON_CORRECTION) * 4 ) );

    Real result( beta1over4SQRTchi * (_W1 - _W2)  );
    result += Celerity0(indz);
    result *= rhoover2SQRTchi;

    if( ii == 1 ) { //! dP/dW1
		return result;
    }
    if( ii == 2 ) { //! dP/dW2
        return -result;
    }
    ERROR_MSG("P(W1,W2)'s differential function has only 2 components.");
    return -1.;
}


/*! compute W1 or W2 given the pressure:
  W1 - W2 = (4 * sqrt(chi) / beta1) * sqrt( beta0 * beta1 / rho ) ) * ( sqrt( P / beta0 + 1 ) - 1 )
*/
Real OneDNonLinModelParam::
W_from_P(const Real& _P, const Real& _W, const ID& ii, const UInt& indz) const
{
    Real SQRTbeta0beta1overrho( _M_PressBeta0[indz] * _M_PressBeta1[indz] / _M_DensityRho );
    SQRTbeta0beta1overrho = std::sqrt( SQRTbeta0beta1overrho );

    Real SQRTchi4overbeta1( std::sqrt(M_ROBERTSON_CORRECTION) * 4 / _M_PressBeta1[indz] );

	Real add( SQRTchi4overbeta1 * SQRTbeta0beta1overrho
              * ( pow( ( _P / _M_PressBeta0[indz] + 1 ), 0.5 ) - 1 ) );

	Debug(6320) << "[OneDNonLinModelParam::W_from_P] "
				<< "SQRTchi4overbeta1 = " << SQRTchi4overbeta1
				<< ", beta0beta1overrho = " << SQRTbeta0beta1overrho
				<< ", pow( ( _P / _M_PressBeta0[indz] + 1 ), 0.5 ) = " << pow( ( _P / _M_PressBeta0[indz] + 1 ), 0.5 ) << "\n";
	Debug(6320) << "[OneDNonLinModelParam::W_from_P] add term = " << add << "\n";

	if( ii == 1 )
		return _W - add;
    if( ii == 2 )
		return _W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}


/*! compute W1 or W2 given the flux: fixed point problem
  ( W1 - W2 + celerity0/K0 )^(2/beta1) * ( W1 + W2 ) = Q / K1

  where
  K0 = beta1 / ( 4 * sqrt(chi) )
  K1 = A0 / 2 * ( rho / (beta0*beta1) )^(1/beta1) * K0^(2/beta1)
*/
Real OneDNonLinModelParam::
W_from_Q(const Real& _Q, const Real& _W_n, const Real& _W, const ID& ii, const UInt& indz) const
{
	Real K0( _M_PressBeta1[indz] / ( std::sqrt(M_ROBERTSON_CORRECTION) * 4 ) );

	Real K1( (_M_Area0[indz] / 2) );
	K1 *= pow( _M_DensityRho / (_M_PressBeta0[indz] * _M_PressBeta1[indz]), 1/_M_PressBeta1[indz] );
	K1 *= pow( K0, 2/_M_PressBeta1[indz] );

	Real w_k = _W_n;
	Real f_k, df_k, tau_k;

	if( ii == 1 ){ // W1 given
		f_k = pow( _W - w_k + Celerity0(indz) / K0, 2/_M_PressBeta1[indz] );
		tau_k = pow( _W - w_k + Celerity0(indz) / K0, 2/_M_PressBeta1[indz] );
		df_k = (-2 / _M_PressBeta1[indz])
            * pow( _W - w_k + Celerity0(indz) / K0, 2/_M_PressBeta1[indz] - 1 );
	}
	if( ii == 2 ){ // W2 given
		f_k = pow( w_k - _W + Celerity0(indz) / K0, 2/_M_PressBeta1[indz] );
		tau_k = pow( w_k - _W + Celerity0(indz) / K0, 2/_M_PressBeta1[indz] );
		df_k = (-2 / _M_PressBeta1[indz]) * pow( w_k - _W + Celerity0(indz) / K0, 2/_M_PressBeta1[indz] - 1 );
	}
	f_k *= (_W + w_k);
	f_k += - _Q / K1;
	df_k *= (_W + w_k);
	df_k += f_k;

	Debug(6320) << "[OneDNonLinModelParam::W_from_Q] "
				<< "K0 = " << K0
				<< ", K1 = " << K1
				<< ", tau_k = " << tau_k << "\n";

	Real w_kp1 = _Q / (K1 * tau_k) - _W;

        /* for debugging purposes
	std::ofstream ofile;
	ofile.open( "imposed_W_from_Q.m", std::ios::app );
	ofile << _Q << "\t"
          << w_kp1 << "\n";
	ofile.close();
        */

	return w_kp1;

    //    ERROR_MSG("You can only find W1 or W2 as function of P");
    //    return -1.;
}


void OneDNonLinModelParam::showMe(std::ostream& c) const
{
    //! parameters
    BloodFlowParam::showMe(c);
    c << "alphaCor = " << _M_AlphaCoriolis << "\n";
    c << "dalphaCordz = " << _M_dAlphaCoriolisdz << "\n";
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++


//+++++++++++++++++++++++++++++++++++++++++++++++++++

LinearSimpleParam::LinearSimpleParam(const GetPot& dfile) :
    //! we suppose so far that all params are constant along the vessel
    //! otherwise replace _M_paramSize(1) by _M_paramSize(_M_dimDof)
    _M_paramSize(dfile("discretization/nb_elem", 0 ) + 1),
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
    return _M_Flux11[ii];
}
Real LinearSimpleParam::Flux12(const UInt& ii) const {
    return _M_Flux12[ii];
}
Real LinearSimpleParam::Flux21(const UInt& ii) const {
    return _M_Flux21[ii];
}
Real LinearSimpleParam::Flux22(const UInt& ii) const {
    return _M_Flux22[ii];
}

Real LinearSimpleParam::Celerity1(const UInt& ii) const {
    return _M_celerity1[ii];
}

Real LinearSimpleParam::Celerity2(const UInt& ii) const {
    return _M_celerity2[ii];
}

Real LinearSimpleParam::LeftEigenVector11(const UInt& ii) const {
    return _M_celer1_left_eigVector1[ii];
}

Real LinearSimpleParam::LeftEigenVector12(const UInt& ii) const {
    return _M_celer1_left_eigVector2[ii];
}

Real LinearSimpleParam::LeftEigenVector21(const UInt& ii) const {
    return _M_celer2_left_eigVector1[ii];
}

Real LinearSimpleParam::LeftEigenVector22(const UInt& ii) const {
    return _M_celer2_left_eigVector2[ii];
}

Real LinearSimpleParam::Source10(const UInt& ii) const {
    return _M_Source10[ii];
}
Real LinearSimpleParam::Source20(const UInt& ii) const {
    return _M_Source20[ii];
}
Real LinearSimpleParam::Source11(const UInt& ii) const {
    return _M_Source11[ii];
}
Real LinearSimpleParam::Source12(const UInt& ii) const {
    return _M_Source12[ii];
}
Real LinearSimpleParam::Source21(const UInt& ii) const {
    return _M_Source21[ii];
}
Real LinearSimpleParam::Source22(const UInt& ii) const {
    return _M_Source22[ii];
}

void LinearSimpleParam::showMe(std::ostream& c) const
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


LinearizedParam::LinearizedParam(const GetPot& dfile) :
    BloodFlowParam( dfile ), LinearSimpleParam( dfile )
{
    /*
      The linearization of Euler model yields

      F = [ Q; A * (c_L)^2];

      B = [ 0; k_R / A0];

      c_L = sqrt( beta0 * beta1 / rho );
    */
    for( UInt indz=0; indz < LinearSimpleParam::_M_paramSize; ++indz )
        {
            _M_celerity1[indz] = std::sqrt( _M_PressBeta0[indz] * _M_PressBeta1[indz] / _M_DensityRho );
            _M_Flux21[indz] = std::pow( _M_celerity1[indz], 2 );
            _M_Source22[indz] = _M_FrictionKr[indz] / Area0(indz);
        }
    _M_celerity2 = - _M_celerity1;

    _M_celer1_left_eigVector1 = _M_celerity1;
    _M_celer1_left_eigVector2 = ScalarVector(LinearSimpleParam::_M_paramSize, 1.);
    _M_celer2_left_eigVector1 = _M_celerity2;
    _M_celer2_left_eigVector2 = ScalarVector(LinearSimpleParam::_M_paramSize, 1.);

    _M_Flux11 = ZeroVector(LinearSimpleParam::_M_paramSize);
    _M_Flux12 = ScalarVector(LinearSimpleParam::_M_paramSize, 1.);
    _M_Flux22 = ZeroVector(LinearSimpleParam::_M_paramSize);

    _M_Source10 = ZeroVector(LinearSimpleParam::_M_paramSize);
    _M_Source11 = ZeroVector(LinearSimpleParam::_M_paramSize);
    _M_Source12 = ZeroVector(LinearSimpleParam::_M_paramSize);
    //
    _M_Source20 = ZeroVector(LinearSimpleParam::_M_paramSize);
    _M_Source21 = ZeroVector(LinearSimpleParam::_M_paramSize);

}

// Reimann Invariants corresponding to data (Q, A) at node indz
/*
  W1,2 = Q +- celerity * ( A - A0 )
*/
void
LinearizedParam::W_from_U( Real& _W1, Real& _W2,
                           const Real& _U1, const Real& _U2,
                           const UInt& indz ) const
{
    _W1 = _U2 + Celerity0(indz) * ( _U1 - Area0(indz) );

    _W2 = _U2 - Celerity0(indz) * ( _U1 - Area0(indz) );

    Debug( 6320 ) << "[LinearizedParam::W_from_U] Q " << _U2 << "\n";
    Debug( 6320 ) << "[LinearizedParam::W_from_U] W1 " << _W1 << "\n";
    Debug( 6320 ) << "[LinearizedParam::W_from_U] W2 " << _W2 << "\n";
    Debug( 6320 ) << "[LinearizedParam::W_from_U] Celerity " << Celerity0(indz) << "\n";
    Debug( 6320 ) << "[LinearizedParam::W_from_U] ( _U1 - Area0(indz) ) " << ( _U1 - Area0(indz) ) << "\n";
}

// Physical variables corresponding to (W1, W2) at node indz
/*
  A = A0 + (W1 - W2) / (2 * celerity)

  Q = (W1 + W2) / 2
*/
void
LinearizedParam::U_from_W( Real& _U1, Real& _U2,
                           const Real& _W1, const Real& _W2,
                           const UInt& indz ) const
{
    _U1 = Area0(indz) + (_W1 - _W2) / ( 2 * Celerity0(indz) );

    _U2 = ( _W1 + _W2 ) / 2;

}


//! compute the pressure as a function of W1, W2:
/*
  P = beta0 * ( ( 1 / Area0 )^(beta1) * ( (W1 - W2) / (2 * celerity0) + Area0 )^(beta1) - 1 )
*/
Real LinearizedParam::
pressure_W(const Real& _W1, const Real& _W2,
           const UInt& indz) const
{
    return ( _M_PressBeta0[indz]
             * ( std::pow( 1 / Area0(indz), _M_PressBeta1[indz] )
                 * std::pow( (_W1 - _W2) / ( 2*Celerity0(indz) ) + Area0(indz), _M_PressBeta1[indz] )
                 - 1 )
             );
}

/*! Derivative of pressure as a function of (W1, W2)
  dP(W1,W2)/dW_1 = beta0 * beta1 / ( 2 * celerity0 * Area0^(beta1) ) * ( (W1 - W2) / ( 2 * celerity0 ) + Area0 )^(beta1-1)
  dP(W1,W2)/dW_2 = - dP(W1,W2)/dW_1
*/
Real LinearizedParam::
pressure_WDiff( const Real& _W1, const Real& _W2,
                const ID& ii,
                const UInt& indz) const
{
    Real beta0beta1overA0beta1 ( _M_PressBeta0[indz] * _M_PressBeta1[indz] / std::pow( Area0(indz), _M_PressBeta1[indz] ) );

    Real oneover2celerity( 1 / ( 2 * Celerity0(indz) ) );

    Real result( beta0beta1overA0beta1 * oneover2celerity );
    result *= ( ( _W1 - _W2 ) * oneover2celerity + Area0(indz) );

    if( ii == 1 ) { //! dP/dW1
        return result;
    }
    if( ii == 2 ) { //! dP/dW2
        return -result;
    }
    ERROR_MSG("P(W1,W2)'s differential function has only 2 components.");
    return -1.;
}


/*! compute W1 or W2 given the pressure:
  W1 - W2 = (2 * celerity * A0) * ( ( P / beta0 + 1 )^(1/beta1) - 1 )
*/
Real LinearizedParam::
W_from_P(const Real& _P, const Real& _W, const ID& ii, const UInt& indz) const
{
    Real add( 2 * Celerity0(indz) * Area0(indz)
              * ( pow( ( _P / _M_PressBeta0[indz] + 1 ), 1 / _M_PressBeta1[indz] ) - 1 ) );

    Debug(6320) << "[W_from_P] "
                << "2 * Celerity0(indz) * Area0(indz) = " << 2 * Celerity0(indz) * Area0(indz)
                << ", pow( ( _P / _M_PressBeta0[indz] + 1 ), 1 / _M_PressBeta1[indz] ) = "
                << pow( ( _P / _M_PressBeta0[indz] + 1 ), 1 / _M_PressBeta1[indz]  ) << "\n";
    Debug(6320) << "[W_from_P] add term = " << add << "\n";

    if( ii == 1 )
        return _W - add;
    if( ii == 2 )
        return _W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}


/*! compute W1 or W2 given the flux

  W1 + W2 = 2 * Q
*/
Real LinearizedParam::
W_from_Q(const Real& _Q, const Real& /*_W_n*/, const Real& _W, const ID& ii, const UInt& /*indz*/) const
{
    Real add( 2 * _Q );
    if( ii == 1 ){ // W1 given
        return add - _W;
    }
    if( ii == 2 ){ // W2 given
        return add - _W;
    }
    ERROR_MSG("You can only find W1 or W2 as function of Q");
    return -1.;
}


//! output
void LinearizedParam::
showMe(std::ostream& c) const
{
    LinearSimpleParam::showMe(c);
    BloodFlowParam::showMe(c);
}

}
