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
  \file ud_functions.hpp
  \author Tiziano Passerini
  \author Lucia Mirabella *queen*

  \date 04/2007

  \brief This file contains functions to impose boundary conditions for 1D tubes.

*/

#ifndef _ONED_FUNCTIONS_1D_
#define _ONED_FUNCTIONS_1D_
#include <cmath>
#include <life/lifecore/life.hpp>

#include <life/lifefem/oneDBCFunctions.hpp>

namespace LifeV
{

  /*!
    \brief A sinusoidal wave
  */
  class Sin : public OneDBCFunctionBase
  {
  public:
    /*!
      \brief The constructor

      \param[in] mean the time average of the sinus \f$ B \f$
      \param[in] scale the sinus amplitude \f$ A \f$
      \param[in] period the sinus period \f$ T \f$
      \param[in] phase the sinus phase \f$ \phi \f$
    */
    Sin(const Real mean = 0, const Real scale = 10,
        const Real period = 1., const Real phase = 0. ) :
      _M_mean(mean), _M_scale(scale), _M_period(period), _M_phase(phase) {}

    /*!
      \brief Compute the sinus at the specified time

      \param[in] time the time
      \return \f$ B + A sin( \frac{2 \pi t}{T} + \phi ) \f$
    */
    Real evaluate( const Real& time )
    {
      return _M_mean + _M_scale*std::sin(_M_phase+2*M_PI*time/_M_period);
    };

    ~Sin() {}
  private:
    Real _M_mean;
    Real _M_scale;
    Real _M_period;
    Real _M_phase;
  };


  /*!
    \brief A superimposition of a sinusoidal and a cosinusoidal waves,
      whose amplitude is damped by exponential terms
  */
  class Cos_min_Sin : public OneDBCFunctionBase
  {
  public:
    /*!
      \brief The constructor

      \param[in] coeff_exp_t_cos the coefficient multiplying time in the
        exponential damping the cosinus amplitude \f$ \lambda_c \f$
      \param[in] mean_cos the time average of the cosinus \f$ B_c \f$
      \param[in] amplitude_cos the cosinus amplitude \f$ A_c \f$
      \param[in] frequency_cos the cosinus period \f$ T_c \f$
      \param[in] phase_cos the cosinus phase \f$ \phi_c \f$
      \param[in] coeff_exp_t_sin the coefficient multiplying time in the
        exponential damping the sinus amplitude \f$ \lambda_s \f$
      \param[in] mean_sin the time average of the sinus \f$ B_s \f$
      \param[in] amplitude_sin the sinus amplitude \f$ A_s \f$
      \param[in] frequency_sin the sinus period \f$ T_s \f$
      \param[in] phase_sin the sinus phase \f$ \phi_s \f$
    */
    Cos_min_Sin(const Real coeff_exp_t_cos = 0, const Real mean_cos = 0,
        const Real amplitude_cos = 10, const Real frequency_cos = 8.*atan(1.),
        const Real phase_cos = 0.,
        const Real coeff_exp_t_sin = 0, const Real mean_sin = 0,
        const Real amplitude_sin = 10, const Real frequency_sin = 8.*atan(1.),
        const Real phase_sin = 0. ) :
    	_M_coeff_exp_t_cos(coeff_exp_t_cos), _M_mean_cos(mean_cos),
    	_M_amplitude_cos(amplitude_cos), _M_frequency_cos(frequency_cos), _M_phase_cos(phase_cos),
    	_M_coeff_exp_t_sin(coeff_exp_t_sin), _M_mean_sin(mean_sin),
    	_M_amplitude_sin(amplitude_sin), _M_frequency_sin(frequency_sin), _M_phase_sin(phase_sin) {}

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
        ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    Real evaluate( const Real& time )
    {
      Real result = _M_mean_cos
      + _M_amplitude_cos * std::cos(_M_phase_cos+time*_M_frequency_cos)
      - ( _M_mean_sin
      + _M_amplitude_sin * std::sin(_M_phase_sin+time*_M_frequency_sin) );

      return result;
    };

    /*! \name Getters
     */
    //@{
    Real& coeff_exp_t_cos(){ return _M_coeff_exp_t_cos;}
    Real& mean_cos(){ return _M_mean_cos;}
    Real& amplitude_cos(){ return _M_amplitude_cos;}
    Real& frequency_cos(){ return _M_frequency_cos;}
    Real& phase_cos(){ return _M_phase_cos;}
    Real& coeff_exp_t_sin(){ return _M_coeff_exp_t_sin;}
    Real& mean_sin(){ return _M_mean_sin;}
    Real& amplitude_sin(){ return _M_amplitude_sin;}
    Real& frequency_sin(){ return _M_frequency_sin;}
    Real& phase_sin(){ return _M_phase_sin;}
    //@}

    ~Cos_min_Sin() {}
  private:
    Real _M_coeff_exp_t_cos;
    Real _M_mean_cos;
    Real _M_amplitude_cos;
    Real _M_frequency_cos;
    Real _M_phase_cos;
    Real _M_coeff_exp_t_sin;
    Real _M_mean_sin;
    Real _M_amplitude_sin;
    Real _M_frequency_sin;
    Real _M_phase_sin;
  };


  /*!
    \brief A particular case of Cos_min_Sin

    This analytical solution is in the form

    U = ( Re(U) + i Im(U) ) * exp( Im(k) x ) *
        exp[ i ( \omega t - Re(k) x ) ] .

    Its real part is
    Re(U) = [ Re(U) * exp( Im(k) x ) ] * cos( \omega t - Re(k) x ) +
     - [ Im(U) * exp( Im(k) x ) ] * sin( \omega t - Re(k) x ) +
  */
  class Analytical_Solution : public Cos_min_Sin
  {
  public:
    /*!
      \brief The constructor

      \param[in] sol_amplitude_Re the real part of the wave amplitude
      \param[in] sol_amplitude_Im the imaginary part of the wave amplitude
      \param[in] kappa_Re the real part of the wave number \f$ k \f$
      \param[in] kappa_Im the imaginary part of the wave number \f$ k \f$
      \param[in] omega the wave time frequency \f$ \omega \f$ (supposed to be real)
    */
    Analytical_Solution( Real const& sol_amplitude_Re, Real const& sol_amplitude_Im,
        Real const& kappa_Re, Real const& kappa_Im,
        Real const& omega ) :
          Cos_min_Sin(0., 0., sol_amplitude_Re, omega, 0.,
              0., 0., sol_amplitude_Im, omega, 0.),
          _M_sol_amplitude_Re(sol_amplitude_Re),
          _M_sol_amplitude_Im(sol_amplitude_Im),
          _M_kappa_Re(kappa_Re),
          _M_kappa_Im(kappa_Im)
          {}

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
        ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    Real evaluate( const Real& time )
    {
      return Cos_min_Sin::evaluate(time);
    };

    /*!
      \brief Compute the wave at the specified time

      \param[in] time the time
      \return \f$ B_c + A_c cos( \frac{2 \pi t}{T_c} + \phi_c ) -
        ( B_s + A_s sin( \frac{2 \pi t}{T_s} + \phi_s ) ) \f$
    */
    void update_x( const Real& _x )
    {
      this->amplitude_cos() = _M_sol_amplitude_Re * std::exp( _M_kappa_Im * _x );
      this->phase_cos() = - _M_kappa_Re * _x;
      this->amplitude_sin() = _M_sol_amplitude_Im * std::exp( _M_kappa_Im * _x );
      this->phase_sin() = - _M_kappa_Re * _x;
    };

    ~Analytical_Solution() {}

  private:
    Real _M_sol_amplitude_Re;
    Real _M_sol_amplitude_Im;
    Real _M_kappa_Re;
    Real _M_kappa_Im;
  };


  class PhysiologicalFlux : public OneDBCFunctionBase
  {
  public:
    PhysiologicalFlux( GetPot const& data_file );

    Real evaluate( const Real& time );

    ~PhysiologicalFlux() {}
  private:
	Real _M_rampT;
	Real _M_time_step;
    Real _M_scale;
  };



  template<class FLUX, class SOURCE, class PARAM>
  class PressureRamp : public Compatibility<FLUX, SOURCE>
  {
  public:
    PressureRamp(const BasicOneDMesh& mesh,
        const FLUX& fluxFun,
        const SOURCE& sourceFun,
        const OneDModelHandler::ScalVec_vector& U_thistime,
        const Real& dt, const std::string& border, const std::string & var,
        const PARAM& onedparam,
        const Real& startT = 0., const Real& duration = 0.05, const Real& endvalue = 106400 );

    Real evaluate( const Real& time );

    ~PressureRamp() {}
  private:
    //! Reference to the solver parameters
    const PARAM& _M_onedparam;

    Real _M_startT;
    Real _M_duration;
    Real _M_endvalue;
  };


  template<class FLUX, class SOURCE, class PARAM>
  PressureRamp<FLUX, SOURCE, PARAM>::PressureRamp (const BasicOneDMesh& mesh,
      const FLUX& fluxFun,
      const SOURCE& sourceFun,
      const OneDModelHandler::ScalVec_vector& U_thistime,
      const Real& dt, const std::string& border, const std::string & var,
      const PARAM& onedparam,
      const Real& startT, const Real& duration, const Real& endvalue )
  :
    Compatibility<FLUX, SOURCE>( mesh, fluxFun, sourceFun, U_thistime,/* W_thistime,*/ dt, border, var),
    _M_onedparam( onedparam ),
    _M_startT( startT ),
    _M_duration( duration ),
    _M_endvalue( endvalue )
    {}



  template<class FLUX, class SOURCE, class PARAM>
  Real PressureRamp<FLUX, SOURCE, PARAM>::evaluate( const Real& time )
  {
    Real W_out, result(0.);

    Real _P = ( time < (_M_startT + _M_duration) ) ?
        ( ( time - _M_startT ) / _M_duration ) : 1;
        _P *= _M_endvalue;

        Debug( 6030 ) << "[PressureRamp::evaluate] imposed pressure = " << _P << "\n";
        Debug( 6030 ) << "[PressureRamp::evaluate] target pressure = " << _M_endvalue << "\n";

        switch( this->_M_oneDBCFunctionsMapStringValues[this->_M_var] )
        {
        case OneDBCW1:
          W_out = this->extrapolate_W( OneDBCW2 );
          result = _M_onedparam.W_from_P( _P, W_out, 2, this->_M_boundaryDof);
          break;
        case OneDBCW2:
          W_out = this->extrapolate_W( OneDBCW1 );
          result = _M_onedparam.W_from_P( _P, W_out, 1, this->_M_boundaryDof);
          break;
        default:
          std::cout << "\n[PressureRamp::evaluate] incorrect variable identifier: " << this->_M_var << std::endl;
        }

        Debug( 6030 ) << "[PressureRamp::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

        return result;

  }



  template<class FLUX, class SOURCE, class PARAM>
  class Heart : public Compatibility<FLUX, SOURCE>
  {
  public:
    Heart ( const GetPot& data_file,
        const PARAM& onedparam,
        const BasicOneDMesh& mesh,
        const FLUX& fluxFun,
        const SOURCE& sourceFun,
        const OneDModelHandler::ScalVec_vector& U_thistime,
        const Real& dt, const std::string& border, const std::string & var,
        const bool type, const Real& startT = 0. );

    Real evaluate( const Real& time );

    ~Heart() {}
  private:
    //! Compute Pv function
    Real PvFunc(const Real& T_reset) const;
    //! Compute Pv from the Elastance value and the current volume
    Real PvCalc(const Real& T_reset);
    //! Compute elastance function
    Real Elastance(const Real& T_reset) const;
    //! Compute compliance function
    Real Compliance(const Real& T_reset) const;

    //! Data file
    GetPot const& _M_data_file;
    //! Reference to the solver parameters
    const PARAM& _M_onedparam;

    //! how to compute ventricular pressure
    bool _M_type;

    Real _M_startT;
    //! Boolean for the aortic valve: it is 1 if the valve is open
    bool _M_OV;
    //! Reference value for the ventricule volume
    Real V0;
    //! Systole period
    Real _M_periodosis;
    //! Whole cycle period
    Real _M_periodotot;
    //! Elastance
    Real El;
    //! Compliance
    Real Co;
    //! Ventricular pressure
    Real Pv;
    //! Aortic pressure
    Real Pa;
    //! Previous step ventricular volume
    Real _M_Vol_old;
    //! Current step ventricular volume
    Real _M_Vol_new;
    //! Previous step Q
    Real _M_Q_old;
    //! current step Q
    Real _M_Q_new;
    //! exiting characteristic variable
    Real _M_charact;
    //! Value of U at the boundary, at the neighboring internal node
    //  OneDModelHandler::Vec2D U_boundary, U_internalBd;
    //*************************************
  };



  template<class FLUX, class SOURCE, class PARAM>
  Heart<FLUX, SOURCE, PARAM>::Heart ( const GetPot& data_file, const PARAM& onedparam,
      const BasicOneDMesh& mesh,
      const FLUX& fluxFun,
      const SOURCE& sourceFun,
      const OneDModelHandler::ScalVec_vector& U_thistime,
      const Real& dt, const std::string& border, const std::string & var,
      const bool type, const Real& startT )
  :
    Compatibility<FLUX, SOURCE>( mesh, fluxFun, sourceFun, U_thistime, /*W_thistime,*/ dt, border, var),
    _M_data_file(data_file),
    _M_onedparam( onedparam ),
    _M_type(type),
    _M_startT(startT),
    _M_OV(0),
    V0(50),
    _M_periodosis(data_file( "time/sysperiod", 0.5) ),
    _M_periodotot(data_file( "time/cycleperiod", 1.) ),
    _M_Vol_old(120)
    {
      // wrong if Heart is istantiated before solver initialization
      _M_Q_old=this->_M_U_thistime[1](this->_M_boundaryDof);
      Debug( 6030 ) << "[Heart::Heart] U1_thistime(this->_M_boundaryDof) = "
      << this->_M_U_thistime[0](this->_M_boundaryDof) << "\n";

    }



  template<class FLUX, class SOURCE, class PARAM>
  Real Heart<FLUX, SOURCE, PARAM>::evaluate( const Real& time )
  {
    Real W_out, result;
    UInt W_out_id;

    Debug( 6030 ) << "[Heart::evaluate] U1_thistime(this->_M_boundaryDof) = " << this->_M_U_thistime[0](this->_M_boundaryDof) << "\n";
    Debug( 6030 ) << "[Heart::evaluate] _M_U_boundary[0] = " << this->_M_U_boundary[0] << "\n";

    Real T_reset( ( time - _M_startT ) );
    T_reset -= static_cast<int>( std::floor( (T_reset + this->_M_time_step/2) /_M_periodotot ) ) * _M_periodotot;

    if (_M_type==0)
      Pv=PvFunc(T_reset);
    else
      Pv=PvCalc(T_reset);

    switch( this->_M_oneDBCFunctionsMapStringValues[this->_M_var] )
    {
    case OneDBCW1:
      W_out = this->extrapolate_W( OneDBCW2 );
      W_out_id = 2;
      break;
    case OneDBCW2:
      W_out = this->extrapolate_W( OneDBCW1 );
      W_out_id = 1;
      break;
    default:
      std::cout << "\n[Heart::evaluate] incorrect variable identifier: " << this->_M_var << std::endl;
    }

    Debug( 6030 ) << "[Heart::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

    Pa=_M_onedparam.pressure(this->_M_U_boundary[0], this->_M_boundaryDof);

    Debug(6030) << "[Heart::evaluate] imposed heart pressure = " << Pv
    << ",\taorta pressure = " << Pa
    << ",\ttime = " << time << "\n";

    if (_M_OV==0)
      {
        Debug(6030)<<"[Heart::evaluate] CLOSED valve\n";
        result = - W_out;
        if (Pv > Pa) _M_OV=1;
      }
    else
      {
        Debug(6030)<<"[Heart::evaluate] OPEN valve\n";
        result = _M_onedparam.W_from_P( Pv, W_out, W_out_id, this->_M_boundaryDof);
        if (this->_M_U_boundary[1]<0) _M_OV=0;
      }

    return result;

  }


  template<class FLUX, class SOURCE, class PARAM>
  Real Heart<FLUX, SOURCE, PARAM>::PvFunc(const Real& T_reset) const
  {
    Debug(6030)<<"[Heart::PvFunc] T_reset" << T_reset << "\n";
    Debug(6030)<<"[Heart::PvFunc] _M_periodosis" << _M_periodosis << "\n";
    if (T_reset < _M_periodosis)
      return 160000 * std::sin(M_PI/_M_periodosis * T_reset);
    else
      return 0;

  }

  template<class FLUX, class SOURCE, class PARAM>
  Real Heart<FLUX, SOURCE, PARAM>::PvCalc(const Real& T_reset)
  {

    El=Elastance(T_reset);
    Co=Compliance(T_reset);
    Debug(6030)<<"[Heart::PvCalc] Co = "<<Co<<"\n";
    if (T_reset == 0)
      {
        _M_Vol_old=120;
      }
    _M_Q_new=this->_M_U_boundary[1];
    _M_Vol_new=_M_Vol_old - (_M_Q_new+_M_Q_old)/2 * this->_M_time_step;
    _M_Vol_old=_M_Vol_new;
    _M_Q_old=_M_Q_new;

    Debug(6030)<< "[Heart::PvCalc] heart volume = " << _M_Vol_new << "\n";
    //return 1/Co*(_M_Vol_new-V0);
    return El*(_M_Vol_new-V0);

  }


  template<class FLUX, class SOURCE, class PARAM>
  Real Heart<FLUX, SOURCE, PARAM>::Elastance(const Real& T_reset) const
  {

    if (T_reset< _M_periodosis)
      return 1333* ( 1 + 0.7 * std::sin(  3*M_PI/(2*_M_periodosis) * T_reset ) );//spiegare i numeri (LM 2006)
    else
      return 0.3 + 0.7 * (T_reset - _M_periodosis)/(_M_periodotot-_M_periodosis);
  }

  template<class FLUX, class SOURCE, class PARAM>
  Real Heart<FLUX, SOURCE, PARAM>::Compliance(const Real& T_reset) const
  {

    if (T_reset< _M_periodosis)
      return 0.11*pow((0.00003/0.0146),((1.-exp(-T_reset/(0.0025*60.)))/(1.-exp(-_M_periodosis/(0.0025*60.)))));//spiegare i numeri
    else
      return 2.25e-4*pow((0.0146/0.00003),((1.-exp(-(T_reset-_M_periodosis)/(0.0075*60.)))/(1-exp(-(_M_periodotot-_M_periodosis)/(0.0075*60.)))));

  }



  template<class FLUX, class SOURCE, class PARAM>
  class Resistance : public Compatibility<FLUX, SOURCE>
  {
  public:
    Resistance(  const Real & resistance, // const GetPot& data_file,
        const PARAM& onedparam,
        const BasicOneDMesh& mesh,
        const FLUX& fluxFun,
        const SOURCE& sourceFun,
        const OneDModelHandler::ScalVec_vector& U_thistime,
        const Real& dt, const std::string& border, const std::string & var, const bool& absorbing = false);

    Real evaluate( const Real& time );

    ~Resistance() {}
  private:
    //! Resistance value
    Real _M_resistance;
    //! Reference to the solver parameters
    const PARAM& _M_onedparam;
    //! to impose absorbing conditions
    bool _M_absorbing;
  };


  template<class FLUX, class SOURCE, class PARAM>
  Resistance<FLUX, SOURCE, PARAM>::Resistance(  const Real & resistance,
      const PARAM& onedparam, // const GetPot& data_file,
      const BasicOneDMesh& mesh,
      const FLUX& fluxFun,
      const SOURCE& sourceFun,
      const OneDModelHandler::ScalVec_vector& U_thistime,
      const Real& dt, const std::string& border, const std::string & var, const bool& absorbing )
  :
    Compatibility<FLUX, SOURCE>( mesh, fluxFun, sourceFun, U_thistime, /*W_thistime,*/ dt, border, var),
    _M_resistance(resistance),
    _M_onedparam( onedparam ),
    _M_absorbing( absorbing )
    {
      Debug( 6030 ) << "[Resistance::Resistance] resistance = " << _M_resistance << "\n";
    }



  template<class FLUX, class SOURCE, class PARAM>
  Real Resistance<FLUX, SOURCE, PARAM>::evaluate( const Real& /*time*/ )
  {
    //! Coefficients
    Real W_out(0.), result;
    Real a1, a2, a11, a22, b1, b2, c1, c2;

    this->update_U_boundary();
    this->update_U_internalBd();

    Debug( 6030 ) << "[Resistance::Resistance] at node " << this->_M_boundaryDof
    << ", A = " << this->_M_U_boundary[0] << "( " << this->_M_U_thistime[0][0] << " ) "
    << ", Q = " << this->_M_U_boundary[1]
    << ", W1 = " << this->_M_W_boundary[0]
    << ", W2 = " << this->_M_W_boundary[1]
    << "\n";

    this->computeEigenValuesVectors();

    a1 = _M_onedparam.pressure(this->_M_U_boundary[0], this->_M_boundaryDof); // pressure at previous time step

    a2 = this->_M_U_boundary[1]; // flux at previous time step

    b1 = _M_onedparam.pressure_WDiff( this->_M_W_boundary[0], this->_M_W_boundary[1], 1, this->_M_boundaryDof);  // dP / dW1

    b2 = this->_M_U_boundary[0] / 2; // dQ / dW1

    c1 = _M_onedparam.pressure_WDiff( this->_M_W_boundary[0], this->_M_W_boundary[1], 2, this->_M_boundaryDof);  // dP / dW2

    c2 = b2; // dQ / dW2

    Debug( 6030 ) << "[Resistance::evaluate] P(A) = " << a1 << "\n";
    Debug( 6030 ) << "[Resistance::evaluate] P(W1,W2) = "
    << _M_onedparam.pressure_W(this->_M_W_boundary[0], this->_M_W_boundary[1], this->_M_boundaryDof) << "\n";

    a11 = a1 - b1*this->_M_W_boundary[0] - c1*this->_M_W_boundary[1];
    a22 = a2 - b2*this->_M_W_boundary[0] - c2*this->_M_W_boundary[1];

    switch( this->_M_oneDBCFunctionsMapStringValues[this->_M_var] )
    {
    case OneDBCW1:
      W_out = this->extrapolate_L_dot_U(this->_M_eigval2, this->_M_left_eigvec2)
      - dot( this->_M_left_eigvec2, this->_M_U_boundary ) + this->_M_W_boundary[1];

      break;
    case OneDBCW2:
      W_out = this->extrapolate_L_dot_U(this->_M_eigval1, this->_M_left_eigvec1)
      - dot( this->_M_left_eigvec1, this->_M_U_boundary ) + this->_M_W_boundary[0];

      break;
    default:
      std::cout << "\n[Resistance::evaluate] incorrect variable identifier: " << this->_M_border << std::endl;
    }

    Debug( 6030 ) << "[Resistance::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

    if( _M_absorbing ) {
      _M_resistance = b1 / b2;
      Debug( 6030 ) << "[Resistance::evaluate] imposing absorbing condition, R = " << _M_resistance << "\n";
    }

    result = W_out * ((b2*_M_resistance-b1)/(c1-c2*_M_resistance))
    + ((a22*_M_resistance-a11)/(c1-c2*_M_resistance));

    Debug( 6030 ) << "[Resistance::evaluate] a1 = " << a1 << "\n";
    Debug( 6030 ) << "[Resistance::evaluate] b1 = " << b1 << "\n";
    Debug( 6030 ) << "[Resistance::evaluate] c1 = " << c1 << "\n";
    Debug( 6030 ) << "[Resistance::evaluate] a2 = " << a2 << "\n";
    Debug( 6030 ) << "[Resistance::evaluate] b2 = " << b2 << "\n";
    Debug( 6030 ) << "[Resistance::evaluate] c2 = " << c2 << "\n";

    return result;
  }

}

#endif
