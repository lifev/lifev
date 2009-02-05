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
  \file ud_functions.cpp
  \author Lucia Mirabella
  \date 04/2007
  \version 1.0

  \brief This file contains functions to impose boundary conditions for 1D tubes.

*/

#include "ud_functions.hpp"

namespace LifeV
{

  PressureRamp::PressureRamp (const BasicOneDMesh& mesh,
  				const NonLinearFluxFun1D& fluxFun,
  				const NonLinearSourceFun1D& sourceFun,
  				const ScalVec& U1_thistime, const ScalVec& U2_thistime,
  				const ScalVec& W1_thistime, const ScalVec& W2_thistime,
				const Real& dt, const std::string& border, const std::string & var,
				const OneDNonLinModelParam& onedparam,
				const Real& startT, const Real& duration, const Real& endvalue ) :
    Compatibility( mesh, fluxFun, sourceFun, U1_thistime, U2_thistime, W1_thistime, W2_thistime, dt, border, var),
    _M_onedparam( onedparam ),
	_M_startT( startT ),
	_M_duration( duration ),
	_M_endvalue( endvalue )
	{}

  Real PressureRamp::evaluate( const Real& time )
  {
    Real W_out, result;

	Real _P = ( time < (_M_startT + _M_duration) ) ?
				( ( time - _M_startT ) / _M_duration ) : 1;
	_P *= _M_endvalue;

	Debug( 6030 ) << "[PressureRamp::evaluate] imposed pressure = " << _P << "\n";
	Debug( 6030 ) << "[PressureRamp::evaluate] target pressure = " << _M_endvalue << "\n";

    switch( _M_oneDBCFunctionsMapStringValues[_M_var] )
      {
      case OneDBCW1:
	W_out = extrapolate_W( OneDBCW2 );
	result = _M_onedparam.W_from_P( _P, W_out, 2, _M_boundaryDof);
	break;
      case OneDBCW2:
	W_out = extrapolate_W( OneDBCW1 );
	result = _M_onedparam.W_from_P( _P, W_out, 1, _M_boundaryDof);
	break;
      default:
	std::cout << "\n[PressureRamp::evaluate] incorrect variable identifier: " << _M_var << std::endl;
      }

	Debug( 6030 ) << "[PressureRamp::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

    return result;

  }


  Heart::Heart ( const GetPot& data_file, const OneDNonLinModelParam& onedparam,
  				const BasicOneDMesh& mesh,
  				const NonLinearFluxFun1D& fluxFun,
  				const NonLinearSourceFun1D& sourceFun,
  				const ScalVec& U1_thistime, const ScalVec& U2_thistime,
  				const ScalVec& W1_thistime, const ScalVec& W2_thistime,
				const Real& dt, const std::string& border, const std::string & var,
				const bool type, const Real& startT ):
    Compatibility( mesh, fluxFun, sourceFun, U1_thistime, U2_thistime, W1_thistime, W2_thistime, dt, border, var),
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
	    _M_Q_old=_M_U2_thistime(_M_boundaryDof);
	Debug( 6030 ) << "[Heart::Heart] U1_thistime(_M_boundaryDof) = " << _M_U1_thistime(_M_boundaryDof) << "\n";

    }



  Real Heart::evaluate( const Real& time )
  {
    Real W_out, result;
    UInt W_out_id;

	Debug( 6030 ) << "[Heart::evaluate] U1_thistime(_M_boundaryDof) = " << _M_U1_thistime(_M_boundaryDof) << "\n";
	Debug( 6030 ) << "[Heart::evaluate] _M_U_boundary.first = " << _M_U_boundary.first << "\n";

    Real T_reset( ( time - _M_startT ) );
    T_reset -= static_cast<int>( std::floor( (T_reset + _M_time_step/2) /_M_periodotot ) ) * _M_periodotot;

    if (_M_type==0)
      Pv=PvFunc(T_reset);
    else
      Pv=PvCalc(T_reset);

    switch( _M_oneDBCFunctionsMapStringValues[_M_var] )
      {
      case OneDBCW1:
	W_out = extrapolate_W( OneDBCW2 );
	W_out_id = 2;
		break;
      case OneDBCW2:
	W_out = extrapolate_W( OneDBCW1 );
	W_out_id = 1;
	break;
      default:
	std::cout << "\n[Heart::evaluate] incorrect variable identifier: " << _M_var << std::endl;
      }

	Debug( 6030 ) << "[Heart::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

    Pa=_M_onedparam.pressure(_M_U_boundary.first, _M_boundaryDof);

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
	result = _M_onedparam.W_from_P( Pv, W_out, W_out_id, _M_boundaryDof);
	if (_M_U_boundary.second<0) _M_OV=0;
      }

    return result;

  }


  Real Heart::PvFunc(const Real& T_reset) const
  {
	Debug(6030)<<"[Heart::PvFunc] T_reset" << T_reset << "\n";
	Debug(6030)<<"[Heart::PvFunc] _M_periodosis" << _M_periodosis << "\n";
   if (T_reset < _M_periodosis)
      return 160000 * std::sin(M_PI/_M_periodosis * T_reset);
    else
      return 0;

  }

  Real Heart::PvCalc(const Real& T_reset)
  {

    El=Elastance(T_reset);
    Co=Compliance(T_reset);
    Debug(6030)<<"[Heart::PvCalc] Co = "<<Co<<"\n";
    if (T_reset == 0)
      {
	_M_Vol_old=120;
      }
    _M_Q_new=_M_U_boundary.second;
    _M_Vol_new=_M_Vol_old - (_M_Q_new+_M_Q_old)/2 * _M_time_step;
    _M_Vol_old=_M_Vol_new;
    _M_Q_old=_M_Q_new;

    Debug(6030)<< "[Heart::PvCalc] heart volume = " << _M_Vol_new << "\n";
    //return 1/Co*(_M_Vol_new-V0);
    return El*(_M_Vol_new-V0);

  }


  Real Heart::Elastance(const Real& T_reset) const
  {

    if (T_reset< _M_periodosis)
      return 1333* ( 1 + 0.7 * std::sin(  3*M_PI/(2*_M_periodosis) * T_reset ) );//spiegare i numeri (LM 2006)
    else
      return 0.3 + 0.7 * (T_reset - _M_periodosis)/(_M_periodotot-_M_periodosis);
  }

  Real Heart::Compliance(const Real& T_reset) const
  {

     if (T_reset< _M_periodosis)
      return 0.11*pow((0.00003/0.0146),((1.-exp(-T_reset/(0.0025*60.)))/(1.-exp(-_M_periodosis/(0.0025*60.)))));//spiegare i numeri
    else
      return 2.25e-4*pow((0.0146/0.00003),((1.-exp(-(T_reset-_M_periodosis)/(0.0075*60.)))/(1-exp(-(_M_periodotot-_M_periodosis)/(0.0075*60.)))));

  }


PhysiologicalFlux::PhysiologicalFlux( GetPot const& data_file )
{
    _M_rampT = data_file("rampT",0.05);
    _M_time_step = data_file("timestep",0.001);
    _M_scale = data_file("scale",10);
}

    Real PhysiologicalFlux::evaluate( const Real& time )
    {
    Real newtime;

    Real strokes=72.0;
    Real percar=60.0/strokes;
    Real Tfin=percar;
    Real pigreco2=6.2831853;
    Real coeff01=0.65;
    Real coeff02=-0.35;
    Real coeff11=0.35;
    Real coeff12=-0.05;
    Real coeff21=0.3;
    Real coeff31=0.32;
    Real coeff41=0.36;
    Real coeff42=-0.04;
    Real prefirst=0.15*Tfin;
    Real first=0.2*Tfin;
    Real presecond=0.3*Tfin;
    Real second=0.51*Tfin;
    Real a,b1,b2,a22,a12,a11,a21,det,ddtt,coeff22,coeff23,coeff32,coeff33;
    Real flux = 0;
    Real Tcorr;
    Real Taux=Tfin;

    if (time<_M_rampT)
      newtime=_M_time_step;
    else
      newtime=time+_M_time_step-_M_rampT;

    while (Taux<newtime) Taux=Taux+Tfin;

    Tcorr=newtime-Taux+Tfin;
    if (Tcorr==Tfin) Tcorr=0;

    if (Tcorr<=prefirst)
      {
	a=pigreco2*Tcorr/first;
	flux = coeff01+coeff02*cos(a);
      }
    else if ((Tcorr>prefirst)&&(Tcorr<=first))
      {
	b1=coeff01-coeff31;
	b2=coeff02*pigreco2/first;
	a22=prefirst-first;
	a12=a22*a22;
	a11=a12*a12;
	a21=4*a12*a22;
	a22=2*a22;
	det=a22*a11-a12*a21;
	coeff32=(a22*b1-a12*b2)/det;
	coeff33=(a11*b2-a21*b1)/det;
	ddtt=Tcorr-first;
	flux=coeff32*ddtt*ddtt*ddtt*ddtt+coeff33*ddtt*ddtt+coeff31;
      }
    else if ((Tcorr>first)&&(Tcorr<=presecond))
      {
	a=pigreco2*(Tcorr)/first;
	flux = coeff41+coeff42*cos(a);
      }
    else if ((Tcorr>presecond)&&(Tcorr<=second))
      {
	a=pigreco2*(Tcorr-first)/first;
	flux = coeff11+coeff12*cos(a);
      }
    else if (Tcorr>second)
      {
	a=pigreco2*(second-first)/first;
	b1=coeff11+coeff12*cos(a)-coeff21;
	b2=-coeff12*pigreco2*sin(a)/first;
	a22=Tfin-second;
	a12=a22*a22;
	a11=a12*a12;
	a21=-4*a12*a22;
	a22=-2*a22;
	det=a22*a11-a12*a21;
	coeff22=(a22*b1-a12*b2)/det;
	coeff23=(a11*b2-a21*b1)/det;
	ddtt=Tcorr-Tfin;
	flux=coeff22*ddtt*ddtt*ddtt*ddtt+coeff23*ddtt*ddtt+coeff21;
      }

    if (time<_M_rampT) flux=(time/_M_rampT)*flux;

    return _M_scale*flux;
    }



Resistence::Resistence(  const Real & resistence,
			const OneDNonLinModelParam& onedparam, // const GetPot& data_file,
  			const BasicOneDMesh& mesh,
    		const NonLinearFluxFun1D& fluxFun,
    		const NonLinearSourceFun1D& sourceFun,
    		const ScalVec& U1_thistime, const ScalVec& U2_thistime,
    		const ScalVec& W1_thistime, const ScalVec& W2_thistime,
			const Real& dt, const std::string& border, const std::string & var ) :
    Compatibility( mesh, fluxFun, sourceFun, U1_thistime, U2_thistime, W1_thistime, W2_thistime, dt, border, var),
   _M_resistence(resistence),
   _M_onedparam( onedparam )
{
	Debug( 6030 ) << "[Resistence::Resistence] resistence = " << _M_resistence << "\n";
}



  Real Resistence::evaluate( const Real& time )
  {
    //! Coefficients
    Real W_out, result;
    Real a1, a2, a11, a22, b1, b2, c1, c2;

	update_U_boundary();
	update_U_internalBd();

	Debug( 6030 ) << "[Resistence::Resistence] at node " << _M_boundaryDof
	<< ", A = " << _M_U_boundary.first << "( " << _M_U1_thistime[0] << " ) "
	<< ", Q = " << _M_U_boundary.second
	<< ", W1 = " << _M_W_boundary.first
	<< ", W2 = " << _M_W_boundary.second
	<< "\n";

	computeEigenValuesVectors();

    a1 = _M_onedparam.pressure(_M_U_boundary.first, _M_boundaryDof); // pressure at previou time step

    a2 = _M_U_boundary.second; // flux at previous time step

    b1 = _M_onedparam.pressure_WDiff( _M_W_boundary.first, _M_W_boundary.second, 1, _M_boundaryDof);  // dP / dW1

    b2 = _M_U_boundary.first / 2; // dQ / dW1

    c1 = _M_onedparam.pressure_WDiff( _M_W_boundary.first, _M_W_boundary.second, 2, _M_boundaryDof);  // dP / dW2

    c2 = b2; // dQ / dW2

	Debug( 6030 ) << "[Resistence::evaluate] P(A) = " << a1 << "\n";
	Debug( 6030 ) << "[Resistence::evaluate] P(W1,W2) = "
				<< _M_onedparam.pressure_W(_M_W_boundary.first, _M_W_boundary.second, _M_boundaryDof) << "\n";

    a11 = a1 - b1*_M_W_boundary.first - c1*_M_W_boundary.second;
    a22 = a2 - b2*_M_W_boundary.first - c2*_M_W_boundary.second;

    switch( _M_oneDBCFunctionsMapStringValues[_M_var] )
      {
      case OneDBCW1:
	W_out = extrapolate_L_dot_U(_M_eigval2, _M_left_eigvec2)
		- dot( _M_left_eigvec2, _M_U_boundary ) + _M_W_boundary.second;

	break;
      case OneDBCW2:
	W_out = extrapolate_L_dot_U(_M_eigval1, _M_left_eigvec1)
		- dot( _M_left_eigvec1, _M_U_boundary ) + _M_W_boundary.first;

	break;
      default:
	std::cout << "\n[Resistence::evaluate] incorrect variable identifier: " << _M_border << std::endl;
      }

	Debug( 6030 ) << "[Resistence::evaluate] extrapolated exiting characteristic = " << W_out << "\n";

    result = W_out * ((b2*_M_resistence-b1)/(c1-c2*_M_resistence))
    + ((a22*_M_resistence-a11)/(c1-c2*_M_resistence));

	Debug( 6030 ) << "[Resistence::evaluate] c1-c2*_M_resistence = " << c1-c2*_M_resistence << "\n";

	Debug( 6030 ) << "[Resistence::evaluate] c1 = " << c1 << "\n";

	Debug( 6030 ) << "[Resistence::evaluate] c2 = " << c2 << "\n";

    return result;
  }

}
