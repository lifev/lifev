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
#include <life/lifesolver/oneDNonLinModelParam.hpp>

#include <cmath>

namespace LifeV
{
 
OneDNonLinModelParam::OneDNonLinModelParam(const GetPot& dfile) :
    //! we suppose so far that all params are constant along the vessel
    //! otherwise replace _M_paramSize(1) by _M_paramSize(_M_dimDof)
    _M_paramSize(dfile("discretization/nb_elem", 0 ) + 1),
    _M_length(0.),
    _M_xleft(dfile("discretization/x_left", 0. )),
    _M_xright(dfile("discretization/x_right", 0. )),
    _M_terminalResistence(dfile("parameters/terminalResistence", 0. )),
    _M_Area0(_M_paramSize),
    _M_dArea0dz(_M_paramSize),
    _M_AlphaCoriolis(_M_paramSize),
    _M_dAlphaCoriolisdz(_M_paramSize),
    _M_PressBeta0(_M_paramSize),
    _M_dPressBeta0dz(_M_paramSize),
    _M_PressBeta1(_M_paramSize),
    _M_dPressBeta1dz(_M_paramSize),
    _M_FrictionKr(_M_paramSize)
{
  _M_length = _M_xright - _M_xleft;

    Real _A0 = dfile("parameters/Area0",M_PI);
    Real _alpha = dfile("parameters/alphaCor",1./M_ROBERTSON_CORRECTION);
    Real _beta0 = dfile("parameters/beta0",1.e6);
    Real _beta1 = dfile("parameters/beta1",0.5);
    Real _kr    = dfile("parameters/Kr",1.);
    Real _dA0dz = 0.;
    Real _dalphadz = 0.;
    Real _dbeta0dz = 0.;
    Real _dbeta1dz = 0.;
    //-------------------------------------------
    //! Initialisation of the parameter variables
    //-------------------------------------------
    for (UInt iz = 0; iz < _M_paramSize ; iz ++ ) {
        _M_Area0[iz] = _A0;
        _M_AlphaCoriolis[iz] = _alpha;
        _M_PressBeta0[iz] = _beta0;
        _M_PressBeta1[iz] = _beta1;
        _M_FrictionKr[iz] = _kr;
        _M_dArea0dz[iz] = _dA0dz;
        _M_dAlphaCoriolisdz[iz] = _dalphadz;
        _M_dPressBeta0dz[iz] = _dbeta0dz;
        _M_dPressBeta1dz[iz] = _dbeta1dz;
    }
    _M_DensityRho = dfile("parameters/rho",1.);

	if( dfile("parameters/use_physical_values",false) )
	{
		Debug( 6030 ) << "[OneDNonLinModelParam::OneDNonLinModelParam] initializing from physical values\n";
		initParam(dfile);
	}
	
	Debug( 6030 ) << "[OneDNonLinModelParam::OneDNonLinModelParam] Robertson correction coefficient = "
		<< M_ROBERTSON_CORRECTION << "\n";
}

//! so far we assume that all params are constant along the vessel
//! otherwise, replace [0] by [ii]
Real OneDNonLinModelParam::Area0(const UInt& ii) const {
    return _M_Area0[ii]; 
}

Real OneDNonLinModelParam::AlphaCor(const UInt& ii) const {
    return _M_AlphaCoriolis[ii];
}

Real OneDNonLinModelParam::Beta0(const UInt& ii) const {
    return _M_PressBeta0[ii];
}

Real OneDNonLinModelParam::Beta1(const UInt& ii) const {
    return _M_PressBeta1[ii];
}

Real OneDNonLinModelParam::dArea0dz(const UInt& ii) const {
    return _M_dArea0dz[ii]; 
}

Real OneDNonLinModelParam::dAlphaCordz(const UInt& ii) const {
    return _M_dAlphaCoriolisdz[ii];
}

Real OneDNonLinModelParam::dBeta0dz(const UInt& ii) const {
    return _M_dPressBeta0dz[ii];
}

Real OneDNonLinModelParam::dBeta1dz(const UInt& ii) const {
    return _M_dPressBeta1dz[ii];
}

Real OneDNonLinModelParam::FrictionKr(const UInt& ii) const {
    return _M_FrictionKr[ii];
}

Real OneDNonLinModelParam::Celerity0(const UInt& ii) const {
    Real celerity0( _M_PressBeta0[ii] * _M_PressBeta1[ii] / _M_DensityRho );
    return std::sqrt( celerity0 );
}

Real OneDNonLinModelParam::DensityRho() const {
    return _M_DensityRho;
}

Real OneDNonLinModelParam::XLeft() const {
    return _M_xleft;
}

Real OneDNonLinModelParam::XRight() const {
    return _M_xright;
}

Real OneDNonLinModelParam::Length() const {
    return _M_length;
}

Real OneDNonLinModelParam::TerminalResistence() const {
    return _M_terminalResistence;
}

UInt OneDNonLinModelParam::ParamSize() const {
    return _M_paramSize;
}

//! initialisation from physical values
void OneDNonLinModelParam::initParam(const GetPot& dfile) 
{  
  Real Young_modulus = dfile("1d_physics/young",5.e6);
  Real thickness  = dfile("1d_physics/thickness",0.05);
  Real reference_radius = dfile("1d_physics/radius",1.);
  Real density   = dfile("1d_physics/density",1.);  //???
  Real viscosity = dfile("1d_physics/viscosity",0.035);  //???
  Real ksi       = dfile("1d_physics/ksi",0.5);  //???
  bool thick_vessel = dfile("1d_physics/thick_vessel",0);
  
//  Real Coriolis_coeff = dfile("parameters/alphaCor",1);

    //-------------------------------------------
    //! Initialisation of the parameter variables
    //-------------------------------------------
    for (UInt iz = 0; iz < _M_paramSize ; iz ++ ) {
        //! A0
        _M_Area0[iz] = M_PI*reference_radius * reference_radius;
        //! alpha
//        _M_AlphaCoriolis[iz] = Coriolis_coeff;
		if( thick_vessel ){
	        //! beta0
	        _M_PressBeta0[iz] = - thickness * Young_modulus * sqrt(M_PI) / 
            ( density * sqrt(_M_Area0[iz]) *
            	( (1 - ksi * ksi)
            	  + ksi * (1 + ksi) * (thickness * sqrt(M_PI)
            	  					  / sqrt(_M_Area0[iz]) )
            	  )
           	);
	        //! beta1
			_M_PressBeta1[iz] = - 0.5;
		}
		else{
	        //! beta0
	        _M_PressBeta0[iz] = thickness * Young_modulus * sqrt(M_PI) / 
            ( density * sqrt(_M_Area0[iz]) * (1 - ksi * ksi) );
        	//! beta1
			_M_PressBeta1[iz] = 0.5;
		}			
        //! Kr
        _M_FrictionKr[iz] = 8 * M_PI * viscosity;
    }

    //! rho
    _M_DensityRho = density;

}

//! Pseudo-characteristic variables Z1, Z2 corresponding to data (Q, A) 
//! hypothesis alphaCor = 1, beta1 = 0.5 :
/*
 * void
OneDModelSolver::pseudo_char_from_U(const Real& _A, const Real& _Q,
				       const Vec2D _lefteigvec1, const Vec2D _lefteigvec2,
				       Real& Z1, Real& Z2,
				       const UInt& indz ) const
{

  Z1 = _lefteigvec1.first * _A + _lefteigvec1.second * _Q;

  Z2 = _lefteigvec2.first * _A + _lefteigvec2.second * _Q;

}
*/

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
		* ( beta1 / (4 * sqrt(chi) ) * (W1 - W2) + celerity0 )(^2/beta1)

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


// compute the pressure : beta0 * ( ( _A / Area0 )^beta1 - 1 )
Real OneDNonLinModelParam::
pressure(const Real& _A, const UInt& indz) const
{
    return ( _M_PressBeta0[indz]
    		 * ( std::pow( _A / _M_Area0[indz], _M_PressBeta1[indz] ) - 1 ) );
}

/*! Derivative of pressure as a function of A
  dP(A)/dA = beta1 * beta0 * ( _A / Area0 )^beta1 / A
*/
Real OneDNonLinModelParam::
pressureDiff(const Real& _A, const UInt& indz) const
{
	Real AoverA0POWbeta1( std::pow( _A / _M_Area0[indz], _M_PressBeta1[indz] ) );

	return _M_PressBeta0[indz] * _M_PressBeta1[indz] * AoverA0POWbeta1 / _A;
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

//! compute area given the pressure: A = A0 * ( P / beta0 + 1 )^(1/beta1)
Real OneDNonLinModelParam::
A_from_P(const Real& _P, const UInt& indz) const
{
    return ( _M_Area0[indz]
    		 * std::pow( _P / _M_PressBeta0[indz] + 1, 1/_M_PressBeta1[indz] )  );
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
				
	Debug(6030) << "[OneDNonLinModelParam::W_from_P] "
				<< "SQRTchi4overbeta1 = " << SQRTchi4overbeta1
				<< ", beta0beta1overrho = " << SQRTbeta0beta1overrho
				<< ", pow( ( _P / _M_PressBeta0[indz] + 1 ), 0.5 ) = " << pow( ( _P / _M_PressBeta0[indz] + 1 ), 0.5 ) << "\n";
	Debug(6030) << "[OneDNonLinModelParam::W_from_P] add term = " << add << "\n";
	
	if( ii == 1 )
		return _W - add;
    if( ii == 2 )
		return _W + add;

    ERROR_MSG("You can only find W1 or W2 as function of P");
    return -1.;
}


/*! Total pressure (used for interface conditions)
  Pt = P + rho/2 * (Q/A)^2
*/
Real OneDNonLinModelParam::
totalPressure(const Real& _A, const Real& _Q,
              const UInt& indz) const
{
    Real vel( _Q / _A );

    return ( pressure( _A, indz ) + (_M_DensityRho/2) * vel * vel );
}

/*! Derivative of Total pressure (used for interface conditions)
  dPt/dU_ii = dP/dU_ii + rho/2 * d(Q/A)^2/dU_ii
*/
Real OneDNonLinModelParam::
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

void OneDNonLinModelParam::stiffenVesselLeft( const Real& xl, const Real& xr,
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
	  //	      _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] * ( -factor * (n / delta) * 
	  //			     ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
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
	  //	( 1 + factor );
	}
      else 
	std::cout << "[OneDNonLinModelParam::stiffenVesselRight] error! out of left boundary" << std::endl;
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

void OneDNonLinModelParam::stiffenVesselRight( const Real& xl, const Real& xr,
					       const Real& factor, const Real& alpha,
					       const Real& delta, const Real& n,
					       const Real& min_deltax, const UInt& yesAdaptive )
{
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
	  //	      _M_dPressBeta0dz[iz] = _M_PressBeta0[iz] * ( -factor * (n / delta) * 
	  //			     ( pow(2,(n-1)) * pow(ratio,(n-1)) ) );
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
	  //	( 1 + factor );
	}
      else 
	std::cout << "\n[OneDNonLinModelParam::stiffenVesselRight] error! out of right boundary" << std::endl;
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
    c << "dalphaCordz = " << _M_dAlphaCoriolisdz << "\n";
    c << "dbeta0dz = " << _M_dPressBeta0dz << "\n";
    c << "dbeta1dz = " << _M_dPressBeta1dz << "\n";
    c << "dArea0dz    = " << _M_dArea0dz << "\n";
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

Real LinearSimpleParam::Flux11(const UInt& /*ii*/) const {
    return _M_Flux11[0];
}
Real LinearSimpleParam::Flux12(const UInt& /*ii*/) const {
    return _M_Flux12[0];
}
Real LinearSimpleParam::Flux21(const UInt& /*ii*/) const {
    return _M_Flux21[0];
}
Real LinearSimpleParam::Flux22(const UInt& /*ii*/) const {
    return _M_Flux22[0];
}

Real LinearSimpleParam::Celerity1(const UInt& /*ii*/) const {
    return _M_celerity1[0];
}

Real LinearSimpleParam::Celerity2(const UInt& /*ii*/) const {
    return _M_celerity2[0];
}

Real LinearSimpleParam::LeftEigenVector11(const UInt& /*ii*/) const {
    return _M_celer1_left_eigVector1[0];
}

Real LinearSimpleParam::LeftEigenVector12(const UInt& /*ii*/) const {
    return _M_celer1_left_eigVector2[0];
}

Real LinearSimpleParam::LeftEigenVector21(const UInt& /*ii*/) const {
    return _M_celer2_left_eigVector1[0];
}

Real LinearSimpleParam::LeftEigenVector22(const UInt& /*ii*/) const {
    return _M_celer2_left_eigVector2[0];
}

Real LinearSimpleParam::Source10(const UInt& /*ii*/) const {
    return _M_Source10[0];
}
Real LinearSimpleParam::Source20(const UInt& /*ii*/) const {
    return _M_Source20[0];
}
Real LinearSimpleParam::Source11(const UInt& /*ii*/) const {
    return _M_Source11[0];
}
Real LinearSimpleParam::Source12(const UInt& /*ii*/) const {
    return _M_Source12[0];
}
Real LinearSimpleParam::Source21(const UInt& /*ii*/) const {
    return _M_Source21[0];
}
Real LinearSimpleParam::Source22(const UInt& /*ii*/) const {
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
