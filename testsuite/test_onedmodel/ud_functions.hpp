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

  class Sin : public OneDBCFunctionBase
  {
  public:
    Sin(const Real mean = 0, const Real scale = 10, const Real period = 1 ) :
      _M_mean(mean), _M_scale(scale), _M_period(period) {}
      
    Real evaluate( const Real& time ) {return _M_mean+_M_scale*std::sin(2*M_PI*time/_M_period);};
    
    ~Sin() {}
  private:
    Real _M_mean;
    Real _M_scale;
    Real _M_period;
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



  class PressureRamp : public Compatibility
  {
  public:
    PressureRamp(const BasicOneDMesh& mesh,
  				const NonLinearFluxFun1D& fluxFun,
  				const NonLinearSourceFun1D& sourceFun,	
  				const ScalVec& U1_thistime, const ScalVec& U2_thistime,  
  				const ScalVec& W1_thistime, const ScalVec& W2_thistime,  
				const Real& dt, const std::string& border, const std::string & var,
				const OneDNonLinModelParam& onedparam,
				const Real& startT = 0., const Real& duration = 0.05, const Real& endvalue = 106400 );
      
    Real evaluate( const Real& time );
    
    ~PressureRamp() {}
  private:
	//! Reference to the solver parameters
    const OneDNonLinModelParam& _M_onedparam;    

    Real _M_startT;
    Real _M_duration;
    Real _M_endvalue;
  };



  class Heart : public Compatibility
  {
  public:
    Heart ( const GetPot& data_file,
    		const OneDNonLinModelParam& onedparam,
  			const BasicOneDMesh& mesh,
    		const NonLinearFluxFun1D& fluxFun,
    		const NonLinearSourceFun1D& sourceFun,
    		const ScalVec& U1_thistime, const ScalVec& U2_thistime,  
    		const ScalVec& W1_thistime, const ScalVec& W2_thistime,  
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
    const OneDNonLinModelParam& _M_onedparam;    

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
//  Vec2D U_boundary, U_internalBd;   
    //*************************************    
  };



class Resistence : public Compatibility
{
public:
	Resistence(  const Real & resistence, // const GetPot& data_file,
			const OneDNonLinModelParam& onedparam,
  			const BasicOneDMesh& mesh,
    		const NonLinearFluxFun1D& fluxFun,
    		const NonLinearSourceFun1D& sourceFun,
    		const ScalVec& U1_thistime, const ScalVec& U2_thistime,  
    		const ScalVec& W1_thistime, const ScalVec& W2_thistime,  
			const Real& dt, const std::string& border, const std::string & var );   
	
	Real evaluate( const Real& time );
	
	~Resistence() {}
	private:
	    //! Resistence value   
    Real _M_resistence;
	//! Reference to the solver parameters
    const OneDNonLinModelParam& _M_onedparam;    
};

}

#endif
