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
#include "user_fct.hpp"
#include "dataFhN.hpp"
#define _TGV_ 1.e9

namespace LifeV
{
  double zero(const double& t, const double& x, const double& y, const double& z, const ID& i) {
    return 0.;
  }
  double one(const double& t, const double& x, const double& y, const double& z, const ID& i) {
    return 1.;
  }

  //======================================================================
  // Stimulation coefficient (Robin b.c.)
  //======================================================================
  
  double stim_coef(const double& t,
		   const double& x, const double& y, const double& z,
		   const ID& i) {
    /*
      N.B. : this term is in the matrix -> independant of time (in
      the present implementation)
    */
    double& x1 = DataFhN::stim_center_1(0);
    double& y1 = DataFhN::stim_center_1(1);
    double& z1 = DataFhN::stim_center_1(2);
    double& r1 = DataFhN::stim_radius_1;
    if(fabs(x-x1) < r1 &&
       fabs(y-y1) < r1 &&
       fabs(z-z1) < r1 ){
      return _TGV_;
    }
    // two sites of stimulation
    double& x2 = DataFhN::stim_center_2(0);
    double& y2 = DataFhN::stim_center_2(1);
    double& z2 = DataFhN::stim_center_2(2);
    double& r2 = DataFhN::stim_radius_2;
    if(fabs(x-x2) < r2 &&
       fabs(y-y2) < r2 &&
       fabs(z-z2) < r2 ){
      return _TGV_;
    }
    return 0.;
  }
  double stim_g(const double& t,
		const double& x, const double& y, const double& z,
		const ID& i) {
    /*
      N.B. : this term is in the r.h.s which is updated at each time step
    */
    double& x1 = DataFhN::stim_center_1(0);
    double& y1 = DataFhN::stim_center_1(1);
    double& z1 = DataFhN::stim_center_1(2);
    double& r1 = DataFhN::stim_radius_1;
    double& v1 = DataFhN::stim_value_1;
    double& t1_start = DataFhN::stim_start_1;
    double& t1_stop = DataFhN::stim_stop_1;
    if(fabs(x-x1) < r1 &&
       fabs(y-y1) < r1 &&
       fabs(z-z1) < r1 &&
       t < t1_stop &&
       t >= t1_start ){
      return v1*_TGV_;
    }
    // two sites of stimulation
    double& x2 = DataFhN::stim_center_2(0);
    double& y2 = DataFhN::stim_center_2(1);
    double& z2 = DataFhN::stim_center_2(2);
    double& r2 = DataFhN::stim_radius_2;
    double& v2 = DataFhN::stim_value_2;
    double& t2_start = DataFhN::stim_start_2;
    double& t2_stop = DataFhN::stim_stop_2;
    if(fabs(x-x2) < r2 &&
       fabs(y-y2) < r2 &&
       fabs(z-z2) < r2 &&
       t < t2_stop &&
       t >= t2_start ){
      return v2*_TGV_;
    }
    return 0.;
  }
  
  

  //======================================================================
  //   Initial data
  //======================================================================
  double u_init1(const double& x,const double& y,const double& z)
  {
    if(z>=1 && z<=20){
      //      return 0.9;
      return 0.;
    } else {
      return 0.;
    }
  }

  double u_init2(const double& x,const double& y,const double& z)
  {
    if(x>=110){
      return 0.;
    } else {
      return 0.;
    }
  }


  //=====================================================================
  // diffusion
  //======================================================================
  double  diff_fct_1(double x,double y,double z)
  {
    //    cout << x << " " << y << " " << z << endl;
    if( ( sqrt((94.-x)*(94.-x) + (98.-y)*(98.-y) + (114.-z)*(114.-z))) < 30. )
      {
	return .0000001;
      } else {
	return 1.;
      }
  }
  
  //=====================================================================
  // selection of points for postprocessing
  //======================================================================
bool select_points1(double x, double y,double z) {
    if(fabs(x-5)<1e-3 && fabs(y-5)<1e-3) return true;
    else return false;
  }

}
