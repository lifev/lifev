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
#include "user_functor.hpp"

namespace LifeV
{

  //======================================================================
  //   Volumic source
  //======================================================================

  Vol_source::Vol_source(const FhNHandler& _fhn):
    fhn(_fhn)
  {}
  
  double Vol_source::operator()(const double& x, const double& y,
                const double& z, const double& t,
                const ID& i) const
  {
    if(t>=fhn.t_start_LV && t<=fhn.t_stop_LV){
      /*
      if(((x-fhn.source_center_LV(0))*(x-fhn.source_center_LV(0))
      +(y-fhn.source_center_LV(1))*(y-fhn.source_center_LV(1))
      +(z-fhn.source_center_LV(2))*(z-fhn.source_center_LV(2)))
     <= (fhn.source_radius_LV*fhn.source_radius_LV)){
    return fhn.source_value_LV;
      }
      */
      if( (x>100) && (-0.13*x-0.99*y-0.041*z + 99 >0)){
    //cout << fhn.source_value_LV << " ";
    return fhn.source_value_LV;
      }
    }
    if(t>=fhn.t_start_RV && t<=fhn.t_stop_RV){
      /*      if(((x-fhn.source_center_RV(0))*(x-fhn.source_center_LV(0))
      +(y-fhn.source_center_RV(1))*(y-fhn.source_center_LV(1))
      +(z-fhn.source_center_RV(2))*(z-fhn.source_center_LV(2)))
     <= (fhn.source_radius_RV*fhn.source_radius_RV)){
    return fhn.source_value_RV;
    }
      */
      if( (x>100) && (-0.13*x-0.99*y-0.041*z + 99 <=0)){
    //cout << fhn.source_value_RV << " ";
    return fhn.source_value_RV;
      }
    }
    //cout << 0. << " ";
    return 0.;
  }
}
