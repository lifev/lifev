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
namespace LifeV
{

  Real fZero(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0.0;
  }

  Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0.;
  }
  
  Real u1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    Real v=1.0;
    
    switch(i){
    case 1:
      return -(2.0*v)/(1+2*v*t)*(x-1.0);
      break;
    case 2: 
      return (2.0*v)/(1+2*v*t)*(y-0.5);
    case 3:
      return 0.0;
      break;
    }
    return 0.0;
  }
  
  // Initial velocity 
  Real u0(const Real& t, const Real& x, const Real& y, const Real& z,
	  const ID& i)
  {
    return 0;
  }
  
  
  Real g(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    switch(i) {
    case 1:
    case 2:
      return 0.0;
      break;
    case 3:
      if ( t < 0.005)
	return 1.3332e4;
      else 
	return 0;
      break; 
    }  
    return 0;
  }

  Real bdDisp(const Real& t, const Real& x, const Real& y, const Real& z,
	      const ID& i) {
    
    Real R=sqrt(x*x+y*y);
    Real omega = acos(-1.0)*40;
    
    switch(i){
    case 1:
      return sin(t*omega)*0.02*x/R; 
      break;
    case 2:
      return sin(t*omega)*0.02*y/R;
    case 3:
      return 0.0; 
      break;
      
    }
    return 0.0;
  }
  
}
