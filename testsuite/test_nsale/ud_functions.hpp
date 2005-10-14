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

#include <math.h>

using namespace std;

namespace LifeV
{
  
   Real g_aneurism(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
   {

     
     // inlet normal aneurism
     Real n1 =  0.1040890;
     Real n2 = -0.9298549;
     Real n3 =  0.3528957;
     Real w=1;
     
     if ( t > 0.005 ) w=0;
     
     switch(i) {
     case 1:
       return -1.3320e4*w*n1;
       break;
     case 2:
       return -1.3320e4*w*n2;
       break;
     case 3:
       return -1.3320e4*w*n3;
       break;
     }
     
     return 0;
   }

    Real g_caronew(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
   {

     
     // inlet normal caronew
     Real n1 =   0.1845595;
     Real n2 =   0.9784991;
     Real n3 = - 0.0920717;

     Real w=1;
     
     if ( t > 0.005 ) w=0;
     
     switch(i) {
     case 1:
       return -1.3320e4*w*n1;
       break;
     case 2:
       return -1.3320e4*w*n2;
       break;
     case 3:
       return -1.3320e4*w*n3;
       break;
     }
     
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
	return 13320.0;
      else
	return 0;
      break;
    }
    return 0;
  }


 Real g_branch(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  { 

    double t_etoile = 0.8;
    double val_p_in=0.;
    double modtps = fmod(t,t_etoile);
    double p1 = 200;
    double p2 = 0.;
    

    double t1 = t_etoile*0.3;
    double t2 = t_etoile*0.35;
    double t3 = t_etoile*0.95;
    double t4 = t_etoile;

    
    if (modtps>=0. && modtps <= t1){
      val_p_in = p1;
    }
    if (modtps>=t1 && modtps <= t2){
      val_p_in = ( (modtps - t1)*p2 + (t2-modtps)*p1 )
	/ (t2-t1);     }
    if (modtps>=t2 && modtps <= t3){
      val_p_in = p2;
    }
    if (modtps>=t3 && modtps <= t4){
      val_p_in = ( (modtps - t3)*p1 + (t4-modtps)*p2 )
	/ (t4-t3);
    }
    
    switch(i) {
    case 1:
    case 3:
      return 0.0;
      break;
    case 2:
      return val_p_in;
      break;
    }
    return 0;
  }








  
  
  Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0.;
  }

  

  
  // Initial velocity 
  Real u0(const Real& t, const Real& x, const Real& y, const Real& z,
      const ID& i)
  {
    return 0.0;

  }
  
  
  Real fZero(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    return 0.;
  }


  Real qND(const Real& t)
  {
  if (t<1.0/7.0) {
    return
      10.9506162830641998*t-193.580197870147003*t*t*t;
  }
  else if ((t>=1.0/7.0)&&(t<2.0/7.0)) {
     return
       1.12874750944692-.901232566128460033*t-82.9629419443489980*(t-1./7.)*(t-1./7.)+281.900989350737973*(t-1./7.)*(t-1./7.)*(t-1./7.);
  }
  else if ((t>=2.0/7.0)&&(t<3.0/7.0)) {
    return
      2.09876743387152-7.34568601855032987*t+37.8517677773960983*(t-2./7.)*(t-2./7.)-21.3007595328068007*(t-2./7.)*(t-2./7.)*(t-2./7.);
  }
  else if ((t>=3.0/7.0)&&(t<4.0/7.0)) {
    return
     -1.26684713156991+2.16497664032978986*t+28.7228708347645991*(t-3./7.)*(t-3./7.)-139.416951219511986*(t-3./7.)*(t-3./7.)*(t-3./7.);
  }
  else if ((t>=4.0/7.0)&&(t<5.0/7.0)) {
     return
       -.899016832703540+1.83577945723119007*t-31.0272511164547993*(t-4./7.)*(t-4./7.)+94.3095644108550033*(t-4./7.)*(t-4./7.)*(t-4./7.);
  } 
  else if ((t>=5.0/7.0)&&(t<6.0/7.0)) {
    return
     .950496049467536-1.25509446925455004*t+9.39113363105458988*(t-5./7.)*(t-5./7.)-22.7603064239091992*(t-5./7.)*(t-5./7.)*(t-5./7.);
  }
  else if ((t>=6.0/7.0)&&(t<1.)) {
     return
       -0.296557883888689e-1+0.345984197870137026e-1*t-.363283407763643984*(t-6./7.)*(t-6./7.)+.847661284781833002*(t-6./7.)*(t-6./7.)*(t-6./7.);
  }
  }
  
  
  Real u_in(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
  {
    Real R=0.8; 		//Tube radius (cm)
    Real T=0.84;			//Period (s)
    Real pi=acos(-1.0);
    Real t_loc=t/T-int(t/T);
    
    switch(i) {
    case 1:
    case 2:
      return 0.0;
      break;
    case 3:  	// profile in tanh
    //			std::cout << t_loc << "  " <<qND(t_loc)  << std::endl;
      return 
	126.18*qND(t_loc)/pi/(R*R)*20/18*(1-(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R*(x*x+y*y)/R/R);
    break;
    }
  return 0.0;
  
  }
  
  //u9 profile: 126.18*qND(t_loc)/pi/(R*R)*11/9*(1-sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R*sqrt(x*x+y*y)/R);
  //profile tanh: 126.18*qND(t_loc)/2.0165035/(R*R)*tanh(3.141592654*(1-sqrt(x*x+y*y)/R));
  //Flow rate (cm3/s)(126.18 for Rep=2700)
  //(m3/s)(2.51E-5 for Rem=500, 1.759e-5 for 300)  
}
