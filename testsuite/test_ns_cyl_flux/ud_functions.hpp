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

Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return 0.;
}

Real u1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.0;
}

// Initial velocity
Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  Real pi = 3.14159265358979;
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
    
    return 32/pi*(1.0/4-(x*x+y*y));// initial condition for flux=1 for the cylinder
    //return 0.325*32/pi*(1.0/4-(x*x+y*y));// initial condition for physiological flux for the cylinder
     //return 2.5774*(1-16/25*(x*x+y*y));
     //return 0.0;
     break;
  }
  return 0;
}

// Initial velocity for stationary NS (NSo)
Real u0o(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
     return 0.0;
     break;
  }
  return 0;
}

// Initial pressure
Real p0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
 return 1.0;
}

// Neumann conditions for NSo
Real uo(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  //Real pi = 3.14159265358979;
  switch(i) {
  case 1:
  case 2:
    return 0.0;
    break;
  case 3:
    return -1.0;
    break;
  }
  return 0.0;
}

// constant flux
Real my_flux_cost( Real time ){
   return -1.;
  }

// sinusoidal flux
Real my_flux_cos( Real time ){
  Real pi = 3.14159265358979;
  return -1.*cos(2*pi*time);
}

// physiological flux
Real my_flux_physio( Real time ){
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
  Real a,b1,b2,a22,a12,a11,a21,det,dt,coeff22,coeff23,coeff32,coeff33;
  Real flux;
  Real Tcorr;
  Real Taux=Tfin;
  while (Taux<time) {
    Taux=Taux+Tfin;
  }
  Tcorr=time-Taux+Tfin;
  if (Tcorr==Tfin) {
    Tcorr=0;
  }
  if (Tcorr<=prefirst) {
         a=pigreco2*Tcorr/first;
         flux = coeff01+coeff02*cos(a);
  }
  else if ((Tcorr>prefirst)&&(Tcorr<=first)) {
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
         dt=Tcorr-first;
         flux=coeff32*dt*dt*dt*dt+coeff33*dt*dt+coeff31;
  } 
  else if ((Tcorr>first)&&(Tcorr<=presecond)) {
         a=pigreco2*(Tcorr)/first;
         flux = coeff41+coeff42*cos(a);
  }
  else if ((Tcorr>presecond)&&(Tcorr<=second)) {
         a=pigreco2*(Tcorr-first)/first;
         flux = coeff11+coeff12*cos(a);
  }
  else if (Tcorr>second) {
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
         dt=Tcorr-Tfin;
         flux=coeff22*dt*dt*dt*dt+coeff23*dt*dt+coeff21;
  }
  return -flux;
}

}
