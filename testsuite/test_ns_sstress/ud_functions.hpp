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
// Classical waveform of the cardiovascular input (normalized: max=1,min=0.3)
// see calc_g.m
// By A. Veneziani, 1996; in C++: 2002
namespace LifeV
{
Real calc_g(const Real& time)
{
  Real  strokes=72.0;
  Real percar=60.0/strokes;
  Real tfin=percar;
  Real pigreco2=6.283185307179586;
  Real coeff01=0.650,
    coeff02=-0.350,
    coeff11=0.350,
    coeff12=-0.050,
    coeff21=0.30,
    coeff31=0.320,
    coeff41=0.360,
    coeff42=-0.040,
    preprimo,
    primo,
    presecondo,
    secondo,
    coeff32,coeff33,coeff23,coeff22,
    a,a11,a12,a21,a22,b1,b2,det,dt;

    preprimo=0.150*tfin;
  primo=0.20*tfin;
  presecondo=0.30*tfin;
  secondo=0.510*tfin;
  double dumm;
  Real tc=time/percar;
  Real tcorr = modf(tc,&dumm);
       if (tcorr <= preprimo){
         a=pigreco2*tcorr/primo;
         return coeff01+coeff02*cos(a);
       }
       else if ((tcorr  > preprimo) & (tcorr <= primo)){
         b1=coeff01-coeff31;
         b2=coeff02*pigreco2/primo;
         a22=preprimo-primo;
         a12=a22*a22;
         a11=a12*a12;
         a21=4*a12*a22;
         a22=2*a22;
         det=a22*a11-a12*a21;
         coeff32=(a22*b1-a12*b2)/det;
         coeff33=(a11*b2-a21*b1)/det;
         dt=tcorr-primo;
         return coeff32*dt*dt*dt*dt+coeff33*dt*dt+coeff31;
       }
       else if (( tcorr > primo) & (tcorr <= presecondo)) {
         a=pigreco2*tcorr/primo;
         return coeff41+coeff42*cos(a);
       }
       else if ((tcorr > presecondo) & (tcorr <= secondo)){
         a=pigreco2*(tcorr-primo)/primo;
         return coeff11+coeff12*cos(a);
       }
       else if (tcorr > secondo) {
         a=pigreco2*(secondo-primo)/primo;
         b1=coeff11+coeff12*cos(a)-coeff21;
         b2=-coeff12*pigreco2*sin(a)/primo;
         a22=tfin-secondo;
         a12=a22*a22;
         a11=a12*a12;
         a21=-4*a12*a22;
         a22=-2*a22;
         det=a22*a11-a12*a21;
         coeff22=(a22*b1-a12*b2)/det;
         coeff23=(a11*b2-a21*b1)/det;
         dt=tcorr - tfin;
         return coeff22*dt*dt*dt*dt+coeff23*dt*dt+coeff21;
      }
       return 0.;
}

// velocity imposed on the walls
Real u_in(const Real& t, const Real& x, const Real& y, const Real& z,
    const ID& i)
{
  Real coef = calc_g(t);
  Real peak = 150; // mm/s
  switch(i){
  case 3:
    //    cout << "Velocity " << peak*coef << endl;
     return peak*coef;
     break;
  case 1: case 2:
    return 0.;
    break;
  }
  return 0.;
}

//
Real u_w(const Real& t, const Real& x, const Real& y, const Real& z,
    const ID& i)
{
  return 0.;
}
//
Real u_out(const Real& t, const Real& x, const Real& y, const Real& z,
    const ID& i)
{
  return 0.;
}


Real g_in(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  Real coef = calc_g(t);
  return coef;
}

Real g_out(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  return 0;
}

Real g_wall(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  return 0;
}

Real alpha(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  return 0.01; //mm^2*s/g
}

Real pdiri(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  return 0;
}
Real pneum(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i) {
  return 0;
}

Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    return 0.;
}

//Real u1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
//{
//  return 0.0;
//}

// Initial velocity
Real u0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
 switch(i) {
     case 1:
           return 0.0;
           break;
     case 2:
           return 0.0;
           break;
     case 3:
      return (0.25-(x*x+y*y))/0.14;
         break;
        }
      return 0.0;
}

// Initial pressure
Real p0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  return 0.;
}
}
