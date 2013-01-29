//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the boundary conditions for the Monolithic Test
 *
 *  @date 2009-04-09
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @maintainer Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#include "ud_functions.hpp"

#define ANEURISM100170
//#define AORTA

namespace LifeV
{

Real uInterpolated(const Real& time, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    //GetPot data_file( getpot::StringVector data_file_name );
    //    GetPot data_file(data_file_name.c_str());

    //    Real scale_factor = data_file("fluid/miscellaneous/my_flux_physio_scale_factor", .1);
    Real scale_factor = -50.;//data_file("fluid/miscellaneous/my_flux_physio_scale_factor", .1);

    Real newtime;
    //     Real intervallorampa = data_file("fluid/miscellaneous/timeramp",0.05);
    //     Real deltat = data_file("fluid/discretization/timestep",0.01);
    Real intervallorampa = 0.05; //data_file("fluid/miscellaneous/timeramp",0.05);
    Real deltat = 0.001; //data_file("fluid/discretization/timestep",0.01);

    Real strokes   =  72.0;
    Real percar    =  60.0/strokes;
    Real Tfin      =  percar;
    Real pigreco2  =  6.2831853;
    Real coeff01   =  0.65;
    Real coeff02   =- 0.35;
    Real coeff11   =  0.35;
    Real coeff12   =- 0.05;
    Real coeff21   =  0.3;
    Real coeff31   =  0.32;
    Real coeff41   =  0.36;
    Real coeff42   =- 0.04;
    Real prefirst  =  0.15*Tfin;
    Real first     =  0.2 *Tfin;
    Real presecond =  0.3 *Tfin;
    Real second    =  0.51*Tfin;
    Real a,b1,b2,a22,a12,a11,a21,det,dt,coeff22,coeff23,coeff32,coeff33;
    Real flux(0);
    Real Tcorr;
    Real Taux      =  Tfin;

    if (time < intervallorampa)
        newtime = deltat;
    else newtime = time + deltat - intervallorampa;

    while (Taux < newtime)
    {
        Taux = Taux + Tfin;
    }
    Tcorr = newtime - Taux + Tfin;

    if (Tcorr == Tfin)
    {
        Tcorr=0;
    }

    if (Tcorr <= prefirst)
    {
        a    = pigreco2*Tcorr/first;
        flux = coeff01 + coeff02*cos(a);
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
        dt=Tcorr-first;
        flux=coeff32*dt*dt*dt*dt+coeff33*dt*dt+coeff31;
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
        a       =  pigreco2*(second - first)/first;
        b1      =  coeff11 + coeff12*cos(a) - coeff21;
        b2      =- coeff12*pigreco2*sin(a)/first;
        a22     =  Tfin-second;
        a12     =  a22*a22;
        a11     =  a12*a12;
        a21     =- 4*a12*a22;
        a22     =- 2*a22;
        det     =  a22*a11-a12*a21;
        coeff22 =  (a22*b1-a12*b2)/det;
        coeff23 =  (a11*b2-a21*b1)/det;
        dt      =  Tcorr-Tfin;
        flux    =  coeff22*dt*dt*dt*dt+coeff23*dt*dt+coeff21;
    }
    if (time < intervallorampa)
    {
        flux = ( time/intervallorampa )*scale_factor*flux;
    }
    else
        flux = scale_factor * flux;

    //  Real pi = 3.14159265358979;
#ifdef AORTA       // ifdef ANEURISM100170
    if ( i == 2 )
        return flux; // *1.42?
    if ( i == 1 )
        return flux; // *1.42?
#endif
#ifdef ANEURISM100170       // ifdef ANEURISM100170
    if ( i == 2 )
        return -flux; // *1.42?
#else
    if ( i == 3 )
        return -flux;
#endif
    return 0;

}



// fct_type getUInterpolated()
// {
//     fct_type f;
//     f = boost::bind(&Cylinder::Private::uInterpolated, this, _1, _2, _3, _4, _5);
//     return f;
// }

Real aortaPhisPress(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    /*switch(i) {
      case 1:
      return 0.0;
      break;*/
    //  case 2:
    if (t<=0.00) return -110170;
    if (t<=0.01) return -109540;
    if (t<=0.02) return -108930;
    if (t<=0.03) return -108320;
    if (t<=0.04) return -107710;
    if (t<=0.05) return -107120;
    if (t<=0.06) return -106530;
    if (t<=0.07) return -111130;
    if (t<=0.08) return -115440;
    if (t<=0.09) return -118690;
    if (t<=0.10) return -121460;
    if (t<=0.11) return -123940;
    if (t<=0.12) return -126350;
    if (t<=0.13) return -128890;
    if (t<=0.14) return -131510;
    if (t<=0.15) return -133980;
    if (t<=0.16) return -136200;
    if (t<=0.17) return -138330;
    if (t<=0.18) return -140350;
    if (t<=0.19) return -142290;
    if (t<=0.20) return -144360;
    if (t<=0.21) return -146130;
    if (t<=0.22) return -147530;
    if (t<=0.23) return -148780;
    if (t<=0.24) return -149740;
    if (t<=0.25) return -150320;
    if (t<=0.26) return -150470;
    if (t<=0.27) return -150250;
    if (t<=0.28) return -149750;
    if (t<=0.29) return -148990;
    if (t<=0.30) return -148220;
    if (t<=0.31) return -147210;
    if (t<=0.32) return -145940;
    if (t<=0.33) return -144960;
    if (t<=0.34) return -143750;
    if (t<=0.35) return -141980;
    if (t<=0.36) return -139900;
    if (t<=0.37) return -137260;
    if (t<=0.38) return -133970;
    if (t<=0.39) return -131670;
    if (t<=0.40) return -131320;
    if (t<=0.41) return -133150;
    if (t<=0.42) return -132710;
    if (t<=0.43) return -131570;
    if (t<=0.44) return -130280;
    if (t<=0.45) return -129750;
    if (t<=0.46) return -129330;
    if (t<=0.47) return -128910;
    if (t<=0.48) return -128360;
    if (t<=0.49) return -127680;
    if (t<=0.50) return -127000;
    if (t<=0.51) return -126410;
    if (t<=0.52) return -125920;
    if (t<=0.53) return -125480;
    if (t<=0.54) return -125040;
    if (t<=0.55) return -124560;
    if (t<=0.56) return -124050;
    if (t<=0.57) return -123530;
    if (t<=0.58) return -123000;
    if (t<=0.59) return -122440;
    if (t<=0.60) return -121840;
    if (t<=0.61) return -121220;
    if (t<=0.62) return -120580;
    if (t<=0.63) return -119950;
    if (t<=0.64) return -119330;
    if (t<=0.65) return -118710;
    if (t<=0.66) return -118100;
    if (t<=0.67) return -117470;
    if (t<=0.68) return -116840;
    if (t<=0.69) return -116200;
    if (t<=0.70) return -115560;
    if (t<=0.71) return -114920;
    if (t<=0.72) return -114280;
    if (t<=0.73) return -113650;
    if (t<=0.74) return -113020;
    if (t<=0.75) return -112400;
    if (t<=0.76) return -111790;
    if (t<=0.77) return -111200;
    if (t<=0.78) return -110620;
    if (t<=0.79) return -110060;
    //    break;
    /*  case 3:
        return 0.0;
        break;}
        return 0.;*/
    return 0.;
}

Real f(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.;
}

Real u1(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real fZero(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

// Initial velocity
Real u0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real p0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}


Real E(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -10000;  // (see paper by Liu, Dang, etc.. about the sourrounding tissue effect on arteries)
}


Real hydro(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -1.33*1e5;
}

Real u2(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
    case 0:
        return 0.0;
        break;
    case 2:
        if ( t <= 0.003 )
            return 1.5;
        //      return 0.01;
        return 0.0;
        break;
    case 1:
        return 0.0;
        //      return 1.3332e4;
        //    else
        //      return 0.0;
        break;
    }
    return 0;
}


Real u2normal(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    //   if (t<=0.003)
        return -1.3332e4;

    return 0.;
}


// Initial displacement and velocity
Real d0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
    case 0:
        return 0.;
        break;
    case 1:
        return 0.;
        break;
    case 2:
        return 0.;
        break;
    default:
        ERROR_MSG("This entry is not allowed: ud_functions.hpp");
        return 0.;
        break;
    }
}

Real w0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch (i)
    {
    case 0:
        return 0.0;
        break;
    case 1:
        return 0.0;
        break;
    case 2:
        return 10.0;
        break;
    default:
        ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
        return 0.;
        break;
    }
}

Real fluxFunction(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -100;
}

Real fluxFunctionAneurysm(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{

  Real fluxFinal;
  Real rampAmpl(0.4);
  Real dt(0.001);

  if ( t <= 0.4 )
    fluxFinal = (1.8334/0.4) * t;
  else
    {


      // We change the flux for our geometry
      const Real pi         = 3.141592653589793;
      //Simone's area
      //const Real area       = 0.0907122;
      //const Real area       = 0.0777195; //FluidSmooth
      const Real area = 0.193529; // BigMesh

      const Real areaFactor = area/((0.6/2)*(0.6/2)*pi);
      const Real Average = (48.21 * pow(area,1.84)) * 60; //Mean Cebral's Flux per minut

      // Unit conversion from ml/min to cm^3/s
      const Real unitFactor = 1./ 60.;

      // T is the period of the cardiac cycle
      const Real T          = 0.8;

      // a0 is the average VFR (the value is taken from Karniadakis p970)
      const Real a0         = 255;
      const Real volumetric = Average/a0; //VolumetricFactor per minut

      // Fourrier
      const Int M(7);
      const Real a[M] = {-0.152001,-0.111619, 0.043304,0.028871,0.002098,-0.027237,-0.000557};
      const Real b[M] = { 0.129013,-0.031435,-0.086106,0.028263,0.010177, 0.012160,-0.026303};

      Real flux(0);
      //      const Real xi(2*pi*t/T);
      const Real xi(2*pi*(t-rampAmpl+dt)/T);

      flux = a0;
      Int k(1);
      for (; k<=M ; ++k)
        flux += a0*(a[k-1]*cos(k*xi) + b[k-1]*sin(k*xi));

      //return - (flux * areaFactor * unitFactor);
      fluxFinal =  (flux * volumetric * unitFactor);

    }

  return fluxFinal;

}

Real aneurismFluxInVectorial(const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  //Components for Simone's mesh
  /*
    Real n1(-0.000019768882940);
    Real n2(-0.978289616345544);
    Real n3( 0.207242433298975);
  */
  /*
    Real n1(-0.67460);
    Real n2(-0.19888);
    Real n3(-0.17089);
  */
    Real n1(-0.671731456);
    Real n2(-0.204172511);
    Real n3(-0.712102827);


    Real flux(fluxFunctionAneurysm(t,x,y,z,i));
    //Simone's area
    //Real area(0.0907122); // Computed with the triangle on the INLET boundary
    
    //Real area(0.0777195); //fluidsmooth
    Real area(0.193529); //fluidBig
    switch(i) {
    case 0:
        return n1*flux/area;
    case 1:
        return n2*flux/area;
    case 2:
        return n3*flux/area;
    default:
        return 0.0;
    }
}



Real squareSinusoidalFluxFunction(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -(t<(0.005/2))*std::sin(2*M_PI*t/0.005)*std::sin(2*M_PI*t/0.005);
}

Real benchmarkP(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
if(t < 0.0025)
    return -50000*std::sin(2*3.1415*t/0.005)*std::sin(2*3.1415*t/0.005);
else return 0;
}

Real LifeV::aortaVelIn::S_timestep;

}


