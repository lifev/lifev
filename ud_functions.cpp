/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

   <short description here>

   Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file ud_functions.cpp
 */

#include "ud_functions.hpp"

//#define ANEURISM100170
#define AORTA
namespace LifeV
{



Real uInterpolated(const Real& time, const Real& x, const Real& y, const Real& z, const ID& i)
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
    Real flux;
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

Real aortaPhisPress(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    /*switch(i) {
      case 1:
      return 0.0;
      break;*/
    //  case 2:
    if(t<=0.00) return -110170;
    if(t<=0.01) return -109540;
    if(t<=0.02) return -108930;
    if(t<=0.03) return -108320;
    if(t<=0.04) return -107710;
    if(t<=0.05) return -107120;
    if(t<=0.06) return -106530;
    if(t<=0.07) return -111130;
    if(t<=0.08) return -115440;
    if(t<=0.09) return -118690;
    if(t<=0.10) return -121460;
    if(t<=0.11) return -123940;
    if(t<=0.12) return -126350;
    if(t<=0.13) return -128890;
    if(t<=0.14) return -131510;
    if(t<=0.15) return -133980;
    if(t<=0.16) return -136200;
    if(t<=0.17) return -138330;
    if(t<=0.18) return -140350;
    if(t<=0.19) return -142290;
    if(t<=0.20) return -144360;
    if(t<=0.21) return -146130;
    if(t<=0.22) return -147530;
    if(t<=0.23) return -148780;
    if(t<=0.24) return -149740;
    if(t<=0.25) return -150320;
    if(t<=0.26) return -150470;
    if(t<=0.27) return -150250;
    if(t<=0.28) return -149750;
    if(t<=0.29) return -148990;
    if(t<=0.30) return -148220;
    if(t<=0.31) return -147210;
    if(t<=0.32) return -145940;
    if(t<=0.33) return -144960;
    if(t<=0.34) return -143750;
    if(t<=0.35) return -141980;
    if(t<=0.36) return -139900;
    if(t<=0.37) return -137260;
    if(t<=0.38) return -133970;
    if(t<=0.39) return -131670;
    if(t<=0.40) return -131320;
    if(t<=0.41) return -133150;
    if(t<=0.42) return -132710;
    if(t<=0.43) return -131570;
    if(t<=0.44) return -130280;
    if(t<=0.45) return -129750;
    if(t<=0.46) return -129330;
    if(t<=0.47) return -128910;
    if(t<=0.48) return -128360;
    if(t<=0.49) return -127680;
    if(t<=0.50) return -127000;
    if(t<=0.51) return -126410;
    if(t<=0.52) return -125920;
    if(t<=0.53) return -125480;
    if(t<=0.54) return -125040;
    if(t<=0.55) return -124560;
    if(t<=0.56) return -124050;
    if(t<=0.57) return -123530;
    if(t<=0.58) return -123000;
    if(t<=0.59) return -122440;
    if(t<=0.60) return -121840;
    if(t<=0.61) return -121220;
    if(t<=0.62) return -120580;
    if(t<=0.63) return -119950;
    if(t<=0.64) return -119330;
    if(t<=0.65) return -118710;
    if(t<=0.66) return -118100;
    if(t<=0.67) return -117470;
    if(t<=0.68) return -116840;
    if(t<=0.69) return -116200;
    if(t<=0.70) return -115560;
    if(t<=0.71) return -114920;
    if(t<=0.72) return -114280;
    if(t<=0.73) return -113650;
    if(t<=0.74) return -113020;
    if(t<=0.75) return -112400;
    if(t<=0.76) return -111790;
    if(t<=0.77) return -111200;
    if(t<=0.78) return -110620;
    if(t<=0.79) return -110060;
    //    break;
    /*  case 3:
        return 0.0;
        break;}
        return 0.;*/
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
    return -29;//*e-1*5; // circa [(110-60)*(133.332*10)]/[10*(2.08-1.85)] (il 10 x via dei mm, il 133... x via dei mmHg)
    // (vedi paper di Liu, Dang, etc). Nel loro grafico sono invertite x e y. 5e1 invece di e5 xche' d e'
    // riscalato * il timestep, che e' 5e-4
}


Real hydro(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    return -1.33*1e5;
}

Real u2(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch(i) {
    case 1:
        return 0.0;
        break;
    case 3:
        if ( t <= 0.003 )
            return 1.3332e5;
        //      return 0.01;
        return 0.0;
        break;
    case 2:
        return 0.0;
        //      return 1.3332e4;
        //    else
        //      return 0.0;
        break;
    }
    return 0;
}


Real u2normal(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if(t<=0.003)
        return -1.3332e4;//1.3332e5;
}


// Initial displacement and velocity
Real d0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch(i) {
    case 1:
        return 0.;
        break;
    case 2:
        return 0.;
        break;
    case 3:
        return 0.;
        break;
    default:
        ERROR_MSG("This entry is not allowed: ud_functions.hpp");
        break;
    }
}

Real w0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch(i) {
    case 1:
        return 0.0;
        break;
    case 2:
        return 0.0;
        break;
    case 3:
        return 0.0;
        break;
    default:
        ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
        break;
    }
}

Real fluxFunction(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t<=0.003)
        return -100.;
    else
        return 0.;
}

}


