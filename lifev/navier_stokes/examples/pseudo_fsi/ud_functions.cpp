/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

   <short description here>

   Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file ud_functions.cpp
 */

#include "ud_functions.hpp"
//#include "flowConditions.hpp"
#include "lifev/core/array/VectorEpetra.hpp"


//#define ANEURISM100170
#define AORTA
namespace LifeV
{




Real linearPressAorticValve(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*1000)/1000;
    Real tii=ti+0.01;
    return (aortaPhisPress(tii)-aortaPhisPress(ti))/(0.01)*(t-(ti))+aortaPhisPress(ti) + 115000;
}


Real linearFluxThoracicAorta(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;
    return (aortaFlux3(tii)-aortaFlux3(ti))/(0.01)*(t-(ti))+aortaFlux3(ti);
}


Real linearFluxIn(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;

    if (t<0.001) 
    {
      return -(linearFluxThoracicAorta(t,0,0,0,0)+ linearFluxLeftSubclavian(t,0,0,0,0)
               +linearFluxRightCarotid(t,0,0,0,0)+ linearFluxRightVertebral(t,0,0,0,0) 
                   +linearFluxRightSubclavian(t,0,0,0,0) +linearFluxLeftVertebral(t,0,0,0,0)+linearFluxCommonCarotid(t,0,0,0,0));
    }

    else 
    {  
      return (aortaFluxIn(tii)-aortaFluxIn(ti))/(0.01)*(t-(ti))+aortaFluxIn(ti);
    }
}


Real linearFluxLeftSubclavian(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;
    return (aortaFlux4(tii)-aortaFlux4(ti))/(0.01)*(t-(ti))+aortaFlux4(ti);
}



Real linearFluxRightCarotid(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;
    return (aortaFlux5(tii)-aortaFlux5(ti))/(0.01)*(t-(ti))+aortaFlux5(ti);
}



Real linearFluxRightVertebral(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;
    return (aortaFlux6(tii)-aortaFlux6(ti))/(0.01)*(t-(ti))+aortaFlux6(ti);
}


Real linearFluxRightSubclavian(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;
    return (aortaFlux7(tii)-aortaFlux7(ti))/(0.01)*(t-(ti))+aortaFlux7(ti);
}


Real linearFluxCommonCarotid(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;
    return (aortaFlux8(tii)-aortaFlux8(ti))/(0.01)*(t-(ti))+aortaFlux8(ti);
}

Real linearFluxLeftVertebral(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real ti=floor(t*100)/100;
    Real tii=ti+0.01;
    return (aortaFlux9(tii)-aortaFlux9(ti))/(0.01)*(t-(ti))+aortaFlux9(ti);
}


// Real linearFlux3_(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
// {
//     Real ti=floor(t*1000)/1000;
//     Real tii=ti+0.001;
//     return (aortaFlux3_(tii)-aortaFlux3_(ti))/(0.001)*(t-(ti))+aortaFlux3_(ti);
// }


// Real linearFlux6_(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
// {
//     Real ti=floor(t*1000)/1000;
//     Real tii=ti+0.001;
//     return (aortaFlux6_(tii)-aortaFlux6_(ti))/(0.001)*(t-(ti))+aortaFlux6_(ti);
// }

Real aortaPhisPress(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    /*switch(i) {
      case 1:
      return 0.0;
      break;*/
    //  case 2:
    Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    if(newt<=0.00+0.01) return -1.e1*11017;
    if(newt<=0.01+0.01) return -1.e1*10954;
    if(newt<=0.02+0.01) return -1.e1*10893;
    if(newt<=0.03+0.01) return -1.e1*10832;
    if(newt<=0.04+0.01) return -1.e1*10771;
    if(newt<=0.05+0.01) return -1.e1*10712;
    if(newt<=0.06+0.01) return -1.e1*10653;
    if(newt<=0.07+0.01) return -1.e1*11113;
    if(newt<=0.08+0.01) return -1.e1*11544;
    if(newt<=0.09+0.01) return -1.e1*11869;
    if(newt<=0.10+0.01) return -1.e1*12146;
    if(newt<=0.11+0.01) return -1.e1*12394;
    if(newt<=0.12+0.01) return -1.e1*12635;
    if(newt<=0.13+0.01) return -1.e1*12889;
    if(newt<=0.14+0.01) return -1.e1*13151;
    if(newt<=0.15+0.01) return -1.e1*13398;
    if(newt<=0.16+0.01) return -1.e1*13620;
    if(newt<=0.17+0.01) return -1.e1*13833;
    if(newt<=0.18+0.01) return -1.e1*14035;
    if(newt<=0.19+0.01) return -1.e1*14229;
    if(newt<=0.20+0.01) return -1.e1*14436;
    if(newt<=0.21+0.01) return -1.e1*14613;
    if(newt<=0.22+0.01) return -1.e1*14753;
    if(newt<=0.23+0.01) return -1.e1*14878;
    if(newt<=0.24+0.01) return -1.e1*14974;
    if(newt<=0.25+0.01) return -1.e1*15032;
    if(newt<=0.26+0.01) return -1.e1*15047;
    if(newt<=0.27+0.01) return -1.e1*15025;
    if(newt<=0.28+0.01) return -1.e1*14975;
    if(newt<=0.29+0.01) return -1.e1*14899;
    if(newt<=0.30+0.01) return -1.e1*14822;
    if(newt<=0.31+0.01) return -1.e1*14721;
    if(newt<=0.32+0.01) return -1.e1*14594;
    if(newt<=0.33+0.01) return -1.e1*14496;
    if(newt<=0.34+0.01) return -1.e1*14375;
    if(newt<=0.35+0.01) return -1.e1*14198;
    if(newt<=0.36+0.01) return -1.e1*13990;
    if(newt<=0.37+0.01) return -1.e1*13726;
    if(newt<=0.38+0.01) return -1.e1*13397;
    if(newt<=0.39+0.01) return -1.e1*13167;
    if(newt<=0.40+0.01) return -1.e1*13132;
    if(newt<=0.41+0.01) return -1.e1*13315;
    if(newt<=0.42+0.01) return -1.e1*13271;
    if(newt<=0.43+0.01) return -1.e1*13157;
    if(newt<=0.44+0.01) return -1.e1*13028;
    if(newt<=0.45+0.01) return -1.e1*12975;
    if(newt<=0.46+0.01) return -1.e1*12933;
    if(newt<=0.47+0.01) return -1.e1*12891;
    if(newt<=0.48+0.01) return -1.e1*12836;
    if(newt<=0.49+0.01) return -1.e1*12768;
    if(newt<=0.50+0.01) return -1.e1*12700;
    if(newt<=0.51+0.01) return -1.e1*12641;
    if(newt<=0.52+0.01) return -1.e1*12592;
    if(newt<=0.53+0.01) return -1.e1*12548;
    if(newt<=0.54+0.01) return -1.e1*12504;
    if(newt<=0.55+0.01) return -1.e1*12456;
    if(newt<=0.56+0.01) return -1.e1*12405;
    if(newt<=0.57+0.01) return -1.e1*12353;
    if(newt<=0.58+0.01) return -1.e1*12300;
    if(newt<=0.59+0.01) return -1.e1*12244;
    if(newt<=0.60+0.01) return -1.e1*12184;
    if(newt<=0.61+0.01) return -1.e1*12122;
    if(newt<=0.62+0.01) return -1.e1*12058;
    if(newt<=0.63+0.01) return -1.e1*11995;
    if(newt<=0.64+0.01) return -1.e1*11933;
    if(newt<=0.65+0.01) return -1.e1*11871;
    if(newt<=0.66+0.01) return -1.e1*11810;
    if(newt<=0.67+0.01) return -1.e1*11747;
    if(newt<=0.68+0.01) return -1.e1*11684;
    if(newt<=0.69+0.01) return -1.e1*11620;
    if(newt<=0.70+0.01) return -1.e1*11556;
    if(newt<=0.71+0.01) return -1.e1*11492;
    if(newt<=0.72+0.01) return -1.e1*11428;
    if(newt<=0.73+0.01) return -1.e1*11365;
    if(newt<=0.74+0.01) return -1.e1*11302;
    if(newt<=0.75+0.01) return -1.e1*11240;
    if(newt<=0.76+0.01) return -1.e1*11179;
    if(newt<=0.77+0.01) return -1.e1*11120;
    if(newt<=0.78+0.01) return -1.e1*11062;
    if(newt<=0.79+0.01) return -1.e1*11006;
    //    break;       
    /*  case 3:
        return 0.0;
        break;}
        return 0.;*/
    else return 0;
}


// Real aortaFlux3_(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
// {
//     Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
//     if(newt<=0.00+0.01) return  33.6             ;
//     if(newt<=0.01+0.01) return  33.83            ;
//     if(newt<=0.02+0.01) return  34.12            ;
//     if(newt<=0.03+0.01) return  34.44            ;
//     if(newt<=0.04+0.01) return  34.76            ;
//     if(newt<=0.05+0.01) return  35.05            ;
//     if(newt<=0.06+0.01) return  35.29            ;
//     if(newt<=0.07+0.01) return  35.49            ;
//     if(newt<=0.08+0.01) return  35.74            ;
//     if(newt<=0.09+0.01) return  36.27            ;
//     if(newt<=0.10+0.01) return  37.83            ;
//     if(newt<=0.11+0.01) return  41.86            ;
//     if(newt<=0.12+0.01) return  50.47            ;
//     if(newt<=0.13+0.01) return  65.55999999999999;
//     if(newt<=0.14+0.01) return  87.10999999999999;
//     if(newt<=0.15+0.01) return  111.8            ;
//     if(newt<=0.16+0.01) return  134.6            ;
//     if(newt<=0.17+0.01) return  152              ;
//     if(newt<=0.18+0.01) return  164.3            ;
//     if(newt<=0.19+0.01) return  173.8            ;
//     if(newt<=0.20+0.01) return  182.2            ;
//     if(newt<=0.21+0.01) return  189.8            ;
//     if(newt<=0.22+0.01) return  195.9            ;
//     if(newt<=0.23+0.01) return  200              ;
//     if(newt<=0.24+0.01) return  202              ;
//     if(newt<=0.25+0.01) return  202.1            ;
//     if(newt<=0.26+0.01) return  200.5            ;
//     if(newt<=0.27+0.01) return  197.3            ;
//     if(newt<=0.28+0.01) return  192.7            ;
//     if(newt<=0.29+0.01) return  186.9            ;
//     if(newt<=0.30+0.01) return  179.9            ;
//     if(newt<=0.31+0.01) return  172              ;
//     if(newt<=0.32+0.01) return  163.3            ;
//     if(newt<=0.33+0.01) return  154.1            ;
//     if(newt<=0.34+0.01) return  144.7            ;
//     if(newt<=0.35+0.01) return  135.1            ;
//     if(newt<=0.36+0.01) return  125.5            ;
//     if(newt<=0.37+0.01) return  115.9            ;
//     if(newt<=0.38+0.01) return  106.3            ;
//     if(newt<=0.39+0.01) return  96.60999999999999;
//     if(newt<=0.40+0.01) return  86.44            ;
//     if(newt<=0.41+0.01) return  75.67999999999999;
//     if(newt<=0.42+0.01) return  64.56999999999999;
//     if(newt<=0.43+0.01) return  53.92            ;
//     if(newt<=0.44+0.01) return  44.94            ;
//     if(newt<=0.45+0.01) return  38.75            ;
//     if(newt<=0.46+0.01) return  35.66            ;
//     if(newt<=0.47+0.01) return  34.83000000000001;
//     if(newt<=0.48+0.01) return  34.67            ;
//     if(newt<=0.49+0.01) return  33.83            ;
//     if(newt<=0.50+0.01) return  32.01            ;
//     if(newt<=0.51+0.01) return  29.85            ;
//     if(newt<=0.52+0.01) return  28.19            ;
//     if(newt<=0.53+0.01) return  27.4             ;
//     if(newt<=0.54+0.01) return  27.28            ;
//     if(newt<=0.55+0.01) return  27.42            ;
//     if(newt<=0.56+0.01) return  27.54            ;
//     if(newt<=0.57+0.01) return  27.63            ;
//     if(newt<=0.58+0.01) return  27.82            ;
//     if(newt<=0.59+0.01) return  28.21            ;
//     if(newt<=0.60+0.01) return  28.76            ;
//     if(newt<=0.61+0.01) return  29.39            ;
//     if(newt<=0.62+0.01) return  30               ;
//     if(newt<=0.63+0.01) return  30.54            ;
//     if(newt<=0.64+0.01) return  31.01            ;
//     if(newt<=0.65+0.01) return  31.42            ;
//     if(newt<=0.66+0.01) return  31.78            ;
//     if(newt<=0.67+0.01) return  32.09            ;
//     if(newt<=0.68+0.01) return  32.34            ;
//     if(newt<=0.69+0.01) return  32.54            ;
//     if(newt<=0.70+0.01) return  32.69            ;
//     if(newt<=0.71+0.01) return  32.8             ;
//     if(newt<=0.72+0.01) return  32.87            ;
//     if(newt<=0.73+0.01) return  32.91            ;
//     if(newt<=0.74+0.01) return  32.92            ;
//     if(newt<=0.75+0.01) return  32.93000000000001;
//     if(newt<=0.76+0.01) return  32.93000000000001;
//     if(newt<=0.77+0.01) return  32.96            ;
//     if(newt<=0.78+0.01) return  33.01000000000001;
//     if(newt<=0.79+0.01) return  33.1             ;
// }


Real aortaFlux3(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
        Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    //if(newt<=0.1) return (417.6)*(t);
    if(newt<=0.00+0.01) return 3.350e-05*1e6;
    if(newt<=0.01+0.01) return 3.373e-05*1e6;
    if(newt<=0.02+0.01) return 3.402e-05*1e6;
    if(newt<=0.03+0.01) return 3.434e-05*1e6;
    if(newt<=0.04+0.01) return 3.466e-05*1e6;
    if(newt<=0.05+0.01) return 3.495e-05*1e6;
    if(newt<=0.06+0.01) return 3.519e-05*1e6;
    if(newt<=0.07+0.01) return 3.539e-05*1e6;
    if(newt<=0.08+0.01) return 3.564e-05*1e6;
    if(newt<=0.09+0.01) return 3.617e-05*1e6;
    if(newt<=0.10+0.01) return 3.773e-05*1e6;
    if(newt<=0.11+0.01) return 4.176e-05*1e6;
    if(newt<=0.12+0.01) return 5.037e-05*1e6;
    if(newt<=0.13+0.01) return 6.546e-05*1e6;
    if(newt<=0.14+0.01) return 8.701e-05*1e6;
    if(newt<=0.15+0.01) return 1.117e-04*1e6;
    if(newt<=0.16+0.01) return 1.345e-04*1e6;
    if(newt<=0.17+0.01) return 1.519e-04*1e6;
    if(newt<=0.18+0.01) return 1.642e-04*1e6;
    if(newt<=0.19+0.01) return 1.737e-04*1e6;
    if(newt<=0.20+0.01) return 1.821e-04*1e6;
    if(newt<=0.21+0.01) return 1.897e-04*1e6;
    if(newt<=0.22+0.01) return 1.958e-04*1e6;
    if(newt<=0.23+0.01) return 1.999e-04*1e6;
    if(newt<=0.24+0.01) return 2.019e-04*1e6;
    if(newt<=0.25+0.01) return 2.020e-04*1e6;
    if(newt<=0.26+0.01) return 2.004e-04*1e6;
    if(newt<=0.27+0.01) return 1.972e-04*1e6;
    if(newt<=0.28+0.01) return 1.926e-04*1e6;
    if(newt<=0.29+0.01) return 1.868e-04*1e6;
    if(newt<=0.30+0.01) return 1.798e-04*1e6;
    if(newt<=0.31+0.01) return 1.719e-04*1e6;
    if(newt<=0.32+0.01) return 1.632e-04*1e6;
    if(newt<=0.33+0.01) return 1.540e-04*1e6;
    if(newt<=0.34+0.01) return 1.446e-04*1e6;
    if(newt<=0.35+0.01) return 1.350E-04*1e6;
    if(newt<=0.36+0.01) return 1.254e-04*1e6;
    if(newt<=0.37+0.01) return 1.158e-04*1e6;
    if(newt<=0.38+0.01) return 1.062e-04*1e6;
    if(newt<=0.39+0.01) return 9.651e-05*1e6;
    if(newt<=0.40+0.01) return 8.634e-05*1e6;
    if(newt<=0.41+0.01) return 7.558e-05*1e6;
    if(newt<=0.42+0.01) return 6.447e-05*1e6;
    if(newt<=0.43+0.01) return 5.382e-05*1e6;
    if(newt<=0.44+0.01) return 4.484e-05*1e6;
    if(newt<=0.45+0.01) return 3.865e-05*1e6;
    if(newt<=0.46+0.01) return 3.556e-05*1e6;
    if(newt<=0.47+0.01) return 3.473e-05*1e6;
    if(newt<=0.48+0.01) return 3.457e-05*1e6;
    if(newt<=0.49+0.01) return 3.373e-05*1e6;
    if(newt<=0.50+0.01) return 3.191e-05*1e6;
    if(newt<=0.51+0.01) return 2.975e-05*1e6;
    if(newt<=0.52+0.01) return 2.809e-05*1e6;
    if(newt<=0.53+0.01) return 2.730e-05*1e6;
    if(newt<=0.54+0.01) return 2.718e-05*1e6;
    if(newt<=0.55+0.01) return 2.732e-05*1e6;
    if(newt<=0.56+0.01) return 2.744e-05*1e6;
    if(newt<=0.57+0.01) return 2.753e-05*1e6;
    if(newt<=0.58+0.01) return 2.772e-05*1e6;
    if(newt<=0.59+0.01) return 2.811e-05*1e6;
    if(newt<=0.60+0.01) return 2.866e-05*1e6;
    if(newt<=0.61+0.01) return 2.929e-05*1e6;
    if(newt<=0.62+0.01) return 2.990e-05*1e6;
    if(newt<=0.63+0.01) return 3.044e-05*1e6;
    if(newt<=0.64+0.01) return 3.091e-05*1e6;
    if(newt<=0.65+0.01) return 3.132e-05*1e6;
    if(newt<=0.66+0.01) return 3.168e-05*1e6;
    if(newt<=0.67+0.01) return 3.199e-05*1e6;
    if(newt<=0.68+0.01) return 3.224e-05*1e6;
    if(newt<=0.69+0.01) return 3.244e-05*1e6;
    if(newt<=0.70+0.01) return 3.259e-05*1e6;
    if(newt<=0.71+0.01) return 3.270e-05*1e6;
    if(newt<=0.72+0.01) return 3.277e-05*1e6;
    if(newt<=0.73+0.01) return 3.281e-05*1e6;
    if(newt<=0.74+0.01) return 3.282e-05*1e6;
    if(newt<=0.75+0.01) return 3.283e-05*1e6;
    if(newt<=0.76+0.01) return 3.283e-05*1e6;
    if(newt<=0.77+0.01) return 3.286e-05*1e6;
    if(newt<=0.78+0.01) return 3.291e-05*1e6;
    if(newt<=0.79+0.01) return 3.300e-05*1e6;
    else return 0;
}//thoracic aorta,

Real aortaFlux5(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
        Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    //  if(newt<=0.1) return (73.87)*(t);
    if(newt<=0.00+0.01) return  3.033e-06*1e6;
    if(newt<=0.01+0.01) return  3.041e-06*1e6;
    if(newt<=0.02+0.01) return  3.062e-06*1e6;
    if(newt<=0.03+0.01) return  3.094e-06*1e6;
    if(newt<=0.04+0.01) return  3.127e-06*1e6;
    if(newt<=0.05+0.01) return  3.150e-06*1e6;
    if(newt<=0.06+0.01) return  3.152e-06*1e6;
    if(newt<=0.07+0.01) return  3.141e-06*1e6;
    if(newt<=0.08+0.01) return  3.196e-06*1e6;
    if(newt<=0.09+0.01) return  3.574e-06*1e6;
    if(newt<=0.10+0.01) return  4.778e-06*1e6;
    if(newt<=0.11+0.01) return  7.387e-06*1e6;
    if(newt<=0.12+0.01) return  1.150e-05*1e6;
    if(newt<=0.13+0.01) return  1.609e-05*1e6;
    if(newt<=0.14+0.01) return  1.933e-05*1e6;
    if(newt<=0.15+0.01) return  2.007e-05*1e6;
    if(newt<=0.16+0.01) return  1.885e-05*1e6;
    if(newt<=0.17+0.01) return  1.706e-05*1e6;
    if(newt<=0.18+0.01) return  1.569e-05*1e6;
    if(newt<=0.19+0.01) return  1.481e-05*1e6;
    if(newt<=0.20+0.01) return  1.401e-05*1e6;
    if(newt<=0.21+0.01) return  1.294e-05*1e6;
    if(newt<=0.22+0.01) return  1.160e-05*1e6;
    if(newt<=0.23+0.01) return  1.018e-05*1e6;
    if(newt<=0.24+0.01) return  8.832e-06*1e6;
    if(newt<=0.25+0.01) return  7.609e-06*1e6;
    if(newt<=0.26+0.01) return  6.578e-06*1e6;
    if(newt<=0.27+0.01) return  5.843e-06*1e6;
    if(newt<=0.28+0.01) return  5.472e-06*1e6;
    if(newt<=0.29+0.01) return  5.412e-06*1e6;
    if(newt<=0.30+0.01) return  5.491e-06*1e6;
    if(newt<=0.31+0.01) return  5.527e-06*1e6;
    if(newt<=0.32+0.01) return  5.420e-06*1e6;
    if(newt<=0.33+0.01) return  5.169e-06*1e6;
    if(newt<=0.34+0.01) return  4.829e-06*1e6;
    if(newt<=0.35+0.01) return  4.465e-06*1e6;
    if(newt<=0.36+0.01) return  4.111e-06*1e6;
    if(newt<=0.37+0.01) return  3.750e-06*1e6;
    if(newt<=0.38+0.01) return  3.304e-06*1e6;
    if(newt<=0.39+0.01) return  2.668e-06*1e6;
    if(newt<=0.40+0.01) return  1.800e-06*1e6;
    if(newt<=0.41+0.01) return  8.269e-07*1e6;
    if(newt<=0.42+0.01) return  9.760e-08*1e6;
    if(newt<=0.43+0.01) return  7.311e-08*1e6;
    if(newt<=0.44+0.01) return  1.041e-06*1e6;
    if(newt<=0.45+0.01) return  2.783e-06*1e6;
    if(newt<=0.46+0.01) return  4.537e-06*1e6;
    if(newt<=0.47+0.01) return  5.488e-06*1e6;
    if(newt<=0.48+0.01) return  5.431e-06*1e6;
    if(newt<=0.49+0.01) return  4.863e-06*1e6;
    if(newt<=0.50+0.01) return  4.452e-06*1e6;
    if(newt<=0.51+0.01) return  4.499e-06*1e6;
    if(newt<=0.52+0.01) return  4.824e-06*1e6;
    if(newt<=0.53+0.01) return  5.059e-06*1e6;
    if(newt<=0.54+0.01) return  4.989e-06*1e6;
    if(newt<=0.55+0.01) return  4.671e-06*1e6;
    if(newt<=0.56+0.01) return  4.292e-06*1e6;
    if(newt<=0.57+0.01) return  3.981e-06*1e6;
    if(newt<=0.58+0.01) return  3.749e-06*1e6;
    if(newt<=0.59+0.01) return  3.553e-06*1e6;
    if(newt<=0.60+0.01) return  3.377e-06*1e6;
    if(newt<=0.61+0.01) return  3.255e-06*1e6;
    if(newt<=0.62+0.01) return  3.224e-06*1e6;
    if(newt<=0.63+0.01) return  3.281e-06*1e6;
    if(newt<=0.64+0.01) return  3.377e-06*1e6;
    if(newt<=0.65+0.01) return  3.452e-06*1e6;
    if(newt<=0.66+0.01) return  3.472e-06*1e6;
    if(newt<=0.67+0.01) return  3.441e-06*1e6;
    if(newt<=0.68+0.01) return  3.389e-06*1e6;
    if(newt<=0.69+0.01) return  3.343e-06*1e6;
    if(newt<=0.70+0.01) return  3.312e-06*1e6;
    if(newt<=0.71+0.01) return  3.289e-06*1e6;
    if(newt<=0.72+0.01) return  3.262e-06*1e6;
    if(newt<=0.73+0.01) return  3.223e-06*1e6;
    if(newt<=0.74+0.01) return  3.177e-06*1e6;
    if(newt<=0.75+0.01) return  3.132e-06*1e6;
    if(newt<=0.76+0.01) return  3.094e-06*1e6;
    if(newt<=0.77+0.01) return  3.065e-06*1e6;
    if(newt<=0.78+0.01) return  3.040e-06*1e6;
    if(newt<=0.79+0.01) return  3.016e-06*1e6;
    else return 0;
}//first branch_1,


// Real aortaFlux6_(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
// {
//         Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
//     if(newt<=0.00+0.01) return    0.6737071250000001;
//     if(newt<=0.01+0.01) return    0.6799071250000001;
//     if(newt<=0.02+0.01) return    0.6897071250000001;
//     if(newt<=0.03+0.01) return    0.6997071250000001;
//     if(newt<=0.04+0.01) return    0.706407125       ;
//     if(newt<=0.05+0.01) return    0.7084071250000001;
//     if(newt<=0.06+0.01) return    0.7056071250000001;
//     if(newt<=0.07+0.01) return    0.7049071250000001;
//     if(newt<=0.08+0.01) return    0.7382071250000001;
//     if(newt<=0.09+0.01) return    0.8778071249999999;
//     if(newt<=0.10+0.01) return    1.209007125       ;
//     if(newt<=0.11+0.01) return    1.746007125       ;
//     if(newt<=0.12+0.01) return    2.368007125       ;
//     if(newt<=0.13+0.01) return    2.897007125       ;
//     if(newt<=0.14+0.01) return    3.258007125       ;
//     if(newt<=0.15+0.01) return    3.518007125       ;
//     if(newt<=0.16+0.01) return    3.752007125       ;
//     if(newt<=0.17+0.01) return    3.936007124999999 ;
//     if(newt<=0.18+0.01) return    3.957007125       ;
//     if(newt<=0.19+0.01) return    3.716007125       ;
//     if(newt<=0.20+0.01) return    3.212007125       ;
//     if(newt<=0.21+0.01) return    2.551007125       ;
//     if(newt<=0.22+0.01) return    1.898007125       ;
//     if(newt<=0.23+0.01) return    1.395007125       ;
//     if(newt<=0.24+0.01) return    1.107007125       ;
//     if(newt<=0.25+0.01) return    1.009007125       ;
//     if(newt<=0.26+0.01) return    1.028007125       ;
//     if(newt<=0.27+0.01) return    1.082007125       ;
//     if(newt<=0.28+0.01) return    1.114007125       ;
//     if(newt<=0.29+0.01) return    1.111007125       ;
//     if(newt<=0.30+0.01) return    1.089007125       ;
//     if(newt<=0.31+0.01) return    1.071007125       ;
//     if(newt<=0.32+0.01) return    1.067007125       ;
//     if(newt<=0.33+0.01) return    1.068007125       ;
//     if(newt<=0.34+0.01) return    1.056007125       ;
//     if(newt<=0.35+0.01) return    1.021007125       ;
//     if(newt<=0.36+0.01) return    0.955007125       ;
//     if(newt<=0.37+0.01) return    0.856707125       ;
//     if(newt<=0.38+0.01) return    0.7230071250000001;
//     if(newt<=0.39+0.01) return    0.5555071250000001;
//     if(newt<=0.40+0.01) return    0.3694071249999999;
//     if(newt<=0.41+0.01) return    0.203607125       ;
//     if(newt<=0.42+0.01) return    0.117107125       ;
//     if(newt<=0.43+0.01) return    0.154707125       ;
//     if(newt<=0.44+0.01) return    0.3019071249999999;
//     if(newt<=0.45+0.01) return    0.4833071250000001;
//     if(newt<=0.46+0.01) return    0.6279071250000001;
//     if(newt<=0.47+0.01) return    0.7323071250000001;
//     if(newt<=0.48+0.01) return    0.8435071250000001;
//     if(newt<=0.49+0.01) return    0.989007125       ;
//     if(newt<=0.50+0.01) return    1.137007125       ;
//     if(newt<=0.51+0.01) return    1.223007125       ;
//     if(newt<=0.52+0.01) return    1.206007125       ;
//     if(newt<=0.53+0.01) return    1.096007125       ;
//     if(newt<=0.54+0.01) return    0.9400071249999998;
//     if(newt<=0.55+0.01) return    0.7924071250000001;
//     if(newt<=0.56+0.01) return    0.6857071250000001;
//     if(newt<=0.57+0.01) return    0.6278071250000001;
//     if(newt<=0.58+0.01) return    0.6083071250000001;
//     if(newt<=0.59+0.01) return    0.610607125       ;
//     if(newt<=0.60+0.01) return    0.620007125       ;
//     if(newt<=0.61+0.01) return    0.6279071250000001;
//     if(newt<=0.62+0.01) return    0.6320071250000001;
//     if(newt<=0.63+0.01) return    0.6343071250000001;
//     if(newt<=0.64+0.01) return    0.6385071250000001;
//     if(newt<=0.65+0.01) return    0.6476071250000001;
//     if(newt<=0.66+0.01) return    0.6620071250000001;
//     if(newt<=0.67+0.01) return    0.679107125       ;
//     if(newt<=0.68+0.01) return    0.6948071250000001;
//     if(newt<=0.69+0.01) return    0.7049071250000001;
//     if(newt<=0.70+0.01) return    0.7071071250000001;
//     if(newt<=0.71+0.01) return    0.7021071250000001;
//     if(newt<=0.72+0.01) return    0.6926071250000001;
//     if(newt<=0.73+0.01) return    0.6825071250000001;
//     if(newt<=0.74+0.01) return    0.675007125       ;
//     if(newt<=0.75+0.01) return    0.671207125       ;
//     if(newt<=0.76+0.01) return    0.6702071250000001;
//     if(newt<=0.77+0.01) return    0.6702071250000001;
//     if(newt<=0.78+0.01) return    0.669507125       ;
//     if(newt<=0.79+0.01) return    0.6674071250000001;
// }

// Real aortaFlux6(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
// {
//         Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
//     //  if(newt<=0.1) return (18.54)*(t);
//     if(newt<=0.00+0.01) return   0.7737071250000001;
//     if(newt<=0.01+0.01) return   0.779907125       ;
//     if(newt<=0.02+0.01) return   0.7897071250000001;
//     if(newt<=0.03+0.01) return   0.7997071250000001;
//     if(newt<=0.04+0.01) return   0.806407125       ;
//     if(newt<=0.05+0.01) return   0.8084071250000001;
//     if(newt<=0.06+0.01) return   0.8056071250000001;
//     if(newt<=0.07+0.01) return   0.8049071250000001;
//     if(newt<=0.08+0.01) return   0.8382071250000001;
//     if(newt<=0.09+0.01) return   0.9778071249999999;
//     if(newt<=0.10+0.01) return   1.309007125       ;
//     if(newt<=0.11+0.01) return   1.846007125       ;
//     if(newt<=0.12+0.01) return   2.468007125       ;
//     if(newt<=0.13+0.01) return   2.997007125       ;
//     if(newt<=0.14+0.01) return   3.358007125       ;
//     if(newt<=0.15+0.01) return   3.618007125       ;
//     if(newt<=0.16+0.01) return   3.852007125       ;
//     if(newt<=0.17+0.01) return   4.036007124999999 ;
//     if(newt<=0.18+0.01) return   4.057007125       ;
//     if(newt<=0.19+0.01) return   3.816007125       ;
//     if(newt<=0.20+0.01) return   3.312007125       ;
//     if(newt<=0.21+0.01) return   2.651007125       ;
//     if(newt<=0.22+0.01) return   1.998007125       ;
//     if(newt<=0.23+0.01) return   1.495007125       ;
//     if(newt<=0.24+0.01) return   1.207007125       ;
//     if(newt<=0.25+0.01) return   1.109007125       ;
//     if(newt<=0.26+0.01) return   1.128007125       ;
//     if(newt<=0.27+0.01) return   1.182007125       ;
//     if(newt<=0.28+0.01) return   1.214007125       ;
//     if(newt<=0.29+0.01) return   1.211007125       ;
//     if(newt<=0.30+0.01) return   1.189007125       ;
//     if(newt<=0.31+0.01) return   1.171007125       ;
//     if(newt<=0.32+0.01) return   1.167007125       ;
//     if(newt<=0.33+0.01) return   1.168007125       ;
//     if(newt<=0.34+0.01) return   1.156007125       ;
//     if(newt<=0.35+0.01) return   1.121007125       ;
//     if(newt<=0.36+0.01) return   1.055007125       ;
//     if(newt<=0.37+0.01) return   0.956707125       ;
//     if(newt<=0.38+0.01) return   0.8230071250000001;
//     if(newt<=0.39+0.01) return   0.6555071250000001;
//     if(newt<=0.40+0.01) return   0.469407125       ;
//     if(newt<=0.41+0.01) return   0.303607125       ;
//     if(newt<=0.42+0.01) return   0.217107125       ;
//     if(newt<=0.43+0.01) return   0.254707125       ;
//     if(newt<=0.44+0.01) return   0.401907125       ;
//     if(newt<=0.45+0.01) return   0.583307125       ;
//     if(newt<=0.46+0.01) return   0.7279071250000001;
//     if(newt<=0.47+0.01) return   0.832307125       ;
//     if(newt<=0.48+0.01) return   0.9435071250000001;
//     if(newt<=0.49+0.01) return   1.089007125       ;
//     if(newt<=0.50+0.01) return   1.237007125       ;
//     if(newt<=0.51+0.01) return   1.323007125       ;
//     if(newt<=0.52+0.01) return   1.306007125       ;
//     if(newt<=0.53+0.01) return   1.196007125       ;
//     if(newt<=0.54+0.01) return   1.040007125       ;
//     if(newt<=0.55+0.01) return   0.8924071250000001;
//     if(newt<=0.56+0.01) return   0.7857071250000001;
//     if(newt<=0.57+0.01) return   0.7278071250000001;
//     if(newt<=0.58+0.01) return   0.708307125       ;
//     if(newt<=0.59+0.01) return   0.710607125       ;
//     if(newt<=0.60+0.01) return   0.720007125       ;
//     if(newt<=0.61+0.01) return   0.7279071250000001;
//     if(newt<=0.62+0.01) return   0.7320071250000001;
//     if(newt<=0.63+0.01) return   0.7343071250000001;
//     if(newt<=0.64+0.01) return   0.738507125       ;
//     if(newt<=0.65+0.01) return   0.747607125       ;
//     if(newt<=0.66+0.01) return   0.7620071250000001;
//     if(newt<=0.67+0.01) return   0.779107125       ;
//     if(newt<=0.68+0.01) return   0.7948071250000001;
//     if(newt<=0.69+0.01) return   0.8049071250000001;
//     if(newt<=0.70+0.01) return   0.807107125       ;
//     if(newt<=0.71+0.01) return   0.802107125       ;
//     if(newt<=0.72+0.01) return   0.7926071250000001;
//     if(newt<=0.73+0.01) return   0.7825071250000001;
//     if(newt<=0.74+0.01) return   0.775007125       ;
//     if(newt<=0.75+0.01) return   0.771207125       ;
//     if(newt<=0.76+0.01) return   0.7702071250000001;
//     if(newt<=0.77+0.01) return   0.7702071250000001;
//     if(newt<=0.78+0.01) return   0.769507125       ;
//     if(newt<=0.79+0.01) return   0.7674071250000001;
// }//branch 1_2 smallest

Real aortaFlux6(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    //  if(newt<=0.1) return (18.54)*(t);
if(newt<=0.00+0.01) return  7.817e-07*1e6;
if(newt<=0.01+0.01) return  7.879e-07*1e6;
if(newt<=0.02+0.01) return  7.977e-07*1e6;
if(newt<=0.03+0.01) return  8.077e-07*1e6;
if(newt<=0.04+0.01) return  8.144e-07*1e6;
if(newt<=0.05+0.01) return  8.164e-07*1e6;
if(newt<=0.06+0.01) return  8.136e-07*1e6;
if(newt<=0.07+0.01) return  8.129e-07*1e6;
if(newt<=0.08+0.01) return  8.462e-07*1e6;
if(newt<=0.09+0.01) return  9.858e-07*1e6;
if(newt<=0.10+0.01) return  1.317e-06*1e6;
if(newt<=0.11+0.01) return  1.854e-06*1e6;
if(newt<=0.12+0.01) return  2.476e-06*1e6;
if(newt<=0.13+0.01) return  3.005e-06*1e6;
if(newt<=0.14+0.01) return  3.366e-06*1e6;
if(newt<=0.15+0.01) return  3.626e-06*1e6;
if(newt<=0.16+0.01) return  3.860e-06*1e6;
if(newt<=0.17+0.01) return  4.044e-06*1e6;
if(newt<=0.18+0.01) return  4.065e-06*1e6;
if(newt<=0.19+0.01) return  3.824e-06*1e6;
if(newt<=0.20+0.01) return  3.320e-06*1e6;
if(newt<=0.21+0.01) return  2.659e-06*1e6;
if(newt<=0.22+0.01) return  2.006e-06*1e6;
if(newt<=0.23+0.01) return  1.503e-06*1e6;
if(newt<=0.24+0.01) return  1.215e-06*1e6;
if(newt<=0.25+0.01) return  1.117e-06*1e6;
if(newt<=0.26+0.01) return  1.136e-06*1e6;
if(newt<=0.27+0.01) return  1.190e-06*1e6;
if(newt<=0.28+0.01) return  1.222e-06*1e6;
if(newt<=0.29+0.01) return  1.219e-06*1e6;
if(newt<=0.30+0.01) return  1.197e-06*1e6;
if(newt<=0.31+0.01) return  1.179e-06*1e6;
if(newt<=0.32+0.01) return  1.175e-06*1e6;
if(newt<=0.33+0.01) return  1.176e-06*1e6;
if(newt<=0.34+0.01) return  1.164e-06*1e6;
if(newt<=0.35+0.01) return  1.129e-06*1e6;
if(newt<=0.36+0.01) return  1.063e-06*1e6;
if(newt<=0.37+0.01) return  9.647e-07*1e6;
if(newt<=0.38+0.01) return  8.310e-07*1e6;
if(newt<=0.39+0.01) return  6.635e-07*1e6;
if(newt<=0.40+0.01) return  4.774e-07*1e6;
if(newt<=0.41+0.01) return  3.116e-07*1e6;
if(newt<=0.42+0.01) return  2.251e-07*1e6;
if(newt<=0.43+0.01) return  2.627e-07*1e6;
if(newt<=0.44+0.01) return  4.099e-07*1e6;
if(newt<=0.45+0.01) return  5.913e-07*1e6;
if(newt<=0.46+0.01) return  7.359e-07*1e6;
if(newt<=0.47+0.01) return  8.403e-07*1e6;
if(newt<=0.48+0.01) return  9.515e-07*1e6;
if(newt<=0.49+0.01) return  1.097e-06*1e6;
if(newt<=0.50+0.01) return  1.245e-06*1e6;
if(newt<=0.51+0.01) return  1.331e-06*1e6;
if(newt<=0.52+0.01) return  1.314e-06*1e6;
if(newt<=0.53+0.01) return  1.204e-06*1e6;
if(newt<=0.54+0.01) return  1.048e-06*1e6;
if(newt<=0.55+0.01) return  9.004e-07*1e6;
if(newt<=0.56+0.01) return  7.937e-07*1e6;
if(newt<=0.57+0.01) return  7.358e-07*1e6;
if(newt<=0.58+0.01) return  7.163e-07*1e6;
if(newt<=0.59+0.01) return  7.186e-07*1e6;
if(newt<=0.60+0.01) return  7.280e-07*1e6;
if(newt<=0.61+0.01) return  7.359e-07*1e6;
if(newt<=0.62+0.01) return  7.400e-07*1e6;
if(newt<=0.63+0.01) return  7.423e-07*1e6;
if(newt<=0.64+0.01) return  7.465e-07*1e6;
if(newt<=0.65+0.01) return  7.556e-07*1e6;
if(newt<=0.66+0.01) return  7.700e-07*1e6;
if(newt<=0.67+0.01) return  7.871e-07*1e6;
if(newt<=0.68+0.01) return  8.028e-07*1e6;
if(newt<=0.69+0.01) return  8.129e-07*1e6;
if(newt<=0.70+0.01) return  8.151e-07*1e6;
if(newt<=0.71+0.01) return  8.101e-07*1e6;
if(newt<=0.72+0.01) return  8.006e-07*1e6;
if(newt<=0.73+0.01) return  7.905e-07*1e6;
if(newt<=0.74+0.01) return  7.830e-07*1e6;
if(newt<=0.75+0.01) return  7.792e-07*1e6;
if(newt<=0.76+0.01) return  7.782e-07*1e6;
if(newt<=0.77+0.01) return  7.782e-07*1e6;
if(newt<=0.78+0.01) return  7.775e-07*1e6;
if(newt<=0.79+0.01) return  7.754e-07*1e6;
 else return 0;
}//branch 1_2 smallest

Real aortaFlux7(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
        Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    //  if(newt<=0.1) return (30.32)*(t);
    if(newt<=0.00+0.01) return  1.930e-06 *1e6;
    if(newt<=0.01+0.01) return  1.710e-06 *1e6;
    if(newt<=0.02+0.01) return  1.495e-06 *1e6;
    if(newt<=0.03+0.01) return  1.289e-06 *1e6;
    if(newt<=0.04+0.01) return  1.096e-06 *1e6;
    if(newt<=0.05+0.01) return  9.184e-07 *1e6;
    if(newt<=0.06+0.01) return  7.568e-07 *1e6;
    if(newt<=0.07+0.01) return  6.240e-07 *1e6;
    if(newt<=0.08+0.01) return  5.854e-07 *1e6;
    if(newt<=0.09+0.01) return  8.113e-07 *1e6;
    if(newt<=0.10+0.01) return  1.559e-06 *1e6;
    if(newt<=0.11+0.01) return  3.032e-06 *1e6;
    if(newt<=0.12+0.01) return  5.221e-06 *1e6;
    if(newt<=0.13+0.01) return  7.903e-06 *1e6;
    if(newt<=0.14+0.01) return  1.081e-05 *1e6;
    if(newt<=0.15+0.01) return  1.376e-05 *1e6;
    if(newt<=0.16+0.01) return  1.661e-05 *1e6;
    if(newt<=0.17+0.01) return  1.919e-05 *1e6;
    if(newt<=0.18+0.01) return  2.134e-05 *1e6;
    if(newt<=0.19+0.01) return  2.292e-05 *1e6;
    if(newt<=0.20+0.01) return  2.389e-05 *1e6;
    if(newt<=0.21+0.01) return  2.433e-05 *1e6;
    if(newt<=0.22+0.01) return  2.435e-05 *1e6;
    if(newt<=0.23+0.01) return  2.406e-05 *1e6;
    if(newt<=0.24+0.01) return  2.353e-05 *1e6;
    if(newt<=0.25+0.01) return  2.280e-05 *1e6;
    if(newt<=0.26+0.01) return  2.188e-05 *1e6;
    if(newt<=0.27+0.01) return  2.078e-05 *1e6;
    if(newt<=0.28+0.01) return  1.947e-05 *1e6;
    if(newt<=0.29+0.01) return  1.797e-05 *1e6;
    if(newt<=0.30+0.01) return  1.628e-05 *1e6;
    if(newt<=0.31+0.01) return  1.445e-05 *1e6;
    if(newt<=0.32+0.01) return  1.254e-05 *1e6;
    if(newt<=0.33+0.01) return  1.060e-05 *1e6;
    if(newt<=0.34+0.01) return  8.684e-06 *1e6;
    if(newt<=0.35+0.01) return  6.838e-06 *1e6;
    if(newt<=0.36+0.01) return  5.084e-06 *1e6;
    if(newt<=0.37+0.01) return  3.412e-06 *1e6;
    if(newt<=0.38+0.01) return  1.784e-06 *1e6;
    if(newt<=0.39+0.01) return  1.534e-07 *1e6;
    if(newt<=0.40+0.01) return  -1.494e-06*1e6;
    if(newt<=0.41+0.01) return  -3.093e-06*1e6;
    if(newt<=0.42+0.01) return  -4.495e-06*1e6;
    if(newt<=0.43+0.01) return  -5.521e-06*1e6;
    if(newt<=0.44+0.01) return  -6.086e-06*1e6;
    if(newt<=0.45+0.01) return  -6.252e-06*1e6;
    if(newt<=0.46+0.01) return  -6.160e-06*1e6;
    if(newt<=0.47+0.01) return  -5.908e-06*1e6;
    if(newt<=0.48+0.01) return  -5.508e-06*1e6;
    if(newt<=0.49+0.01) return  -4.944e-06*1e6;
    if(newt<=0.50+0.01) return  -4.234e-06*1e6;
    if(newt<=0.51+0.01) return  -3.441e-06*1e6;
    if(newt<=0.52+0.01) return  -2.639e-06*1e6;
    if(newt<=0.53+0.01) return  -1.873e-06*1e6;
    if(newt<=0.54+0.01) return  -1.160e-06*1e6;
    if(newt<=0.55+0.01) return  -5.019e-07*1e6;
    if(newt<=0.56+0.01) return  1.024e-07 *1e6;
    if(newt<=0.57+0.01) return  6.543e-07 *1e6;
    if(newt<=0.58+0.01) return  1.156e-06 *1e6;
    if(newt<=0.59+0.01) return  1.614e-06 *1e6;
    if(newt<=0.60+0.01) return  2.034e-06 *1e6;
    if(newt<=0.61+0.01) return  2.419e-06 *1e6;
    if(newt<=0.62+0.01) return  2.765e-06 *1e6;
    if(newt<=0.63+0.01) return  3.065e-06 *1e6;
    if(newt<=0.64+0.01) return  3.307e-06 *1e6;
    if(newt<=0.65+0.01) return  3.486e-06 *1e6;
    if(newt<=0.66+0.01) return  3.600e-06 *1e6;
    if(newt<=0.67+0.01) return  3.654e-06 *1e6;
    if(newt<=0.68+0.01) return  3.657e-06 *1e6;
    if(newt<=0.69+0.01) return  3.621e-06 *1e6;
    if(newt<=0.70+0.01) return  3.554e-06 *1e6;
    if(newt<=0.71+0.01) return  3.465e-06 *1e6;
    if(newt<=0.72+0.01) return  3.357e-06 *1e6;
    if(newt<=0.73+0.01) return  3.233e-06 *1e6;
    if(newt<=0.74+0.01) return  3.094e-06 *1e6;
    if(newt<=0.75+0.01) return  2.941e-06 *1e6;
    if(newt<=0.76+0.01) return  2.774e-06 *1e6;
    if(newt<=0.77+0.01) return  2.591e-06 *1e6;
    if(newt<=0.78+0.01) return  2.395e-06 *1e6;
    if(newt<=0.79+0.01) return  2.185e-06 *1e6;
    else return 0;
}//R. Brachia, branch 1_3

Real aortaFlux8(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
        Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    //  if(newt<=0.1) return (69.56)*(t);
    if(newt<=0.00+0.01) return  2.445e-06*1e6;
    if(newt<=0.01+0.01) return  2.470e-06*1e6;
    if(newt<=0.02+0.01) return  2.519e-06*1e6;
    if(newt<=0.03+0.01) return  2.578e-06*1e6;
    if(newt<=0.04+0.01) return  2.631e-06*1e6;
    if(newt<=0.05+0.01) return  2.665e-06*1e6;
    if(newt<=0.06+0.01) return  2.674e-06*1e6;
    if(newt<=0.07+0.01) return  2.674e-06*1e6;
    if(newt<=0.08+0.01) return  2.777e-06*1e6;
    if(newt<=0.09+0.01) return  3.260e-06*1e6;
    if(newt<=0.10+0.01) return  4.539e-06*1e6;
    if(newt<=0.11+0.01) return  6.956e-06*1e6;
    if(newt<=0.12+0.01) return  1.045e-05*1e6;
    if(newt<=0.13+0.01) return  1.438e-05*1e6;
    if(newt<=0.14+0.01) return  1.769e-05*1e6;
    if(newt<=0.15+0.01) return  1.958e-05*1e6;
    if(newt<=0.16+0.01) return  1.991e-05*1e6;
    if(newt<=0.17+0.01) return  1.910e-05*1e6;
    if(newt<=0.18+0.01) return  1.775e-05*1e6;
    if(newt<=0.19+0.01) return  1.628e-05*1e6;
    if(newt<=0.20+0.01) return  1.479e-05*1e6;
    if(newt<=0.21+0.01) return  1.326e-05*1e6;
    if(newt<=0.22+0.01) return  1.165e-05*1e6;
    if(newt<=0.23+0.01) return  1.001e-05*1e6;
    if(newt<=0.24+0.01) return  8.448e-06*1e6;
    if(newt<=0.25+0.01) return  7.065e-06*1e6;
    if(newt<=0.26+0.01) return  5.930e-06*1e6;
    if(newt<=0.27+0.01) return  5.089e-06*1e6;
    if(newt<=0.28+0.01) return  4.566e-06*1e6;
    if(newt<=0.29+0.01) return  4.350e-06*1e6;
    if(newt<=0.30+0.01) return  4.382e-06*1e6;
    if(newt<=0.31+0.01) return  4.552e-06*1e6;
    if(newt<=0.32+0.01) return  4.729e-06*1e6;
    if(newt<=0.33+0.01) return  4.798e-06*1e6;
    if(newt<=0.34+0.01) return  4.698e-06*1e6;
    if(newt<=0.35+0.01) return  4.428e-06*1e6;
    if(newt<=0.36+0.01) return  4.026e-06*1e6;
    if(newt<=0.37+0.01) return  3.513e-06*1e6;
    if(newt<=0.38+0.01) return  2.875e-06*1e6;
    if(newt<=0.39+0.01) return  2.076e-06*1e6;
    if(newt<=0.40+0.01) return  1.131e-06*1e6;
    if(newt<=0.41+0.01) return  1.861e-07*1e6;
    if(newt<=0.42+0.01) return  -4.577e-07*1e6;
    if(newt<=0.43+0.01) return  -4.657e-07*1e6;
    if(newt<=0.44+0.01) return  3.124e-07*1e6;
    if(newt<=0.45+0.01) return  1.684e-06*1e6;
    if(newt<=0.46+0.01) return  3.174e-06*1e6;
    if(newt<=0.47+0.01) return  4.306e-06*1e6;
    if(newt<=0.48+0.01) return  4.873e-06*1e6;
    if(newt<=0.49+0.01) return  4.980e-06*1e6;
    if(newt<=0.50+0.01) return  4.876e-06*1e6;
    if(newt<=0.51+0.01) return  4.757e-06*1e6;
    if(newt<=0.52+0.01) return  4.681e-06*1e6;
    if(newt<=0.53+0.01) return  4.603e-06*1e6;
    if(newt<=0.54+0.01) return  4.458e-06*1e6;
    if(newt<=0.55+0.01) return  4.224e-06*1e6;
    if(newt<=0.56+0.01) return  3.928e-06*1e6;
    if(newt<=0.57+0.01) return  3.618e-06*1e6;
    if(newt<=0.58+0.01) return  3.335e-06*1e6;
    if(newt<=0.59+0.01) return  3.098e-06*1e6;
    if(newt<=0.60+0.01) return  2.913e-06*1e6;
    if(newt<=0.61+0.01) return  2.778e-06*1e6;
    if(newt<=0.62+0.01) return  2.697e-06*1e6;
    if(newt<=0.63+0.01) return  2.670e-06*1e6;
    if(newt<=0.64+0.01) return  2.692e-06*1e6;
    if(newt<=0.65+0.01) return  2.746e-06*1e6;
    if(newt<=0.66+0.01) return  2.811e-06*1e6;
    if(newt<=0.67+0.01) return  2.868e-06*1e6;
    if(newt<=0.68+0.01) return  2.904e-06*1e6;
    if(newt<=0.69+0.01) return  2.915e-06*1e6;
    if(newt<=0.70+0.01) return  2.903e-06*1e6;
    if(newt<=0.71+0.01) return  2.872e-06*1e6;
    if(newt<=0.72+0.01) return  2.826e-06*1e6;
    if(newt<=0.73+0.01) return  2.769e-06*1e6;
    if(newt<=0.74+0.01) return  2.704e-06*1e6;
    if(newt<=0.75+0.01) return  2.636e-06*1e6;
    if(newt<=0.76+0.01) return  2.569e-06*1e6;
    if(newt<=0.77+0.01) return  2.511e-06*1e6;
    if(newt<=0.78+0.01) return  2.465e-06*1e6;
    if(newt<=0.79+0.01) return  2.433e-06*1e6;
    else return 0;
}// 15, LCCA, branch 2

Real aortaFlux9(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
        Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    //  if(newt<=0.1) return (10.55)*(t);
    if(newt<=0.00+0.01) return  7.854e-07*1e6;
    if(newt<=0.01+0.01) return  7.900e-07*1e6;
    if(newt<=0.02+0.01) return  7.988e-07*1e6;
    if(newt<=0.03+0.01) return  8.077e-07*1e6;
    if(newt<=0.04+0.01) return  8.131e-07*1e6;
    if(newt<=0.05+0.01) return  8.141e-07*1e6;
    if(newt<=0.06+0.01) return  8.120e-07*1e6;
    if(newt<=0.07+0.01) return  8.085e-07*1e6;
    if(newt<=0.08+0.01) return  8.082e-07*1e6;
    if(newt<=0.09+0.01) return  8.255e-07*1e6;
    if(newt<=0.10+0.01) return  8.927e-07*1e6;
    if(newt<=0.11+0.01) return  1.055e-06*1e6;
    if(newt<=0.12+0.01) return  1.343e-06*1e6;
    if(newt<=0.13+0.01) return  1.728e-06*1e6;
    if(newt<=0.14+0.01) return  2.118e-06*1e6;
    if(newt<=0.15+0.01) return  2.419e-06*1e6;
    if(newt<=0.16+0.01) return  2.598e-06*1e6;
    if(newt<=0.17+0.01) return  2.673e-06*1e6;
    if(newt<=0.18+0.01) return  2.642e-06*1e6;
    if(newt<=0.19+0.01) return  2.465e-06*1e6;
    if(newt<=0.20+0.01) return  2.118e-06*1e6;
    if(newt<=0.21+0.01) return  1.660e-06*1e6;
    if(newt<=0.22+0.01) return  1.210e-06*1e6;
    if(newt<=0.23+0.01) return  8.723e-07*1e6;
    if(newt<=0.24+0.01) return  6.812e-07*1e6;
    if(newt<=0.25+0.01) return  6.155e-07*1e6;
    if(newt<=0.26+0.01) return  6.323e-07*1e6;
    if(newt<=0.27+0.01) return  6.927e-07*1e6;
    if(newt<=0.28+0.01) return  7.708e-07*1e6;
    if(newt<=0.29+0.01) return  8.538e-07*1e6;
    if(newt<=0.30+0.01) return  9.397e-07*1e6;
    if(newt<=0.31+0.01) return  1.030e-06*1e6;
    if(newt<=0.32+0.01) return  1.120e-06*1e6;
    if(newt<=0.33+0.01) return  1.203e-06*1e6;
    if(newt<=0.34+0.01) return  1.263e-06*1e6;
    if(newt<=0.35+0.01) return  1.292e-06*1e6;
    if(newt<=0.36+0.01) return  1.288e-06*1e6;
    if(newt<=0.37+0.01) return  1.258e-06*1e6;
    if(newt<=0.38+0.01) return  1.207e-06*1e6;
    if(newt<=0.39+0.01) return  1.142e-06*1e6;
    if(newt<=0.40+0.01) return  1.064e-06*1e6;
    if(newt<=0.41+0.01) return  9.816e-07*1e6;
    if(newt<=0.42+0.01) return  9.133e-07*1e6;
    if(newt<=0.43+0.01) return  8.870e-07*1e6;
    if(newt<=0.44+0.01) return  9.268e-07*1e6;
    if(newt<=0.45+0.01) return  1.035e-06*1e6;
    if(newt<=0.46+0.01) return  1.183e-06*1e6;
    if(newt<=0.47+0.01) return  1.329e-06*1e6;
    if(newt<=0.48+0.01) return  1.445e-06*1e6;
    if(newt<=0.49+0.01) return  1.529e-06*1e6;
    if(newt<=0.50+0.01) return  1.587e-06*1e6;
    if(newt<=0.51+0.01) return  1.612e-06*1e6;
    if(newt<=0.52+0.01) return  1.587e-06*1e6;
    if(newt<=0.53+0.01) return  1.501e-06*1e6;
    if(newt<=0.54+0.01) return  1.370e-06*1e6;
    if(newt<=0.55+0.01) return  1.230e-06*1e6;
    if(newt<=0.56+0.01) return  1.112e-06*1e6;
    if(newt<=0.57+0.01) return  1.029e-06*1e6;
    if(newt<=0.58+0.01) return  9.776e-07*1e6;
    if(newt<=0.59+0.01) return  9.460e-07*1e6;
    if(newt<=0.60+0.01) return  9.239e-07*1e6;
    if(newt<=0.61+0.01) return  9.061e-07*1e6;
    if(newt<=0.62+0.01) return  8.917e-07*1e6;
    if(newt<=0.63+0.01) return  8.821e-07*1e6;
    if(newt<=0.64+0.01) return  8.786e-07*1e6;
    if(newt<=0.65+0.01) return  8.809e-07*1e6;
    if(newt<=0.66+0.01) return  8.871e-07*1e6;
    if(newt<=0.67+0.01) return  8.940e-07*1e6;
    if(newt<=0.68+0.01) return  8.978e-07*1e6;
    if(newt<=0.69+0.01) return  8.959e-07*1e6;
    if(newt<=0.70+0.01) return  8.874e-07*1e6;
    if(newt<=0.71+0.01) return  8.734e-07*1e6;
    if(newt<=0.72+0.01) return  8.561e-07*1e6;
    if(newt<=0.73+0.01) return  8.383e-07*1e6;
    if(newt<=0.74+0.01) return  8.221e-07*1e6;
    if(newt<=0.75+0.01) return  8.088e-07*1e6;
    if(newt<=0.76+0.01) return  7.986e-07*1e6;
    if(newt<=0.77+0.01) return  7.909e-07*1e6;
    if(newt<=0.78+0.01) return  7.852e-07*1e6;
    if(newt<=0.79+0.01) return  7.807e-07*1e6;
    else
      return 0;
}// 20 LVA branch 3_1

Real aortaFlux4(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
        Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    //  if(newt<=0.1) return (24.60)*(t);
    if(newt<=0.00+0.01) return  2.581e-06 *1e6;
    if(newt<=0.01+0.01) return  2.499e-06 *1e6;
    if(newt<=0.02+0.01) return  2.399e-06 *1e6;
    if(newt<=0.03+0.01) return  2.281e-06 *1e6;
    if(newt<=0.04+0.01) return  2.147e-06 *1e6;
    if(newt<=0.05+0.01) return  2.002e-06 *1e6;
    if(newt<=0.06+0.01) return  1.849e-06 *1e6;
    if(newt<=0.07+0.01) return  1.695e-06 *1e6;
    if(newt<=0.08+0.01) return  1.565e-06 *1e6;
    if(newt<=0.09+0.01) return  1.537e-06 *1e6;
    if(newt<=0.10+0.01) return  1.766e-06 *1e6;
    if(newt<=0.11+0.01) return  2.460e-06 *1e6;
    if(newt<=0.12+0.01) return  3.773e-06 *1e6;
    if(newt<=0.13+0.01) return  5.709e-06 *1e6;
    if(newt<=0.14+0.01) return  8.131e-06 *1e6;
    if(newt<=0.15+0.01) return  1.087e-05 *1e6;
    if(newt<=0.16+0.01) return  1.379e-05 *1e6;
    if(newt<=0.17+0.01) return  1.675e-05 *1e6;
    if(newt<=0.18+0.01) return  1.957e-05 *1e6;
    if(newt<=0.19+0.01) return  2.205e-05 *1e6;
    if(newt<=0.20+0.01) return  2.403e-05 *1e6;
    if(newt<=0.21+0.01) return  2.538e-05 *1e6;
    if(newt<=0.22+0.01) return  2.606e-05 *1e6;
    if(newt<=0.23+0.01) return  2.611e-05 *1e6;
    if(newt<=0.24+0.01) return  2.562e-05 *1e6;
    if(newt<=0.25+0.01) return  2.474e-05 *1e6;
    if(newt<=0.26+0.01) return  2.362e-05 *1e6;
    if(newt<=0.27+0.01) return  2.238e-05 *1e6;
    if(newt<=0.28+0.01) return  2.111e-05 *1e6;
    if(newt<=0.29+0.01) return  1.981e-05 *1e6;
    if(newt<=0.30+0.01) return  1.850e-05 *1e6;
    if(newt<=0.31+0.01) return  1.715e-05 *1e6;
    if(newt<=0.32+0.01) return  1.576e-05 *1e6;
    if(newt<=0.33+0.01) return  1.432e-05 *1e6;
    if(newt<=0.34+0.01) return  1.284e-05 *1e6;
    if(newt<=0.35+0.01) return  1.132e-05 *1e6;
    if(newt<=0.36+0.01) return  9.768e-06 *1e6;
    if(newt<=0.37+0.01) return  8.180e-06 *1e6;
    if(newt<=0.38+0.01) return  6.543e-06 *1e6;
    if(newt<=0.39+0.01) return  4.831e-06 *1e6;
    if(newt<=0.40+0.01) return  3.030e-06 *1e6;
    if(newt<=0.41+0.01) return  1.163e-06 *1e6;
    if(newt<=0.42+0.01) return  -6.817e-07*1e6;
    if(newt<=0.43+0.01) return  -2.362e-06*1e6;
    if(newt<=0.44+0.01) return  -3.738e-06*1e6;
    if(newt<=0.45+0.01) return  -4.742e-06*1e6;
    if(newt<=0.46+0.01) return  -5.400e-06*1e6;
    if(newt<=0.47+0.01) return  -5.785e-06*1e6;
    if(newt<=0.48+0.01) return  -5.956e-06*1e6;
    if(newt<=0.49+0.01) return  -5.932e-06*1e6;
    if(newt<=0.50+0.01) return  -5.723e-06*1e6;
    if(newt<=0.51+0.01) return  -5.358e-06*1e6;
    if(newt<=0.52+0.01) return  -4.889e-06*1e6;
    if(newt<=0.53+0.01) return  -4.370e-06*1e6;
    if(newt<=0.54+0.01) return  -3.846e-06*1e6;
    if(newt<=0.55+0.01) return  -3.341e-06*1e6;
    if(newt<=0.56+0.01) return  -2.866e-06*1e6;
    if(newt<=0.57+0.01) return  -2.422e-06*1e6;
    if(newt<=0.58+0.01) return  -2.004e-06*1e6;
    if(newt<=0.59+0.01) return  -1.601e-06*1e6;
    if(newt<=0.60+0.01) return  -1.206e-06*1e6;
    if(newt<=0.61+0.01) return  -8.123e-07*1e6;
    if(newt<=0.62+0.01) return  -4.195e-07*1e6;
    if(newt<=0.63+0.01) return  -3.198e-08*1e6;
    if(newt<=0.64+0.01) return  3.429e-07 *1e6;
    if(newt<=0.65+0.01) return  6.970e-07 *1e6;
    if(newt<=0.66+0.01) return  1.024e-06 *1e6;
    if(newt<=0.67+0.01) return  1.320e-06 *1e6;
    if(newt<=0.68+0.01) return  1.583e-06 *1e6;
    if(newt<=0.69+0.01) return  1.813e-06 *1e6;
    if(newt<=0.70+0.01) return  2.011e-06 *1e6;
    if(newt<=0.71+0.01) return  2.180e-06 *1e6;
    if(newt<=0.72+0.01) return  2.323e-06 *1e6;
    if(newt<=0.73+0.01) return  2.442e-06 *1e6;
    if(newt<=0.74+0.01) return  2.539e-06 *1e6;
    if(newt<=0.75+0.01) return  2.613e-06 *1e6;
    if(newt<=0.76+0.01) return  2.665e-06 *1e6;
    if(newt<=0.77+0.01) return  2.691e-06 *1e6;
    if(newt<=0.78+0.01) return  2.690e-06 *1e6;
    if(newt<=0.79+0.01) return  2.661e-06 *1e6;
    else return 0.0;
}//21, L. Brachia, bhanch 3_2


// Real aortaPhisPress2(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
// {
//     /*switch(i) {
//       case 1:
//       return 0.0;
//       break;*/
//     //  case 2:
//         Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
//     if(newt<=0.00+0.01) return -11017;
//     if(newt<=0.01+0.01) return -10954;
//     if(newt<=0.02+0.01) return -10893;
//     if(newt<=0.03+0.01) return -10832;
//     if(newt<=0.04+0.01) return -10771;
//     if(newt<=0.05+0.01) return -10712;
//     if(newt<=0.06+0.01) return -10653;
//     if(newt<=0.07+0.01) return -10800;
//     if(newt<=0.08+0.01) return -11000;
//     if(newt<=0.09+0.01) return -11200;
//     if(newt<=0.10+0.01) return -11400;
//     if(newt<=0.11+0.01) return -11600;
//     if(newt<=0.12+0.01) return -11800;
//     if(newt<=0.13+0.01) return -12000;
//     if(newt<=0.14+0.01) return -12200;
//     if(newt<=0.15+0.01) return -12400;
//     if(newt<=0.16+0.01) return -12600;
//     if(newt<=0.17+0.01) return -12800;
//     if(newt<=0.18+0.01) return -13000;
//     if(newt<=0.19+0.01) return -13200;
//     if(newt<=0.20+0.01) return -13436;
//     if(newt<=0.21+0.01) return -13613;
//     if(newt<=0.22+0.01) return -13800;
//     if(newt<=0.23+0.01) return -14000;
//     if(newt<=0.24+0.01) return -14200;
//     if(newt<=0.25+0.01) return -14400;
//     if(newt<=0.26+0.01) return -14600;
//     if(newt<=0.27+0.01) return -14800;
//     if(newt<=0.28+0.01) return -14975;
//     if(newt<=0.29+0.01) return -14899;
//     if(newt<=0.30+0.01) return -14822;
//     if(newt<=0.31+0.01) return -14721;
//     if(newt<=0.32+0.01) return -14594;
//     if(newt<=0.33+0.01) return -14496;
//     if(newt<=0.34+0.01) return -14375;
//     if(newt<=0.35+0.01) return -14198;
//     if(newt<=0.36+0.01) return -13990;
//     if(newt<=0.37+0.01) return -13726;
//     if(newt<=0.38+0.01) return -13397;
//     if(newt<=0.39+0.01) return -13167;
//     if(newt<=0.40+0.01) return -13132;
//     if(newt<=0.41+0.01) return -13315;
//     if(newt<=0.42+0.01) return -13271;
//     if(newt<=0.43+0.01) return -13157;
//     if(newt<=0.44+0.01) return -13028;
//     if(newt<=0.45+0.01) return -12975;
//     if(newt<=0.46+0.01) return -12933;
//     if(newt<=0.47+0.01) return -12891;
//     if(newt<=0.48+0.01) return -12836;
//     if(newt<=0.49+0.01) return -12768;
//     if(newt<=0.50+0.01) return -12700;
//     if(newt<=0.51+0.01) return -12641;
//     if(newt<=0.52+0.01) return -12592;
//     if(newt<=0.53+0.01) return -12548;
//     if(newt<=0.54+0.01) return -12504;
//     if(newt<=0.55+0.01) return -12456;
//     if(newt<=0.56+0.01) return -12405;
//     if(newt<=0.57+0.01) return -12353;
//     if(newt<=0.58+0.01) return -12300;
//     if(newt<=0.59+0.01) return -12244;
//     if(newt<=0.60+0.01) return -12184;
//     if(newt<=0.61+0.01) return -12122;
//     if(newt<=0.62+0.01) return -12058;
//     if(newt<=0.63+0.01) return -11995;
//     if(newt<=0.64+0.01) return -11933;
//     if(newt<=0.65+0.01) return -11871;
//     if(newt<=0.66+0.01) return -11810;
//     if(newt<=0.67+0.01) return -11747;
//     if(newt<=0.68+0.01) return -11684;
//     if(newt<=0.69+0.01) return -11620;
//     if(newt<=0.70+0.01) return -11556;
//     if(newt<=0.71+0.01) return -11492;
//     if(newt<=0.72+0.01) return -11428;
//     if(newt<=0.73+0.01) return -11365;
//     if(newt<=0.74+0.01) return -11302;
//     if(newt<=0.75+0.01) return -11240;
//     if(newt<=0.76+0.01) return -11179;
//     if(newt<=0.77+0.01) return -11120;
//     if(newt<=0.78+0.01) return -11062;
//     if(newt<=0.79+0.01) return -11006;
//     //    break;
//     /*  case 3:
//         return 0.0;
//         break;}
//         return 0.;*/
// }

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


Real E(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    if(z<=5. && z>0.)
        return /*-*/60000.;//-29;//*e-1*5; // circa [(110-60)*(133.332*10)]/[10*(2.08-1.8)] (il 10 x via dei mm, il 133... x via dei mmHg)
    // (vedi paper di Liu, Dang, etc). Nel loro grafico sono invertite x e y. 5e1 invece di e5 xche' d e'
    // riscalato * il timestep, che e' 5e-4
    if(z>5 && z<6)
      return 60000.+30000*(z-5);
    if(z>=6)
        return 90000;
  if ( z<0 && z>=-3 )
    return 60000.;//500-500*(z);
  if(z<-3)
    return 60000.;
  else 
    return 0;}

Real hydro(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 110170;
}

Real u2(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch(i) {
    case 1:
        return 0.0;
        break;
    case 3:
        if ( t <= 0.003 )
            return 1.3332e4;//1.3332e5;
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


Real aortaFluxIn(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
    if(newt<=0.00+0.01) return	0.0000e+00;
    if(newt<=0.01+0.01) return	0.0000e+00;
    if(newt<=0.02+0.01) return	0.0000e+00;
    if(newt<=0.03+0.01) return	0.0000e+00;
    if(newt<=0.04+0.01) return	-9.1759e-06*1e6;
    if(newt<=0.05+0.01) return	-3.0930e-05*1e6;
    if(newt<=0.06+0.01) return	-6.2639e-05*1e6;
    if(newt<=0.07+0.01) return	-1.0212e-04*1e6;
    if(newt<=0.08+0.01) return	-1.4760e-04*1e6;
    if(newt<=0.09+0.01) return	-1.9726e-04*1e6;
    if(newt<=0.10+0.01) return	-2.4980e-04*1e6;
    if(newt<=0.11+0.01) return	-2.9526e-04*1e6;
    if(newt<=0.12+0.01) return	-3.2956e-04*1e6;
    if(newt<=0.13+0.01) return	-3.5469e-04*1e6;
    if(newt<=0.14+0.01) return	-3.7250e-04*1e6;
    if(newt<=0.15+0.01) return	-3.8429e-04*1e6;
    if(newt<=0.16+0.01) return	-3.9123e-04*1e6;
    if(newt<=0.17+0.01) return	-3.9431e-04*1e6;
    if(newt<=0.18+0.01) return	-3.9349e-04*1e6;
    if(newt<=0.19+0.01) return	-3.8858e-04*1e6;
    if(newt<=0.20+0.01) return	-3.7985e-04*1e6;
    if(newt<=0.21+0.01) return	-3.6756e-04*1e6;
    if(newt<=0.22+0.01) return	-3.5207e-04*1e6;
    if(newt<=0.23+0.01) return	-3.3408e-04*1e6;
    if(newt<=0.24+0.01) return	-3.1402e-04*1e6;
    if(newt<=0.25+0.01) return	-2.9288e-04*1e6;
    if(newt<=0.26+0.01) return	-2.7154e-04*1e6;
    if(newt<=0.27+0.01) return	-2.5054e-04*1e6;
    if(newt<=0.28+0.01) return	-2.2979e-04*1e6;
    if(newt<=0.29+0.01) return	-2.0904e-04*1e6;
    if(newt<=0.30+0.01) return	-1.8880e-04*1e6;
    if(newt<=0.31+0.01) return	-1.6899e-04*1e6;
    if(newt<=0.32+0.01) return	-1.4864e-04*1e6;
    if(newt<=0.33+0.01) return	-1.2730e-04*1e6;
    if(newt<=0.34+0.01) return	-1.0400e-04*1e6;
    if(newt<=0.35+0.01) return	-7.9755e-05*1e6;
    if(newt<=0.36+0.01) return	-5.8719e-05*1e6;
    if(newt<=0.37+0.01) return	-4.0345e-05*1e6;
    if(newt<=0.38+0.01) return	-2.4596e-05*1e6;
    if(newt<=0.39+0.01) return	-1.2259e-05*1e6;
    if(newt<=0.40+0.01) return	-3.8110e-06*1e6;
    if(newt<=0.41+0.01) return	0.0000e+00;
    if(newt<=0.42+0.01) return	0.0000e+00;
    if(newt<=0.43+0.01) return	0.0000e+00;
    if(newt<=0.44+0.01) return	0.0000e+00;
    if(newt<=0.45+0.01) return	0.0000e+00;
    if(newt<=0.46+0.01) return	0.0000e+00;
    if(newt<=0.47+0.01) return	0.0000e+00;
    if(newt<=0.48+0.01) return	0.0000e+00;
    if(newt<=0.49+0.01) return	0.0000e+00;
    if(newt<=0.50+0.01) return	0.0000e+00;
    if(newt<=0.51+0.01) return	0.0000e+00;
    if(newt<=0.52+0.01) return	0.0000e+00;
    if(newt<=0.53+0.01) return	0.0000e+00;
    if(newt<=0.54+0.01) return	0.0000e+00;
    if(newt<=0.55+0.01) return	0.0000e+00;
    if(newt<=0.56+0.01) return	0.0000e+00;
    if(newt<=0.57+0.01) return	0.0000e+00;
    if(newt<=0.58+0.01) return	0.0000e+00;
    if(newt<=0.59+0.01) return	0.0000e+00;
    if(newt<=0.60+0.01) return	0.0000e+00;
    if(newt<=0.61+0.01) return	0.0000e+00;
    if(newt<=0.62+0.01) return	0.0000e+00;
    if(newt<=0.63+0.01) return	0.0000e+00;
    if(newt<=0.64+0.01) return	0.0000e+00;
    if(newt<=0.65+0.01) return	0.0000e+00;
    if(newt<=0.66+0.01) return	0.0000e+00;
    if(newt<=0.67+0.01) return	0.0000e+00;
    if(newt<=0.68+0.01) return	0.0000e+00;
    if(newt<=0.69+0.01) return	0.0000e+00;
    if(newt<=0.70+0.01) return	0.0000e+00;
    if(newt<=0.71+0.01) return	0.0000e+00;
    if(newt<=0.72+0.01) return	0.0000e+00;
    if(newt<=0.73+0.01) return	0.0000e+00;
    if(newt<=0.74+0.01) return	0.0000e+00;
    if(newt<=0.75+0.01) return	0.0000e+00;
    if(newt<=0.76+0.01) return	0.0000e+00;
    if(newt<=0.77+0.01) return	0.0000e+00;
    if(newt<=0.78+0.01) return	0.0000e+00;
    if(newt<=0.79+0.01) return	0.0000e+00;
    else return 0;
    //    break;       
}

Real aortaVelocityIn(const Real&  t, const Real&/* x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return linearFluxIn(t,0,0,0,0)/7.23681;  
}

Real aortaVelIn::func(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(i==1)
    {
        if(t<=0.00+0.01) return	0.0000e+00;
        if(t<=0.01+0.01) return	0.0000e+00;
        if(t<=0.02+0.01) return	0.0000e+00;
        if(t<=0.03+0.01) return	0.0000e+00;
        if(t<=0.04+0.01) return	9.1759e-06*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.05+0.01) return	3.0930e-05*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.06+0.01) return	6.2639e-05*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.07+0.01) return	1.0212e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.08+0.01) return	1.4760e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.09+0.01) return	1.9726e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.10+0.01) return	2.4980e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.11+0.01) return	2.9526e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.12+0.01) return	3.2956e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.13+0.01) return	3.5469e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.14+0.01) return	3.7250e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.15+0.01) return	3.8429e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.16+0.01) return	3.9123e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.17+0.01) return	3.9431e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.18+0.01) return	3.9349e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.19+0.01) return	3.8858e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.20+0.01) return	3.7985e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.21+0.01) return	3.6756e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.22+0.01) return	3.5207e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.23+0.01) return	3.3408e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.24+0.01) return	3.1402e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.25+0.01) return	2.9288e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.26+0.01) return	2.7154e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.27+0.01) return	2.5054e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.28+0.01) return	2.2979e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.29+0.01) return	2.0904e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.30+0.01) return	1.8880e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.31+0.01) return	1.6899e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.32+0.01) return	1.4864e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.33+0.01) return	1.2730e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.34+0.01) return	1.0400e-04*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.35+0.01) return	7.9755e-05*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.36+0.01) return	5.8719e-05*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.37+0.01) return	4.0345e-05*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.38+0.01) return	2.4596e-05*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.39+0.01) return	1.2259e-05*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.40+0.01) return	3.8110e-06*1e6/aortaVelIn::A*(-0.07478994824163);
        if(t<=0.41+0.01) return	0.0000e+00;
        if(t<=0.42+0.01) return	0.0000e+00;
        if(t<=0.43+0.01) return	0.0000e+00;
        if(t<=0.44+0.01) return	0.0000e+00;
        if(t<=0.45+0.01) return	0.0000e+00;
        if(t<=0.46+0.01) return	0.0000e+00;
        if(t<=0.47+0.01) return	0.0000e+00;
        if(t<=0.48+0.01) return	0.0000e+00;
        if(t<=0.49+0.01) return	0.0000e+00;
        if(t<=0.50+0.01) return	0.0000e+00;
        if(t<=0.51+0.01) return	0.0000e+00;
        if(t<=0.52+0.01) return	0.0000e+00;
        if(t<=0.53+0.01) return	0.0000e+00;
        if(t<=0.54+0.01) return	0.0000e+00;
        if(t<=0.55+0.01) return	0.0000e+00;
        if(t<=0.56+0.01) return	0.0000e+00;
        if(t<=0.57+0.01) return	0.0000e+00;
        if(t<=0.58+0.01) return	0.0000e+00;
        if(t<=0.59+0.01) return	0.0000e+00;
        if(t<=0.60+0.01) return	0.0000e+00;
        if(t<=0.61+0.01) return	0.0000e+00;
        if(t<=0.62+0.01) return	0.0000e+00;
        if(t<=0.63+0.01) return	0.0000e+00;
        if(t<=0.64+0.01) return	0.0000e+00;
        if(t<=0.65+0.01) return	0.0000e+00;
        if(t<=0.66+0.01) return	0.0000e+00;
        if(t<=0.67+0.01) return	0.0000e+00;
        if(t<=0.68+0.01) return	0.0000e+00;
        if(t<=0.69+0.01) return	0.0000e+00;
        if(t<=0.70+0.01) return	0.0000e+00;
        if(t<=0.71+0.01) return	0.0000e+00;
        if(t<=0.72+0.01) return	0.0000e+00;
        if(t<=0.73+0.01) return	0.0000e+00;
        if(t<=0.74+0.01) return	0.0000e+00;
        if(t<=0.75+0.01) return	0.0000e+00;
        if(t<=0.76+0.01) return	0.0000e+00;
        if(t<=0.77+0.01) return	0.0000e+00;
        if(t<=0.78+0.01) return	0.0000e+00;
        if(t<=0.79+0.01) return	0.0000e+00;
    }
    if(i==3)
    {
        if(t<=0.00+0.01) return	0.0000e+00;
        if(t<=0.01+0.01) return	0.0000e+00;
        if(t<=0.02+0.01) return	0.0000e+00;
        if(t<=0.03+0.01) return	0.0000e+00;
        if(t<=0.04+0.01) return	9.1759e-06*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.05+0.01) return	3.0930e-05*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.06+0.01) return	6.2639e-05*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.07+0.01) return	1.0212e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.08+0.01) return	1.4760e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.09+0.01) return	1.9726e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.10+0.01) return	2.4980e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.11+0.01) return	2.9526e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.12+0.01) return	3.2956e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.13+0.01) return	3.5469e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.14+0.01) return	3.7250e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.15+0.01) return	3.8429e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.16+0.01) return	3.9123e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.17+0.01) return	3.9431e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.18+0.01) return	3.9349e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.19+0.01) return	3.8858e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.20+0.01) return	3.7985e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.21+0.01) return	3.6756e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.22+0.01) return	3.5207e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.23+0.01) return	3.3408e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.24+0.01) return	3.1402e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.25+0.01) return	2.9288e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.26+0.01) return	2.7154e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.27+0.01) return	2.5054e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.28+0.01) return	2.2979e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.29+0.01) return	2.0904e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.30+0.01) return	1.8880e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.31+0.01) return	1.6899e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.32+0.01) return	1.4864e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.33+0.01) return	1.2730e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.34+0.01) return	1.0400e-04*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.35+0.01) return	7.9755e-05*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.36+0.01) return	5.8719e-05*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.37+0.01) return	4.0345e-05*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.38+0.01) return	2.4596e-05*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.39+0.01) return	1.2259e-05*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.40+0.01) return	3.8110e-06*1e6/aortaVelIn::A*(0.997199309888456);
        if(t<=0.41+0.01) return	0.0000e+00;
        if(t<=0.42+0.01) return	0.0000e+00;
        if(t<=0.43+0.01) return	0.0000e+00;
        if(t<=0.44+0.01) return	0.0000e+00;
        if(t<=0.45+0.01) return	0.0000e+00;
        if(t<=0.46+0.01) return	0.0000e+00;
        if(t<=0.47+0.01) return	0.0000e+00;
        if(t<=0.48+0.01) return	0.0000e+00;
        if(t<=0.49+0.01) return	0.0000e+00;
        if(t<=0.50+0.01) return	0.0000e+00;
        if(t<=0.51+0.01) return	0.0000e+00;
        if(t<=0.52+0.01) return	0.0000e+00;
        if(t<=0.53+0.01) return	0.0000e+00;
        if(t<=0.54+0.01) return	0.0000e+00;
        if(t<=0.55+0.01) return	0.0000e+00;
        if(t<=0.56+0.01) return	0.0000e+00;
        if(t<=0.57+0.01) return	0.0000e+00;
        if(t<=0.58+0.01) return	0.0000e+00;
        if(t<=0.59+0.01) return	0.0000e+00;
        if(t<=0.60+0.01) return	0.0000e+00;
        if(t<=0.61+0.01) return	0.0000e+00;
        if(t<=0.62+0.01) return	0.0000e+00;
        if(t<=0.63+0.01) return	0.0000e+00;
        if(t<=0.64+0.01) return	0.0000e+00;
        if(t<=0.65+0.01) return	0.0000e+00;
        if(t<=0.66+0.01) return	0.0000e+00;
        if(t<=0.67+0.01) return	0.0000e+00;
        if(t<=0.68+0.01) return	0.0000e+00;
        if(t<=0.69+0.01) return	0.0000e+00;
        if(t<=0.70+0.01) return	0.0000e+00;
        if(t<=0.71+0.01) return	0.0000e+00;
        if(t<=0.72+0.01) return	0.0000e+00;
        if(t<=0.73+0.01) return	0.0000e+00;
        if(t<=0.74+0.01) return	0.0000e+00;
        if(t<=0.75+0.01) return	0.0000e+00;
        if(t<=0.76+0.01) return	0.0000e+00;
        if(t<=0.77+0.01) return	0.0000e+00;
        if(t<=0.78+0.01) return	0.0000e+00;
        if(t<=0.79+0.01) return	0.0000e+00;
    }
    if(i==2)
        return 0.;
    else
      return 0.;
}
void setArea(const Real& area)
{
    LifeV::aortaVelIn::A=area;
}
Real LifeV::aortaVelIn::A;

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

Real cylinderFlux(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/) //outlet flux from the bifurcation
{
        Real newt = ((Real)(((int)(t*1000)) % 1096))/1000;
if( newt <= 0.    +0.004 ) return     4.5          *(-1.0); else 
if( newt <= 0.004 +0.004 ) return     4.4687       *(-1.0); else 
if( newt <= 0.008 +0.004 ) return     4.4291       *(-1.0); else 
if( newt <= 0.012 +0.004 ) return     4.3816       *(-1.0); else 
if( newt <= 0.016 +0.004 ) return     4.3269       *(-1.0); else 
if( newt <= 0.02  +0.004 ) return     4.266        *(-1.0); else         
if( newt <= 0.024 +0.004 ) return     4.2007       *(-1.0); else 
if( newt <= 0.028 +0.004 ) return     4.1331       *(-1.0); else         
if( newt <= 0.032 +0.004 ) return     4.0656       *(-1.0); else 
if( newt <= 0.036 +0.004 ) return     4.0014       *(-1.0); else 
if( newt <= 0.04  +0.004 ) return     3.944        *(-1.0); else         
if( newt <= 0.044 +0.004 ) return     3.8969       *(-1.0); else 
if( newt <= 0.048 +0.004 ) return     3.8643       *(-1.0); else 
if( newt <= 0.052 +0.004 ) return     3.8503       *(-1.0); else 
if( newt <= 0.056 +0.004 ) return     3.8589       *(-1.0); else 
if( newt <= 0.06  +0.004 ) return     3.8942       *(-1.0); else         
if( newt <= 0.064 +0.004 ) return     3.9601       *(-1.0); else 
if( newt <= 0.068 +0.004 ) return     4.0597       *(-1.0); else 
if( newt <= 0.072 +0.004 ) return     4.1962       *(-1.0); else 
if( newt <= 0.076 +0.004 ) return     4.3716       *(-1.0); else         
if( newt <= 0.08  +0.004 ) return     4.5876       *(-1.0); else         
if( newt <= 0.084 +0.004 ) return     4.8447       *(-1.0); else 
if( newt <= 0.088 +0.004 ) return     5.1428       *(-1.0); else         
if( newt <= 0.092 +0.004 ) return     5.4808       *(-1.0); else 
if( newt <= 0.096 +0.004 ) return     5.8566       *(-1.0); else 
if( newt <= 0.1   +0.004 ) return     6.2672       *(-1.0); else         
if( newt <= 0.104 +0.004 ) return     6.7088       *(-1.0); else 
if( newt <= 0.108 +0.004 ) return     7.1767       *(-1.0); else 
if( newt <= 0.112 +0.004 ) return     7.6657       *(-1.0); else 
if( newt <= 0.116 +0.004 ) return     8.1698       *(-1.0); else 
if( newt <= 0.12  +0.004 ) return     8.6827       *(-1.0); else         
if( newt <= 0.124 +0.004 ) return     9.198        *(-1.0); else 
if( newt <= 0.128 +0.004 ) return     9.7091       *(-1.0); else 
if( newt <= 0.132 +0.004 ) return     10.209       *(-1.0); else 
if( newt <= 0.136 +0.004 ) return     10.693       *(-1.0); else 
if( newt <= 0.14  +0.004 ) return     11.153       *(-1.0); else         
if( newt <= 0.144 +0.004 ) return     11.586       *(-1.0); else 
if( newt <= 0.148 +0.004 ) return     11.985       *(-1.0); else 
if( newt <= 0.152 +0.004 ) return     12.349       *(-1.0); else 
if( newt <= 0.156 +0.004 ) return     12.673       *(-1.0); else 
if( newt <= 0.16  +0.004 ) return     12.957       *(-1.0); else         
if( newt <= 0.164 +0.004 ) return     13.198       *(-1.0); else 
if( newt <= 0.168 +0.004 ) return     13.397       *(-1.0); else 
if( newt <= 0.172 +0.004 ) return     13.555       *(-1.0); else 
if( newt <= 0.176 +0.004 ) return     13.674       *(-1.0); else         
if( newt <= 0.18  +0.004 ) return     13.754       *(-1.0); else         
if( newt <= 0.184 +0.004 ) return     13.801       *(-1.0); else         
if( newt <= 0.188 +0.004 ) return     13.816       *(-1.0); else 
if( newt <= 0.192 +0.004 ) return     13.803       *(-1.0); else         
if( newt <= 0.196 +0.004 ) return     13.767       *(-1.0); else 
if( newt <= 0.2   +0.004 ) return     13.711       *(-1.0); else         
if( newt <= 0.204 +0.004 ) return     13.638       *(-1.0); else 
if( newt <= 0.208 +0.004 ) return     13.553       *(-1.0); else 
if( newt <= 0.212 +0.004 ) return     13.458       *(-1.0); else 
if( newt <= 0.216 +0.004 ) return     13.357       *(-1.0); else 
if( newt <= 0.22  +0.004 ) return     13.251       *(-1.0); else         
if( newt <= 0.224 +0.004 ) return     13.142       *(-1.0); else 
if( newt <= 0.228 +0.004 ) return     13.031       *(-1.0); else 
if( newt <= 0.232 +0.004 ) return     12.92        *(-1.0); else 
if( newt <= 0.236 +0.004 ) return     12.808       *(-1.0); else 
if( newt <= 0.24  +0.004 ) return     12.695       *(-1.0); else         
if( newt <= 0.244 +0.004 ) return     12.581       *(-1.0); else 
if( newt <= 0.248 +0.004 ) return     12.465       *(-1.0); else 
if( newt <= 0.252 +0.004 ) return     12.346       *(-1.0); else         
if( newt <= 0.256 +0.004 ) return     12.223       *(-1.0); else 
if( newt <= 0.26  +0.004 ) return     12.095       *(-1.0); else         
if( newt <= 0.264 +0.004 ) return     11.961       *(-1.0); else 
if( newt <= 0.268 +0.004 ) return     11.82        *(-1.0); else 
if( newt <= 0.272 +0.004 ) return     11.672       *(-1.0); else 
if( newt <= 0.276 +0.004 ) return     11.516       *(-1.0); else 
if( newt <= 0.28  +0.004 ) return     11.352       *(-1.0); else 
if( newt <= 0.284 +0.004 ) return     11.18        *(-1.0); else 
if( newt <= 0.288 +0.004 ) return     11.001       *(-1.0); else 
if( newt <= 0.292 +0.004 ) return     10.816       *(-1.0); else 
if( newt <= 0.296 +0.004 ) return     10.625       *(-1.0); else 
if( newt <= 0.3   +0.004 ) return     10.43        *(-1.0); else         
if( newt <= 0.304 +0.004 ) return     10.232       *(-1.0); else 
if( newt <= 0.308 +0.004 ) return     10.032       *(-1.0); else
if( newt <= 0.312 +0.004 ) return     9.8312       *(-1.0); else
if( newt <= 0.316 +0.004 ) return     9.6316       *(-1.0); else
if( newt <= 0.32  +0.004 ) return     9.434        *(-1.0); else
if( newt <= 0.324 +0.004 ) return     9.2395       *(-1.0); else
if( newt <= 0.328 +0.004 ) return     9.0488       *(-1.0); else
if( newt <= 0.332 +0.004 ) return     8.8625       *(-1.0); else
if( newt <= 0.336 +0.004 ) return     8.6811       *(-1.0); else
if( newt <= 0.34  +0.004 ) return     8.5048       *(-1.0); else
if( newt <= 0.344 +0.004 ) return     8.3336       *(-1.0); else
if( newt <= 0.348 +0.004 ) return     8.1673       *(-1.0); else
if( newt <= 0.352 +0.004 ) return     8.0057       *(-1.0); else
if( newt <= 0.356 +0.004 ) return     7.8485       *(-1.0); else
if( newt <= 0.36  +0.004 ) return     7.6953       *(-1.0); else
if( newt <= 0.364 +0.004 ) return     7.5459       *(-1.0); else
if( newt <= 0.368 +0.004 ) return     7.4          *(-1.0); else
if( newt <= 0.372 +0.004 ) return     7.2575       *(-1.0); else
if( newt <= 0.376 +0.004 ) return     7.1185       *(-1.0); else
if( newt <= 0.38  +0.004 ) return     6.9833       *(-1.0); else
if( newt <= 0.384 +0.004 ) return     6.8522       *(-1.0); else
if( newt <= 0.388 +0.004 ) return     6.7259       *(-1.0); else
if( newt <= 0.392 +0.004 ) return     6.6052       *(-1.0); else
if( newt <= 0.396 +0.004 ) return     6.4911       *(-1.0); else
if( newt <= 0.4   +0.004 ) return     6.3847       *(-1.0); else        
if( newt <= 0.404 +0.004 ) return     6.2871       *(-1.0); else
if( newt <= 0.408 +0.004 ) return     6.1994       *(-1.0); else
if( newt <= 0.412 +0.004 ) return     6.1227       *(-1.0); else
if( newt <= 0.416 +0.004 ) return     6.0581       *(-1.0); else
if( newt <= 0.42  +0.004 ) return     6.0063       *(-1.0); else
if( newt <= 0.424 +0.004 ) return     5.9679       *(-1.0); else
if( newt <= 0.428 +0.004 ) return     5.9431       *(-1.0); else
if( newt <= 0.432 +0.004 ) return     5.932        *(-1.0); else
if( newt <= 0.436 +0.004 ) return     5.9341       *(-1.0); else
if( newt <= 0.44  +0.004 ) return     5.9487       *(-1.0); else
if( newt <= 0.444 +0.004 ) return     5.9746       *(-1.0); else
if( newt <= 0.448 +0.004 ) return     6.0106       *(-1.0); else
if( newt <= 0.452 +0.004 ) return     6.0547       *(-1.0); else
if( newt <= 0.456 +0.004 ) return     6.1052       *(-1.0); else
if( newt <= 0.46  +0.004 ) return     6.1598       *(-1.0); else
if( newt <= 0.464 +0.004 ) return     6.2162       *(-1.0); else
if( newt <= 0.468 +0.004 ) return     6.2722       *(-1.0); else
if( newt <= 0.472 +0.004 ) return     6.3256       *(-1.0); else
if( newt <= 0.476 +0.004 ) return     6.3741       *(-1.0); else
if( newt <= 0.48  +0.004 ) return     6.4158       *(-1.0); else
if( newt <= 0.484 +0.004 ) return     6.4492       *(-1.0); else
if( newt <= 0.488 +0.004 ) return     6.4727       *(-1.0); else
if( newt <= 0.492 +0.004 ) return     6.4853       *(-1.0); else
if( newt <= 0.496 +0.004 ) return     6.4866       *(-1.0); else
if( newt <= 0.5   +0.004 ) return     6.4761       *(-1.0); else
if( newt <= 0.504 +0.004 ) return     6.4541       *(-1.0); else
if( newt <= 0.508 +0.004 ) return     6.4211       *(-1.0); else
if( newt <= 0.512 +0.004 ) return     6.3779       *(-1.0); else
if( newt <= 0.516 +0.004 ) return     6.3257       *(-1.0); else
if( newt <= 0.52  +0.004 ) return     6.2658       *(-1.0); else
if( newt <= 0.524 +0.004 ) return     6.1998       *(-1.0); else
if( newt <= 0.528 +0.004 ) return     6.1292       *(-1.0); else
if( newt <= 0.532 +0.004 ) return     6.0558       *(-1.0); else
if( newt <= 0.536 +0.004 ) return     5.981        *(-1.0); else
if( newt <= 0.54  +0.004 ) return     5.9063       *(-1.0); else
if( newt <= 0.544 +0.004 ) return     5.833        *(-1.0); else
if( newt <= 0.548 +0.004 ) return     5.7621       *(-1.0); else
if( newt <= 0.552 +0.004 ) return     5.6944       *(-1.0); else
if( newt <= 0.556 +0.004 ) return     5.6304       *(-1.0); else
if( newt <= 0.56  +0.004 ) return     5.5704       *(-1.0); else
if( newt <= 0.564 +0.004 ) return     5.5144       *(-1.0); else
if( newt <= 0.568 +0.004 ) return     5.4619       *(-1.0); else
if( newt <= 0.572 +0.004 ) return     5.4127       *(-1.0); else
if( newt <= 0.576 +0.004 ) return     5.3659       *(-1.0); else
if( newt <= 0.58  +0.004 ) return     5.3208       *(-1.0); else
if( newt <= 0.584 +0.004 ) return     5.2767       *(-1.0); else
if( newt <= 0.588 +0.004 ) return     5.2326       *(-1.0); else
if( newt <= 0.592 +0.004 ) return     5.1877       *(-1.0); else
if( newt <= 0.596 +0.004 ) return     5.1415       *(-1.0); else
if( newt <= 0.6   +0.004 ) return     5.0934       *(-1.0); else
if( newt <= 0.604 +0.004 ) return     5.043        *(-1.0); else
if( newt <= 0.608 +0.004 ) return     4.9902       *(-1.0); else
if( newt <= 0.612 +0.004 ) return     4.9351       *(-1.0); else
if( newt <= 0.616 +0.004 ) return     4.8781       *(-1.0); else
if( newt <= 0.62  +0.004 ) return     4.8197       *(-1.0); else
if( newt <= 0.624 +0.004 ) return     4.7605       *(-1.0); else
if( newt <= 0.628 +0.004 ) return     4.7015       *(-1.0); else
if( newt <= 0.632 +0.004 ) return     4.6436       *(-1.0); else
if( newt <= 0.636 +0.004 ) return     4.5879       *(-1.0); else
if( newt <= 0.64  +0.004 ) return     4.5355       *(-1.0); else
if( newt <= 0.644 +0.004 ) return     4.4872       *(-1.0); else
if( newt <= 0.648 +0.004 ) return     4.444        *(-1.0); else
if( newt <= 0.652 +0.004 ) return     4.4067       *(-1.0); else
if( newt <= 0.656 +0.004 ) return     4.3758       *(-1.0); else
if( newt <= 0.66  +0.004 ) return     4.3517       *(-1.0); else
if( newt <= 0.664 +0.004 ) return     4.3346       *(-1.0); else
if( newt <= 0.668 +0.004 ) return     4.3244       *(-1.0); else
if( newt <= 0.672 +0.004 ) return     4.3206       *(-1.0); else
if( newt <= 0.676 +0.004 ) return     4.3228       *(-1.0); else
if( newt <= 0.68  +0.004 ) return     4.3302       *(-1.0); else
if( newt <= 0.684 +0.004 ) return     4.3419       *(-1.0); else
if( newt <= 0.688 +0.004 ) return     4.3568       *(-1.0); else
if( newt <= 0.692 +0.004 ) return     4.3739       *(-1.0); else
if( newt <= 0.696 +0.004 ) return     4.392        *(-1.0); else
if( newt <= 0.7   +0.004 ) return     4.41         *(-1.0); else
if( newt <= 0.704 +0.004 ) return     4.4268       *(-1.0); else
if( newt <= 0.708 +0.004 ) return     4.4416       *(-1.0); else
if( newt <= 0.712 +0.004 ) return     4.4537       *(-1.0); else
if( newt <= 0.716 +0.004 ) return     4.4624       *(-1.0); else
if( newt <= 0.72  +0.004 ) return     4.4673       *(-1.0); else
if( newt <= 0.724 +0.004 ) return     4.4684       *(-1.0); else
if( newt <= 0.728 +0.004 ) return     4.4657       *(-1.0); else
if( newt <= 0.732 +0.004 ) return     4.4594       *(-1.0); else
if( newt <= 0.736 +0.004 ) return     4.45         *(-1.0); else
if( newt <= 0.74  +0.004 ) return     4.4381       *(-1.0); else
if( newt <= 0.744 +0.004 ) return     4.4243       *(-1.0); else
if( newt <= 0.748 +0.004 ) return     4.4095       *(-1.0); else
if( newt <= 0.752 +0.004 ) return     4.3945       *(-1.0); else
if( newt <= 0.756 +0.004 ) return     4.38         *(-1.0); else
if( newt <= 0.76  +0.004 ) return     4.3669       *(-1.0); else
if( newt <= 0.764 +0.004 ) return     4.3557       *(-1.0); else
if( newt <= 0.768 +0.004 ) return     4.347        *(-1.0); else
if( newt <= 0.772 +0.004 ) return     4.3412       *(-1.0); else
if( newt <= 0.776 +0.004 ) return     4.3385       *(-1.0); else
if( newt <= 0.78  +0.004 ) return     4.339        *(-1.0); else
if( newt <= 0.784 +0.004 ) return     4.3425       *(-1.0); else
if( newt <= 0.788 +0.004 ) return     4.3489       *(-1.0); else
if( newt <= 0.792 +0.004 ) return     4.3576       *(-1.0); else  
if( newt <= 0.796 +0.004 ) return     4.3683       *(-1.0); else 
if( newt <= 0.800  +0.004 ) return     4.3803      *(-1.0); else 
if( newt <= 0.804 +0.004 ) return     4.3929       *(-1.0); else 
if( newt <= 0.808 +0.004 ) return     4.4055       *(-1.0); else 
if( newt <= 0.812 +0.004 ) return     4.4174       *(-1.0); else         
if( newt <= 0.816 +0.004 ) return     4.4281       *(-1.0); else 
if( newt <= 0.82  +0.004 ) return     4.4371       *(-1.0); else         
if( newt <= 0.824 +0.004 ) return     4.4439       *(-1.0); else 
if( newt <= 0.828 +0.004 ) return     4.4483       *(-1.0); else 
if( newt <= 0.832 +0.004 ) return     4.4501       *(-1.0); else         
if( newt <= 0.836 +0.004 ) return     4.4494       *(-1.0); else 
if( newt <= 0.84  +0.004 ) return     4.4463       *(-1.0); else 
if( newt <= 0.844 +0.004 ) return     4.441        *(-1.0); else 
if( newt <= 0.848 +0.004 ) return     4.4339       *(-1.0); else 
if( newt <= 0.852 +0.004 ) return     4.4254       *(-1.0); else         
if( newt <= 0.856 +0.004 ) return     4.4161       *(-1.0); else 
if( newt <= 0.86  +0.004 ) return     4.4064       *(-1.0); else 
if( newt <= 0.864 +0.004 ) return     4.3969       *(-1.0); else 
if( newt <= 0.868 +0.004 ) return     4.3881       *(-1.0); else         
if( newt <= 0.872 +0.004 ) return     4.3804       *(-1.0); else         
if( newt <= 0.876 +0.004 ) return     4.3742       *(-1.0); else 
if( newt <= 0.88  +0.004 ) return     4.3698       *(-1.0); else         
if( newt <= 0.884 +0.004 ) return     4.3674       *(-1.0); else 
if( newt <= 0.888 +0.004 ) return     4.3671       *(-1.0); else 
if( newt <= 0.892 +0.004 ) return     4.3687       *(-1.0); else         
if( newt <= 0.896 +0.004 ) return     4.3722       *(-1.0); else 
if( newt <= 0.9   +0.004 ) return     4.3774       *(-1.0); else 
if( newt <= 0.904 +0.004 ) return     4.384        *(-1.0); else 
if( newt <= 0.908 +0.004 ) return     4.3915       *(-1.0); else 
if( newt <= 0.912 +0.004 ) return     4.3996       *(-1.0); else         
if( newt <= 0.916 +0.004 ) return     4.408        *(-1.0); else 
if( newt <= 0.92  +0.004 ) return     4.4162       *(-1.0); else 
if( newt <= 0.924 +0.004 ) return     4.4238       *(-1.0); else 
if( newt <= 0.928 +0.004 ) return     4.4306       *(-1.0); else 
if( newt <= 0.932 +0.004 ) return     4.4362       *(-1.0); else         
if( newt <= 0.936 +0.004 ) return     4.4406       *(-1.0); else 
if( newt <= 0.94  +0.004 ) return     4.4434       *(-1.0); else 
if( newt <= 0.944 +0.004 ) return     4.4449       *(-1.0); else 
if( newt <= 0.948 +0.004 ) return     4.4448       *(-1.0); else 
if( newt <= 0.952 +0.004 ) return     4.4434       *(-1.0); else         
if( newt <= 0.956 +0.004 ) return     4.4408       *(-1.0); else 
if( newt <= 0.96  +0.004 ) return     4.4371       *(-1.0); else 
if( newt <= 0.964 +0.004 ) return     4.4326       *(-1.0); else 
if( newt <= 0.968 +0.004 ) return     4.4276       *(-1.0); else         
if( newt <= 0.972 +0.004 ) return     4.4221       *(-1.0); else         
if( newt <= 0.976 +0.004 ) return     4.4165       *(-1.0); else         
if( newt <= 0.98  +0.004 ) return     4.411        *(-1.0); else 
if( newt <= 0.984 +0.004 ) return     4.4056       *(-1.0); else         
if( newt <= 0.988 +0.004 ) return     4.4007       *(-1.0); else 
if( newt <= 0.992 +0.004 ) return     4.3962       *(-1.0); else         
if( newt <= 0.996 +0.004 ) return     4.3924       *(-1.0); else 
if( newt <= 1.0   +0.004 ) return     4.3892       *(-1.0); else 
if( newt <= 1.004 +0.004 ) return     4.3867       *(-1.0); else 
if( newt <= 1.008 +0.004 ) return     4.385        *(-1.0); else 
if( newt <= 1.012 +0.004 ) return     4.3842       *(-1.0); else         
if( newt <= 1.016 +0.004 ) return     4.3843       *(-1.0); else 
if( newt <= 1.02  +0.004 ) return     4.3854       *(-1.0); else 
if( newt <= 1.024 +0.004 ) return     4.3876       *(-1.0); else 
if( newt <= 1.028 +0.004 ) return     4.391        *(-1.0); else 
if( newt <= 1.032 +0.004 ) return     4.3957       *(-1.0); else         
if( newt <= 1.036 +0.004 ) return     4.4019       *(-1.0); else 
if( newt <= 1.04  +0.004 ) return     4.4096       *(-1.0); else 
if( newt <= 1.044 +0.004 ) return     4.4188       *(-1.0); else         
if( newt <= 1.048 +0.004 ) return     4.4296       *(-1.0); else 
if( newt <= 1.052 +0.004 ) return     4.4418       *(-1.0); else         
if( newt <= 1.056 +0.004 ) return     4.4554       *(-1.0); else 
if( newt <= 1.06  +0.004 ) return     4.4699       *(-1.0); else 
if( newt <= 1.064 +0.004 ) return     4.485        *(-1.0); else 
if( newt <= 1.068 +0.004 ) return     4.5002       *(-1.0); else 
if( newt <= 1.072 +0.004 ) return     4.5147       *(-1.0); else 
if( newt <= 1.076 +0.004 ) return     4.5278       *(-1.0); else 
if( newt <= 1.08  +0.004 ) return     4.5386       *(-1.0); else 
if( newt <= 1.084 +0.004 ) return     4.546        *(-1.0); else 
if( newt <= 1.088 +0.004 ) return     4.5492       *(-1.0); else 
if( newt <= 1.092 +0.004 ) return     4.5471       *(-1.0); else         
if( newt <= 1.096 +0.004 ) return     4.5387       *(-1.0); else
  return 0;

}


Real linearCylinderFlux(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
  
    Real dt = ((Real)(((int)(t*1000)) % 40))/1000;

    Real newt1 = t - dt;

    Real newt2 = t + ( 0.004 - dt );

    Real flux1 = cylinderFlux( newt1,0,0,0,0 );

    Real flux2 = cylinderFlux( newt2,0,0,0,0 );

    Real m = ( flux2 - flux1 ) / ( newt2 - newt1 );

    return ( m * ( t - newt1 ) + flux1 );

}


}

