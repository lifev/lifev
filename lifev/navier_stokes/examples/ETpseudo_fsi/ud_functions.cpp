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
Real aortaVelIn(const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real vBar=linearFluxIn(t,x,y,z,i)/7.64808;
    switch (i)
    {
        case 0:
            return vBar*0.0778;
            break;
        case 1:
            return 0.;
            break;
        case 2:
            return vBar*(-0.9969);
            break;
        default:
            ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
            break;
    }

}

Real bypassVelInlet2(const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real vBar=0.8*linearPontdist(t,x,y,z,i)/0.357434;
    switch (i)
    {
        case 0:
            return vBar*0.113589;
            break;
        case 1:
            return -0.213165*vBar;
            break;
        case 2:
            return vBar*(0.970391);
            break;
        default:
            ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
            break;
    }

}

Real bypassVelInlet4(const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real vBar=0.2*linearPontdist(t,x,y,z,i)/0.357434;
    switch (i)
    {
        case 0:
            return vBar*-0.127337;
            break;
        case 1:
            return 0.3022305*vBar;
            break;
        case 2:
            return vBar*(0.9446915);
            break;
        default:
            ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
            break;
    }

}


Real bypassVelInMag(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    Real vBar=linearPontdist(t,0,0,0,0)/0.35;
    return vBar;
}


Real aortaFluxIn(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if(t<=0.00+0.01) return	0.0000e+00;
    if(t<=0.01+0.01) return	0.0000e+00;
    if(t<=0.02+0.01) return	0.0000e+00;
    if(t<=0.03+0.01) return	0.0000e+00;
    if(t<=0.04+0.01) return	-9.1759e-06*1e6;
    if(t<=0.05+0.01) return	-3.0930e-05*1e6;
    if(t<=0.06+0.01) return	-6.2639e-05*1e6;
    if(t<=0.07+0.01) return	-1.0212e-04*1e6;
    if(t<=0.08+0.01) return	-1.4760e-04*1e6;
    if(t<=0.09+0.01) return	-1.9726e-04*1e6;
    if(t<=0.10+0.01) return	-2.4980e-04*1e6;
    if(t<=0.11+0.01) return	-2.9526e-04*1e6;
    if(t<=0.12+0.01) return	-3.2956e-04*1e6;
    if(t<=0.13+0.01) return	-3.5469e-04*1e6;
    if(t<=0.14+0.01) return	-3.7250e-04*1e6;
    if(t<=0.15+0.01) return	-3.8429e-04*1e6;
    if(t<=0.16+0.01) return	-3.9123e-04*1e6;
    if(t<=0.17+0.01) return	-3.9431e-04*1e6;
    if(t<=0.18+0.01) return	-3.9349e-04*1e6;
    if(t<=0.19+0.01) return	-3.8858e-04*1e6;
    if(t<=0.20+0.01) return	-3.7985e-04*1e6;
    if(t<=0.21+0.01) return	-3.6756e-04*1e6;
    if(t<=0.22+0.01) return	-3.5207e-04*1e6;
    if(t<=0.23+0.01) return	-3.3408e-04*1e6;
    if(t<=0.24+0.01) return	-3.1402e-04*1e6;
    if(t<=0.25+0.01) return	-2.9288e-04*1e6;
    if(t<=0.26+0.01) return	-2.7154e-04*1e6;
    if(t<=0.27+0.01) return	-2.5054e-04*1e6;
    if(t<=0.28+0.01) return	-2.2979e-04*1e6;
    if(t<=0.29+0.01) return	-2.0904e-04*1e6;
    if(t<=0.30+0.01) return	-1.8880e-04*1e6;
    if(t<=0.31+0.01) return	-1.6899e-04*1e6;
    if(t<=0.32+0.01) return	-1.4864e-04*1e6;
    if(t<=0.33+0.01) return	-1.2730e-04*1e6;
    if(t<=0.34+0.01) return	-1.0400e-04*1e6;
    if(t<=0.35+0.01) return	-7.9755e-05*1e6;
    if(t<=0.36+0.01) return	-5.8719e-05*1e6;
    if(t<=0.37+0.01) return	-4.0345e-05*1e6;
    if(t<=0.38+0.01) return	-2.4596e-05*1e6;
    if(t<=0.39+0.01) return	-1.2259e-05*1e6;
    if(t<=0.40+0.01) return	-3.8110e-06*1e6;
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

Real aortaPhisPress(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    /*switch(i) {
    case 1:
    return 0.0;
    break;*/
    //  case 2:
    if (t<=0.00+0.01) return -1.e1*11017;
    if (t<=0.01+0.01) return -1.e1*10954;
    if (t<=0.02+0.01) return -1.e1*10893;
    if (t<=0.03+0.01) return -1.e1*10832;
    if (t<=0.04+0.01) return -1.e1*10771;
    if (t<=0.05+0.01) return -1.e1*10712;
    if (t<=0.06+0.01) return -1.e1*10653;
    if (t<=0.07+0.01) return -1.e1*11113;
    if (t<=0.08+0.01) return -1.e1*11544;
    if (t<=0.09+0.01) return -1.e1*11869;
    if (t<=0.10+0.01) return -1.e1*12146;
    if (t<=0.11+0.01) return -1.e1*12394;
    if (t<=0.12+0.01) return -1.e1*12635;
    if (t<=0.13+0.01) return -1.e1*12889;
    if (t<=0.14+0.01) return -1.e1*13151;
    if (t<=0.15+0.01) return -1.e1*13398;
    if (t<=0.16+0.01) return -1.e1*13620;
    if (t<=0.17+0.01) return -1.e1*13833;
    if (t<=0.18+0.01) return -1.e1*14035;
    if (t<=0.19+0.01) return -1.e1*14229;
    if (t<=0.20+0.01) return -1.e1*14436;
    if (t<=0.21+0.01) return -1.e1*14613;
    if (t<=0.22+0.01) return -1.e1*14753;
    if (t<=0.23+0.01) return -1.e1*14878;
    if (t<=0.24+0.01) return -1.e1*14974;
    if (t<=0.25+0.01) return -1.e1*15032;
    if (t<=0.26+0.01) return -1.e1*15047;
    if (t<=0.27+0.01) return -1.e1*15025;
    if (t<=0.28+0.01) return -1.e1*14975;
    if (t<=0.29+0.01) return -1.e1*14899;
    if (t<=0.30+0.01) return -1.e1*14822;
    if (t<=0.31+0.01) return -1.e1*14721;
    if (t<=0.32+0.01) return -1.e1*14594;
    if (t<=0.33+0.01) return -1.e1*14496;
    if (t<=0.34+0.01) return -1.e1*14375;
    if (t<=0.35+0.01) return -1.e1*14198;
    if (t<=0.36+0.01) return -1.e1*13990;
    if (t<=0.37+0.01) return -1.e1*13726;
    if (t<=0.38+0.01) return -1.e1*13397;
    if (t<=0.39+0.01) return -1.e1*13167;
    if (t<=0.40+0.01) return -1.e1*13132;
    if (t<=0.41+0.01) return -1.e1*13315;
    if (t<=0.42+0.01) return -1.e1*13271;
    if (t<=0.43+0.01) return -1.e1*13157;
    if (t<=0.44+0.01) return -1.e1*13028;
    if (t<=0.45+0.01) return -1.e1*12975;
    if (t<=0.46+0.01) return -1.e1*12933;
    if (t<=0.47+0.01) return -1.e1*12891;
    if (t<=0.48+0.01) return -1.e1*12836;
    if (t<=0.49+0.01) return -1.e1*12768;
    if (t<=0.50+0.01) return -1.e1*12700;
    if (t<=0.51+0.01) return -1.e1*12641;
    if (t<=0.52+0.01) return -1.e1*12592;
    if (t<=0.53+0.01) return -1.e1*12548;
    if (t<=0.54+0.01) return -1.e1*12504;
    if (t<=0.55+0.01) return -1.e1*12456;
    if (t<=0.56+0.01) return -1.e1*12405;
    if (t<=0.57+0.01) return -1.e1*12353;
    if (t<=0.58+0.01) return -1.e1*12300;
    if (t<=0.59+0.01) return -1.e1*12244;
    if (t<=0.60+0.01) return -1.e1*12184;
    if (t<=0.61+0.01) return -1.e1*12122;
    if (t<=0.62+0.01) return -1.e1*12058;
    if (t<=0.63+0.01) return -1.e1*11995;
    if (t<=0.64+0.01) return -1.e1*11933;
    if (t<=0.65+0.01) return -1.e1*11871;
    if (t<=0.66+0.01) return -1.e1*11810;
    if (t<=0.67+0.01) return -1.e1*11747;
    if (t<=0.68+0.01) return -1.e1*11684;
    if (t<=0.69+0.01) return -1.e1*11620;
    if (t<=0.70+0.01) return -1.e1*11556;
    if (t<=0.71+0.01) return -1.e1*11492;
    if (t<=0.72+0.01) return -1.e1*11428;
    if (t<=0.73+0.01) return -1.e1*11365;
    if (t<=0.74+0.01) return -1.e1*11302;
    if (t<=0.75+0.01) return -1.e1*11240;
    if (t<=0.76+0.01) return -1.e1*11179;
    if (t<=0.77+0.01) return -1.e1*11120;
    if (t<=0.78+0.01) return -1.e1*11062;
    if (t<=0.79+0.01) return -1.e1*11006;
    //    break;
    /*  case 3:
    return 0.0;
    break;}
    return 0.;*/
}


Real aortaFlux3(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t<=0.00+0.01) return 3.350E-05*1e6;
    if (t<=0.01+0.01) return 3.373E-05*1e6;
    if (t<=0.02+0.01) return 3.402E-05*1e6;
    if (t<=0.03+0.01) return 3.434E-05*1e6;
    if (t<=0.04+0.01) return 3.466E-05*1e6;
    if (t<=0.05+0.01) return 3.495E-05*1e6;
    if (t<=0.06+0.01) return 3.519E-05*1e6;
    if (t<=0.07+0.01) return 3.539E-05*1e6;
    if (t<=0.08+0.01) return 3.564E-05*1e6;
    if (t<=0.09+0.01) return 3.617E-05*1e6;
    if (t<=0.10+0.01) return 3.773E-05*1e6;
    if (t<=0.11+0.01) return 4.176E-05*1e6;
    if (t<=0.12+0.01) return 5.037E-05*1e6;
    if (t<=0.13+0.01) return 6.546E-05*1e6;
    if (t<=0.14+0.01) return 8.701E-05*1e6;
    if (t<=0.15+0.01) return 1.117E-04*1e6;
    if (t<=0.16+0.01) return 1.345E-04*1e6;
    if (t<=0.17+0.01) return 1.519E-04*1e6;
    if (t<=0.18+0.01) return 1.642E-04*1e6;
    if (t<=0.19+0.01) return 1.737E-04*1e6;
    if (t<=0.20+0.01) return 1.821E-04*1e6;
    if (t<=0.21+0.01) return 1.897E-04*1e6;
    if (t<=0.22+0.01) return 1.958E-04*1e6;
    if (t<=0.23+0.01) return 1.999E-04*1e6;
    if (t<=0.24+0.01) return 2.019E-04*1e6;
    if (t<=0.25+0.01) return 2.020E-04*1e6;
    if (t<=0.26+0.01) return 2.004E-04*1e6;
    if (t<=0.27+0.01) return 1.972E-04*1e6;
    if (t<=0.28+0.01) return 1.926E-04*1e6;
    if (t<=0.29+0.01) return 1.868E-04*1e6;
    if (t<=0.30+0.01) return 1.798E-04*1e6;
    if (t<=0.31+0.01) return 1.719E-04*1e6;
    if (t<=0.32+0.01) return 1.632E-04*1e6;
    if (t<=0.33+0.01) return 1.540E-04*1e6;
    if (t<=0.34+0.01) return 1.446E-04*1e6;
    if (t<=0.35+0.01) return 1.350E-04*1e6;
    if (t<=0.36+0.01) return 1.254E-04*1e6;
    if (t<=0.37+0.01) return 1.158E-04*1e6;
    if (t<=0.38+0.01) return 1.062E-04*1e6;
    if (t<=0.39+0.01) return 9.651E-05*1e6;
    if (t<=0.40+0.01) return 8.634E-05*1e6;
    if (t<=0.41+0.01) return 7.558E-05*1e6;
    if (t<=0.42+0.01) return 6.447E-05*1e6;
    if (t<=0.43+0.01) return 5.382E-05*1e6;
    if (t<=0.44+0.01) return 4.484E-05*1e6;
    if (t<=0.45+0.01) return 3.865E-05*1e6;
    if (t<=0.46+0.01) return 3.556E-05*1e6;
    if (t<=0.47+0.01) return 3.473E-05*1e6;
    if (t<=0.48+0.01) return 3.457E-05*1e6;
    if (t<=0.49+0.01) return 3.373E-05*1e6;
    if (t<=0.50+0.01) return 3.191E-05*1e6;
    if (t<=0.51+0.01) return 2.975E-05*1e6;
    if (t<=0.52+0.01) return 2.809E-05*1e6;
    if (t<=0.53+0.01) return 2.730E-05*1e6;
    if (t<=0.54+0.01) return 2.718E-05*1e6;
    if (t<=0.55+0.01) return 2.732E-05*1e6;
    if (t<=0.56+0.01) return 2.744E-05*1e6;
    if (t<=0.57+0.01) return 2.753E-05*1e6;
    if (t<=0.58+0.01) return 2.772E-05*1e6;
    if (t<=0.59+0.01) return 2.811E-05*1e6;
    if (t<=0.60+0.01) return 2.866E-05*1e6;
    if (t<=0.61+0.01) return 2.929E-05*1e6;
    if (t<=0.62+0.01) return 2.990E-05*1e6;
    if (t<=0.63+0.01) return 3.044E-05*1e6;
    if (t<=0.64+0.01) return 3.091E-05*1e6;
    if (t<=0.65+0.01) return 3.132E-05*1e6;
    if (t<=0.66+0.01) return 3.168E-05*1e6;
    if (t<=0.67+0.01) return 3.199E-05*1e6;
    if (t<=0.68+0.01) return 3.224E-05*1e6;
    if (t<=0.69+0.01) return 3.244E-05*1e6;
    if (t<=0.70+0.01) return 3.259E-05*1e6;
    if (t<=0.71+0.01) return 3.270E-05*1e6;
    if (t<=0.72+0.01) return 3.277E-05*1e6;
    if (t<=0.73+0.01) return 3.281E-05*1e6;
    if (t<=0.74+0.01) return 3.282E-05*1e6;
    if (t<=0.75+0.01) return 3.283E-05*1e6;
    if (t<=0.76+0.01) return 3.283E-05*1e6;
    if (t<=0.77+0.01) return 3.286E-05*1e6;
    if (t<=0.78+0.01) return 3.291E-05*1e6;
    if (t<=0.79+0.01) return 3.300E-05*1e6;
}//thoracic aorta,

Real aortaFlux5(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t<=0.00+0.01) return  3.033E-06*1e6;
    if (t<=0.01+0.01) return  3.041E-06*1e6;
    if (t<=0.02+0.01) return  3.062E-06*1e6;
    if (t<=0.03+0.01) return  3.094E-06*1e6;
    if (t<=0.04+0.01) return  3.127E-06*1e6;
    if (t<=0.05+0.01) return  3.150E-06*1e6;
    if (t<=0.06+0.01) return  3.152E-06*1e6;
    if (t<=0.07+0.01) return  3.141E-06*1e6;
    if (t<=0.08+0.01) return  3.196E-06*1e6;
    if (t<=0.09+0.01) return  3.574E-06*1e6;
    if (t<=0.10+0.01) return  4.778E-06*1e6;
    if (t<=0.11+0.01) return  7.387E-06*1e6;
    if (t<=0.12+0.01) return  1.150E-05*1e6;
    if (t<=0.13+0.01) return  1.609E-05*1e6;
    if (t<=0.14+0.01) return  1.933E-05*1e6;
    if (t<=0.15+0.01) return  2.007E-05*1e6;
    if (t<=0.16+0.01) return  1.885E-05*1e6;
    if (t<=0.17+0.01) return  1.706E-05*1e6;
    if (t<=0.18+0.01) return  1.569E-05*1e6;
    if (t<=0.19+0.01) return  1.481E-05*1e6;
    if (t<=0.20+0.01) return  1.401E-05*1e6;
    if (t<=0.21+0.01) return  1.294E-05*1e6;
    if (t<=0.22+0.01) return  1.160E-05*1e6;
    if (t<=0.23+0.01) return  1.018E-05*1e6;
    if (t<=0.24+0.01) return  8.832E-06*1e6;
    if (t<=0.25+0.01) return  7.609E-06*1e6;
    if (t<=0.26+0.01) return  6.578E-06*1e6;
    if (t<=0.27+0.01) return  5.843E-06*1e6;
    if (t<=0.28+0.01) return  5.472E-06*1e6;
    if (t<=0.29+0.01) return  5.412E-06*1e6;
    if (t<=0.30+0.01) return  5.491E-06*1e6;
    if (t<=0.31+0.01) return  5.527E-06*1e6;
    if (t<=0.32+0.01) return  5.420E-06*1e6;
    if (t<=0.33+0.01) return  5.169E-06*1e6;
    if (t<=0.34+0.01) return  4.829E-06*1e6;
    if (t<=0.35+0.01) return  4.465E-06*1e6;
    if (t<=0.36+0.01) return  4.111E-06*1e6;
    if (t<=0.37+0.01) return  3.750E-06*1e6;
    if (t<=0.38+0.01) return  3.304E-06*1e6;
    if (t<=0.39+0.01) return  2.668E-06*1e6;
    if (t<=0.40+0.01) return  1.800E-06*1e6;
    if (t<=0.41+0.01) return  8.269E-07*1e6;
    if (t<=0.42+0.01) return  9.760E-08*1e6;
    if (t<=0.43+0.01) return  7.311E-08*1e6;
    if (t<=0.44+0.01) return  1.041E-06*1e6;
    if (t<=0.45+0.01) return  2.783E-06*1e6;
    if (t<=0.46+0.01) return  4.537E-06*1e6;
    if (t<=0.47+0.01) return  5.488E-06*1e6;
    if (t<=0.48+0.01) return  5.431E-06*1e6;
    if (t<=0.49+0.01) return  4.863E-06*1e6;
    if (t<=0.50+0.01) return  4.452E-06*1e6;
    if (t<=0.51+0.01) return  4.499E-06*1e6;
    if (t<=0.52+0.01) return  4.824E-06*1e6;
    if (t<=0.53+0.01) return  5.059E-06*1e6;
    if (t<=0.54+0.01) return  4.989E-06*1e6;
    if (t<=0.55+0.01) return  4.671E-06*1e6;
    if (t<=0.56+0.01) return  4.292E-06*1e6;
    if (t<=0.57+0.01) return  3.981E-06*1e6;
    if (t<=0.58+0.01) return  3.749E-06*1e6;
    if (t<=0.59+0.01) return  3.553E-06*1e6;
    if (t<=0.60+0.01) return  3.377E-06*1e6;
    if (t<=0.61+0.01) return  3.255E-06*1e6;
    if (t<=0.62+0.01) return  3.224E-06*1e6;
    if (t<=0.63+0.01) return  3.281E-06*1e6;
    if (t<=0.64+0.01) return  3.377E-06*1e6;
    if (t<=0.65+0.01) return  3.452E-06*1e6;
    if (t<=0.66+0.01) return  3.472E-06*1e6;
    if (t<=0.67+0.01) return  3.441E-06*1e6;
    if (t<=0.68+0.01) return  3.389E-06*1e6;
    if (t<=0.69+0.01) return  3.343E-06*1e6;
    if (t<=0.70+0.01) return  3.312E-06*1e6;
    if (t<=0.71+0.01) return  3.289E-06*1e6;
    if (t<=0.72+0.01) return  3.262E-06*1e6;
    if (t<=0.73+0.01) return  3.223E-06*1e6;
    if (t<=0.74+0.01) return  3.177E-06*1e6;
    if (t<=0.75+0.01) return  3.132E-06*1e6;
    if (t<=0.76+0.01) return  3.094E-06*1e6;
    if (t<=0.77+0.01) return  3.065E-06*1e6;
    if (t<=0.78+0.01) return  3.040E-06*1e6;
    if (t<=0.79+0.01) return  3.016E-06*1e6;
}//first branch_1,

Real aortaFlux6(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t<=0.00+0.01) return  7.817E-07*1e6;
    if (t<=0.01+0.01) return  7.879E-07*1e6;
    if (t<=0.02+0.01) return  7.977E-07*1e6;
    if (t<=0.03+0.01) return  8.077E-07*1e6;
    if (t<=0.04+0.01) return  8.144E-07*1e6;
    if (t<=0.05+0.01) return  8.164E-07*1e6;
    if (t<=0.06+0.01) return  8.136E-07*1e6;
    if (t<=0.07+0.01) return  8.129E-07*1e6;
    if (t<=0.08+0.01) return  8.462E-07*1e6;
    if (t<=0.09+0.01) return  9.858E-07*1e6;
    if (t<=0.10+0.01) return  1.317E-06*1e6;
    if (t<=0.11+0.01) return  1.854E-06*1e6;
    if (t<=0.12+0.01) return  2.476E-06*1e6;
    if (t<=0.13+0.01) return  3.005E-06*1e6;
    if (t<=0.14+0.01) return  3.366E-06*1e6;
    if (t<=0.15+0.01) return  3.626E-06*1e6;
    if (t<=0.16+0.01) return  3.860E-06*1e6;
    if (t<=0.17+0.01) return  4.044E-06*1e6;
    if (t<=0.18+0.01) return  4.065E-06*1e6;
    if (t<=0.19+0.01) return  3.824E-06*1e6;
    if (t<=0.20+0.01) return  3.320E-06*1e6;
    if (t<=0.21+0.01) return  2.659E-06*1e6;
    if (t<=0.22+0.01) return  2.006E-06*1e6;
    if (t<=0.23+0.01) return  1.503E-06*1e6;
    if (t<=0.24+0.01) return  1.215E-06*1e6;
    if (t<=0.25+0.01) return  1.117E-06*1e6;
    if (t<=0.26+0.01) return  1.136E-06*1e6;
    if (t<=0.27+0.01) return  1.190E-06*1e6;
    if (t<=0.28+0.01) return  1.222E-06*1e6;
    if (t<=0.29+0.01) return  1.219E-06*1e6;
    if (t<=0.30+0.01) return  1.197E-06*1e6;
    if (t<=0.31+0.01) return  1.179E-06*1e6;
    if (t<=0.32+0.01) return  1.175E-06*1e6;
    if (t<=0.33+0.01) return  1.176E-06*1e6;
    if (t<=0.34+0.01) return  1.164E-06*1e6;
    if (t<=0.35+0.01) return  1.129E-06*1e6;
    if (t<=0.36+0.01) return  1.063E-06*1e6;
    if (t<=0.37+0.01) return  9.647E-07*1e6;
    if (t<=0.38+0.01) return  8.310E-07*1e6;
    if (t<=0.39+0.01) return  6.635E-07*1e6;
    if (t<=0.40+0.01) return  4.774E-07*1e6;
    if (t<=0.41+0.01) return  3.116E-07*1e6;
    if (t<=0.42+0.01) return  2.251E-07*1e6;
    if (t<=0.43+0.01) return  2.627E-07*1e6;
    if (t<=0.44+0.01) return  4.099E-07*1e6;
    if (t<=0.45+0.01) return  5.913E-07*1e6;
    if (t<=0.46+0.01) return  7.359E-07*1e6;
    if (t<=0.47+0.01) return  8.403E-07*1e6;
    if (t<=0.48+0.01) return  9.515E-07*1e6;
    if (t<=0.49+0.01) return  1.097E-06*1e6;
    if (t<=0.50+0.01) return  1.245E-06*1e6;
    if (t<=0.51+0.01) return  1.331E-06*1e6;
    if (t<=0.52+0.01) return  1.314E-06*1e6;
    if (t<=0.53+0.01) return  1.204E-06*1e6;
    if (t<=0.54+0.01) return  1.048E-06*1e6;
    if (t<=0.55+0.01) return  9.004E-07*1e6;
    if (t<=0.56+0.01) return  7.937E-07*1e6;
    if (t<=0.57+0.01) return  7.358E-07*1e6;
    if (t<=0.58+0.01) return  7.163E-07*1e6;
    if (t<=0.59+0.01) return  7.186E-07*1e6;
    if (t<=0.60+0.01) return  7.280E-07*1e6;
    if (t<=0.61+0.01) return  7.359E-07*1e6;
    if (t<=0.62+0.01) return  7.400E-07*1e6;
    if (t<=0.63+0.01) return  7.423E-07*1e6;
    if (t<=0.64+0.01) return  7.465E-07*1e6;
    if (t<=0.65+0.01) return  7.556E-07*1e6;
    if (t<=0.66+0.01) return  7.700E-07*1e6;
    if (t<=0.67+0.01) return  7.871E-07*1e6;
    if (t<=0.68+0.01) return  8.028E-07*1e6;
    if (t<=0.69+0.01) return  8.129E-07*1e6;
    if (t<=0.70+0.01) return  8.151E-07*1e6;
    if (t<=0.71+0.01) return  8.101E-07*1e6;
    if (t<=0.72+0.01) return  8.006E-07*1e6;
    if (t<=0.73+0.01) return  7.905E-07*1e6;
    if (t<=0.74+0.01) return  7.830E-07*1e6;
    if (t<=0.75+0.01) return  7.792E-07*1e6;
    if (t<=0.76+0.01) return  7.782E-07*1e6;
    if (t<=0.77+0.01) return  7.782E-07*1e6;
    if (t<=0.78+0.01) return  7.775E-07*1e6;
    if (t<=0.79+0.01) return  7.754E-07*1e6;
}//branch 1_2 smallest

Real aortaFlux7(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t<=0.00+0.01) return  1.930E-06 *1e6;
    if (t<=0.01+0.01) return  1.710E-06 *1e6;
    if (t<=0.02+0.01) return  1.495E-06 *1e6;
    if (t<=0.03+0.01) return  1.289E-06 *1e6;
    if (t<=0.04+0.01) return  1.096E-06 *1e6;
    if (t<=0.05+0.01) return  9.184E-07 *1e6;
    if (t<=0.06+0.01) return  7.568E-07 *1e6;
    if (t<=0.07+0.01) return  6.240E-07 *1e6;
    if (t<=0.08+0.01) return  5.854E-07 *1e6;
    if (t<=0.09+0.01) return  8.113E-07 *1e6;
    if (t<=0.10+0.01) return  1.559E-06 *1e6;
    if (t<=0.11+0.01) return  3.032E-06 *1e6;
    if (t<=0.12+0.01) return  5.221E-06 *1e6;
    if (t<=0.13+0.01) return  7.903E-06 *1e6;
    if (t<=0.14+0.01) return  1.081E-05 *1e6;
    if (t<=0.15+0.01) return  1.376E-05 *1e6;
    if (t<=0.16+0.01) return  1.661E-05 *1e6;
    if (t<=0.17+0.01) return  1.919E-05 *1e6;
    if (t<=0.18+0.01) return  2.134E-05 *1e6;
    if (t<=0.19+0.01) return  2.292E-05 *1e6;
    if (t<=0.20+0.01) return  2.389E-05 *1e6;
    if (t<=0.21+0.01) return  2.433E-05 *1e6;
    if (t<=0.22+0.01) return  2.435E-05 *1e6;
    if (t<=0.23+0.01) return  2.406E-05 *1e6;
    if (t<=0.24+0.01) return  2.353E-05 *1e6;
    if (t<=0.25+0.01) return  2.280E-05 *1e6;
    if (t<=0.26+0.01) return  2.188E-05 *1e6;
    if (t<=0.27+0.01) return  2.078E-05 *1e6;
    if (t<=0.28+0.01) return  1.947E-05 *1e6;
    if (t<=0.29+0.01) return  1.797E-05 *1e6;
    if (t<=0.30+0.01) return  1.628E-05 *1e6;
    if (t<=0.31+0.01) return  1.445E-05 *1e6;
    if (t<=0.32+0.01) return  1.254E-05 *1e6;
    if (t<=0.33+0.01) return  1.060E-05 *1e6;
    if (t<=0.34+0.01) return  8.684E-06 *1e6;
    if (t<=0.35+0.01) return  6.838E-06 *1e6;
    if (t<=0.36+0.01) return  5.084E-06 *1e6;
    if (t<=0.37+0.01) return  3.412E-06 *1e6;
    if (t<=0.38+0.01) return  1.784E-06 *1e6;
    if (t<=0.39+0.01) return  1.534E-07 *1e6;
    if (t<=0.40+0.01) return  -1.494E-06*1e6;
    if (t<=0.41+0.01) return  -3.093E-06*1e6;
    if (t<=0.42+0.01) return  -4.495E-06*1e6;
    if (t<=0.43+0.01) return  -5.521E-06*1e6;
    if (t<=0.44+0.01) return  -6.086E-06*1e6;
    if (t<=0.45+0.01) return  -6.252E-06*1e6;
    if (t<=0.46+0.01) return  -6.160E-06*1e6;
    if (t<=0.47+0.01) return  -5.908E-06*1e6;
    if (t<=0.48+0.01) return  -5.508E-06*1e6;
    if (t<=0.49+0.01) return  -4.944E-06*1e6;
    if (t<=0.50+0.01) return  -4.234E-06*1e6;
    if (t<=0.51+0.01) return  -3.441E-06*1e6;
    if (t<=0.52+0.01) return  -2.639E-06*1e6;
    if (t<=0.53+0.01) return  -1.873E-06*1e6;
    if (t<=0.54+0.01) return  -1.160E-06*1e6;
    if (t<=0.55+0.01) return  -5.019E-07*1e6;
    if (t<=0.56+0.01) return  1.024E-07 *1e6;
    if (t<=0.57+0.01) return  6.543E-07 *1e6;
    if (t<=0.58+0.01) return  1.156E-06 *1e6;
    if (t<=0.59+0.01) return  1.614E-06 *1e6;
    if (t<=0.60+0.01) return  2.034E-06 *1e6;
    if (t<=0.61+0.01) return  2.419E-06 *1e6;
    if (t<=0.62+0.01) return  2.765E-06 *1e6;
    if (t<=0.63+0.01) return  3.065E-06 *1e6;
    if (t<=0.64+0.01) return  3.307E-06 *1e6;
    if (t<=0.65+0.01) return  3.486E-06 *1e6;
    if (t<=0.66+0.01) return  3.600E-06 *1e6;
    if (t<=0.67+0.01) return  3.654E-06 *1e6;
    if (t<=0.68+0.01) return  3.657E-06 *1e6;
    if (t<=0.69+0.01) return  3.621E-06 *1e6;
    if (t<=0.70+0.01) return  3.554E-06 *1e6;
    if (t<=0.71+0.01) return  3.465E-06 *1e6;
    if (t<=0.72+0.01) return  3.357E-06 *1e6;
    if (t<=0.73+0.01) return  3.233E-06 *1e6;
    if (t<=0.74+0.01) return  3.094E-06 *1e6;
    if (t<=0.75+0.01) return  2.941E-06 *1e6;
    if (t<=0.76+0.01) return  2.774E-06 *1e6;
    if (t<=0.77+0.01) return  2.591E-06 *1e6;
    if (t<=0.78+0.01) return  2.395E-06 *1e6;
    if (t<=0.79+0.01) return  2.185E-06 *1e6;
}//R. Brachia, branch 1_3

Real aortaFlux8(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t<=0.00+0.01) return  2.445E-06*1e6;
    if (t<=0.01+0.01) return  2.470E-06*1e6;
    if (t<=0.02+0.01) return  2.519E-06*1e6;
    if (t<=0.03+0.01) return  2.578E-06*1e6;
    if (t<=0.04+0.01) return  2.631E-06*1e6;
    if (t<=0.05+0.01) return  2.665E-06*1e6;
    if (t<=0.06+0.01) return  2.674E-06*1e6;
    if (t<=0.07+0.01) return  2.674E-06*1e6;
    if (t<=0.08+0.01) return  2.777E-06*1e6;
    if (t<=0.09+0.01) return  3.260E-06*1e6;
    if (t<=0.10+0.01) return  4.539E-06*1e6;
    if (t<=0.11+0.01) return  6.956E-06*1e6;
    if (t<=0.12+0.01) return  1.045E-05*1e6;
    if (t<=0.13+0.01) return  1.438E-05*1e6;
    if (t<=0.14+0.01) return  1.769E-05*1e6;
    if (t<=0.15+0.01) return  1.958E-05*1e6;
    if (t<=0.16+0.01) return  1.991E-05*1e6;
    if (t<=0.17+0.01) return  1.910E-05*1e6;
    if (t<=0.18+0.01) return  1.775E-05*1e6;
    if (t<=0.19+0.01) return  1.628E-05*1e6;
    if (t<=0.20+0.01) return  1.479E-05*1e6;
    if (t<=0.21+0.01) return  1.326E-05*1e6;
    if (t<=0.22+0.01) return  1.165E-05*1e6;
    if (t<=0.23+0.01) return  1.001E-05*1e6;
    if (t<=0.24+0.01) return  8.448E-06*1e6;
    if (t<=0.25+0.01) return  7.065E-06*1e6;
    if (t<=0.26+0.01) return  5.930E-06*1e6;
    if (t<=0.27+0.01) return  5.089E-06*1e6;
    if (t<=0.28+0.01) return  4.566E-06*1e6;
    if (t<=0.29+0.01) return  4.350E-06*1e6;
    if (t<=0.30+0.01) return  4.382E-06*1e6;
    if (t<=0.31+0.01) return  4.552E-06*1e6;
    if (t<=0.32+0.01) return  4.729E-06*1e6;
    if (t<=0.33+0.01) return  4.798E-06*1e6;
    if (t<=0.34+0.01) return  4.698E-06*1e6;
    if (t<=0.35+0.01) return  4.428E-06*1e6;
    if (t<=0.36+0.01) return  4.026E-06*1e6;
    if (t<=0.37+0.01) return  3.513E-06*1e6;
    if (t<=0.38+0.01) return  2.875E-06*1e6;
    if (t<=0.39+0.01) return  2.076E-06*1e6;
    if (t<=0.40+0.01) return  1.131E-06*1e6;
    if (t<=0.41+0.01) return  1.861E-07*1e6;
    if (t<=0.42+0.01) return  -4.577E-07*1e6;
    if (t<=0.43+0.01) return  -4.657E-07*1e6;
    if (t<=0.44+0.01) return  3.124E-07*1e6;
    if (t<=0.45+0.01) return  1.684E-06*1e6;
    if (t<=0.46+0.01) return  3.174E-06*1e6;
    if (t<=0.47+0.01) return  4.306E-06*1e6;
    if (t<=0.48+0.01) return  4.873E-06*1e6;
    if (t<=0.49+0.01) return  4.980E-06*1e6;
    if (t<=0.50+0.01) return  4.876E-06*1e6;
    if (t<=0.51+0.01) return  4.757E-06*1e6;
    if (t<=0.52+0.01) return  4.681E-06*1e6;
    if (t<=0.53+0.01) return  4.603E-06*1e6;
    if (t<=0.54+0.01) return  4.458E-06*1e6;
    if (t<=0.55+0.01) return  4.224E-06*1e6;
    if (t<=0.56+0.01) return  3.928E-06*1e6;
    if (t<=0.57+0.01) return  3.618E-06*1e6;
    if (t<=0.58+0.01) return  3.335E-06*1e6;
    if (t<=0.59+0.01) return  3.098E-06*1e6;
    if (t<=0.60+0.01) return  2.913E-06*1e6;
    if (t<=0.61+0.01) return  2.778E-06*1e6;
    if (t<=0.62+0.01) return  2.697E-06*1e6;
    if (t<=0.63+0.01) return  2.670E-06*1e6;
    if (t<=0.64+0.01) return  2.692E-06*1e6;
    if (t<=0.65+0.01) return  2.746E-06*1e6;
    if (t<=0.66+0.01) return  2.811E-06*1e6;
    if (t<=0.67+0.01) return  2.868E-06*1e6;
    if (t<=0.68+0.01) return  2.904E-06*1e6;
    if (t<=0.69+0.01) return  2.915E-06*1e6;
    if (t<=0.70+0.01) return  2.903E-06*1e6;
    if (t<=0.71+0.01) return  2.872E-06*1e6;
    if (t<=0.72+0.01) return  2.826E-06*1e6;
    if (t<=0.73+0.01) return  2.769E-06*1e6;
    if (t<=0.74+0.01) return  2.704E-06*1e6;
    if (t<=0.75+0.01) return  2.636E-06*1e6;
    if (t<=0.76+0.01) return  2.569E-06*1e6;
    if (t<=0.77+0.01) return  2.511E-06*1e6;
    if (t<=0.78+0.01) return  2.465E-06*1e6;
    if (t<=0.79+0.01) return  2.433E-06*1e6;
}// 15, LCCA, branch 2

Real aortaFlux9(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t<=0.00+0.01) return  7.854E-07*1e6;
    if (t<=0.01+0.01) return  7.900E-07*1e6;
    if (t<=0.02+0.01) return  7.988E-07*1e6;
    if (t<=0.03+0.01) return  8.077E-07*1e6;
    if (t<=0.04+0.01) return  8.131E-07*1e6;
    if (t<=0.05+0.01) return  8.141E-07*1e6;
    if (t<=0.06+0.01) return  8.120E-07*1e6;
    if (t<=0.07+0.01) return  8.085E-07*1e6;
    if (t<=0.08+0.01) return  8.082E-07*1e6;
    if (t<=0.09+0.01) return  8.255E-07*1e6;
    if (t<=0.10+0.01) return  8.927E-07*1e6;
    if (t<=0.11+0.01) return  1.055E-06*1e6;
    if (t<=0.12+0.01) return  1.343E-06*1e6;
    if (t<=0.13+0.01) return  1.728E-06*1e6;
    if (t<=0.14+0.01) return  2.118E-06*1e6;
    if (t<=0.15+0.01) return  2.419E-06*1e6;
    if (t<=0.16+0.01) return  2.598E-06*1e6;
    if (t<=0.17+0.01) return  2.673E-06*1e6;
    if (t<=0.18+0.01) return  2.642E-06*1e6;
    if (t<=0.19+0.01) return  2.465E-06*1e6;
    if (t<=0.20+0.01) return  2.118E-06*1e6;
    if (t<=0.21+0.01) return  1.660E-06*1e6;
    if (t<=0.22+0.01) return  1.210E-06*1e6;
    if (t<=0.23+0.01) return  8.723E-07*1e6;
    if (t<=0.24+0.01) return  6.812E-07*1e6;
    if (t<=0.25+0.01) return  6.155E-07*1e6;
    if (t<=0.26+0.01) return  6.323E-07*1e6;
    if (t<=0.27+0.01) return  6.927E-07*1e6;
    if (t<=0.28+0.01) return  7.708E-07*1e6;
    if (t<=0.29+0.01) return  8.538E-07*1e6;
    if (t<=0.30+0.01) return  9.397E-07*1e6;
    if (t<=0.31+0.01) return  1.030E-06*1e6;
    if (t<=0.32+0.01) return  1.120E-06*1e6;
    if (t<=0.33+0.01) return  1.203E-06*1e6;
    if (t<=0.34+0.01) return  1.263E-06*1e6;
    if (t<=0.35+0.01) return  1.292E-06*1e6;
    if (t<=0.36+0.01) return  1.288E-06*1e6;
    if (t<=0.37+0.01) return  1.258E-06*1e6;
    if (t<=0.38+0.01) return  1.207E-06*1e6;
    if (t<=0.39+0.01) return  1.142E-06*1e6;
    if (t<=0.40+0.01) return  1.064E-06*1e6;
    if (t<=0.41+0.01) return  9.816E-07*1e6;
    if (t<=0.42+0.01) return  9.133E-07*1e6;
    if (t<=0.43+0.01) return  8.870E-07*1e6;
    if (t<=0.44+0.01) return  9.268E-07*1e6;
    if (t<=0.45+0.01) return  1.035E-06*1e6;
    if (t<=0.46+0.01) return  1.183E-06*1e6;
    if (t<=0.47+0.01) return  1.329E-06*1e6;
    if (t<=0.48+0.01) return  1.445E-06*1e6;
    if (t<=0.49+0.01) return  1.529E-06*1e6;
    if (t<=0.50+0.01) return  1.587E-06*1e6;
    if (t<=0.51+0.01) return  1.612E-06*1e6;
    if (t<=0.52+0.01) return  1.587E-06*1e6;
    if (t<=0.53+0.01) return  1.501E-06*1e6;
    if (t<=0.54+0.01) return  1.370E-06*1e6;
    if (t<=0.55+0.01) return  1.230E-06*1e6;
    if (t<=0.56+0.01) return  1.112E-06*1e6;
    if (t<=0.57+0.01) return  1.029E-06*1e6;
    if (t<=0.58+0.01) return  9.776E-07*1e6;
    if (t<=0.59+0.01) return  9.460E-07*1e6;
    if (t<=0.60+0.01) return  9.239E-07*1e6;
    if (t<=0.61+0.01) return  9.061E-07*1e6;
    if (t<=0.62+0.01) return  8.917E-07*1e6;
    if (t<=0.63+0.01) return  8.821E-07*1e6;
    if (t<=0.64+0.01) return  8.786E-07*1e6;
    if (t<=0.65+0.01) return  8.809E-07*1e6;
    if (t<=0.66+0.01) return  8.871E-07*1e6;
    if (t<=0.67+0.01) return  8.940E-07*1e6;
    if (t<=0.68+0.01) return  8.978E-07*1e6;
    if (t<=0.69+0.01) return  8.959E-07*1e6;
    if (t<=0.70+0.01) return  8.874E-07*1e6;
    if (t<=0.71+0.01) return  8.734E-07*1e6;
    if (t<=0.72+0.01) return  8.561E-07*1e6;
    if (t<=0.73+0.01) return  8.383E-07*1e6;
    if (t<=0.74+0.01) return  8.221E-07*1e6;
    if (t<=0.75+0.01) return  8.088E-07*1e6;
    if (t<=0.76+0.01) return  7.986E-07*1e6;
    if (t<=0.77+0.01) return  7.909E-07*1e6;
    if (t<=0.78+0.01) return  7.852E-07*1e6;
    if (t<=0.79+0.01) return  7.807E-07*1e6;
}// 20 LVA branch 3_1

Real aortaFlux4(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t<=0.00+0.01) return  2.581E-06 *1e6;
    if (t<=0.01+0.01) return  2.499E-06 *1e6;
    if (t<=0.02+0.01) return  2.399E-06 *1e6;
    if (t<=0.03+0.01) return  2.281E-06 *1e6;
    if (t<=0.04+0.01) return  2.147E-06 *1e6;
    if (t<=0.05+0.01) return  2.002E-06 *1e6;
    if (t<=0.06+0.01) return  1.849E-06 *1e6;
    if (t<=0.07+0.01) return  1.695E-06 *1e6;
    if (t<=0.08+0.01) return  1.565E-06 *1e6;
    if (t<=0.09+0.01) return  1.537E-06 *1e6;
    if (t<=0.10+0.01) return  1.766E-06 *1e6;
    if (t<=0.11+0.01) return  2.460E-06 *1e6;
    if (t<=0.12+0.01) return  3.773E-06 *1e6;
    if (t<=0.13+0.01) return  5.709E-06 *1e6;
    if (t<=0.14+0.01) return  8.131E-06 *1e6;
    if (t<=0.15+0.01) return  1.087E-05 *1e6;
    if (t<=0.16+0.01) return  1.379E-05 *1e6;
    if (t<=0.17+0.01) return  1.675E-05 *1e6;
    if (t<=0.18+0.01) return  1.957E-05 *1e6;
    if (t<=0.19+0.01) return  2.205E-05 *1e6;
    if (t<=0.20+0.01) return  2.403E-05 *1e6;
    if (t<=0.21+0.01) return  2.538E-05 *1e6;
    if (t<=0.22+0.01) return  2.606E-05 *1e6;
    if (t<=0.23+0.01) return  2.611E-05 *1e6;
    if (t<=0.24+0.01) return  2.562E-05 *1e6;
    if (t<=0.25+0.01) return  2.474E-05 *1e6;
    if (t<=0.26+0.01) return  2.362E-05 *1e6;
    if (t<=0.27+0.01) return  2.238E-05 *1e6;
    if (t<=0.28+0.01) return  2.111E-05 *1e6;
    if (t<=0.29+0.01) return  1.981E-05 *1e6;
    if (t<=0.30+0.01) return  1.850E-05 *1e6;
    if (t<=0.31+0.01) return  1.715E-05 *1e6;
    if (t<=0.32+0.01) return  1.576E-05 *1e6;
    if (t<=0.33+0.01) return  1.432E-05 *1e6;
    if (t<=0.34+0.01) return  1.284E-05 *1e6;
    if (t<=0.35+0.01) return  1.132E-05 *1e6;
    if (t<=0.36+0.01) return  9.768E-06 *1e6;
    if (t<=0.37+0.01) return  8.180E-06 *1e6;
    if (t<=0.38+0.01) return  6.543E-06 *1e6;
    if (t<=0.39+0.01) return  4.831E-06 *1e6;
    if (t<=0.40+0.01) return  3.030E-06 *1e6;
    if (t<=0.41+0.01) return  1.163E-06 *1e6;
    if (t<=0.42+0.01) return  -6.817E-07*1e6;
    if (t<=0.43+0.01) return  -2.362E-06*1e6;
    if (t<=0.44+0.01) return  -3.738E-06*1e6;
    if (t<=0.45+0.01) return  -4.742E-06*1e6;
    if (t<=0.46+0.01) return  -5.400E-06*1e6;
    if (t<=0.47+0.01) return  -5.785E-06*1e6;
    if (t<=0.48+0.01) return  -5.956E-06*1e6;
    if (t<=0.49+0.01) return  -5.932E-06*1e6;
    if (t<=0.50+0.01) return  -5.723E-06*1e6;
    if (t<=0.51+0.01) return  -5.358E-06*1e6;
    if (t<=0.52+0.01) return  -4.889E-06*1e6;
    if (t<=0.53+0.01) return  -4.370E-06*1e6;
    if (t<=0.54+0.01) return  -3.846E-06*1e6;
    if (t<=0.55+0.01) return  -3.341E-06*1e6;
    if (t<=0.56+0.01) return  -2.866E-06*1e6;
    if (t<=0.57+0.01) return  -2.422E-06*1e6;
    if (t<=0.58+0.01) return  -2.004E-06*1e6;
    if (t<=0.59+0.01) return  -1.601E-06*1e6;
    if (t<=0.60+0.01) return  -1.206E-06*1e6;
    if (t<=0.61+0.01) return  -8.123E-07*1e6;
    if (t<=0.62+0.01) return  -4.195E-07*1e6;
    if (t<=0.63+0.01) return  -3.198E-08*1e6;
    if (t<=0.64+0.01) return  3.429E-07 *1e6;
    if (t<=0.65+0.01) return  6.970E-07 *1e6;
    if (t<=0.66+0.01) return  1.024E-06 *1e6;
    if (t<=0.67+0.01) return  1.320E-06 *1e6;
    if (t<=0.68+0.01) return  1.583E-06 *1e6;
    if (t<=0.69+0.01) return  1.813E-06 *1e6;
    if (t<=0.70+0.01) return  2.011E-06 *1e6;
    if (t<=0.71+0.01) return  2.180E-06 *1e6;
    if (t<=0.72+0.01) return  2.323E-06 *1e6;
    if (t<=0.73+0.01) return  2.442E-06 *1e6;
    if (t<=0.74+0.01) return  2.539E-06 *1e6;
    if (t<=0.75+0.01) return  2.613E-06 *1e6;
    if (t<=0.76+0.01) return  2.665E-06 *1e6;
    if (t<=0.77+0.01) return  2.691E-06 *1e6;
    if (t<=0.78+0.01) return  2.690E-06 *1e6;
    if (t<=0.79+0.01) return  2.661E-06 *1e6;
}//21, L. Brachia, bhanch 3_2







Real f(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 30.;
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


Real E(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    if(z<5 && z>=0)
        return  60*1e3;// 1.0 ;//-29;//*e-1*5; // about [(110-60)*(133.332*10)]/[10*(2.08-1.8)] ( 10 because of mm, 133... because of mmHg)
    //    (see paper by Liu, Dang, etc). in their plot x and y are inverted.
    if(z>5 && z<6)
        return  60*1e3+30*(z-5)*1e3; //20.0*(z-5) + 10.0;//
    if(z>=6)
        return 90*1e3;  //30.0; //
    // if ( x < -4.0 && x >= -4.5 )
    //     return 1.0 - 4.0 * (x + 4.0);
    // if (  x < -4.5 )
    //     return 3.0;

    //    return 1.0;
    if ( z<0 && z>=-3 )
        return 60*1e3;  //(10.0/3.0)*(-z) + 10.0; //60*1e3; ////60*1e3;
    if(z<-3)
        return 60*1e3; //20.0; //60*1e3; //10.0;//60*1e3;
}


// Initial displacement and velocity
Real d0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
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
        ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
        break;
    }
}

Real w0(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch (i)
    {
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



Real linearPress2( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100.0)/100.0;
    Real tii=ti+0.01;
    return (aortaPhisPress(tii)-aortaPhisPress(ti))/(0.01)*(t-(ti))+aortaPhisPress(ti) + 115000.0;
}


Real linearFlux3( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000.0))%800)/1000.0;

    Real ti=floor(t*100.0)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux3(tii)-aortaFlux3(ti))/(0.01)*(t-(ti))+aortaFlux3(ti);
}


Real linearFluxIn(Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100.0)/100.0;
    Real tii=ti+0.01;

    return (aortaFluxIn(tii)-aortaFluxIn(ti))/(0.01)*(t-(ti))+aortaFluxIn(ti);
}


Real linearFlux4(Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100.0)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux4(tii)-aortaFlux4(ti))/(0.01)*(t-(ti))+aortaFlux4(ti);
}



Real linearFlux5( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux5(tii)-aortaFlux5(ti))/(0.01)*(t-(ti))+aortaFlux5(ti);
}



Real linearFlux6( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux6(tii)-aortaFlux6(ti))/(0.01)*(t-(ti))+aortaFlux6(ti);
}


Real linearFlux7( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux7(tii)-aortaFlux7(ti))/(0.01)*(t-(ti))+aortaFlux7(ti);
}


Real linearFlux8(Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux8(tii)-aortaFlux8(ti))/(0.01)*(t-(ti))+aortaFlux8(ti);
}

Real linearFlux9( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux9(tii)-aortaFlux9(ti))/(0.01)*(t-(ti))+aortaFlux9(ti);
}


Real linearFlux3_(Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux3_(tii)-aortaFlux3_(ti))/(0.01)*(t-(ti))+aortaFlux3_(ti);
}


Real linearFlux6_( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t>0.8)
        t=(((int)floor(t*1000))%800)/1000.0;

    Real ti=floor(t*100)/100.0;
    Real tii=ti+0.01;
    return (aortaFlux6_(tii)-aortaFlux6_(ti))/(0.01)*(t-(ti))+aortaFlux6_(ti);
}

Real aortaFlux3_(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t<=0.00+0.01) return  33.6             ;
    if(t<=0.01+0.01) return  33.83            ;
    if(t<=0.02+0.01) return  34.12            ;
    if(t<=0.03+0.01) return  34.44            ;
    if(t<=0.04+0.01) return  34.76            ;
    if(t<=0.05+0.01) return  35.05            ;
    if(t<=0.06+0.01) return  35.29            ;
    if(t<=0.07+0.01) return  35.49            ;
    if(t<=0.08+0.01) return  35.74            ;
    if(t<=0.09+0.01) return  36.27            ;
    if(t<=0.10+0.01) return  37.83            ;
    if(t<=0.11+0.01) return  41.86            ;
    if(t<=0.12+0.01) return  50.47            ;
    if(t<=0.13+0.01) return  65.55999999999999;
    if(t<=0.14+0.01) return  87.10999999999999;
    if(t<=0.15+0.01) return  111.8            ;
    if(t<=0.16+0.01) return  134.6            ;
    if(t<=0.17+0.01) return  152              ;
    if(t<=0.18+0.01) return  164.3            ;
    if(t<=0.19+0.01) return  173.8            ;
    if(t<=0.20+0.01) return  182.2            ;
    if(t<=0.21+0.01) return  189.8            ;
    if(t<=0.22+0.01) return  195.9            ;
    if(t<=0.23+0.01) return  200              ;
    if(t<=0.24+0.01) return  202              ;
    if(t<=0.25+0.01) return  202.1            ;
    if(t<=0.26+0.01) return  200.5            ;
    if(t<=0.27+0.01) return  197.3            ;
    if(t<=0.28+0.01) return  192.7            ;
    if(t<=0.29+0.01) return  186.9            ;
    if(t<=0.30+0.01) return  179.9            ;
    if(t<=0.31+0.01) return  172              ;
    if(t<=0.32+0.01) return  163.3            ;
    if(t<=0.33+0.01) return  154.1            ;
    if(t<=0.34+0.01) return  144.7            ;
    if(t<=0.35+0.01) return  135.1            ;
    if(t<=0.36+0.01) return  125.5            ;
    if(t<=0.37+0.01) return  115.9            ;
    if(t<=0.38+0.01) return  106.3            ;
    if(t<=0.39+0.01) return  96.60999999999999;
    if(t<=0.40+0.01) return  86.44            ;
    if(t<=0.41+0.01) return  75.67999999999999;
    if(t<=0.42+0.01) return  64.56999999999999;
    if(t<=0.43+0.01) return  53.92            ;
    if(t<=0.44+0.01) return  44.94            ;
    if(t<=0.45+0.01) return  38.75            ;
    if(t<=0.46+0.01) return  35.66            ;
    if(t<=0.47+0.01) return  34.83000000000001;
    if(t<=0.48+0.01) return  34.67            ;
    if(t<=0.49+0.01) return  33.83            ;
    if(t<=0.50+0.01) return  32.01            ;
    if(t<=0.51+0.01) return  29.85            ;
    if(t<=0.52+0.01) return  28.19            ;
    if(t<=0.53+0.01) return  27.4             ;
    if(t<=0.54+0.01) return  27.28            ;
    if(t<=0.55+0.01) return  27.42            ;
    if(t<=0.56+0.01) return  27.54            ;
    if(t<=0.57+0.01) return  27.63            ;
    if(t<=0.58+0.01) return  27.82            ;
    if(t<=0.59+0.01) return  28.21            ;
    if(t<=0.60+0.01) return  28.76            ;
    if(t<=0.61+0.01) return  29.39            ;
    if(t<=0.62+0.01) return  30               ;
    if(t<=0.63+0.01) return  30.54            ;
    if(t<=0.64+0.01) return  31.01            ;
    if(t<=0.65+0.01) return  31.42            ;
    if(t<=0.66+0.01) return  31.78            ;
    if(t<=0.67+0.01) return  32.09            ;
    if(t<=0.68+0.01) return  32.34            ;
    if(t<=0.69+0.01) return  32.54            ;
    if(t<=0.70+0.01) return  32.69            ;
    if(t<=0.71+0.01) return  32.8             ;
    if(t<=0.72+0.01) return  32.87            ;
    if(t<=0.73+0.01) return  32.91            ;
    if(t<=0.74+0.01) return  32.92            ;
    if(t<=0.75+0.01) return  32.93000000000001;
    if(t<=0.76+0.01) return  32.93000000000001;
    if(t<=0.77+0.01) return  32.96            ;
    if(t<=0.78+0.01) return  33.01000000000001;
    if(t<=0.79+0.01) return  33.1             ;
}

Real aortaFlux6_(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    if(t<=0.00+0.01) return    0.6737071250000001;
    if(t<=0.01+0.01) return    0.6799071250000001;
    if(t<=0.02+0.01) return    0.6897071250000001;
    if(t<=0.03+0.01) return    0.6997071250000001;
    if(t<=0.04+0.01) return    0.706407125       ;
    if(t<=0.05+0.01) return    0.7084071250000001;
    if(t<=0.06+0.01) return    0.7056071250000001;
    if(t<=0.07+0.01) return    0.7049071250000001;
    if(t<=0.08+0.01) return    0.7382071250000001;
    if(t<=0.09+0.01) return    0.8778071249999999;
    if(t<=0.10+0.01) return    1.209007125       ;
    if(t<=0.11+0.01) return    1.746007125       ;
    if(t<=0.12+0.01) return    2.368007125       ;
    if(t<=0.13+0.01) return    2.897007125       ;
    if(t<=0.14+0.01) return    3.258007125       ;
    if(t<=0.15+0.01) return    3.518007125       ;
    if(t<=0.16+0.01) return    3.752007125       ;
    if(t<=0.17+0.01) return    3.936007124999999 ;
    if(t<=0.18+0.01) return    3.957007125       ;
    if(t<=0.19+0.01) return    3.716007125       ;
    if(t<=0.20+0.01) return    3.212007125       ;
    if(t<=0.21+0.01) return    2.551007125       ;
    if(t<=0.22+0.01) return    1.898007125       ;
    if(t<=0.23+0.01) return    1.395007125       ;
    if(t<=0.24+0.01) return    1.107007125       ;
    if(t<=0.25+0.01) return    1.009007125       ;
    if(t<=0.26+0.01) return    1.028007125       ;
    if(t<=0.27+0.01) return    1.082007125       ;
    if(t<=0.28+0.01) return    1.114007125       ;
    if(t<=0.29+0.01) return    1.111007125       ;
    if(t<=0.30+0.01) return    1.089007125       ;
    if(t<=0.31+0.01) return    1.071007125       ;
    if(t<=0.32+0.01) return    1.067007125       ;
    if(t<=0.33+0.01) return    1.068007125       ;
    if(t<=0.34+0.01) return    1.056007125       ;
    if(t<=0.35+0.01) return    1.021007125       ;
    if(t<=0.36+0.01) return    0.955007125       ;
    if(t<=0.37+0.01) return    0.856707125       ;
    if(t<=0.38+0.01) return    0.7230071250000001;
    if(t<=0.39+0.01) return    0.5555071250000001;
    if(t<=0.40+0.01) return    0.3694071249999999;
    if(t<=0.41+0.01) return    0.203607125       ;
    if(t<=0.42+0.01) return    0.117107125       ;
    if(t<=0.43+0.01) return    0.154707125       ;
    if(t<=0.44+0.01) return    0.3019071249999999;
    if(t<=0.45+0.01) return    0.4833071250000001;
    if(t<=0.46+0.01) return    0.6279071250000001;
    if(t<=0.47+0.01) return    0.7323071250000001;
    if(t<=0.48+0.01) return    0.8435071250000001;
    if(t<=0.49+0.01) return    0.989007125       ;
    if(t<=0.50+0.01) return    1.137007125       ;
    if(t<=0.51+0.01) return    1.223007125       ;
    if(t<=0.52+0.01) return    1.206007125       ;
    if(t<=0.53+0.01) return    1.096007125       ;
    if(t<=0.54+0.01) return    0.9400071249999998;
    if(t<=0.55+0.01) return    0.7924071250000001;
    if(t<=0.56+0.01) return    0.6857071250000001;
    if(t<=0.57+0.01) return    0.6278071250000001;
    if(t<=0.58+0.01) return    0.6083071250000001;
    if(t<=0.59+0.01) return    0.610607125       ;
    if(t<=0.60+0.01) return    0.620007125       ;
    if(t<=0.61+0.01) return    0.6279071250000001;
    if(t<=0.62+0.01) return    0.6320071250000001;
    if(t<=0.63+0.01) return    0.6343071250000001;
    if(t<=0.64+0.01) return    0.6385071250000001;
    if(t<=0.65+0.01) return    0.6476071250000001;
    if(t<=0.66+0.01) return    0.6620071250000001;
    if(t<=0.67+0.01) return    0.679107125       ;
    if(t<=0.68+0.01) return    0.6948071250000001;
    if(t<=0.69+0.01) return    0.7049071250000001;
    if(t<=0.70+0.01) return    0.7071071250000001;
    if(t<=0.71+0.01) return    0.7021071250000001;
    if(t<=0.72+0.01) return    0.6926071250000001;
    if(t<=0.73+0.01) return    0.6825071250000001;
    if(t<=0.74+0.01) return    0.675007125       ;
    if(t<=0.75+0.01) return    0.671207125       ;
    if(t<=0.76+0.01) return    0.6702071250000001;
    if(t<=0.77+0.01) return    0.6702071250000001;
    if(t<=0.78+0.01) return    0.669507125       ;
    if(t<=0.79+0.01) return    0.6674071250000001;
}

Real u2(Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -1e4;
}

Real inletCylinder(Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)  {
if(t<=       0+0.005) return 4.5724       ;
if(t<=  0.0050+0.005) return    4.5402    ;
if(t<=  0.0100+0.005) return    4.5010    ;
if(t<=  0.0150+0.005) return    4.4491    ;
if(t<=  0.0200+0.005) return    4.3803    ;
if(t<=  0.0250+0.005) return    4.2893    ;
if(t<=  0.0300+0.005) return    4.1854    ;
if(t<=  0.0350+0.005) return    4.0762    ;
if(t<=  0.0400+0.005) return    3.9711    ;
if(t<=  0.0450+0.005) return    3.8780    ;
if(t<=  0.0500+0.005) return    3.8102    ;
if(t<=  0.0550+0.005) return    3.7778    ;
if(t<=  0.0600+0.005) return    3.7894    ;
if(t<=  0.0650+0.005) return    3.8468    ;
if(t<=  0.0700+0.005) return    3.9632    ;
if(t<=  0.0750+0.005) return    4.1447    ;
if(t<=  0.0800+0.005) return    4.3945    ;
if(t<=  0.0850+0.005) return    4.7092    ;
if(t<=  0.0900+0.005) return    5.0923    ;
if(t<=  0.0950+0.005) return    5.5407    ;
if(t<=  0.1000+0.005) return    6.0483    ;
if(t<=  0.1050+0.005) return    6.6022    ;
if(t<=  0.1100+0.005) return    7.1970    ;
if(t<=  0.1150+0.005) return    7.8215    ;
if(t<=  0.1200+0.005) return    8.4631    ;
if(t<=  0.1250+0.005) return    9.1040    ;
if(t<=  0.1300+0.005) return    9.7359    ;
if(t<=  0.1350+0.005) return   10.3469    ;
if(t<=  0.1400+0.005) return   10.9254    ;
if(t<=  0.1450+0.005) return   11.4571    ;
if(t<=  0.1500+0.005) return   11.9384    ;
if(t<=  0.1550+0.005) return   12.3639    ;
if(t<=  0.1600+0.005) return   12.7294    ;
if(t<=  0.1650+0.005) return   13.0291    ;
if(t<=  0.1700+0.005) return   13.2675    ;
if(t<=  0.1750+0.005) return   13.4473    ;
if(t<=  0.1800+0.005) return   13.5722    ;
if(t<=  0.1850+0.005) return   13.6435    ;
if(t<=  0.1900+0.005) return   13.6715    ;
if(t<=  0.1950+0.005) return   13.6635    ;
if(t<=  0.2000+0.005) return   13.6264    ;
if(t<=  0.2050+0.005) return   13.5633    ;
if(t<=  0.2100+0.005) return   13.4843    ;
if(t<=  0.2150+0.005) return   13.3953    ;
if(t<=  0.2200+0.005) return   13.3006    ;
if(t<=  0.2250+0.005) return   13.1996    ;
if(t<=  0.2300+0.005) return   13.0986    ;
if(t<=  0.2350+0.005) return   12.9990    ;
if(t<=  0.2400+0.005) return   12.9005    ;
if(t<=  0.2450+0.005) return   12.7979    ;
if(t<=  0.2500+0.005) return   12.6940    ;
if(t<=  0.2550+0.005) return   12.5867    ;
if(t<=  0.2600+0.005) return   12.4734    ;
if(t<=  0.2650+0.005) return   12.3465    ;
if(t<=  0.2700+0.005) return   12.2089    ;
if(t<=  0.2750+0.005) return   12.0592    ;
if(t<=  0.2800+0.005) return   11.8960    ;
if(t<=  0.2850+0.005) return    11.7134   ;
if(t<=  0.2900+0.005) return   11.5169    ;
if(t<=  0.2950+0.005) return   11.3076    ;
if(t<=  0.3000+0.005) return   11.0867    ;
if(t<=  0.3050+0.005) return   10.8509    ;
if(t<=  0.3100+0.005) return   10.6082    ;
if(t<=  0.3150+0.005) return   10.3619    ;
if(t<=  0.3200+0.005) return   10.1146    ;
if(t<=  0.3250+0.005) return    9.8641    ;
if(t<=  0.3300+0.005) return    9.6182    ;
if(t<=  0.3350+0.005) return    9.3794    ;
if(t<=  0.3400+0.005) return    9.1490    ;
if(t<=  0.3450+0.005) return    8.9231    ;
if(t<=  0.3500+0.005) return    8.7074    ;
if(t<=  0.3550+0.005) return    8.5024    ;
if(t<=  0.3600+0.005) return    8.3071    ;
if(t<=  0.3650+0.005) return    8.1159    ;
if(t<=  0.3700+0.005) return    7.9326    ;
if(t<=  0.3750+0.005) return    7.7563    ;
if(t<=  0.3800+0.005) return    7.5854    ;
if(t<=  0.3850+0.005) return    7.4141    ;
if(t<=  0.3900+0.005) return    7.2466    ;
if(t<=  0.3950+0.005) return    7.0834    ;
if(t<=  0.4000+0.005) return    6.9244    ;
if(t<=  0.4050+0.005) return    6.7657    ;
if(t<=  0.4100+0.005) return    6.6137    ;
if(t<=  0.4150+0.005) return    6.4706    ;
if(t<=  0.4200+0.005) return    6.3383    ;
if(t<=  0.4250+0.005) return    6.2147    ;
if(t<=  0.4300+0.005) return    6.1069    ;
if(t<=  0.4350+0.005) return    6.0177    ;
if(t<=  0.4400+0.005) return    5.9488    ;
if(t<=  0.4450+0.005) return    5.8975    ;
if(t<=  0.4500+0.005) return    5.8690    ;
if(t<=  0.4550+0.005) return    5.8637    ;
if(t<=  0.4600+0.005) return    5.8803    ;
if(t<=  0.4650+0.005) return    5.9135    ;
if(t<=  0.4700+0.005) return    5.9645    ;
if(t<=  0.4750+0.005) return    6.0304    ;
if(t<=  0.4800+0.005) return    6.1070    ;
if(t<=  0.4850+0.005) return    6.1862    ;
if(t<=  0.4900+0.005) return    6.2669    ;
if(t<=  0.4950+0.005) return    6.3445    ;
if(t<=  0.5000+0.005) return    6.4144    ;
if(t<=  0.5050+0.005) return    6.4689    ;
if(t<=  0.5100+0.005) return    6.5080    ;
if(t<=  0.5150+0.005) return    6.5293    ;
if(t<=  0.5200+0.005) return    6.5309    ;
if(t<=  0.5250+0.005) return    6.5090    ;
if(t<=  0.5300+0.005) return    6.4671    ;
if(t<=  0.5350+0.005) return    6.4067    ;
if(t<=  0.5400+0.005) return    6.3297    ;
if(t<=  0.5450+0.005) return    6.2363    ;
if(t<=  0.5500+0.005) return    6.1327    ;
if(t<=  0.5550+0.005) return    6.0232    ;
if(t<=  0.5600+0.005) return    5.9114    ;
if(t<=  0.5650+0.005) return    5.7985    ;
if(t<=  0.5700+0.005) return    5.6906    ;
if(t<=  0.5750+0.005) return    5.5909    ;
if(t<=  0.5800+0.005) return    5.5014    ;
if(t<=  0.5850+0.005) return    5.4210    ;
if(t<=  0.5900+0.005) return    5.3530    ;
if(t<=  0.5950+0.005) return    5.2974    ;
if(t<=  0.6000+0.005) return    5.2531    ;
if(t<=  0.6050+0.005) return    5.2156    ;
if(t<=  0.6100+0.005) return    5.1857    ;
if(t<=  0.6150+0.005) return    5.1609    ;
if(t<=  0.6200+0.005) return    5.1384    ;
if(t<=  0.6250+0.005) return    5.1127    ;
if(t<=  0.6300+0.005) return    5.0844    ;
if(t<=  0.6350+0.005) return    5.0519    ;
if(t<=  0.6400+0.005) return    5.0136    ;
if(t<=  0.6450+0.005) return    4.9658    ;
if(t<=  0.6500+0.005) return    4.9116    ;
if(t<=  0.6550+0.005) return    4.8517    ;
if(t<=  0.6600+0.005) return    4.7873    ;
if(t<=  0.6650+0.005) return    4.7177    ;
if(t<=  0.6700+0.005) return    4.6474    ;
if(t<=  0.6750+0.005) return    4.5791    ;
if(t<=  0.6800+0.005) return    4.5152    ;
if(t<=  0.6850+0.005) return    4.4540    ;
if(t<=  0.6900+0.005) return    4.4031    ;
if(t<=  0.6950+0.005) return    4.3654    ;
if(t<=  0.7000+0.005) return    4.3416    ;
if(t<=  0.7050+0.005) return    4.3318    ;
if(t<=  0.7100+0.005) return    4.3315    ;
if(t<=  0.7150+0.005) return    4.3363    ;
if(t<=  0.7200+0.005) return    4.3418    ;
if(t<=  0.7250+0.005) return    4.3392    ;
if(t<=  0.7300+0.005) return    4.3368    ;
if(t<=  0.7350+0.005) return    4.3380    ;
if(t<=  0.7400+0.005) return    4.3462    ;
if(t<=  0.7450+0.005) return    4.3651    ;
if(t<=  0.7500+0.005) return    4.3925    ;
if(t<=  0.7550+0.005) return    4.4259    ;
if(t<=  0.7600+0.005) return    4.4609    ;
if(t<=  0.7650+0.005) return    4.4888    ;
if(t<=  0.7700+0.005) return    4.5116    ;
if(t<=  0.7750+0.005) return    4.5276    ;
if(t<=  0.7800+0.005) return    4.5358    ;
if(t<=  0.7850+0.005) return    4.5343    ;
if(t<=  0.7900+0.005) return    4.5247    ;
if(t<=  0.7950+0.005) return    4.5081    ;
if(t<=  0.8000+0.005) return    4.4854    ;
if(t<=  0.8050+0.005) return    4.4561    ;
if(t<=  0.8100+0.005) return    4.4246    ;
if(t<=  0.8150+0.005) return    4.3934    ;
if(t<=  0.8200+0.005) return    4.3647    ;
if(t<=  0.8250+0.005) return    4.3382    ;
if(t<=  0.8300+0.005) return    4.3182    ;
if(t<=  0.8350+0.005) return    4.3062    ;
if(t<=  0.8400+0.005) return    4.3028    ;
if(t<=  0.8450+0.005) return    4.3061    ;
if(t<=  0.8500+0.005) return    4.3179    ;
if(t<=  0.8550+0.005) return    4.3374    ;
if(t<=  0.8600+0.005) return    4.3631    ;
if(t<=  0.8650+0.005) return    4.3906    ;
if(t<=  0.8700+0.005) return    4.4202    ;
if(t<=  0.8750+0.005) return    4.4498    ;
if(t<=  0.8800+0.005) return    4.4769    ;
if(t<=  0.8850+0.005) return    4.4972    ;
if(t<=  0.8900+0.005) return    4.5114    ;
if(t<=  0.8950+0.005) return    4.5186    ;
if(t<=  0.9000+0.005) return    4.5180    ;
if(t<=  0.9050+0.005) return    4.5072    ;
if(t<=  0.9100+0.005) return    4.4894    ;
if(t<=  0.9150+0.005) return    4.4659    ;
if(t<=  0.9200+0.005) return    4.4382    ;
if(t<=  0.9250+0.005) return    4.4061    ;
if(t<=  0.9300+0.005) return    4.3740    ;
if(t<=  0.9350+0.005) return    4.3445    ;
if(t<=  0.9400+0.005) return    4.3196    ;
if(t<=  0.9450+0.005) return    4.2988    ;
if(t<=  0.9500+0.005) return    4.2859    ;
if(t<=  0.9550+0.005) return    4.2818    ;
if(t<=  0.9600+0.005) return    4.2866    ;
if(t<=  0.9650+0.005) return    4.2976    ;
if(t<=  0.9700+0.005) return    4.3160    ;
if(t<=  0.9750+0.005) return    4.3403    ;
if(t<=  0.9800+0.005) return    4.3684    ;
if(t<=  0.9850+0.005) return    4.3956    ;
if(t<=  0.9900+0.005) return    4.4216    ;
if(t<=  0.9950+0.005) return    4.4444    ;
if(t<=  1.0000+0.005) return    4.4616    ;
if(t<=  1.0050+0.005) return    4.4690    ;
if(t<=  1.0100+0.005) return    4.4680    ;
if(t<=  1.0150+0.005) return    4.4583    ;
if(t<=  1.0200+0.005) return    4.4400    ;
if(t<=  1.0250+0.005) return    4.4118    ;
if(t<=  1.0300+0.005) return    4.3779    ;
if(t<=  1.0350+0.005) return    4.3406    ;
if(t<=  1.0400+0.005) return    4.3027    ;
if(t<=  1.0450+0.005) return    4.2645    ;
if(t<=  1.0500+0.005) return    4.2318    ;
if(t<=  1.0550+0.005) return    4.2079    ;
if(t<=  1.0600+0.005) return    4.1949    ;
if(t<=  1.0650+0.005) return    4.1924    ;
if(t<=  1.0700+0.005) return    4.2030    ;
if(t<=  1.0750+0.005) return    4.2261    ;
if(t<=  1.0800+0.005) return    4.2606    ;
if(t<=  1.0850+0.005) return    4.3021    ;
if(t<=  1.0900+0.005) return    4.3510    ;
if(t<=  1.0950+0.005) return    4.4055    ;
if(t<=  1.1000+0.005) return    4.4630    ;
if(t<=  1.1050+0.005) return    4.5187    ;
if(t<=  1.1100+0.005) return    4.5724    ;
                          }


Real linearInletCylinder( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{

    Real tNew = t;
    if(t>1.102)
        tNew=(((int)floor(t*10000))%11000)/10000.0;

    Real dt =  static_cast<int>( std::floor( std::ceil( tNew * 1000000.0 ) / 1000.0 ) ) % 5 / 1000.0;
    Real ti=floor(tNew*1000)/1000.0 - dt;
    Real tii=ti+0.005;
    return -((inletCylinder(tii,0,0,0,0)-inletCylinder(ti,0,0,0,0))/(0.005)*(tNew-(ti))+inletCylinder(ti,0,0,0,0));
}


Real linearVelInletCylinder( Real  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  //Components for Simone's mesh

  Real n1(0.0);
  Real n2(0.0);
  Real n3(-1.0);

  Real flux( -7.0 /*linearInletCylinder(t,x,y,z,i)*/);
    //Center for the aneurym with one IO
    Real x0(0);
    Real y0(0);
    Real area(0.5*0.5*3.141592653589793); //fluidBig

    Real radiusSquared(0);
    radiusSquared = 0.5*0.5; //area / 3.1415962;

    Real peak(0);
    peak = ( 2.0 * flux ) / ( area );

    switch(i) {
    case 0:
      return n1*( peak * std::max( 0.0, ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ));
    case 1:
      return n2*( peak * std::max( 0.0, ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ));
    case 2:
      return n3*( peak * std::max( 0.0, ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ));
    default:
        return 0.0;
    }
}


Real oneVelInletCylinder( Real  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  //Components for Simone's mesh

  Real n1(0.0);
  Real n2(0.0);
  Real n3(-1.0);

    Real flux(1.0);
    //Center for the aneurym with one IO
    Real x0(0);
    Real y0(0);
    Real area(0.5*0.5*3.141592653589793); //fluidBig

    Real radiusSquared(0);
    radiusSquared = 0.5*0.5; //area / 3.1415962;

    Real peak(0);
    peak = ( 2.0 * flux ) / ( area );

    switch(i) {
    case 0:
      return n1*( peak * std::max( 0.0, ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ));
    case 1:
      return n2*( peak * std::max( 0.0, ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ));
    case 2:
      return n3*( peak * std::max( 0.0, ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ));
    default:
        return 0.0;
    }
}

Real popliteal(const Real t, const Real& , const Real& , const Real& , const ID& i) //outlet flux from the bifurcation
{
        Real newt = ((Real)(((int)(t*1000)) % 792))/1000;
if( newt <= 0.    +0.004 ) return      0.627074	     ; else
if( newt <= 0.004 +0.004 ) return          0.822842	     ; else
if( newt <= 0.008 +0.004 ) return          1.15228	     ; else
if( newt <= 0.012 +0.004 ) return          1.46407	     ; else
if( newt <= 0.016 +0.004 ) return          1.70808	     ; else
if( newt <= 0.02  +0.004 ) return         1.96562	     ; else
if( newt <= 0.024 +0.004 ) return          2.23675	     ; else
if( newt <= 0.028 +0.004 ) return          2.4672	     ; else
if( newt <= 0.032 +0.004 ) return          2.68409	     ; else
if( newt <= 0.036 +0.004 ) return          2.88745	     ; else
if( newt <= 0.04  +0.004 ) return         3.18553	     ; else
if( newt <= 0.044 +0.004 ) return          3.45675	     ; else
if( newt <= 0.048 +0.004 ) return          3.80916	     ; else
if( newt <= 0.052 +0.004 ) return          4.13451	     ; else
if( newt <= 0.056 +0.004 ) return          4.41926	     ; else
if( newt <= 0.06  +0.004 ) return         4.609	     ; else
if( newt <= 0.064 +0.004 ) return          4.69048	     ; else
if( newt <= 0.068 +0.004 ) return          4.81226	     ; else
if( newt <= 0.072 +0.004 ) return          4.97494	     ; else
if( newt <= 0.076 +0.004 ) return          5.1785	     ; else
if( newt <= 0.08  +0.004 ) return         5.34118	     ; else
if( newt <= 0.084 +0.004 ) return          5.46296	     ; else
if( newt <= 0.088 +0.004 ) return          5.6121	     ; else
if( newt <= 0.092 +0.004 ) return          5.70712	     ; else
if( newt <= 0.096 +0.004 ) return          5.74771	     ; else
if( newt <= 0.1   +0.004 ) return        5.78831	     ; else
if( newt <= 0.104 +0.004 ) return          5.91039	     ; else
if( newt <= 0.108 +0.004 ) return          6.03247	     ; else
if( newt <= 0.112 +0.004 ) return          6.16808	     ; else
if( newt <= 0.116 +0.004 ) return          6.12719	     ; else
if( newt <= 0.12  +0.004 ) return         6.046	     ; else
if( newt <= 0.124 +0.004 ) return          5.99187	     ; else
if( newt <= 0.128 +0.004 ) return          6.08659	     ; else
if( newt <= 0.132 +0.004 ) return          6.19514	     ; else
if( newt <= 0.136 +0.004 ) return          6.20867	     ; else
if( newt <= 0.14  +0.004 ) return         6.12719	     ; else
if( newt <= 0.144 +0.004 ) return          6.05953	     ; else
if( newt <= 0.148 +0.004 ) return          6.16808	     ; else
if( newt <= 0.152 +0.004 ) return          6.50696	     ; else
if( newt <= 0.156 +0.004 ) return          7.15765	     ; else
if( newt <= 0.16  +0.004 ) return         8.14723	     ; else
if( newt <= 0.164 +0.004 ) return          9.05532	     ; else
if( newt <= 0.168 +0.004 ) return          9.31301	     ; else
if( newt <= 0.172 +0.004 ) return          8.48611	     ; else
if( newt <= 0.176 +0.004 ) return          7.483	     ; else
if( newt <= 0.18  +0.004 ) return         6.52049	     ; else
if( newt <= 0.184 +0.004 ) return          5.8292	     ; else
if( newt <= 0.188 +0.004 ) return          5.40883	     ; else
if( newt <= 0.192 +0.004 ) return          5.1785	     ; else
if( newt <= 0.196 +0.004 ) return          5.12408	     ; else
if( newt <= 0.2   +0.004 ) return        5.13761	     ; else
if( newt <= 0.204 +0.004 ) return          5.15144	     ; else
if( newt <= 0.208 +0.004 ) return          5.15144	     ; else
if( newt <= 0.212 +0.004 ) return          4.94787	     ; else
if( newt <= 0.216 +0.004 ) return          4.67695	     ; else
if( newt <= 0.22  +0.004 ) return         4.3516	     ; else
if( newt <= 0.224 +0.004 ) return          4.03978	     ; else
if( newt <= 0.228 +0.004 ) return          3.80916	     ; else
if( newt <= 0.232 +0.004 ) return          3.56529	     ; else
if( newt <= 0.236 +0.004 ) return          3.29407	     ; else
if( newt <= 0.24  +0.004 ) return         2.9281	     ; else
if( newt <= 0.244 +0.004 ) return          2.52142	     ; else
if( newt <= 0.248 +0.004 ) return          1.89784	     ; else
if( newt <= 0.252 +0.004 ) return          1.4505	     ; else
if( newt <= 0.256 +0.004 ) return          1.24715	     ; else
if( newt <= 0.26  +0.004 ) return         1.05738	     ; else
if( newt <= 0.264 +0.004 ) return          0.745594	     ; else
if( newt <= 0.268 +0.004 ) return          0.474461	     ; else
if( newt <= 0.272 +0.004 ) return          0.0135561	     ; else
if( newt <= 0.276 +0.004 ) return          -0.515144	     ; else
if( newt <= 0.28  +0.004 ) return         -0.826931	     ; else
if( newt <= 0.284 +0.004 ) return          -1.12516	     ; else
if( newt <= 0.288 +0.004 ) return          -1.59962	     ; else
if( newt <= 0.292 +0.004 ) return          -1.74873	     ; else
if( newt <= 0.296 +0.004 ) return          -1.77585	     ; else
if( newt <= 0.3   +0.004 ) return        -1.72161	     ; else
if( newt <= 0.304 +0.004 ) return          -1.76229	     ; else
if( newt <= 0.308 +0.004 ) return          -1.80294	     ; else
if( newt <= 0.312 +0.004 ) return          -1.76229	     ; else
if( newt <= 0.316 +0.004 ) return          -1.73517	     ; else
if( newt <= 0.32  +0.004 ) return         -1.70808	     ; else
if( newt <= 0.324 +0.004 ) return          -1.65383	     ; else
if( newt <= 0.328 +0.004 ) return          -1.62674	     ; else
if( newt <= 0.332 +0.004 ) return          -1.58606	     ; else
if( newt <= 0.336 +0.004 ) return          -1.58606	     ; else
if( newt <= 0.34  +0.004 ) return         -1.55893	     ; else
if( newt <= 0.344 +0.004 ) return          -1.55893	     ; else
if( newt <= 0.348 +0.004 ) return          -1.58606	     ; else
if( newt <= 0.352 +0.004 ) return          -1.55893	     ; else
if( newt <= 0.356 +0.004 ) return          -1.53184	     ; else
if( newt <= 0.36  +0.004 ) return         -1.58606	     ; else
if( newt <= 0.364 +0.004 ) return          -1.65383	     ; else
if( newt <= 0.368 +0.004 ) return          -1.66739	     ; else
if( newt <= 0.372 +0.004 ) return          -1.55893	     ; else
if( newt <= 0.376 +0.004 ) return          -1.4505	     ; else
if( newt <= 0.38  +0.004 ) return         -1.36917	     ; else
if( newt <= 0.384 +0.004 ) return          -1.28783	     ; else
if( newt <= 0.388 +0.004 ) return          -1.16581	     ; else
if( newt <= 0.392 +0.004 ) return          -1.04382	     ; else
if( newt <= 0.396 +0.004 ) return          -1.00314	     ; else
if( newt <= 0.4   +0.004 ) return        -1.04382	     ; else
if( newt <= 0.404 +0.004 ) return          -1.07094	     ; else
if( newt <= 0.408 +0.004 ) return          -1.0167	     ; else
if( newt <= 0.412 +0.004 ) return          -0.921799	     ; else
if( newt <= 0.416 +0.004 ) return          -0.799808	     ; else
if( newt <= 0.42  +0.004 ) return         -0.664257	     ; else
if( newt <= 0.424 +0.004 ) return          -0.491405	     ; else
if( newt <= 0.428 +0.004 ) return          -0.332114	     ; else
if( newt <= 0.432 +0.004 ) return          -0.227064	     ; else
if( newt <= 0.436 +0.004 ) return          -0.244008	     ; else
if( newt <= 0.44  +0.004 ) return         -0.220286	     ; else
if( newt <= 0.444 +0.004 ) return          -0.20108	     ; else
if( newt <= 0.448 +0.004 ) return          -0.141209	     ; else
if( newt <= 0.452 +0.004 ) return          -0.0813369     ; else
if( newt <= 0.456 +0.004 ) return          0.023723	     ; else
if( newt <= 0.46  +0.004 ) return         0.128783	     ; else
if( newt <= 0.464 +0.004 ) return          0.260955	     ; else
if( newt <= 0.468 +0.004 ) return          0.420246	     ; else
if( newt <= 0.472 +0.004 ) return          0.593068	     ; else
if( newt <= 0.476 +0.004 ) return          0.684584	     ; else
if( newt <= 0.48  +0.004 ) return         0.721854	     ; else
if( newt <= 0.484 +0.004 ) return          0.70491	     ; else
if( newt <= 0.488 +0.004 ) return          0.745594	     ; else
if( newt <= 0.492 +0.004 ) return          0.813369	     ; else
if( newt <= 0.496 +0.004 ) return          0.962482	     ; else
if( newt <= 0.5   +0.004 ) return        1.04382	     ; else
if( newt <= 0.504 +0.004 ) return          1.07094	     ; else
if( newt <= 0.508 +0.004 ) return          1.09803	     ; else
if( newt <= 0.512 +0.004 ) return          1.15228	     ; else
if( newt <= 0.516 +0.004 ) return          1.12516	     ; else
if( newt <= 0.52  +0.004 ) return         1.07094	     ; else
if( newt <= 0.524 +0.004 ) return          1.0167	     ; else
if( newt <= 0.528 +0.004 ) return          0.989604	     ; else
if( newt <= 0.532 +0.004 ) return          0.894706	     ; else
if( newt <= 0.536 +0.004 ) return          0.921799	     ; else
if( newt <= 0.54  +0.004 ) return         1.04382	     ; else
if( newt <= 0.544 +0.004 ) return          0.959099	     ; else
if( newt <= 0.548 +0.004 ) return          0.847257	     ; else
if( newt <= 0.552 +0.004 ) return          0.870967	     ; else
if( newt <= 0.556 +0.004 ) return          0.948921	     ; else
if( newt <= 0.56  +0.004 ) return         1.15228	     ; else
if( newt <= 0.564 +0.004 ) return          1.13872	     ; else
if( newt <= 0.568 +0.004 ) return          1.1116	     ; else
if( newt <= 0.572 +0.004 ) return          1.0167	     ; else
if( newt <= 0.576 +0.004 ) return          1.04382	     ; else
if( newt <= 0.58  +0.004 ) return         0.989604	     ; else
if( newt <= 0.584 +0.004 ) return          1.05738	     ; else
if( newt <= 0.588 +0.004 ) return          1.17937	     ; else
if( newt <= 0.592 +0.004 ) return          1.27427	     ; else
if( newt <= 0.596 +0.004 ) return          1.28783	     ; else
if( newt <= 0.6   +0.004 ) return        1.31492	     ; else
if( newt <= 0.604 +0.004 ) return          1.30139	     ; else
if( newt <= 0.608 +0.004 ) return          1.30139	     ; else
if( newt <= 0.612 +0.004 ) return          1.31492	     ; else
if( newt <= 0.616 +0.004 ) return          1.34204	     ; else
if( newt <= 0.62  +0.004 ) return         1.36917	     ; else
if( newt <= 0.624 +0.004 ) return          1.40982	     ; else
if( newt <= 0.628 +0.004 ) return          1.36917	     ; else
if( newt <= 0.632 +0.004 ) return          1.35561	     ; else
if( newt <= 0.636 +0.004 ) return          1.32848	     ; else
if( newt <= 0.64  +0.004 ) return         1.32848	     ; else
if( newt <= 0.644 +0.004 ) return          1.31492	     ; else
if( newt <= 0.648 +0.004 ) return          1.31492	     ; else
if( newt <= 0.652 +0.004 ) return          1.27427	     ; else
if( newt <= 0.656 +0.004 ) return          1.22005	     ; else
if( newt <= 0.66  +0.004 ) return         1.15228	     ; else
if( newt <= 0.664 +0.004 ) return          1.12516	     ; else
if( newt <= 0.668 +0.004 ) return          1.12516	     ; else
if( newt <= 0.672 +0.004 ) return          1.16581	     ; else
if( newt <= 0.676 +0.004 ) return          1.17937	     ; else
if( newt <= 0.68  +0.004 ) return         1.15228	     ; else
if( newt <= 0.684 +0.004 ) return          1.15228	     ; else
if( newt <= 0.688 +0.004 ) return          1.03026	     ; else
if( newt <= 0.692 +0.004 ) return          0.908267	     ; else
if( newt <= 0.696 +0.004 ) return          0.244008	     ; else
if( newt <= 0.7   +0.004 ) return        0.0948921	     ; else
if( newt <= 0.704 +0.004 ) return          0.284676	     ; else
if( newt <= 0.708 +0.004 ) return          0.433807	     ; else
if( newt <= 0.712 +0.004 ) return          0.58292	     ; else
if( newt <= 0.716 +0.004 ) return          0.650696	     ; else
if( newt <= 0.72  +0.004 ) return         0.623573	     ; else
if( newt <= 0.724 +0.004 ) return          0.515144	     ; else
if( newt <= 0.728 +0.004 ) return          0.176229	     ; else
if( newt <= 0.732 +0.004 ) return          0		     ; else
if( newt <= 0.736 +0.004 ) return          -0.216898	     ; else
if( newt <= 0.74  +0.004 ) return         -0.691349	     ; else
if( newt <= 0.744 +0.004 ) return          -0.93536	     ; else
if( newt <= 0.748 +0.004 ) return          -0.987663	     ; else
if( newt <= 0.752 +0.004 ) return          -0.985721	     ; else
if( newt <= 0.756 +0.004 ) return          -0.916004	     ; else
if( newt <= 0.76  +0.004 ) return         -0.832726	     ; else
if( newt <= 0.764 +0.004 ) return          -0.708793	     ; else
if( newt <= 0.768 +0.004 ) return          -0.655343	     ; else
if( newt <= 0.772 +0.004 ) return          -0.601893	     ; else
if( newt <= 0.776 +0.004 ) return          -0.453545	     ; else
if( newt <= 0.78  +0.004 ) return         -0.305197	     ; else
if( newt <= 0.784 +0.004 ) return          -0.183974	     ; else
if( newt <= 0.788 +0.004 ) return          -0.144083	     ; else
if( newt <= 0.792 +0.004 ) return          -0.104188      ;
}


Real pont_dist(const Real t, const Real& , const Real& , const Real& , const ID& i) //inlet flux from the bypass
{
        Real newt = ((Real)(((int)(t*1000)) % 792))/1000;
if( newt <= 0.    +0.004 ) return  1.96970000e-02   * (-10); else
if( newt <= 0.004 +0.004 ) return  1.93150000e-02   * (-10); else
if( newt <= 0.008 +0.004 ) return  2.12310000e-02   * (-10); else
if( newt <= 0.012 +0.004 ) return  3.04300000e-02   * (-10); else
if( newt <= 0.016 +0.004 ) return  4.50000000e-02   * (-10); else
if( newt <= 0.02  +0.004 ) return  6.18670000e-02   * (-10); else
if( newt <= 0.024 +0.004 ) return  8.48730000e-02   * (-10); else
if( newt <= 0.028 +0.004 ) return  1.07874000e-01   * (-10); else
if( newt <= 0.032 +0.004 ) return  1.32410000e-01   * (-10); else
if( newt <= 0.036 +0.004 ) return  1.57714000e-01   * (-10); else
if( newt <= 0.04  +0.004 ) return  1.86849000e-01   * (-10); else
if( newt <= 0.044 +0.004 ) return  2.19821000e-01   * (-10); else
if( newt <= 0.048 +0.004 ) return  2.52788000e-01   * (-10); else
if( newt <= 0.052 +0.004 ) return  2.85760000e-01   * (-10); else
if( newt <= 0.056 +0.004 ) return  3.17198000e-01   * (-10); else
if( newt <= 0.06  +0.004 ) return  3.44798000e-01   * (-10); else
if( newt <= 0.064 +0.004 ) return  3.70869000e-01   * (-10); else
if( newt <= 0.068 +0.004 ) return  3.96959000e-01   * (-10); else
if( newt <= 0.072 +0.004 ) return  4.26080000e-01   * (-10); else
if( newt <= 0.076 +0.004 ) return  4.53675000e-01   * (-10); else
if( newt <= 0.08  +0.004 ) return  4.82034000e-01   * (-10); else
if( newt <= 0.084 +0.004 ) return  5.09630000e-01   * (-10); else
if( newt <= 0.088 +0.004 ) return  5.27264000e-01   * (-10); else
if( newt <= 0.092 +0.004 ) return  5.38798000e-01   * (-10); else
if( newt <= 0.096 +0.004 ) return  5.48759000e-01   * (-10); else
if( newt <= 0.1   +0.004 ) return  5.52572000e-01   * (-10); else
if( newt <= 0.104 +0.004 ) return  5.63344000e-01   * (-10); else
if( newt <= 0.108 +0.004 ) return  5.80978000e-01   * (-10); else
if( newt <= 0.112 +0.004 ) return  6.00138000e-01   * (-10); else
if( newt <= 0.116 +0.004 ) return  6.19298000e-01   * (-10); else
if( newt <= 0.12  +0.004 ) return  6.32310000e-01   * (-10); else
if( newt <= 0.124 +0.004 ) return  6.40746000e-01   * (-10); else
if( newt <= 0.128 +0.004 ) return  6.41508000e-01   * (-10); else
if( newt <= 0.132 +0.004 ) return  6.43081000e-01   * (-10); else
if( newt <= 0.136 +0.004 ) return  6.49944000e-01   * (-10); else
if( newt <= 0.14  +0.004 ) return  6.59190000e-01   * (-10); else
if( newt <= 0.144 +0.004 ) return  6.64528000e-01   * (-10); else
if( newt <= 0.148 +0.004 ) return  6.60716000e-01   * (-10); else
if( newt <= 0.152 +0.004 ) return  6.56092000e-01   * (-10); else
if( newt <= 0.156 +0.004 ) return  6.59190000e-01   * (-10); else
if( newt <= 0.16  +0.004 ) return  6.62241000e-01   * (-10); else
if( newt <= 0.164 +0.004 ) return  6.62241000e-01   * (-10); else
if( newt <= 0.168 +0.004 ) return  6.55330000e-01   * (-10); else
if( newt <= 0.172 +0.004 ) return  6.46131000e-01   * (-10); else
if( newt <= 0.176 +0.004 ) return  6.35408000e-01   * (-10); else
if( newt <= 0.18  +0.004 ) return  6.23111000e-01   * (-10); else
if( newt <= 0.184 +0.004 ) return  6.17010000e-01   * (-10); else
if( newt <= 0.188 +0.004 ) return  6.10862000e-01   * (-10); else
if( newt <= 0.192 +0.004 ) return  5.99376000e-01   * (-10); else
if( newt <= 0.196 +0.004 ) return  5.86316000e-01   * (-10); else
if( newt <= 0.2   +0.004 ) return  5.70970000e-01   * (-10); else
if( newt <= 0.204 +0.004 ) return  5.54145000e-01   * (-10); else
if( newt <= 0.208 +0.004 ) return  5.43374000e-01   * (-10); else
if( newt <= 0.212 +0.004 ) return  5.35700000e-01   * (-10); else
if( newt <= 0.216 +0.004 ) return  5.26502000e-01   * (-10); else
if( newt <= 0.22  +0.004 ) return  5.11965000e-01   * (-10); else
if( newt <= 0.224 +0.004 ) return  4.91232000e-01   * (-10); else
if( newt <= 0.228 +0.004 ) return  4.69022000e-01   * (-10); else
if( newt <= 0.232 +0.004 ) return  4.46002000e-01   * (-10); else
if( newt <= 0.236 +0.004 ) return  4.27605000e-01   * (-10); else
if( newt <= 0.24  +0.004 ) return  4.10733000e-01   * (-10); else
if( newt <= 0.244 +0.004 ) return  3.86950000e-01   * (-10); else
if( newt <= 0.248 +0.004 ) return  3.65502000e-01   * (-10); else
if( newt <= 0.252 +0.004 ) return  3.42501000e-01   * (-10); else
if( newt <= 0.256 +0.004 ) return  3.24862000e-01   * (-10); else
if( newt <= 0.26  +0.004 ) return  3.07994000e-01   * (-10); else
if( newt <= 0.264 +0.004 ) return  2.88058000e-01   * (-10); else
if( newt <= 0.268 +0.004 ) return  2.69656000e-01   * (-10); else
if( newt <= 0.272 +0.004 ) return  2.47422000e-01   * (-10); else
if( newt <= 0.276 +0.004 ) return  2.23653000e-01   * (-10); else
if( newt <= 0.28  +0.004 ) return  2.04484000e-01   * (-10); else
if( newt <= 0.284 +0.004 ) return  1.84547000e-01   * (-10); else
if( newt <= 0.288 +0.004 ) return  1.67680000e-01   * (-10); else
if( newt <= 0.292 +0.004 ) return  1.54644000e-01   * (-10); else
if( newt <= 0.296 +0.004 ) return  1.42376000e-01   * (-10); else
if( newt <= 0.3   +0.004 ) return  1.32410000e-01   * (-10); else
if( newt <= 0.304 +0.004 ) return  1.20142000e-01   * (-10); else
if( newt <= 0.308 +0.004 ) return  1.06340000e-01   * (-10); else
if( newt <= 0.312 +0.004 ) return  9.25370000e-02   * (-10); else
if( newt <= 0.316 +0.004 ) return  7.64370000e-02   * (-10); else
if( newt <= 0.32  +0.004 ) return  5.88030000e-02   * (-10); else
if( newt <= 0.324 +0.004 ) return  9.73060000e-03   * (-10); else
if( newt <= 0.328 +0.004 ) return  -3.09101000e-02  * (-10); else
if( newt <= 0.332 +0.004 ) return  -6.84790000e-02  * (-10); else
if( newt <= 0.336 +0.004 ) return  -1.58956600e-01  * (-10); else
if( newt <= 0.34  +0.004 ) return  -2.08029000e-01  * (-10); else
if( newt <= 0.344 +0.004 ) return  -2.37931000e-01  * (-10); else
if( newt <= 0.348 +0.004 ) return  -2.68601000e-01  * (-10); else
if( newt <= 0.352 +0.004 ) return  -2.71666000e-01  * (-10); else
if( newt <= 0.356 +0.004 ) return  -2.74735000e-01  * (-10); else
if( newt <= 0.36  +0.004 ) return  -2.70136000e-01  * (-10); else
if( newt <= 0.364 +0.004 ) return  -2.65532000e-01  * (-10); else
if( newt <= 0.368 +0.004 ) return  -2.59398000e-01  * (-10); else
if( newt <= 0.372 +0.004 ) return  -2.59398000e-01  * (-10); else
if( newt <= 0.376 +0.004 ) return  -2.60165000e-01  * (-10); else
if( newt <= 0.38  +0.004 ) return  -2.60932000e-01  * (-10); else
if( newt <= 0.384 +0.004 ) return  -2.60165000e-01  * (-10); else
if( newt <= 0.388 +0.004 ) return  -2.57101000e-01  * (-10); else
if( newt <= 0.392 +0.004 ) return  -2.58635000e-01  * (-10); else
if( newt <= 0.396 +0.004 ) return  -2.56333000e-01  * (-10); else
if( newt <= 0.4   +0.004 ) return  -2.49432000e-01  * (-10); else
if( newt <= 0.404 +0.004 ) return  -2.43298000e-01  * (-10); else
if( newt <= 0.408 +0.004 ) return  -2.23361000e-01  * (-10); else
if( newt <= 0.412 +0.004 ) return  -2.09563000e-01  * (-10); else
if( newt <= 0.416 +0.004 ) return  -1.71992000e-01  * (-10); else
if( newt <= 0.42  +0.004 ) return  -1.22151700e-01  * (-10); else
if( newt <= 0.424 +0.004 ) return  -9.14815000e-02  * (-10); else
if( newt <= 0.428 +0.004 ) return  -6.08113000e-02  * (-10); else
if( newt <= 0.432 +0.004 ) return  -5.69779000e-02  * (-10); else
if( newt <= 0.436 +0.004 ) return  -5.16103000e-02  * (-10); else
if( newt <= 0.44  +0.004 ) return  -3.93414000e-02  * (-10); else
if( newt <= 0.444 +0.004 ) return  -1.78748000e-02  * (-10); else
if( newt <= 0.448 +0.004 ) return  2.42960000e-02   * (-10); else
if( newt <= 0.452 +0.004 ) return  3.80990000e-02   * (-10); else
if( newt <= 0.456 +0.004 ) return  4.26980000e-02   * (-10); else
if( newt <= 0.46  +0.004 ) return  4.26980000e-02   * (-10); else
if( newt <= 0.464 +0.004 ) return  4.50000000e-02   * (-10); else
if( newt <= 0.468 +0.004 ) return  4.57670000e-02   * (-10); else
if( newt <= 0.472 +0.004 ) return  4.73020000e-02   * (-10); else
if( newt <= 0.476 +0.004 ) return  4.88320000e-02   * (-10); else
if( newt <= 0.48  +0.004 ) return  5.26690000e-02   * (-10); else
if( newt <= 0.484 +0.004 ) return  5.49660000e-02   * (-10); else
if( newt <= 0.488 +0.004 ) return  5.49660000e-02   * (-10); else
if( newt <= 0.492 +0.004 ) return  5.72680000e-02   * (-10); else
if( newt <= 0.496 +0.004 ) return  5.95700000e-02   * (-10); else
if( newt <= 0.5   +0.004 ) return  6.11000000e-02   * (-10); else
if( newt <= 0.504 +0.004 ) return  6.18670000e-02   * (-10); else
if( newt <= 0.508 +0.004 ) return  6.49370000e-02   * (-10); else
if( newt <= 0.512 +0.004 ) return  6.87680000e-02   * (-10); else
if( newt <= 0.516 +0.004 ) return  6.64710000e-02   * (-10); else
if( newt <= 0.52  +0.004 ) return  6.41690000e-02   * (-10); else
if( newt <= 0.524 +0.004 ) return  6.57040000e-02   * (-10); else
if( newt <= 0.528 +0.004 ) return  6.87680000e-02   * (-10); else
if( newt <= 0.532 +0.004 ) return  7.18380000e-02   * (-10); else
if( newt <= 0.536 +0.004 ) return  7.64370000e-02   * (-10); else
if( newt <= 0.54  +0.004 ) return  7.64370000e-02   * (-10); else
if( newt <= 0.544 +0.004 ) return  7.41350000e-02   * (-10); else
if( newt <= 0.548 +0.004 ) return  7.49020000e-02   * (-10); else
if( newt <= 0.552 +0.004 ) return  7.56700000e-02   * (-10); else
if( newt <= 0.556 +0.004 ) return  7.33680000e-02   * (-10); else
if( newt <= 0.56  +0.004 ) return  6.87680000e-02   * (-10); else
if( newt <= 0.564 +0.004 ) return  7.03030000e-02   * (-10); else
if( newt <= 0.568 +0.004 ) return  6.72340000e-02   * (-10); else
if( newt <= 0.572 +0.004 ) return  6.41690000e-02   * (-10); else
if( newt <= 0.576 +0.004 ) return  6.11000000e-02   * (-10); else
if( newt <= 0.58  +0.004 ) return  6.11000000e-02   * (-10); else
if( newt <= 0.584 +0.004 ) return  6.41690000e-02   * (-10); else
if( newt <= 0.588 +0.004 ) return  7.18380000e-02   * (-10); else
if( newt <= 0.592 +0.004 ) return  7.95020000e-02   * (-10); else
if( newt <= 0.596 +0.004 ) return  8.56360000e-02   * (-10); else
if( newt <= 0.6   +0.004 ) return  8.25710000e-02   * (-10); else
if( newt <= 0.604 +0.004 ) return  7.64370000e-02   * (-10); else
if( newt <= 0.608 +0.004 ) return  7.10710000e-02   * (-10); else
if( newt <= 0.612 +0.004 ) return  7.33680000e-02   * (-10); else
if( newt <= 0.616 +0.004 ) return  7.18380000e-02   * (-10); else
if( newt <= 0.62  +0.004 ) return  7.10710000e-02   * (-10); else
if( newt <= 0.624 +0.004 ) return  6.80010000e-02   * (-10); else
if( newt <= 0.628 +0.004 ) return  6.95360000e-02   * (-10); else
if( newt <= 0.632 +0.004 ) return  7.03030000e-02   * (-10); else
if( newt <= 0.636 +0.004 ) return  6.72340000e-02   * (-10); else
if( newt <= 0.64  +0.004 ) return  6.18670000e-02   * (-10); else
if( newt <= 0.644 +0.004 ) return  5.65010000e-02   * (-10); else
if( newt <= 0.648 +0.004 ) return  5.65010000e-02   * (-10); else
if( newt <= 0.652 +0.004 ) return  5.49660000e-02   * (-10); else
if( newt <= 0.656 +0.004 ) return  5.65010000e-02   * (-10); else
if( newt <= 0.66  +0.004 ) return  5.88030000e-02   * (-10); else
if( newt <= 0.664 +0.004 ) return  5.95700000e-02   * (-10); else
if( newt <= 0.668 +0.004 ) return  6.03370000e-02   * (-10); else
if( newt <= 0.672 +0.004 ) return  6.03370000e-02   * (-10); else
if( newt <= 0.676 +0.004 ) return  5.95700000e-02   * (-10); else
if( newt <= 0.68  +0.004 ) return  6.11000000e-02   * (-10); else
if( newt <= 0.684 +0.004 ) return  6.11000000e-02   * (-10); else
if( newt <= 0.688 +0.004 ) return  5.95700000e-02   * (-10); else
if( newt <= 0.692 +0.004 ) return  5.88030000e-02   * (-10); else
if( newt <= 0.696 +0.004 ) return  5.80350000e-02   * (-10); else
if( newt <= 0.7   +0.004 ) return  6.03370000e-02   * (-10); else
if( newt <= 0.704 +0.004 ) return  5.95700000e-02   * (-10); else
if( newt <= 0.708 +0.004 ) return  5.80350000e-02   * (-10); else
if( newt <= 0.712 +0.004 ) return  5.42030000e-02   * (-10); else
if( newt <= 0.716 +0.004 ) return  5.19010000e-02   * (-10); else
if( newt <= 0.72  +0.004 ) return  4.65350000e-02   * (-10); else
if( newt <= 0.724 +0.004 ) return  4.50000000e-02   * (-10); else
if( newt <= 0.728 +0.004 ) return  4.19350000e-02   * (-10); else
if( newt <= 0.732 +0.004 ) return  4.26980000e-02   * (-10); else
if( newt <= 0.736 +0.004 ) return  4.42330000e-02   * (-10); else
if( newt <= 0.74  +0.004 ) return  4.80690000e-02   * (-10); else
if( newt <= 0.744 +0.004 ) return  4.42330000e-02   * (-10); else
if( newt <= 0.748 +0.004 ) return  4.11680000e-02   * (-10); else
if( newt <= 0.752 +0.004 ) return  3.65640000e-02   * (-10); else
if( newt <= 0.756 +0.004 ) return  3.58010000e-02   * (-10); else
if( newt <= 0.76  +0.004 ) return  3.42670000e-02   * (-10); else
if( newt <= 0.764 +0.004 ) return  3.27320000e-02   * (-10); else
if( newt <= 0.768 +0.004 ) return  3.11970000e-02   * (-10); else
if( newt <= 0.772 +0.004 ) return  3.04300000e-02   * (-10); else
if( newt <= 0.776 +0.004 ) return  2.96670000e-02   * (-10); else
if( newt <= 0.78  +0.004 ) return  3.04300000e-02   * (-10); else
if( newt <= 0.784 +0.004 ) return  2.96670000e-02   * (-10); else
if( newt <= 0.788 +0.004 ) return  2.50630000e-02   * (-10); else
if( newt <= 0.792 +0.004 ) return  2.27660000e-02   * (-10); else
if( newt <= 0.796 +0.004 ) return  2.12310000e-02   * (-10); else
if( newt <= 0.8   +0.004 ) return  2.23800000e-02   * (-10);
}

Real linearPopliteal( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if(t>0.8)
        t=(((int)floor(t*1000.0))%800)/1000.0;

    Real dt=((int)floor(t*1000.0))%4/1000.0;
    Real ti=floor(t*1000.0)/1000.0 -dt;
    Real tii=ti+0.004;
    return (popliteal(tii)-popliteal(ti))/(0.004)*(t-(ti))+popliteal(ti);
}

Real linearPontdist( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if(t>0.8)
        t=(((int)floor(t*1000.0))%800)/1000.0;

    Real dt=((int)floor(t*1000.0))%4/1000.0;
    Real ti=floor(t*1000.0)/1000.0 -dt;
    Real tii=ti+0.004;
    return (pont_dist(tii)-pont_dist(ti))/(0.004)*(t-(ti))+pont_dist(ti);
}


Real poplitealPressure( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{   Real newt = ((Real)(((int)(t*1000)) % 800))/1000;
if( newt <=          0 +0.001 ) return          0   * (1.0e4); else
if( newt <=     0.0010 +0.001 ) return    -0.2595   * (1.0e4); else
if( newt <=     0.0020 +0.001 ) return    -0.1949   * (1.0e4); else
if( newt <=     0.0030 +0.001 ) return    -0.1639   * (1.0e4); else
if( newt <=     0.0040 +0.001 ) return    -0.1380   * (1.0e4); else
if( newt <=     0.0050 +0.001 ) return    -0.1340   * (1.0e4); else
if( newt <=     0.0060 +0.001 ) return    -0.1268   * (1.0e4); else
if( newt <=     0.0070 +0.001 ) return    -0.1215   * (1.0e4); else
if( newt <=     0.0080 +0.001 ) return    -0.1210   * (1.0e4); else
if( newt <=     0.0090 +0.001 ) return    -0.1411   * (1.0e4); else
if( newt <=     0.0100 +0.001 ) return    -0.1653   * (1.0e4); else
if( newt <=     0.0110 +0.001 ) return    -0.1953   * (1.0e4); else
if( newt <=     0.0120 +0.001 ) return    -0.2301   * (1.0e4); else
if( newt <=     0.0130 +0.001 ) return    -0.2657   * (1.0e4); else
if( newt <=     0.0140 +0.001 ) return    -0.3023   * (1.0e4); else
if( newt <=     0.0150 +0.001 ) return    -0.3377   * (1.0e4); else
if( newt <=     0.0160 +0.001 ) return    -0.3705   * (1.0e4); else
if( newt <=     0.0170 +0.001 ) return    -0.3932   * (1.0e4); else
if( newt <=     0.0180 +0.001 ) return    -0.4152   * (1.0e4); else
if( newt <=     0.0190 +0.001 ) return    -0.4363   * (1.0e4); else
if( newt <=     0.0200 +0.001 ) return    -0.4580   * (1.0e4); else
if( newt <=     0.0210 +0.001 ) return    -0.4829   * (1.0e4); else
if( newt <=     0.0220 +0.001 ) return    -0.5100   * (1.0e4); else
if( newt <=     0.0230 +0.001 ) return    -0.5399   * (1.0e4); else
if( newt <=     0.0240 +0.001 ) return    -0.5725   * (1.0e4); else
if( newt <=     0.0250 +0.001 ) return    -0.6086   * (1.0e4); else
if( newt <=     0.0260 +0.001 ) return    -0.6456   * (1.0e4); else
if( newt <=     0.0270 +0.001 ) return    -0.6830   * (1.0e4); else
if( newt <=     0.0280 +0.001 ) return    -0.7199   * (1.0e4); else
if( newt <=     0.0290 +0.001 ) return    -0.7515   * (1.0e4); else
if( newt <=     0.0300 +0.001 ) return    -0.7830   * (1.0e4); else
if( newt <=     0.0310 +0.001 ) return    -0.8136   * (1.0e4); else
if( newt <=     0.0320 +0.001 ) return    -0.8436   * (1.0e4); else
if( newt <=     0.0330 +0.001 ) return    -0.8720   * (1.0e4); else
if( newt <=     0.0340 +0.001 ) return    -0.9012   * (1.0e4); else
if( newt <=     0.0350 +0.001 ) return    -0.9312   * (1.0e4); else
if( newt <=     0.0360 +0.001 ) return    -0.9621   * (1.0e4); else
if( newt <=     0.0370 +0.001 ) return    -0.9926   * (1.0e4); else
if( newt <=     0.0380 +0.001 ) return    -1.0240   * (1.0e4); else
if( newt <=     0.0390 +0.001 ) return    -1.0559   * (1.0e4); else
if( newt <=     0.0400 +0.001 ) return    -1.2131   * (1.0e4); else
if( newt <=     0.0410 +0.001 ) return    -1.1907   * (1.0e4); else
if( newt <=     0.0420 +0.001 ) return    -1.1903   * (1.0e4); else
if( newt <=     0.0430 +0.001 ) return    -1.1610   * (1.0e4); else
if( newt <=     0.0440 +0.001 ) return    -1.2035   * (1.0e4); else
if( newt <=     0.0450 +0.001 ) return    -1.2322   * (1.0e4); else
if( newt <=     0.0460 +0.001 ) return    -1.2613   * (1.0e4); else
if( newt <=     0.0470 +0.001 ) return    -1.2932   * (1.0e4); else
if( newt <=     0.0480 +0.001 ) return    -1.3292   * (1.0e4); else
if( newt <=     0.0490 +0.001 ) return    -1.3769   * (1.0e4); else
if( newt <=     0.0500 +0.001 ) return    -1.4237   * (1.0e4); else
if( newt <=     0.0510 +0.001 ) return    -1.4686   * (1.0e4); else
if( newt <=     0.0520 +0.001 ) return    -1.5097   * (1.0e4); else
if( newt <=     0.0530 +0.001 ) return    -1.5436   * (1.0e4); else
if( newt <=     0.0540 +0.001 ) return    -1.5738   * (1.0e4); else
if( newt <=     0.0550 +0.001 ) return    -1.6007   * (1.0e4); else
if( newt <=     0.0560 +0.001 ) return    -1.7456   * (1.0e4); else
if( newt <=     0.0570 +0.001 ) return    -1.7042   * (1.0e4); else
if( newt <=     0.0580 +0.001 ) return    -1.6879   * (1.0e4); else
if( newt <=     0.0590 +0.001 ) return    -1.6474   * (1.0e4); else
if( newt <=     0.0600 +0.001 ) return    -1.6796   * (1.0e4); else
if( newt <=     0.0610 +0.001 ) return    -1.6954   * (1.0e4); else
if( newt <=     0.0620 +0.001 ) return    -1.7171   * (1.0e4); else
if( newt <=     0.0630 +0.001 ) return    -1.7454   * (1.0e4); else
if( newt <=     0.0640 +0.001 ) return    -1.7803   * (1.0e4); else
if( newt <=     0.0650 +0.001 ) return    -1.8083   * (1.0e4); else
if( newt <=     0.0660 +0.001 ) return    -1.8407   * (1.0e4); else
if( newt <=     0.0670 +0.001 ) return    -1.8726   * (1.0e4); else
if( newt <=     0.0680 +0.001 ) return    -1.9017   * (1.0e4); else
if( newt <=     0.0690 +0.001 ) return    -1.9309   * (1.0e4); else
if( newt <=     0.0700 +0.001 ) return    -1.9547   * (1.0e4); else
if( newt <=     0.0710 +0.001 ) return    -1.9744   * (1.0e4); else
if( newt <=     0.0720 +0.001 ) return    -1.9910   * (1.0e4); else
if( newt <=     0.0730 +0.001 ) return    -2.0102   * (1.0e4); else
if( newt <=     0.0740 +0.001 ) return    -2.0277   * (1.0e4); else
if( newt <=     0.0750 +0.001 ) return    -2.0449   * (1.0e4); else
if( newt <=     0.0760 +0.001 ) return    -2.0622   * (1.0e4); else
if( newt <=     0.0770 +0.001 ) return    -2.0840   * (1.0e4); else
if( newt <=     0.0780 +0.001 ) return    -2.1044   * (1.0e4); else
if( newt <=     0.0790 +0.001 ) return    -2.1239   * (1.0e4); else
if( newt <=     0.0800 +0.001 ) return    -2.1421   * (1.0e4); else
if( newt <=     0.0810 +0.001 ) return    -2.1545   * (1.0e4); else
if( newt <=     0.0820 +0.001 ) return    -2.16662  * (1.0e4); else
if( newt <=     0.0830 +0.001 ) return    -2.17782  * (1.0e4); else
if( newt <=     0.0840 +0.001 ) return    -2.18851  * (1.0e4); else
if( newt <=     0.0850 +0.001 ) return    -2.19451  * (1.0e4); else
if( newt <=     0.0860 +0.001 ) return    -2.20121  * (1.0e4); else
if( newt <=     0.0870 +0.001 ) return    -2.20831  * (1.0e4); else
if( newt <=     0.0880 +0.001 ) return    -2.21571  * (1.0e4); else
if( newt <=     0.0890 +0.001 ) return    -2.22621  * (1.0e4); else
if( newt <=     0.0900 +0.001 ) return    -2.23631  * (1.0e4); else
if( newt <=     0.0910 +0.001 ) return    -2.24611  * (1.0e4); else
if( newt <=     0.0920 +0.001 ) return    -2.25551  * (1.0e4); else
if( newt <=     0.0930 +0.001 ) return    -2.25891  * (1.0e4); else
if( newt <=     0.0940 +0.001 ) return    -2.26351  * (1.0e4); else
if( newt <=     0.0950 +0.001 ) return    -2.26861  * (1.0e4); else
if( newt <=     0.0960 +0.001 ) return    -2.27411  * (1.0e4); else
if( newt <=     0.0970 +0.001 ) return    -2.27451  * (1.0e4); else
if( newt <=     0.0980 +0.001 ) return    -2.27681  * (1.0e4); else
if( newt <=     0.0990 +0.001 ) return    -2.28041  * (1.0e4); else
if( newt <=     0.1000 +0.001 ) return    -2.28481  * (1.0e4); else
if( newt <=     0.1010 +0.001 ) return    -2.29021  * (1.0e4); else
if( newt <=     0.1020 +0.001 ) return    -2.29621  * (1.0e4); else
if( newt <=     0.1030 +0.001 ) return    -2.30261  * (1.0e4); else
if( newt <=     0.1040 +0.001 ) return    -2.30941  * (1.0e4); else
if( newt <=     0.1050 +0.001 ) return    -2.32491  * (1.0e4); else
if( newt <=     0.1060 +0.001 ) return    -2.33822  * (1.0e4); else
if( newt <=     0.1070 +0.001 ) return    -2.35052  * (1.0e4); else
if( newt <=     0.1080 +0.001 ) return    -2.36182  * (1.0e4); else
if( newt <=     0.1090 +0.001 ) return    -2.37172  * (1.0e4); else
if( newt <=     0.1100 +0.001 ) return    -2.38012  * (1.0e4); else
if( newt <=     0.1110 +0.001 ) return    -2.38712  * (1.0e4); else
if( newt <=     0.1120 +0.001 ) return    -2.3926   * (1.0e4); else
if( newt <=     0.1130 +0.001 ) return    -2.3984   * (1.0e4); else
if( newt <=     0.1140 +0.001 ) return    -2.4031   * (1.0e4); else
if( newt <=     0.1150 +0.001 ) return    -2.4070   * (1.0e4); else
if( newt <=     0.1160 +0.001 ) return    -2.4106   * (1.0e4); else
if( newt <=     0.1170 +0.001 ) return    -2.3952   * (1.0e4); else
if( newt <=     0.1180 +0.001 ) return    -2.3841   * (1.0e4); else
if( newt <=     0.1190 +0.001 ) return    -2.3751   * (1.0e4); else
if( newt <=     0.1200 +0.001 ) return    -2.3675   * (1.0e4); else
if( newt <=     0.1210 +0.001 ) return    -2.3569   * (1.0e4); else
if( newt <=     0.1220 +0.001 ) return    -2.3482   * (1.0e4); else
if( newt <=     0.1230 +0.001 ) return    -2.3403   * (1.0e4); else
if( newt <=     0.1240 +0.001 ) return    -2.3329   * (1.0e4); else
if( newt <=     0.1250 +0.001 ) return    -2.3284   * (1.0e4); else
if( newt <=     0.1260 +0.001 ) return    -2.3228   * (1.0e4); else
if( newt <=     0.1270 +0.001 ) return    -2.3163   * (1.0e4); else
if( newt <=     0.1280 +0.001 ) return    -2.3086   * (1.0e4); else
if( newt <=     0.1290 +0.001 ) return    -2.3155   * (1.0e4); else
if( newt <=     0.1300 +0.001 ) return    -2.3172   * (1.0e4); else
if( newt <=     0.1310 +0.001 ) return    -2.3160   * (1.0e4); else
if( newt <=     0.1320 +0.001 ) return    -2.3125   * (1.0e4); else
if( newt <=     0.1330 +0.001 ) return    -2.3087   * (1.0e4); else
if( newt <=     0.1340 +0.001 ) return    -2.3033   * (1.0e4); else
if( newt <=     0.1350 +0.001 ) return    -2.2973   * (1.0e4); else
if( newt <=     0.1360 +0.001 ) return    -2.2912   * (1.0e4); else
if( newt <=     0.1370 +0.001 ) return    -2.2752   * (1.0e4); else
if( newt <=     0.1380 +0.001 ) return    -2.2620   * (1.0e4); else
if( newt <=     0.1390 +0.001 ) return    -2.2505   * (1.0e4); else
if( newt <=     0.1400 +0.001 ) return    -2.2405   * (1.0e4); else
if( newt <=     0.1410 +0.001 ) return    -2.2216   * (1.0e4); else
if( newt <=     0.1420 +0.001 ) return    -2.2060   * (1.0e4); else
if( newt <=     0.1430 +0.001 ) return    -2.1920   * (1.0e4); else
if( newt <=     0.1440 +0.001 ) return    -2.1790   * (1.0e4); else
if( newt <=     0.1450 +0.001 ) return    -2.1681   * (1.0e4); else
if( newt <=     0.1460 +0.001 ) return    -2.1569   * (1.0e4); else
if( newt <=     0.1470 +0.001 ) return    -2.1452   * (1.0e4); else
if( newt <=     0.1480 +0.001 ) return    -2.1326   * (1.0e4); else
if( newt <=     0.1490 +0.001 ) return    -2.1381   * (1.0e4); else
if( newt <=     0.1500 +0.001 ) return    -2.1382   * (1.0e4); else
if( newt <=     0.1510 +0.001 ) return    -2.1357   * (1.0e4); else
if( newt <=     0.1520 +0.001 ) return    -2.1310   * (1.0e4); else
if( newt <=     0.1530 +0.001 ) return    -2.1493   * (1.0e4); else
if( newt <=     0.1540 +0.001 ) return    -2.1608   * (1.0e4); else
if( newt <=     0.1550 +0.001 ) return    -2.1695   * (1.0e4); else
if( newt <=     0.1560 +0.001 ) return    -2.1763   * (1.0e4); else
if( newt <=     0.1570 +0.001 ) return    -2.2154   * (1.0e4); else
if( newt <=     0.1580 +0.001 ) return    -2.2463   * (1.0e4); else
if( newt <=     0.1590 +0.001 ) return    -2.2742   * (1.0e4); else
if( newt <=     0.1600 +0.001 ) return    -2.3002   * (1.0e4); else
if( newt <=     0.1610 +0.001 ) return    -2.3619   * (1.0e4); else
if( newt <=     0.1620 +0.001 ) return    -2.4153   * (1.0e4); else
if( newt <=     0.1630 +0.001 ) return    -2.4664   * (1.0e4); else
if( newt <=     0.1640 +0.001 ) return    -2.5167   * (1.0e4); else
if( newt <=     0.1650 +0.001 ) return    -2.5589   * (1.0e4); else
if( newt <=     0.1660 +0.001 ) return    -2.6054   * (1.0e4); else
if( newt <=     0.1670 +0.001 ) return    -2.6563   * (1.0e4); else
if( newt <=     0.1680 +0.001 ) return    -2.7130   * (1.0e4); else
if( newt <=     0.1690 +0.001 ) return    -2.7059   * (1.0e4); else
if( newt <=     0.1700 +0.001 ) return    -2.7218   * (1.0e4); else
if( newt <=     0.1710 +0.001 ) return    -2.7521   * (1.0e4); else
if( newt <=     0.1720 +0.001 ) return    -2.7960   * (1.0e4); else
if( newt <=     0.1730 +0.001 ) return    -2.7344   * (1.0e4); else
if( newt <=     0.1740 +0.001 ) return    -2.7099   * (1.0e4); else
if( newt <=     0.1750 +0.001 ) return    -2.7057   * (1.0e4); else
if( newt <=     0.1760 +0.001 ) return    -2.7178   * (1.0e4); else
if( newt <=     0.1770 +0.001 ) return    -2.7231   * (1.0e4); else
if( newt <=     0.1780 +0.001 ) return    -2.7414   * (1.0e4); else
if( newt <=     0.1790 +0.001 ) return    -2.7661   * (1.0e4); else
if( newt <=     0.1800 +0.001 ) return    -2.7930   * (1.0e4); else
if( newt <=     0.1810 +0.001 ) return    -2.8229   * (1.0e4); else
if( newt <=     0.1820 +0.001 ) return    -2.8472   * (1.0e4); else
if( newt <=     0.1830 +0.001 ) return    -2.8642   * (1.0e4); else
if( newt <=     0.1840 +0.001 ) return    -2.8726   * (1.0e4); else
if( newt <=     0.1850 +0.001 ) return    -2.9007   * (1.0e4); else
if( newt <=     0.1860 +0.001 ) return    -2.9128   * (1.0e4); else
if( newt <=     0.1870 +0.001 ) return    -2.9136   * (1.0e4); else
if( newt <=     0.1880 +0.001 ) return    -2.9048   * (1.0e4); else
if( newt <=     0.1890 +0.001 ) return    -2.9170   * (1.0e4); else
if( newt <=     0.1900 +0.001 ) return    -2.9159   * (1.0e4); else
if( newt <=     0.1910 +0.001 ) return    -2.9069   * (1.0e4); else
if( newt <=     0.1920 +0.001 ) return    -2.8916   * (1.0e4); else
if( newt <=     0.1930 +0.001 ) return    -2.8918   * (1.0e4); else
if( newt <=     0.1940 +0.001 ) return    -2.8832   * (1.0e4); else
if( newt <=     0.1950 +0.001 ) return    -2.8695   * (1.0e4); else
if( newt <=     0.1960 +0.001 ) return    -2.8520   * (1.0e4); else
if( newt <=     0.1970 +0.001 ) return    -2.8502   * (1.0e4); else
if( newt <=     0.1980 +0.001 ) return    -2.8418   * (1.0e4); else
if( newt <=     0.1990 +0.001 ) return    -2.8301   * (1.0e4); else
if( newt <=     0.2000 +0.001 ) return    -2.8160   * (1.0e4); else
if( newt <=     0.2010 +0.001 ) return    -2.8076   * (1.0e4); else
if( newt <=     0.2020 +0.001 ) return    -2.7965   * (1.0e4); else
if( newt <=     0.2030 +0.001 ) return    -2.7843   * (1.0e4); else
if( newt <=     0.2040 +0.001 ) return    -2.7715   * (1.0e4); else
if( newt <=     0.2050 +0.001 ) return    -2.7587   * (1.0e4); else
if( newt <=     0.2060 +0.001 ) return    -2.7465   * (1.0e4); else
if( newt <=     0.2070 +0.001 ) return    -2.7351   * (1.0e4); else
if( newt <=     0.2080 +0.001 ) return    -2.7250   * (1.0e4); else
if( newt <=     0.2090 +0.001 ) return    -2.7147   * (1.0e4); else
if( newt <=     0.2100 +0.001 ) return    -2.7061   * (1.0e4); else
if( newt <=     0.2110 +0.001 ) return    -2.6990   * (1.0e4); else
if( newt <=     0.2120 +0.001 ) return    -2.6933   * (1.0e4); else
if( newt <=     0.2130 +0.001 ) return    -2.6668   * (1.0e4); else
if( newt <=     0.2140 +0.001 ) return    -2.6461   * (1.0e4); else
if( newt <=     0.2150 +0.001 ) return    -2.6279   * (1.0e4); else
if( newt <=     0.2160 +0.001 ) return    -2.6115   * (1.0e4); else
if( newt <=     0.2170 +0.001 ) return    -2.5895   * (1.0e4); else
if( newt <=     0.2180 +0.001 ) return    -2.5705   * (1.0e4); else
if( newt <=     0.2190 +0.001 ) return    -2.5531   * (1.0e4); else
if( newt <=     0.2200 +0.001 ) return    -2.5369   * (1.0e4); else
if( newt <=     0.2210 +0.001 ) return    -2.5157   * (1.0e4); else
if( newt <=     0.2220 +0.001 ) return    -2.4961   * (1.0e4); else
if( newt <=     0.2230 +0.001 ) return    -2.4770   * (1.0e4); else
if( newt <=     0.2240 +0.001 ) return    -2.4576   * (1.0e4); else
if( newt <=     0.2250 +0.001 ) return    -2.4391   * (1.0e4); else
if( newt <=     0.2260 +0.001 ) return    -2.4195   * (1.0e4); else
if( newt <=     0.2270 +0.001 ) return    -2.3990   * (1.0e4); else
if( newt <=     0.2280 +0.001 ) return    -2.3776   * (1.0e4); else
if( newt <=     0.2290 +0.001 ) return    -2.3639   * (1.0e4); else
if( newt <=     0.2300 +0.001 ) return    -2.3474   * (1.0e4); else
if( newt <=     0.2310 +0.001 ) return    -2.3297   * (1.0e4); else
if( newt <=     0.2320 +0.001 ) return    -2.3107   * (1.0e4); else
if( newt <=     0.2330 +0.001 ) return    -2.2894   * (1.0e4); else
if( newt <=     0.2340 +0.001 ) return    -2.2676   * (1.0e4); else
if( newt <=     0.2350 +0.001 ) return    -2.2453   * (1.0e4); else
if( newt <=     0.2360 +0.001 ) return    -2.2225   * (1.0e4); else
if( newt <=     0.2370 +0.001 ) return    -2.1966   * (1.0e4); else
if( newt <=     0.2380 +0.001 ) return    -2.1713   * (1.0e4); else
if( newt <=     0.2390 +0.001 ) return    -2.1463   * (1.0e4); else
if( newt <=     0.2400 +0.001 ) return    -2.1216   * (1.0e4); else
if( newt <=     0.2410 +0.001 ) return    -2.0871   * (1.0e4); else
if( newt <=     0.2420 +0.001 ) return    -2.0551   * (1.0e4); else
if( newt <=     0.2430 +0.001 ) return    -2.0239   * (1.0e4); else
if( newt <=     0.2440 +0.001 ) return    -1.9930   * (1.0e4); else
if( newt <=     0.2450 +0.001 ) return    -1.9580   * (1.0e4); else
if( newt <=     0.2460 +0.001 ) return    -1.9240   * (1.0e4); else
if( newt <=     0.2470 +0.001 ) return    -1.8902   * (1.0e4); else
if( newt <=     0.2480 +0.001 ) return    -1.8563   * (1.0e4); else
if( newt <=     0.2490 +0.001 ) return    -1.7995   * (1.0e4); else
if( newt <=     0.2500 +0.001 ) return    -1.7480   * (1.0e4); else
if( newt <=     0.2510 +0.001 ) return    -1.6984   * (1.0e4); else
if( newt <=     0.2520 +0.001 ) return    -1.6502   * (1.0e4); else
if( newt <=     0.2530 +0.001 ) return    -1.6212   * (1.0e4); else
if( newt <=     0.2540 +0.001 ) return    -1.5877   * (1.0e4); else
if( newt <=     0.2550 +0.001 ) return    -1.5517   * (1.0e4); else
if( newt <=     0.2560 +0.001 ) return    -1.5129   * (1.0e4); else
if( newt <=     0.2570 +0.001 ) return    -1.4965   * (1.0e4); else
if( newt <=     0.2580 +0.001 ) return    -1.4704   * (1.0e4); else
if( newt <=     0.2590 +0.001 ) return    -1.4382   * (1.0e4); else
if( newt <=     0.2600 +0.001 ) return    -1.4005   * (1.0e4); else
if( newt <=     0.2610 +0.001 ) return    -1.3592   * (1.0e4); else
if( newt <=     0.2620 +0.001 ) return    -1.3134   * (1.0e4); else
if( newt <=     0.2630 +0.001 ) return    -1.2641   * (1.0e4); else
if( newt <=     0.2640 +0.001 ) return    -1.2126   * (1.0e4); else
if( newt <=     0.2650 +0.001 ) return    -1.1473   * (1.0e4); else
if( newt <=     0.2660 +0.001 ) return    -1.0855   * (1.0e4); else
if( newt <=     0.2670 +0.001 ) return    -1.0260   * (1.0e4); else
if( newt <=     0.2680 +0.001 ) return    -0.9694   * (1.0e4); else
if( newt <=     0.2690 +0.001 ) return    -0.9197   * (1.0e4); else
if( newt <=     0.2700 +0.001 ) return    -0.8715   * (1.0e4); else
if( newt <=     0.2710 +0.001 ) return    -0.8246   * (1.0e4); else
if( newt <=     0.2720 +0.001 ) return    -0.7783   * (1.0e4); else
if( newt <=     0.2730 +0.001 ) return    -0.7117   * (1.0e4); else
if( newt <=     0.2740 +0.001 ) return    -0.6489   * (1.0e4); else
if( newt <=     0.2750 +0.001 ) return    -0.5866   * (1.0e4); else
if( newt <=     0.2760 +0.001 ) return    -0.5240   * (1.0e4); else
if( newt <=     0.2770 +0.001 ) return    -0.4540   * (1.0e4); else
if( newt <=     0.2780 +0.001 ) return    -0.3857   * (1.0e4); else
if( newt <=     0.2790 +0.001 ) return    -0.3182   * (1.0e4); else
if( newt <=     0.2800 +0.001 ) return    -0.2513   * (1.0e4); else
if( newt <=     0.2810 +0.001 ) return    -0.2075   * (1.0e4); else
if( newt <=     0.2820 +0.001 ) return    -0.15802  * (1.0e4); else
if( newt <=     0.2830 +0.001 ) return    -0.10562  * (1.0e4); else
if( newt <=     0.2840 +0.001 ) return    -0.05011  * (1.0e4); else
if( newt <=     0.2850 +0.001 ) return     0.00771  * (1.0e4); else
if( newt <=     0.2860 +0.001 ) return     0.06971  * (1.0e4); else
if( newt <=     0.2870 +0.001 ) return     0.13571  * (1.0e4); else
if( newt <=     0.2880 +0.001 ) return     0.20521  * (1.0e4); else
if( newt <=     0.2890 +0.001 ) return     0.29581  * (1.0e4); else
if( newt <=     0.2900 +0.001 ) return     0.38341  * (1.0e4); else
if( newt <=     0.2910 +0.001 ) return     0.46981  * (1.0e4); else
if( newt <=     0.2920 +0.001 ) return     0.55461  * (1.0e4); else
if( newt <=     0.2930 +0.001 ) return     0.60371  * (1.0e4); else
if( newt <=     0.2940 +0.001 ) return     0.66001  * (1.0e4); else
if( newt <=     0.2950 +0.001 ) return     0.71861  * (1.0e4); else
if( newt <=     0.2960 +0.001 ) return     0.77981  * (1.0e4); else
if( newt <=     0.2970 +0.001 ) return     0.83161  * (1.0e4); else
if( newt <=     0.2980 +0.001 ) return     0.89061  * (1.0e4); else
if( newt <=     0.2990 +0.001 ) return     0.95491  * (1.0e4); else
if( newt <=     0.3000 +0.001 ) return     1.02421  * (1.0e4); else
if( newt <=     0.3010 +0.001 ) return     1.08901  * (1.0e4); else
if( newt <=     0.3020 +0.001 ) return     1.15901  * (1.0e4); else
if( newt <=     0.3030 +0.001 ) return     1.23141  * (1.0e4); else
if( newt <=     0.3040 +0.001 ) return     1.30471  * (1.0e4); else
if( newt <=     0.3050 +0.001 ) return     1.38741  * (1.0e4); else
if( newt <=     0.3060 +0.001 ) return     1.46602  * (1.0e4); else
if( newt <=     0.3070 +0.001 ) return     1.54152  * (1.0e4); else
if( newt <=     0.3080 +0.001 ) return     1.61382  * (1.0e4); else
if( newt <=     0.3090 +0.001 ) return     1.68312  * (1.0e4); else
if( newt <=     0.3100 +0.001 ) return     1.74982  * (1.0e4); else
if( newt <=     0.3110 +0.001 ) return     1.81452  * (1.0e4); else
if( newt <=     0.3120 +0.001 ) return     1.8778   * (1.0e4); else
if( newt <=     0.3130 +0.001 ) return     1.9319   * (1.0e4); else
if( newt <=     0.3140 +0.001 ) return     1.9880   * (1.0e4); else
if( newt <=     0.3150 +0.001 ) return     2.0449   * (1.0e4); else
if( newt <=     0.3160 +0.001 ) return     2.1028   * (1.0e4); else
if( newt <=     0.3170 +0.001 ) return     2.1629   * (1.0e4); else
if( newt <=     0.3180 +0.001 ) return     2.2232   * (1.0e4); else
if( newt <=     0.3190 +0.001 ) return     2.2834   * (1.0e4); else
if( newt <=     0.3200 +0.001 ) return     2.3434   * (1.0e4); else
if( newt <=     0.3210 +0.001 ) return     2.4023   * (1.0e4); else
if( newt <=     0.3220 +0.001 ) return     2.4598   * (1.0e4); else
if( newt <=     0.3230 +0.001 ) return     2.5156   * (1.0e4); else
if( newt <=     0.3240 +0.001 ) return     2.5693   * (1.0e4); else
if( newt <=     0.3250 +0.001 ) return     2.6180   * (1.0e4); else
if( newt <=     0.3260 +0.001 ) return     2.6651   * (1.0e4); else
if( newt <=     0.3270 +0.001 ) return     2.7099   * (1.0e4); else
if( newt <=     0.3280 +0.001 ) return     2.7522   * (1.0e4); else
if( newt <=     0.3290 +0.001 ) return     2.7945   * (1.0e4); else
if( newt <=     0.3300 +0.001 ) return     2.8325   * (1.0e4); else
if( newt <=     0.3310 +0.001 ) return     2.8664   * (1.0e4); else
if( newt <=     0.3320 +0.001 ) return     2.8961   * (1.0e4); else
if( newt <=     0.3330 +0.001 ) return     2.9207   * (1.0e4); else
if( newt <=     0.3340 +0.001 ) return     2.9428   * (1.0e4); else
if( newt <=     0.3350 +0.001 ) return     2.9629   * (1.0e4); else
if( newt <=     0.3360 +0.001 ) return     2.9814   * (1.0e4); else
if( newt <=     0.3370 +0.001 ) return     3.0024   * (1.0e4); else
if( newt <=     0.3380 +0.001 ) return     3.0201   * (1.0e4); else
if( newt <=     0.3390 +0.001 ) return     3.0343   * (1.0e4); else
if( newt <=     0.3400 +0.001 ) return     3.0438   * (1.0e4); else
if( newt <=     0.3410 +0.001 ) return     3.0450   * (1.0e4); else
if( newt <=     0.3420 +0.001 ) return     3.0401   * (1.0e4); else
if( newt <=     0.3430 +0.001 ) return     3.0280   * (1.0e4); else
if( newt <=     0.3440 +0.001 ) return     3.0092   * (1.0e4); else
if( newt <=     0.3450 +0.001 ) return     2.9883   * (1.0e4); else
if( newt <=     0.3460 +0.001 ) return     2.9642   * (1.0e4); else
if( newt <=     0.3470 +0.001 ) return     2.9395   * (1.0e4); else
if( newt <=     0.3480 +0.001 ) return     2.9157   * (1.0e4); else
if( newt <=     0.3490 +0.001 ) return     2.8961   * (1.0e4); else
if( newt <=     0.3500 +0.001 ) return     2.8770   * (1.0e4); else
if( newt <=     0.3510 +0.001 ) return     2.8581   * (1.0e4); else
if( newt <=     0.3520 +0.001 ) return     2.8384   * (1.0e4); else
if( newt <=     0.3530 +0.001 ) return     2.8114   * (1.0e4); else
if( newt <=     0.3540 +0.001 ) return     2.7836   * (1.0e4); else
if( newt <=     0.3550 +0.001 ) return     2.7540   * (1.0e4); else
if( newt <=     0.3560 +0.001 ) return     2.7225   * (1.0e4); else
if( newt <=     0.3570 +0.001 ) return     2.6901   * (1.0e4); else
if( newt <=     0.3580 +0.001 ) return     2.6578   * (1.0e4); else
if( newt <=     0.3590 +0.001 ) return     2.6268   * (1.0e4); else
if( newt <=     0.3600 +0.001 ) return     2.5974   * (1.0e4); else
if( newt <=     0.3610 +0.001 ) return     2.5780   * (1.0e4); else
if( newt <=     0.3620 +0.001 ) return     2.5576   * (1.0e4); else
if( newt <=     0.3630 +0.001 ) return     2.5370   * (1.0e4); else
if( newt <=     0.3640 +0.001 ) return     2.5156   * (1.0e4); else
if( newt <=     0.3650 +0.001 ) return     2.4947   * (1.0e4); else
if( newt <=     0.3660 +0.001 ) return     2.4724   * (1.0e4); else
if( newt <=     0.3670 +0.001 ) return     2.4493   * (1.0e4); else
if( newt <=     0.3680 +0.001 ) return     2.4258   * (1.0e4); else
if( newt <=     0.3690 +0.001 ) return     2.3965   * (1.0e4); else
if( newt <=     0.3700 +0.001 ) return     2.3691   * (1.0e4); else
if( newt <=     0.3710 +0.001 ) return     2.3429   * (1.0e4); else
if( newt <=     0.3720 +0.001 ) return     2.3181   * (1.0e4); else
if( newt <=     0.3730 +0.001 ) return     2.2823   * (1.0e4); else
if( newt <=     0.3740 +0.001 ) return     2.2514   * (1.0e4); else
if( newt <=     0.3750 +0.001 ) return     2.2231   * (1.0e4); else
if( newt <=     0.3760 +0.001 ) return     2.1969   * (1.0e4); else
if( newt <=     0.3770 +0.001 ) return     2.1722   * (1.0e4); else
if( newt <=     0.3780 +0.001 ) return     2.1484   * (1.0e4); else
if( newt <=     0.3790 +0.001 ) return     2.1245   * (1.0e4); else
if( newt <=     0.3800 +0.001 ) return     2.1001   * (1.0e4); else
if( newt <=     0.3810 +0.001 ) return     2.0776   * (1.0e4); else
if( newt <=     0.3820 +0.001 ) return     2.0531   * (1.0e4); else
if( newt <=     0.3830 +0.001 ) return     2.0271   * (1.0e4); else
if( newt <=     0.3840 +0.001 ) return     1.9996   * (1.0e4); else
if( newt <=     0.3850 +0.001 ) return     1.9708   * (1.0e4); else
if( newt <=     0.3860 +0.001 ) return     1.9408   * (1.0e4); else
if( newt <=     0.3870 +0.001 ) return     1.9099   * (1.0e4); else
if( newt <=     0.3880 +0.001 ) return     1.8783   * (1.0e4); else
if( newt <=     0.3890 +0.001 ) return     1.8422   * (1.0e4); else
if( newt <=     0.3900 +0.001 ) return     1.8070   * (1.0e4); else
if( newt <=     0.3910 +0.001 ) return     1.7723   * (1.0e4); else
if( newt <=     0.3920 +0.001 ) return     1.7380   * (1.0e4); else
if( newt <=     0.3930 +0.001 ) return     1.7041   * (1.0e4); else
if( newt <=     0.3940 +0.001 ) return     1.6704   * (1.0e4); else
if( newt <=     0.3950 +0.001 ) return     1.6368   * (1.0e4); else
if( newt <=     0.3960 +0.001 ) return     1.6027   * (1.0e4); else
if( newt <=     0.3970 +0.001 ) return     1.5765   * (1.0e4); else
if( newt <=     0.3980 +0.001 ) return     1.5471   * (1.0e4); else
if( newt <=     0.3990 +0.001 ) return     1.5157   * (1.0e4); else
if( newt <=     0.4000 +0.001 ) return     1.4824   * (1.0e4); else
if( newt <=     0.4010 +0.001 ) return     1.4560   * (1.0e4); else
if( newt <=     0.4020 +0.001 ) return     1.4260   * (1.0e4); else
if( newt <=     0.4030 +0.001 ) return     1.3944   * (1.0e4); else
if( newt <=     0.4040 +0.001 ) return     1.3617   * (1.0e4); else
if( newt <=     0.4050 +0.001 ) return     1.3271   * (1.0e4); else
if( newt <=     0.4060 +0.001 ) return     1.2930   * (1.0e4); else
if( newt <=     0.4070 +0.001 ) return     1.2593   * (1.0e4); else
if( newt <=     0.4080 +0.001 ) return     1.2264   * (1.0e4); else
if( newt <=     0.4090 +0.001 ) return     1.1862   * (1.0e4); else
if( newt <=     0.4100 +0.001 ) return     1.1493   * (1.0e4); else
if( newt <=     0.4110 +0.001 ) return     1.1147   * (1.0e4); else
if( newt <=     0.4120 +0.001 ) return     1.0823   * (1.0e4); else
if( newt <=     0.4130 +0.001 ) return     1.0478   * (1.0e4); else
if( newt <=     0.4140 +0.001 ) return     1.0168   * (1.0e4); else
if( newt <=     0.4150 +0.001 ) return     0.9883   * (1.0e4); else
if( newt <=     0.4160 +0.001 ) return     0.9621   * (1.0e4); else
if( newt <=     0.4170 +0.001 ) return     0.9347   * (1.0e4); else
if( newt <=     0.4180 +0.001 ) return     0.9092   * (1.0e4); else
if( newt <=     0.4190 +0.001 ) return     0.8848   * (1.0e4); else
if( newt <=     0.4200 +0.001 ) return     0.8611   * (1.0e4); else
if( newt <=     0.4210 +0.001 ) return     0.8372   * (1.0e4); else
if( newt <=     0.4220 +0.001 ) return     0.8156   * (1.0e4); else
if( newt <=     0.4230 +0.001 ) return     0.7968   * (1.0e4); else
if( newt <=     0.4240 +0.001 ) return     0.7811   * (1.0e4); else
if( newt <=     0.4250 +0.001 ) return     0.7648   * (1.0e4); else
if( newt <=     0.4260 +0.001 ) return     0.7525   * (1.0e4); else
if( newt <=     0.4270 +0.001 ) return     0.7431   * (1.0e4); else
if( newt <=     0.4280 +0.001 ) return     0.7359   * (1.0e4); else
if( newt <=     0.4290 +0.001 ) return     0.7315   * (1.0e4); else
if( newt <=     0.4300 +0.001 ) return     0.7271   * (1.0e4); else
if( newt <=     0.4310 +0.001 ) return     0.7222   * (1.0e4); else
if( newt <=     0.4320 +0.001 ) return     0.7166   * (1.0e4); else
if( newt <=     0.4330 +0.001 ) return     0.7158   * (1.0e4); else
if( newt <=     0.4340 +0.001 ) return     0.7130   * (1.0e4); else
if( newt <=     0.4350 +0.001 ) return     0.7092   * (1.0e4); else
if( newt <=     0.4360 +0.001 ) return     0.7047   * (1.0e4); else
if( newt <=     0.4370 +0.001 ) return     0.7121   * (1.0e4); else
if( newt <=     0.4380 +0.001 ) return     0.7147   * (1.0e4); else
if( newt <=     0.4390 +0.001 ) return     0.7143   * (1.0e4); else
if( newt <=     0.4400 +0.001 ) return     0.7113   * (1.0e4); else
if( newt <=     0.4410 +0.001 ) return     0.7017   * (1.0e4); else
if( newt <=     0.4420 +0.001 ) return     0.6919   * (1.0e4); else
if( newt <=     0.4430 +0.001 ) return     0.6821   * (1.0e4); else
if( newt <=     0.4440 +0.001 ) return     0.6731   * (1.0e4); else
if( newt <=     0.4450 +0.001 ) return     0.6662   * (1.0e4); else
if( newt <=     0.4460 +0.001 ) return     0.6611   * (1.0e4); else
if( newt <=     0.4470 +0.001 ) return     0.6585   * (1.0e4); else
if( newt <=     0.4480 +0.001 ) return     0.6582   * (1.0e4); else
if( newt <=     0.4490 +0.001 ) return     0.6559   * (1.0e4); else
if( newt <=     0.4500 +0.001 ) return     0.6572   * (1.0e4); else
if( newt <=     0.4510 +0.001 ) return     0.6613   * (1.0e4); else
if( newt <=     0.4520 +0.001 ) return     0.6679   * (1.0e4); else
if( newt <=     0.4530 +0.001 ) return     0.6770   * (1.0e4); else
if( newt <=     0.4540 +0.001 ) return     0.6883   * (1.0e4); else
if( newt <=     0.4550 +0.001 ) return     0.7018   * (1.0e4); else
if( newt <=     0.4560 +0.001 ) return     0.7170   * (1.0e4); else
if( newt <=     0.4570 +0.001 ) return     0.7279   * (1.0e4); else
if( newt <=     0.4580 +0.001 ) return     0.7396   * (1.0e4); else
if( newt <=     0.4590 +0.001 ) return     0.7502   * (1.0e4); else
if( newt <=     0.4600 +0.001 ) return     0.7594   * (1.0e4); else
if( newt <=     0.4610 +0.001 ) return     0.7671   * (1.0e4); else
if( newt <=     0.4620 +0.001 ) return     0.7735   * (1.0e4); else
if( newt <=     0.4630 +0.001 ) return     0.7790   * (1.0e4); else
if( newt <=     0.4640 +0.001 ) return     0.7842   * (1.0e4); else
if( newt <=     0.4650 +0.001 ) return     0.7864   * (1.0e4); else
if( newt <=     0.4660 +0.001 ) return     0.7893   * (1.0e4); else
if( newt <=     0.4670 +0.001 ) return     0.7925   * (1.0e4); else
if( newt <=     0.4680 +0.001 ) return     0.7959   * (1.0e4); else
if( newt <=     0.4690 +0.001 ) return     0.7963   * (1.0e4); else
if( newt <=     0.4700 +0.001 ) return     0.7973   * (1.0e4); else
if( newt <=     0.4710 +0.001 ) return     0.7983   * (1.0e4); else
if( newt <=     0.4720 +0.001 ) return     0.7989   * (1.0e4); else
if( newt <=     0.4730 +0.001 ) return     0.7977   * (1.0e4); else
if( newt <=     0.4740 +0.001 ) return     0.7960   * (1.0e4); else
if( newt <=     0.4750 +0.001 ) return     0.7937   * (1.0e4); else
if( newt <=     0.4760 +0.001 ) return     0.7905   * (1.0e4); else
if( newt <=     0.4770 +0.001 ) return     0.7949   * (1.0e4); else
if( newt <=     0.4780 +0.001 ) return     0.7963   * (1.0e4); else
if( newt <=     0.4790 +0.001 ) return     0.7959   * (1.0e4); else
if( newt <=     0.4800 +0.001 ) return     0.7940   * (1.0e4); else
if( newt <=     0.4810 +0.001 ) return     0.7962   * (1.0e4); else
if( newt <=     0.4820 +0.001 ) return     0.79562  * (1.0e4); else
if( newt <=     0.4830 +0.001 ) return     0.79342  * (1.0e4); else
if( newt <=     0.4840 +0.001 ) return     0.78971  * (1.0e4); else
if( newt <=     0.4850 +0.001 ) return     0.79061  * (1.0e4); else
if( newt <=     0.4860 +0.001 ) return     0.78931  * (1.0e4); else
if( newt <=     0.4870 +0.001 ) return     0.78711  * (1.0e4); else
if( newt <=     0.4880 +0.001 ) return     0.78441  * (1.0e4); else
if( newt <=     0.4890 +0.001 ) return     0.77531  * (1.0e4); else
if( newt <=     0.4900 +0.001 ) return     0.76761  * (1.0e4); else
if( newt <=     0.4910 +0.001 ) return     0.76051  * (1.0e4); else
if( newt <=     0.4920 +0.001 ) return     0.75401  * (1.0e4); else
if( newt <=     0.4930 +0.001 ) return     0.74501  * (1.0e4); else
if( newt <=     0.4940 +0.001 ) return     0.73741  * (1.0e4); else
if( newt <=     0.4950 +0.001 ) return     0.73031  * (1.0e4); else
if( newt <=     0.4960 +0.001 ) return     0.72381  * (1.0e4); else
if( newt <=     0.4970 +0.001 ) return     0.70921  * (1.0e4); else
if( newt <=     0.4980 +0.001 ) return     0.69701  * (1.0e4); else
if( newt <=     0.4990 +0.001 ) return     0.68591  * (1.0e4); else
if( newt <=     0.5000 +0.001 ) return     0.67541  * (1.0e4); else
if( newt <=     0.5010 +0.001 ) return     0.67231  * (1.0e4); else
if( newt <=     0.5020 +0.001 ) return     0.66751  * (1.0e4); else
if( newt <=     0.5030 +0.001 ) return     0.66181  * (1.0e4); else
if( newt <=     0.5040 +0.001 ) return     0.65491  * (1.0e4); else
if( newt <=     0.5050 +0.001 ) return     0.65251  * (1.0e4); else
if( newt <=     0.5060 +0.001 ) return     0.64732  * (1.0e4); else
if( newt <=     0.5070 +0.001 ) return     0.64022  * (1.0e4); else
if( newt <=     0.5080 +0.001 ) return     0.63142  * (1.0e4); else
if( newt <=     0.5090 +0.001 ) return     0.62132  * (1.0e4); else
if( newt <=     0.5100 +0.001 ) return     0.61002  * (1.0e4); else
if( newt <=     0.5110 +0.001 ) return     0.59812  * (1.0e4); else
if( newt <=     0.5120 +0.001 ) return     0.5858   * (1.0e4); else
if( newt <=     0.5130 +0.001 ) return     0.5707   * (1.0e4); else
if( newt <=     0.5140 +0.001 ) return     0.5569   * (1.0e4); else
if( newt <=     0.5150 +0.001 ) return     0.5439   * (1.0e4); else
if( newt <=     0.5160 +0.001 ) return     0.5318   * (1.0e4); else
if( newt <=     0.5170 +0.001 ) return     0.5291   * (1.0e4); else
if( newt <=     0.5180 +0.001 ) return     0.5248   * (1.0e4); else
if( newt <=     0.5190 +0.001 ) return     0.5199   * (1.0e4); else
if( newt <=     0.5200 +0.001 ) return     0.5141   * (1.0e4); else
if( newt <=     0.5210 +0.001 ) return     0.5099   * (1.0e4); else
if( newt <=     0.5220 +0.001 ) return     0.5037   * (1.0e4); else
if( newt <=     0.5230 +0.001 ) return     0.4958   * (1.0e4); else
if( newt <=     0.5240 +0.001 ) return     0.4864   * (1.0e4); else
if( newt <=     0.5250 +0.001 ) return     0.4759   * (1.0e4); else
if( newt <=     0.5260 +0.001 ) return     0.4649   * (1.0e4); else
if( newt <=     0.5270 +0.001 ) return     0.4539   * (1.0e4); else
if( newt <=     0.5280 +0.001 ) return     0.4434   * (1.0e4); else
if( newt <=     0.5290 +0.001 ) return     0.4309   * (1.0e4); else
if( newt <=     0.5300 +0.001 ) return     0.4204   * (1.0e4); else
if( newt <=     0.5310 +0.001 ) return     0.4115   * (1.0e4); else
if( newt <=     0.5320 +0.001 ) return     0.4040   * (1.0e4); else
if( newt <=     0.5330 +0.001 ) return     0.4048   * (1.0e4); else
if( newt <=     0.5340 +0.001 ) return     0.4048   * (1.0e4); else
if( newt <=     0.5350 +0.001 ) return     0.4047   * (1.0e4); else
if( newt <=     0.5360 +0.001 ) return     0.4043   * (1.0e4); else
if( newt <=     0.5370 +0.001 ) return     0.3908   * (1.0e4); else
if( newt <=     0.5380 +0.001 ) return     0.3801   * (1.0e4); else
if( newt <=     0.5390 +0.001 ) return     0.3705   * (1.0e4); else
if( newt <=     0.5400 +0.001 ) return     0.3620   * (1.0e4); else
if( newt <=     0.5410 +0.001 ) return     0.3447   * (1.0e4); else
if( newt <=     0.5420 +0.001 ) return     0.3312   * (1.0e4); else
if( newt <=     0.5430 +0.001 ) return     0.3200   * (1.0e4); else
if( newt <=     0.5440 +0.001 ) return     0.3107   * (1.0e4); else
if( newt <=     0.5450 +0.001 ) return     0.3246   * (1.0e4); else
if( newt <=     0.5460 +0.001 ) return     0.3339   * (1.0e4); else
if( newt <=     0.5470 +0.001 ) return     0.3412   * (1.0e4); else
if( newt <=     0.5480 +0.001 ) return     0.3461   * (1.0e4); else
if( newt <=     0.5490 +0.001 ) return     0.3512   * (1.0e4); else
if( newt <=     0.5500 +0.001 ) return     0.3529   * (1.0e4); else
if( newt <=     0.5510 +0.001 ) return     0.3520   * (1.0e4); else
if( newt <=     0.5520 +0.001 ) return     0.3489   * (1.0e4); else
if( newt <=     0.5530 +0.001 ) return     0.3303   * (1.0e4); else
if( newt <=     0.5540 +0.001 ) return     0.3151   * (1.0e4); else
if( newt <=     0.5550 +0.001 ) return     0.3019   * (1.0e4); else
if( newt <=     0.5560 +0.001 ) return     0.2913   * (1.0e4); else
if( newt <=     0.5570 +0.001 ) return     0.2778   * (1.0e4); else
if( newt <=     0.5580 +0.001 ) return     0.2686   * (1.0e4); else
if( newt <=     0.5590 +0.001 ) return     0.2624   * (1.0e4); else
if( newt <=     0.5600 +0.001 ) return     0.2583   * (1.0e4); else
if( newt <=     0.5610 +0.001 ) return     0.2423   * (1.0e4); else
if( newt <=     0.5620 +0.001 ) return     0.2298   * (1.0e4); else
if( newt <=     0.5630 +0.001 ) return     0.2179   * (1.0e4); else
if( newt <=     0.5640 +0.001 ) return     0.2057   * (1.0e4); else
if( newt <=     0.5650 +0.001 ) return     0.2154   * (1.0e4); else
if( newt <=     0.5660 +0.001 ) return     0.2183   * (1.0e4); else
if( newt <=     0.5670 +0.001 ) return     0.2177   * (1.0e4); else
if( newt <=     0.5680 +0.001 ) return     0.2140   * (1.0e4); else
if( newt <=     0.5690 +0.001 ) return     0.2090   * (1.0e4); else
if( newt <=     0.5700 +0.001 ) return     0.2016   * (1.0e4); else
if( newt <=     0.5710 +0.001 ) return     0.1924   * (1.0e4); else
if( newt <=     0.5720 +0.001 ) return     0.1819   * (1.0e4); else
if( newt <=     0.5730 +0.001 ) return     0.1777   * (1.0e4); else
if( newt <=     0.5740 +0.001 ) return     0.1711   * (1.0e4); else
if( newt <=     0.5750 +0.001 ) return     0.1635   * (1.0e4); else
if( newt <=     0.5760 +0.001 ) return     0.1554   * (1.0e4); else
if( newt <=     0.5770 +0.001 ) return     0.1343   * (1.0e4); else
if( newt <=     0.5780 +0.001 ) return     0.1166   * (1.0e4); else
if( newt <=     0.5790 +0.001 ) return     0.1007   * (1.0e4); else
if( newt <=     0.5800 +0.001 ) return     0.0864   * (1.0e4); else
if( newt <=     0.5810 +0.001 ) return     0.0822   * (1.0e4); else
if( newt <=     0.5820 +0.001 ) return     0.0772   * (1.0e4); else
if( newt <=     0.5830 +0.001 ) return     0.0723   * (1.0e4); else
if( newt <=     0.5840 +0.001 ) return     0.0672   * (1.0e4); else
if( newt <=     0.5850 +0.001 ) return     0.0492   * (1.0e4); else
if( newt <=     0.5860 +0.001 ) return     0.0340   * (1.0e4); else
if( newt <=     0.5870 +0.001 ) return     0.0197   * (1.0e4); else
if( newt <=     0.5880 +0.001 ) return     0.0061   * (1.0e4); else
if( newt <=     0.5890 +0.001 ) return    -0.0124   * (1.0e4); else
if( newt <=     0.5900 +0.001 ) return    -0.0287   * (1.0e4); else
if( newt <=     0.5910 +0.001 ) return    -0.0436   * (1.0e4); else
if( newt <=     0.5920 +0.001 ) return    -0.0571   * (1.0e4); else
if( newt <=     0.5930 +0.001 ) return    -0.0666   * (1.0e4); else
if( newt <=     0.5940 +0.001 ) return    -0.0756   * (1.0e4); else
if( newt <=     0.5950 +0.001 ) return    -0.0842   * (1.0e4); else
if( newt <=     0.5960 +0.001 ) return    -0.0926   * (1.0e4); else
if( newt <=     0.5970 +0.001 ) return    -0.0928   * (1.0e4); else
if( newt <=     0.5980 +0.001 ) return    -0.0958   * (1.0e4); else
if( newt <=     0.5990 +0.001 ) return    -0.1007   * (1.0e4); else
if( newt <=     0.6000 +0.001 ) return    -0.1074   * (1.0e4); else
if( newt <=     0.6010 +0.001 ) return    -0.1172   * (1.0e4); else
if( newt <=     0.6020 +0.001 ) return    -0.1281   * (1.0e4); else
if( newt <=     0.6030 +0.001 ) return    -0.1400   * (1.0e4); else
if( newt <=     0.6040 +0.001 ) return    -0.1526   * (1.0e4); else
if( newt <=     0.6050 +0.001 ) return    -0.1613   * (1.0e4); else
if( newt <=     0.6060 +0.001 ) return    -0.1714   * (1.0e4); else
if( newt <=     0.6070 +0.001 ) return    -0.1822   * (1.0e4); else
if( newt <=     0.6080 +0.001 ) return    -0.1936   * (1.0e4); else
if( newt <=     0.6090 +0.001 ) return    -0.2069   * (1.0e4); else
if( newt <=     0.6100 +0.001 ) return    -0.2203   * (1.0e4); else
if( newt <=     0.6110 +0.001 ) return    -0.2340   * (1.0e4); else
if( newt <=     0.6120 +0.001 ) return    -0.2481   * (1.0e4); else
if( newt <=     0.6130 +0.001 ) return    -0.2638   * (1.0e4); else
if( newt <=     0.6140 +0.001 ) return    -0.2792   * (1.0e4); else
if( newt <=     0.6150 +0.001 ) return    -0.2945   * (1.0e4); else
if( newt <=     0.6160 +0.001 ) return    -0.3096   * (1.0e4); else
if( newt <=     0.6170 +0.001 ) return    -0.3256   * (1.0e4); else
if( newt <=     0.6180 +0.001 ) return    -0.3408   * (1.0e4); else
if( newt <=     0.6190 +0.001 ) return    -0.3552   * (1.0e4); else
if( newt <=     0.6200 +0.001 ) return    -0.3690   * (1.0e4); else
if( newt <=     0.6210 +0.001 ) return    -0.3825   * (1.0e4); else
if( newt <=     0.6220 +0.001 ) return    -0.3961   * (1.0e4); else
if( newt <=     0.6230 +0.001 ) return    -0.4101   * (1.0e4); else
if( newt <=     0.6240 +0.001 ) return    -0.4247   * (1.0e4); else
if( newt <=     0.6250 +0.001 ) return    -0.4414   * (1.0e4); else
if( newt <=     0.6260 +0.001 ) return    -0.4583   * (1.0e4); else
if( newt <=     0.6270 +0.001 ) return    -0.4756   * (1.0e4); else
if( newt <=     0.6280 +0.001 ) return    -0.4932   * (1.0e4); else
if( newt <=     0.6290 +0.001 ) return    -0.5024   * (1.0e4); else
if( newt <=     0.6300 +0.001 ) return    -0.5139   * (1.0e4); else
if( newt <=     0.6310 +0.001 ) return    -0.5264   * (1.0e4); else
if( newt <=     0.6320 +0.001 ) return    -0.5396   * (1.0e4); else
if( newt <=     0.6330 +0.001 ) return    -0.5564   * (1.0e4); else
if( newt <=     0.6340 +0.001 ) return    -0.5729   * (1.0e4); else
if( newt <=     0.6350 +0.001 ) return    -0.5893   * (1.0e4); else
if( newt <=     0.6360 +0.001 ) return    -0.6055   * (1.0e4); else
if( newt <=     0.6370 +0.001 ) return    -0.6200   * (1.0e4); else
if( newt <=     0.6380 +0.001 ) return    -0.6345   * (1.0e4); else
if( newt <=     0.6390 +0.001 ) return    -0.6490   * (1.0e4); else
if( newt <=     0.6400 +0.001 ) return    -0.6633   * (1.0e4); else
if( newt <=     0.6410 +0.001 ) return    -0.6805   * (1.0e4); else
if( newt <=     0.6420 +0.001 ) return    -0.6971   * (1.0e4); else
if( newt <=     0.6430 +0.001 ) return    -0.7137   * (1.0e4); else
if( newt <=     0.6440 +0.001 ) return    -0.7303   * (1.0e4); else
if( newt <=     0.6450 +0.001 ) return    -0.7457   * (1.0e4); else
if( newt <=     0.6460 +0.001 ) return    -0.7617   * (1.0e4); else
if( newt <=     0.6470 +0.001 ) return    -0.7782   * (1.0e4); else
if( newt <=     0.6480 +0.001 ) return    -0.7951   * (1.0e4); else
if( newt <=     0.6490 +0.001 ) return    -0.8137   * (1.0e4); else
if( newt <=     0.6500 +0.001 ) return    -0.8322   * (1.0e4); else
if( newt <=     0.6510 +0.001 ) return    -0.8507   * (1.0e4); else
if( newt <=     0.6520 +0.001 ) return    -0.8692   * (1.0e4); else
if( newt <=     0.6530 +0.001 ) return    -0.8831   * (1.0e4); else
if( newt <=     0.6540 +0.001 ) return    -0.8979   * (1.0e4); else
if( newt <=     0.6550 +0.001 ) return    -0.9127   * (1.0e4); else
if( newt <=     0.6560 +0.001 ) return    -0.9275   * (1.0e4); else
if( newt <=     0.6570 +0.001 ) return    -0.9410   * (1.0e4); else
if( newt <=     0.6580 +0.001 ) return    -0.9551   * (1.0e4); else
if( newt <=     0.6590 +0.001 ) return    -0.9695   * (1.0e4); else
if( newt <=     0.6600 +0.001 ) return    -0.9843   * (1.0e4); else
if( newt <=     0.6610 +0.001 ) return    -0.9978   * (1.0e4); else
if( newt <=     0.6620 +0.001 ) return    -1.0116   * (1.0e4); else
if( newt <=     0.6630 +0.001 ) return    -1.0252   * (1.0e4); else
if( newt <=     0.6640 +0.001 ) return    -1.0384   * (1.0e4); else
if( newt <=     0.6650 +0.001 ) return    -1.0553   * (1.0e4); else
if( newt <=     0.6660 +0.001 ) return    -1.0704   * (1.0e4); else
if( newt <=     0.6670 +0.001 ) return    -1.0846   * (1.0e4); else
if( newt <=     0.6680 +0.001 ) return    -1.0977   * (1.0e4); else
if( newt <=     0.6690 +0.001 ) return    -1.1129   * (1.0e4); else
if( newt <=     0.6700 +0.001 ) return    -1.1267   * (1.0e4); else
if( newt <=     0.6710 +0.001 ) return    -1.1397   * (1.0e4); else
if( newt <=     0.6720 +0.001 ) return    -1.1522   * (1.0e4); else
if( newt <=     0.6730 +0.001 ) return    -1.1684   * (1.0e4); else
if( newt <=     0.6740 +0.001 ) return    -1.1833   * (1.0e4); else
if( newt <=     0.6750 +0.001 ) return    -1.1975   * (1.0e4); else
if( newt <=     0.6760 +0.001 ) return    -1.2113   * (1.0e4); else
if( newt <=     0.6770 +0.001 ) return    -1.2220   * (1.0e4); else
if( newt <=     0.6780 +0.001 ) return    -1.2334   * (1.0e4); else
if( newt <=     0.6790 +0.001 ) return    -1.2450   * (1.0e4); else
if( newt <=     0.6800 +0.001 ) return    -1.2571   * (1.0e4); else
if( newt <=     0.6810 +0.001 ) return    -1.2653   * (1.0e4); else
if( newt <=     0.6820 +0.001 ) return    -1.27522  * (1.0e4); else
if( newt <=     0.6830 +0.001 ) return    -1.28592  * (1.0e4); else
if( newt <=     0.6840 +0.001 ) return    -1.29731  * (1.0e4); else
if( newt <=     0.6850 +0.001 ) return    -1.31221  * (1.0e4); else
if( newt <=     0.6860 +0.001 ) return    -1.32661  * (1.0e4); else
if( newt <=     0.6870 +0.001 ) return    -1.34091  * (1.0e4); else
if( newt <=     0.6880 +0.001 ) return    -1.35481  * (1.0e4); else
if( newt <=     0.6890 +0.001 ) return    -1.35561  * (1.0e4); else
if( newt <=     0.6900 +0.001 ) return    -1.35921  * (1.0e4); else
if( newt <=     0.6910 +0.001 ) return    -1.36381  * (1.0e4); else
if( newt <=     0.6920 +0.001 ) return    -1.36911  * (1.0e4); else
if( newt <=     0.6930 +0.001 ) return    -1.37531  * (1.0e4); else
if( newt <=     0.6940 +0.001 ) return    -1.38211  * (1.0e4); else
if( newt <=     0.6950 +0.001 ) return    -1.38951  * (1.0e4); else
if( newt <=     0.6960 +0.001 ) return    -1.39721  * (1.0e4); else
if( newt <=     0.6970 +0.001 ) return    -1.34821  * (1.0e4); else
if( newt <=     0.6980 +0.001 ) return    -1.31321  * (1.0e4); else
if( newt <=     0.6990 +0.001 ) return    -1.28351  * (1.0e4); else
if( newt <=     0.7000 +0.001 ) return    -1.25771  * (1.0e4); else
if( newt <=     0.7010 +0.001 ) return    -1.28881  * (1.0e4); else
if( newt <=     0.7020 +0.001 ) return    -1.30841  * (1.0e4); else
if( newt <=     0.7030 +0.001 ) return    -1.32321  * (1.0e4); else
if( newt <=     0.7040 +0.001 ) return    -1.33291  * (1.0e4); else
if( newt <=     0.7050 +0.001 ) return    -1.37221  * (1.0e4); else
if( newt <=     0.7060 +0.001 ) return    -1.39612  * (1.0e4); else
if( newt <=     0.7070 +0.001 ) return    -1.40982  * (1.0e4); else
if( newt <=     0.7080 +0.001 ) return    -1.41422  * (1.0e4); else
if( newt <=     0.7090 +0.001 ) return    -1.40682  * (1.0e4); else
if( newt <=     0.7100 +0.001 ) return    -1.39472  * (1.0e4); else
if( newt <=     0.7110 +0.001 ) return    -1.37952  * (1.0e4); else
if( newt <=     0.7120 +0.001 ) return    -1.3640   * (1.0e4); else
if( newt <=     0.7130 +0.001 ) return    -1.3503   * (1.0e4); else
if( newt <=     0.7140 +0.001 ) return    -1.3404   * (1.0e4); else
if( newt <=     0.7150 +0.001 ) return    -1.3350   * (1.0e4); else
if( newt <=     0.7160 +0.001 ) return    -1.3346   * (1.0e4); else
if( newt <=     0.7170 +0.001 ) return    -1.3301   * (1.0e4); else
if( newt <=     0.7180 +0.001 ) return    -1.3313   * (1.0e4); else
if( newt <=     0.7190 +0.001 ) return    -1.3356   * (1.0e4); else
if( newt <=     0.7200 +0.001 ) return    -1.3417   * (1.0e4); else
if( newt <=     0.7210 +0.001 ) return    -1.3385   * (1.0e4); else
if( newt <=     0.7220 +0.001 ) return    -1.3378   * (1.0e4); else
if( newt <=     0.7230 +0.001 ) return    -1.3377   * (1.0e4); else
if( newt <=     0.7240 +0.001 ) return    -1.3380   * (1.0e4); else
if( newt <=     0.7250 +0.001 ) return    -1.3303   * (1.0e4); else
if( newt <=     0.7260 +0.001 ) return    -1.3255   * (1.0e4); else
if( newt <=     0.7270 +0.001 ) return    -1.3226   * (1.0e4); else
if( newt <=     0.7280 +0.001 ) return    -1.3214   * (1.0e4); else
if( newt <=     0.7290 +0.001 ) return    -1.2975   * (1.0e4); else
if( newt <=     0.7300 +0.001 ) return    -1.2807   * (1.0e4); else
if( newt <=     0.7310 +0.001 ) return    -1.2667   * (1.0e4); else
if( newt <=     0.7320 +0.001 ) return    -1.2546   * (1.0e4); else
if( newt <=     0.7330 +0.001 ) return    -1.2608   * (1.0e4); else
if( newt <=     0.7340 +0.001 ) return    -1.2631   * (1.0e4); else
if( newt <=     0.7350 +0.001 ) return    -1.2635   * (1.0e4); else
if( newt <=     0.7360 +0.001 ) return    -1.2615   * (1.0e4); else
if( newt <=     0.7370 +0.001 ) return    -1.2526   * (1.0e4); else
if( newt <=     0.7380 +0.001 ) return    -1.2416   * (1.0e4); else
if( newt <=     0.7390 +0.001 ) return    -1.2278   * (1.0e4); else
if( newt <=     0.7400 +0.001 ) return    -1.2112   * (1.0e4); else
if( newt <=     0.7410 +0.001 ) return    -1.1652   * (1.0e4); else
if( newt <=     0.7420 +0.001 ) return    -1.1240   * (1.0e4); else
if( newt <=     0.7430 +0.001 ) return    -1.0842   * (1.0e4); else
if( newt <=     0.7440 +0.001 ) return    -1.0457   * (1.0e4); else
if( newt <=     0.7450 +0.001 ) return    -1.0325   * (1.0e4); else
if( newt <=     0.7460 +0.001 ) return    -1.0144   * (1.0e4); else
if( newt <=     0.7470 +0.001 ) return    -0.9944   * (1.0e4); else
if( newt <=     0.7480 +0.001 ) return    -0.9724   * (1.0e4); else
if( newt <=     0.7490 +0.001 ) return    -0.9675   * (1.0e4); else
if( newt <=     0.7500 +0.001 ) return    -0.9542   * (1.0e4); else
if( newt <=     0.7510 +0.001 ) return    -0.9354   * (1.0e4); else
if( newt <=     0.7520 +0.001 ) return    -0.9113   * (1.0e4); else
if( newt <=     0.7530 +0.001 ) return    -0.8879   * (1.0e4); else
if( newt <=     0.7540 +0.001 ) return    -0.8588   * (1.0e4); else
if( newt <=     0.7550 +0.001 ) return    -0.8261   * (1.0e4); else
if( newt <=     0.7560 +0.001 ) return    -0.7912   * (1.0e4); else
if( newt <=     0.7570 +0.001 ) return    -0.7625   * (1.0e4); else
if( newt <=     0.7580 +0.001 ) return    -0.7323   * (1.0e4); else
if( newt <=     0.7590 +0.001 ) return    -0.7026   * (1.0e4); else
if( newt <=     0.7600 +0.001 ) return    -0.6739   * (1.0e4); else
if( newt <=     0.7610 +0.001 ) return    -0.6479   * (1.0e4); else
if( newt <=     0.7620 +0.001 ) return    -0.6224   * (1.0e4); else
if( newt <=     0.7630 +0.001 ) return    -0.5974   * (1.0e4); else
if( newt <=     0.7640 +0.001 ) return    -0.5724   * (1.0e4); else
if( newt <=     0.7650 +0.001 ) return    -0.5516   * (1.0e4); else
if( newt <=     0.7660 +0.001 ) return    -0.5295   * (1.0e4); else
if( newt <=     0.7670 +0.001 ) return    -0.5068   * (1.0e4); else
if( newt <=     0.7680 +0.001 ) return    -0.4838   * (1.0e4); else
if( newt <=     0.7690 +0.001 ) return    -0.4535   * (1.0e4); else
if( newt <=     0.7700 +0.001 ) return    -0.4255   * (1.0e4); else
if( newt <=     0.7710 +0.001 ) return    -0.3988   * (1.0e4); else
if( newt <=     0.7720 +0.001 ) return    -0.3736   * (1.0e4); else
if( newt <=     0.7730 +0.001 ) return    -0.3498   * (1.0e4); else
if( newt <=     0.7740 +0.001 ) return    -0.3275   * (1.0e4); else
if( newt <=     0.7750 +0.001 ) return    -0.3065   * (1.0e4); else
if( newt <=     0.7760 +0.001 ) return    -0.2865   * (1.0e4); else
if( newt <=     0.7770 +0.001 ) return    -0.2772   * (1.0e4); else
if( newt <=     0.7780 +0.001 ) return    -0.2659   * (1.0e4); else
if( newt <=     0.7790 +0.001 ) return    -0.2538   * (1.0e4); else
if( newt <=     0.7800 +0.001 ) return    -0.2409   * (1.0e4); else
if( newt <=     0.7810 +0.001 ) return    -0.2271   * (1.0e4); else
if( newt <=     0.7820 +0.001 ) return    -0.2127   * (1.0e4); else
if( newt <=     0.7830 +0.001 ) return    -0.1979   * (1.0e4); else
if( newt <=     0.7840 +0.001 ) return    -0.1829   * (1.0e4); else
if( newt <=     0.7850 +0.001 ) return    -0.1654   * (1.0e4); else
if( newt <=     0.7860 +0.001 ) return    -0.1494   * (1.0e4); else
if( newt <=     0.7870 +0.001 ) return    -0.1347   * (1.0e4); else
if( newt <=     0.7880 +0.001 ) return    -0.1215   * (1.0e4); else
if( newt <=     0.7890 +0.001 ) return    -0.1818   * (1.0e4); else
if( newt <=     0.7900 +0.001 ) return    -0.2258   * (1.0e4); else
if( newt <=     0.7910 +0.001 ) return    -0.2643   * (1.0e4); else
if( newt <=     0.7920 +0.001 ) return    -0.2985   * (1.0e4); else
if( newt <=     0.7930 +0.001 ) return    -0.2442   * (1.0e4); else
if( newt <=     0.7940 +0.001 ) return    -0.2085   * (1.0e4); else
if( newt <=     0.7950 +0.001 ) return    -0.1803   * (1.0e4); else
if( newt <=     0.7960 +0.001 ) return    -0.1595   * (1.0e4); else
if( newt <=     0.7970 +0.001 ) return    -0.1669   * (1.0e4); else
if( newt <=     0.7980 +0.001 ) return    -0.1773   * (1.0e4); else
if( newt <=     0.7990 +0.001 ) return    -0.1934   * (1.0e4);

    }

}
