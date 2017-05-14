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

Real abdominalAorta (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1e2 * (std::cos (6.283 * t) - 1) / 2;
}

Real aortaFluxIn (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.01 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.02 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.03 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.04 + 0.01)
    {
        return    -9.1759e-06 * 1e6;
    }
    if (t <= 0.05 + 0.01)
    {
        return    -3.0930e-05 * 1e6;
    }
    if (t <= 0.06 + 0.01)
    {
        return    -6.2639e-05 * 1e6;
    }
    if (t <= 0.07 + 0.01)
    {
        return    -1.0212e-04 * 1e6;
    }
    if (t <= 0.08 + 0.01)
    {
        return    -1.4760e-04 * 1e6;
    }
    if (t <= 0.09 + 0.01)
    {
        return    -1.9726e-04 * 1e6;
    }
    if (t <= 0.10 + 0.01)
    {
        return    -2.4980e-04 * 1e6;
    }
    if (t <= 0.11 + 0.01)
    {
        return    -2.9526e-04 * 1e6;
    }
    if (t <= 0.12 + 0.01)
    {
        return    -3.2956e-04 * 1e6;
    }
    if (t <= 0.13 + 0.01)
    {
        return    -3.5469e-04 * 1e6;
    }
    if (t <= 0.14 + 0.01)
    {
        return    -3.7250e-04 * 1e6;
    }
    if (t <= 0.15 + 0.01)
    {
        return    -3.8429e-04 * 1e6;
    }
    if (t <= 0.16 + 0.01)
    {
        return    -3.9123e-04 * 1e6;
    }
    if (t <= 0.17 + 0.01)
    {
        return    -3.9431e-04 * 1e6;
    }
    if (t <= 0.18 + 0.01)
    {
        return    -3.9349e-04 * 1e6;
    }
    if (t <= 0.19 + 0.01)
    {
        return    -3.8858e-04 * 1e6;
    }
    if (t <= 0.20 + 0.01)
    {
        return    -3.7985e-04 * 1e6;
    }
    if (t <= 0.21 + 0.01)
    {
        return    -3.6756e-04 * 1e6;
    }
    if (t <= 0.22 + 0.01)
    {
        return    -3.5207e-04 * 1e6;
    }
    if (t <= 0.23 + 0.01)
    {
        return    -3.3408e-04 * 1e6;
    }
    if (t <= 0.24 + 0.01)
    {
        return    -3.1402e-04 * 1e6;
    }
    if (t <= 0.25 + 0.01)
    {
        return    -2.9288e-04 * 1e6;
    }
    if (t <= 0.26 + 0.01)
    {
        return    -2.7154e-04 * 1e6;
    }
    if (t <= 0.27 + 0.01)
    {
        return    -2.5054e-04 * 1e6;
    }
    if (t <= 0.28 + 0.01)
    {
        return    -2.2979e-04 * 1e6;
    }
    if (t <= 0.29 + 0.01)
    {
        return    -2.0904e-04 * 1e6;
    }
    if (t <= 0.30 + 0.01)
    {
        return    -1.8880e-04 * 1e6;
    }
    if (t <= 0.31 + 0.01)
    {
        return    -1.6899e-04 * 1e6;
    }
    if (t <= 0.32 + 0.01)
    {
        return    -1.4864e-04 * 1e6;
    }
    if (t <= 0.33 + 0.01)
    {
        return    -1.2730e-04 * 1e6;
    }
    if (t <= 0.34 + 0.01)
    {
        return    -1.0400e-04 * 1e6;
    }
    if (t <= 0.35 + 0.01)
    {
        return    -7.9755e-05 * 1e6;
    }
    if (t <= 0.36 + 0.01)
    {
        return    -5.8719e-05 * 1e6;
    }
    if (t <= 0.37 + 0.01)
    {
        return    -4.0345e-05 * 1e6;
    }
    if (t <= 0.38 + 0.01)
    {
        return    -2.4596e-05 * 1e6;
    }
    if (t <= 0.39 + 0.01)
    {
        return    -1.2259e-05 * 1e6;
    }
    if (t <= 0.40 + 0.01)
    {
        return    -3.8110e-06 * 1e6;
    }
    if (t <= 0.41 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.42 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.43 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.44 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.45 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.46 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.47 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.48 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.49 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.50 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.51 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.52 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.53 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.54 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.55 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.56 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.57 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.58 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.59 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.60 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.61 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.62 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.63 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.64 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.65 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.66 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.67 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.68 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.69 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.70 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.71 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.72 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.73 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.74 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.75 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.76 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.77 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.78 + 0.01)
    {
        return    0.0000e+00;
    }
    if (t <= 0.79 + 0.01)
    {
        return    0.0000e+00;
    }
}

Real aortaPhisPress (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    /*switch(i) {
    case 1:
    return 0.0;
    break;*/
    //  case 2:
    if (t <= 0.00)
    {
        return -1.e1 * 11017;
    }
    if (t <= 0.01)
    {
        return -1.e1 * 10954;
    }
    if (t <= 0.02)
    {
        return -1.e1 * 10893;
    }
    if (t <= 0.03)
    {
        return -1.e1 * 10832;
    }
    if (t <= 0.04)
    {
        return -1.e1 * 10771;
    }
    if (t <= 0.05)
    {
        return -1.e1 * 10712;
    }
    if (t <= 0.06)
    {
        return -1.e1 * 10653;
    }
    if (t <= 0.07)
    {
        return -1.e1 * 11113;
    }
    if (t <= 0.08)
    {
        return -1.e1 * 11544;
    }
    if (t <= 0.09)
    {
        return -1.e1 * 11869;
    }
    if (t <= 0.10)
    {
        return -1.e1 * 12146;
    }
    if (t <= 0.11)
    {
        return -1.e1 * 12394;
    }
    if (t <= 0.12)
    {
        return -1.e1 * 12635;
    }
    if (t <= 0.13)
    {
        return -1.e1 * 12889;
    }
    if (t <= 0.14)
    {
        return -1.e1 * 13151;
    }
    if (t <= 0.15)
    {
        return -1.e1 * 13398;
    }
    if (t <= 0.16)
    {
        return -1.e1 * 13620;
    }
    if (t <= 0.17)
    {
        return -1.e1 * 13833;
    }
    if (t <= 0.18)
    {
        return -1.e1 * 14035;
    }
    if (t <= 0.19)
    {
        return -1.e1 * 14229;
    }
    if (t <= 0.20)
    {
        return -1.e1 * 14436;
    }
    if (t <= 0.21)
    {
        return -1.e1 * 14613;
    }
    if (t <= 0.22)
    {
        return -1.e1 * 14753;
    }
    if (t <= 0.23)
    {
        return -1.e1 * 14878;
    }
    if (t <= 0.24)
    {
        return -1.e1 * 14974;
    }
    if (t <= 0.25)
    {
        return -1.e1 * 15032;
    }
    if (t <= 0.26)
    {
        return -1.e1 * 15047;
    }
    if (t <= 0.27)
    {
        return -1.e1 * 15025;
    }
    if (t <= 0.28)
    {
        return -1.e1 * 14975;
    }
    if (t <= 0.29)
    {
        return -1.e1 * 14899;
    }
    if (t <= 0.30)
    {
        return -1.e1 * 14822;
    }
    if (t <= 0.31)
    {
        return -1.e1 * 14721;
    }
    if (t <= 0.32)
    {
        return -1.e1 * 14594;
    }
    if (t <= 0.33)
    {
        return -1.e1 * 14496;
    }
    if (t <= 0.34)
    {
        return -1.e1 * 14375;
    }
    if (t <= 0.35)
    {
        return -1.e1 * 14198;
    }
    if (t <= 0.36)
    {
        return -1.e1 * 13990;
    }
    if (t <= 0.37)
    {
        return -1.e1 * 13726;
    }
    if (t <= 0.38)
    {
        return -1.e1 * 13397;
    }
    if (t <= 0.39)
    {
        return -1.e1 * 13167;
    }
    if (t <= 0.40)
    {
        return -1.e1 * 13132;
    }
    if (t <= 0.41)
    {
        return -1.e1 * 13315;
    }
    if (t <= 0.42)
    {
        return -1.e1 * 13271;
    }
    if (t <= 0.43)
    {
        return -1.e1 * 13157;
    }
    if (t <= 0.44)
    {
        return -1.e1 * 13028;
    }
    if (t <= 0.45)
    {
        return -1.e1 * 12975;
    }
    if (t <= 0.46)
    {
        return -1.e1 * 12933;
    }
    if (t <= 0.47)
    {
        return -1.e1 * 12891;
    }
    if (t <= 0.48)
    {
        return -1.e1 * 12836;
    }
    if (t <= 0.49)
    {
        return -1.e1 * 12768;
    }
    if (t <= 0.50)
    {
        return -1.e1 * 12700;
    }
    if (t <= 0.51)
    {
        return -1.e1 * 12641;
    }
    if (t <= 0.52)
    {
        return -1.e1 * 12592;
    }
    if (t <= 0.53)
    {
        return -1.e1 * 12548;
    }
    if (t <= 0.54)
    {
        return -1.e1 * 12504;
    }
    if (t <= 0.55)
    {
        return -1.e1 * 12456;
    }
    if (t <= 0.56)
    {
        return -1.e1 * 12405;
    }
    if (t <= 0.57)
    {
        return -1.e1 * 12353;
    }
    if (t <= 0.58)
    {
        return -1.e1 * 12300;
    }
    if (t <= 0.59)
    {
        return -1.e1 * 12244;
    }
    if (t <= 0.60)
    {
        return -1.e1 * 12184;
    }
    if (t <= 0.61)
    {
        return -1.e1 * 12122;
    }
    if (t <= 0.62)
    {
        return -1.e1 * 12058;
    }
    if (t <= 0.63)
    {
        return -1.e1 * 11995;
    }
    if (t <= 0.64)
    {
        return -1.e1 * 11933;
    }
    if (t <= 0.65)
    {
        return -1.e1 * 11871;
    }
    if (t <= 0.66)
    {
        return -1.e1 * 11810;
    }
    if (t <= 0.67)
    {
        return -1.e1 * 11747;
    }
    if (t <= 0.68)
    {
        return -1.e1 * 11684;
    }
    if (t <= 0.69)
    {
        return -1.e1 * 11620;
    }
    if (t <= 0.70)
    {
        return -1.e1 * 11556;
    }
    if (t <= 0.71)
    {
        return -1.e1 * 11492;
    }
    if (t <= 0.72)
    {
        return -1.e1 * 11428;
    }
    if (t <= 0.73)
    {
        return -1.e1 * 11365;
    }
    if (t <= 0.74)
    {
        return -1.e1 * 11302;
    }
    if (t <= 0.75)
    {
        return -1.e1 * 11240;
    }
    if (t <= 0.76)
    {
        return -1.e1 * 11179;
    }
    if (t <= 0.77)
    {
        return -1.e1 * 11120;
    }
    if (t <= 0.78)
    {
        return -1.e1 * 11062;
    }
    if (t <= 0.79)
    {
        return -1.e1 * 11006;
    }
    //    break;
    /*  case 3:
    return 0.0;
    break;}
    return 0.;*/
}


Real aortaFlux3 (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00)
    {
        return 3.350E-05 * 1e6;
    }
    if (t <= 0.01)
    {
        return 3.373E-05 * 1e6;
    }
    if (t <= 0.02)
    {
        return 3.402E-05 * 1e6;
    }
    if (t <= 0.03)
    {
        return 3.434E-05 * 1e6;
    }
    if (t <= 0.04)
    {
        return 3.466E-05 * 1e6;
    }
    if (t <= 0.05)
    {
        return 3.495E-05 * 1e6;
    }
    if (t <= 0.06)
    {
        return 3.519E-05 * 1e6;
    }
    if (t <= 0.07)
    {
        return 3.539E-05 * 1e6;
    }
    if (t <= 0.08)
    {
        return 3.564E-05 * 1e6;
    }
    if (t <= 0.09)
    {
        return 3.617E-05 * 1e6;
    }
    if (t <= 0.10)
    {
        return 3.773E-05 * 1e6;
    }
    if (t <= 0.11)
    {
        return 4.176E-05 * 1e6;
    }
    if (t <= 0.12)
    {
        return 5.037E-05 * 1e6;
    }
    if (t <= 0.13)
    {
        return 6.546E-05 * 1e6;
    }
    if (t <= 0.14)
    {
        return 8.701E-05 * 1e6;
    }
    if (t <= 0.15)
    {
        return 1.117E-04 * 1e6;
    }
    if (t <= 0.16)
    {
        return 1.345E-04 * 1e6;
    }
    if (t <= 0.17)
    {
        return 1.519E-04 * 1e6;
    }
    if (t <= 0.18)
    {
        return 1.642E-04 * 1e6;
    }
    if (t <= 0.19)
    {
        return 1.737E-04 * 1e6;
    }
    if (t <= 0.20)
    {
        return 1.821E-04 * 1e6;
    }
    if (t <= 0.21)
    {
        return 1.897E-04 * 1e6;
    }
    if (t <= 0.22)
    {
        return 1.958E-04 * 1e6;
    }
    if (t <= 0.23)
    {
        return 1.999E-04 * 1e6;
    }
    if (t <= 0.24)
    {
        return 2.019E-04 * 1e6;
    }
    if (t <= 0.25)
    {
        return 2.020E-04 * 1e6;
    }
    if (t <= 0.26)
    {
        return 2.004E-04 * 1e6;
    }
    if (t <= 0.27)
    {
        return 1.972E-04 * 1e6;
    }
    if (t <= 0.28)
    {
        return 1.926E-04 * 1e6;
    }
    if (t <= 0.29)
    {
        return 1.868E-04 * 1e6;
    }
    if (t <= 0.30)
    {
        return 1.798E-04 * 1e6;
    }
    if (t <= 0.31)
    {
        return 1.719E-04 * 1e6;
    }
    if (t <= 0.32)
    {
        return 1.632E-04 * 1e6;
    }
    if (t <= 0.33)
    {
        return 1.540E-04 * 1e6;
    }
    if (t <= 0.34)
    {
        return 1.446E-04 * 1e6;
    }
    if (t <= 0.35)
    {
        return 1.350E-04 * 1e6;
    }
    if (t <= 0.36)
    {
        return 1.254E-04 * 1e6;
    }
    if (t <= 0.37)
    {
        return 1.158E-04 * 1e6;
    }
    if (t <= 0.38)
    {
        return 1.062E-04 * 1e6;
    }
    if (t <= 0.39)
    {
        return 9.651E-05 * 1e6;
    }
    if (t <= 0.40)
    {
        return 8.634E-05 * 1e6;
    }
    if (t <= 0.41)
    {
        return 7.558E-05 * 1e6;
    }
    if (t <= 0.42)
    {
        return 6.447E-05 * 1e6;
    }
    if (t <= 0.43)
    {
        return 5.382E-05 * 1e6;
    }
    if (t <= 0.44)
    {
        return 4.484E-05 * 1e6;
    }
    if (t <= 0.45)
    {
        return 3.865E-05 * 1e6;
    }
    if (t <= 0.46)
    {
        return 3.556E-05 * 1e6;
    }
    if (t <= 0.47)
    {
        return 3.473E-05 * 1e6;
    }
    if (t <= 0.48)
    {
        return 3.457E-05 * 1e6;
    }
    if (t <= 0.49)
    {
        return 3.373E-05 * 1e6;
    }
    if (t <= 0.50)
    {
        return 3.191E-05 * 1e6;
    }
    if (t <= 0.51)
    {
        return 2.975E-05 * 1e6;
    }
    if (t <= 0.52)
    {
        return 2.809E-05 * 1e6;
    }
    if (t <= 0.53)
    {
        return 2.730E-05 * 1e6;
    }
    if (t <= 0.54)
    {
        return 2.718E-05 * 1e6;
    }
    if (t <= 0.55)
    {
        return 2.732E-05 * 1e6;
    }
    if (t <= 0.56)
    {
        return 2.744E-05 * 1e6;
    }
    if (t <= 0.57)
    {
        return 2.753E-05 * 1e6;
    }
    if (t <= 0.58)
    {
        return 2.772E-05 * 1e6;
    }
    if (t <= 0.59)
    {
        return 2.811E-05 * 1e6;
    }
    if (t <= 0.60)
    {
        return 2.866E-05 * 1e6;
    }
    if (t <= 0.61)
    {
        return 2.929E-05 * 1e6;
    }
    if (t <= 0.62)
    {
        return 2.990E-05 * 1e6;
    }
    if (t <= 0.63)
    {
        return 3.044E-05 * 1e6;
    }
    if (t <= 0.64)
    {
        return 3.091E-05 * 1e6;
    }
    if (t <= 0.65)
    {
        return 3.132E-05 * 1e6;
    }
    if (t <= 0.66)
    {
        return 3.168E-05 * 1e6;
    }
    if (t <= 0.67)
    {
        return 3.199E-05 * 1e6;
    }
    if (t <= 0.68)
    {
        return 3.224E-05 * 1e6;
    }
    if (t <= 0.69)
    {
        return 3.244E-05 * 1e6;
    }
    if (t <= 0.70)
    {
        return 3.259E-05 * 1e6;
    }
    if (t <= 0.71)
    {
        return 3.270E-05 * 1e6;
    }
    if (t <= 0.72)
    {
        return 3.277E-05 * 1e6;
    }
    if (t <= 0.73)
    {
        return 3.281E-05 * 1e6;
    }
    if (t <= 0.74)
    {
        return 3.282E-05 * 1e6;
    }
    if (t <= 0.75)
    {
        return 3.283E-05 * 1e6;
    }
    if (t <= 0.76)
    {
        return 3.283E-05 * 1e6;
    }
    if (t <= 0.77)
    {
        return 3.286E-05 * 1e6;
    }
    if (t <= 0.78)
    {
        return 3.291E-05 * 1e6;
    }
    if (t <= 0.79)
    {
        return 3.300E-05 * 1e6;
    }
}//thoracic aorta,

Real aortaFlux5 (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00)
    {
        return  3.033E-06 * 1e6;
    }
    if (t <= 0.01)
    {
        return  3.041E-06 * 1e6;
    }
    if (t <= 0.02)
    {
        return  3.062E-06 * 1e6;
    }
    if (t <= 0.03)
    {
        return  3.094E-06 * 1e6;
    }
    if (t <= 0.04)
    {
        return  3.127E-06 * 1e6;
    }
    if (t <= 0.05)
    {
        return  3.150E-06 * 1e6;
    }
    if (t <= 0.06)
    {
        return  3.152E-06 * 1e6;
    }
    if (t <= 0.07)
    {
        return  3.141E-06 * 1e6;
    }
    if (t <= 0.08)
    {
        return  3.196E-06 * 1e6;
    }
    if (t <= 0.09)
    {
        return  3.574E-06 * 1e6;
    }
    if (t <= 0.10)
    {
        return  4.778E-06 * 1e6;
    }
    if (t <= 0.11)
    {
        return  7.387E-06 * 1e6;
    }
    if (t <= 0.12)
    {
        return  1.150E-05 * 1e6;
    }
    if (t <= 0.13)
    {
        return  1.609E-05 * 1e6;
    }
    if (t <= 0.14)
    {
        return  1.933E-05 * 1e6;
    }
    if (t <= 0.15)
    {
        return  2.007E-05 * 1e6;
    }
    if (t <= 0.16)
    {
        return  1.885E-05 * 1e6;
    }
    if (t <= 0.17)
    {
        return  1.706E-05 * 1e6;
    }
    if (t <= 0.18)
    {
        return  1.569E-05 * 1e6;
    }
    if (t <= 0.19)
    {
        return  1.481E-05 * 1e6;
    }
    if (t <= 0.20)
    {
        return  1.401E-05 * 1e6;
    }
    if (t <= 0.21)
    {
        return  1.294E-05 * 1e6;
    }
    if (t <= 0.22)
    {
        return  1.160E-05 * 1e6;
    }
    if (t <= 0.23)
    {
        return  1.018E-05 * 1e6;
    }
    if (t <= 0.24)
    {
        return  8.832E-06 * 1e6;
    }
    if (t <= 0.25)
    {
        return  7.609E-06 * 1e6;
    }
    if (t <= 0.26)
    {
        return  6.578E-06 * 1e6;
    }
    if (t <= 0.27)
    {
        return  5.843E-06 * 1e6;
    }
    if (t <= 0.28)
    {
        return  5.472E-06 * 1e6;
    }
    if (t <= 0.29)
    {
        return  5.412E-06 * 1e6;
    }
    if (t <= 0.30)
    {
        return  5.491E-06 * 1e6;
    }
    if (t <= 0.31)
    {
        return  5.527E-06 * 1e6;
    }
    if (t <= 0.32)
    {
        return  5.420E-06 * 1e6;
    }
    if (t <= 0.33)
    {
        return  5.169E-06 * 1e6;
    }
    if (t <= 0.34)
    {
        return  4.829E-06 * 1e6;
    }
    if (t <= 0.35)
    {
        return  4.465E-06 * 1e6;
    }
    if (t <= 0.36)
    {
        return  4.111E-06 * 1e6;
    }
    if (t <= 0.37)
    {
        return  3.750E-06 * 1e6;
    }
    if (t <= 0.38)
    {
        return  3.304E-06 * 1e6;
    }
    if (t <= 0.39)
    {
        return  2.668E-06 * 1e6;
    }
    if (t <= 0.40)
    {
        return  1.800E-06 * 1e6;
    }
    if (t <= 0.41)
    {
        return  8.269E-07 * 1e6;
    }
    if (t <= 0.42)
    {
        return  9.760E-08 * 1e6;
    }
    if (t <= 0.43)
    {
        return  7.311E-08 * 1e6;
    }
    if (t <= 0.44)
    {
        return  1.041E-06 * 1e6;
    }
    if (t <= 0.45)
    {
        return  2.783E-06 * 1e6;
    }
    if (t <= 0.46)
    {
        return  4.537E-06 * 1e6;
    }
    if (t <= 0.47)
    {
        return  5.488E-06 * 1e6;
    }
    if (t <= 0.48)
    {
        return  5.431E-06 * 1e6;
    }
    if (t <= 0.49)
    {
        return  4.863E-06 * 1e6;
    }
    if (t <= 0.50)
    {
        return  4.452E-06 * 1e6;
    }
    if (t <= 0.51)
    {
        return  4.499E-06 * 1e6;
    }
    if (t <= 0.52)
    {
        return  4.824E-06 * 1e6;
    }
    if (t <= 0.53)
    {
        return  5.059E-06 * 1e6;
    }
    if (t <= 0.54)
    {
        return  4.989E-06 * 1e6;
    }
    if (t <= 0.55)
    {
        return  4.671E-06 * 1e6;
    }
    if (t <= 0.56)
    {
        return  4.292E-06 * 1e6;
    }
    if (t <= 0.57)
    {
        return  3.981E-06 * 1e6;
    }
    if (t <= 0.58)
    {
        return  3.749E-06 * 1e6;
    }
    if (t <= 0.59)
    {
        return  3.553E-06 * 1e6;
    }
    if (t <= 0.60)
    {
        return  3.377E-06 * 1e6;
    }
    if (t <= 0.61)
    {
        return  3.255E-06 * 1e6;
    }
    if (t <= 0.62)
    {
        return  3.224E-06 * 1e6;
    }
    if (t <= 0.63)
    {
        return  3.281E-06 * 1e6;
    }
    if (t <= 0.64)
    {
        return  3.377E-06 * 1e6;
    }
    if (t <= 0.65)
    {
        return  3.452E-06 * 1e6;
    }
    if (t <= 0.66)
    {
        return  3.472E-06 * 1e6;
    }
    if (t <= 0.67)
    {
        return  3.441E-06 * 1e6;
    }
    if (t <= 0.68)
    {
        return  3.389E-06 * 1e6;
    }
    if (t <= 0.69)
    {
        return  3.343E-06 * 1e6;
    }
    if (t <= 0.70)
    {
        return  3.312E-06 * 1e6;
    }
    if (t <= 0.71)
    {
        return  3.289E-06 * 1e6;
    }
    if (t <= 0.72)
    {
        return  3.262E-06 * 1e6;
    }
    if (t <= 0.73)
    {
        return  3.223E-06 * 1e6;
    }
    if (t <= 0.74)
    {
        return  3.177E-06 * 1e6;
    }
    if (t <= 0.75)
    {
        return  3.132E-06 * 1e6;
    }
    if (t <= 0.76)
    {
        return  3.094E-06 * 1e6;
    }
    if (t <= 0.77)
    {
        return  3.065E-06 * 1e6;
    }
    if (t <= 0.78)
    {
        return  3.040E-06 * 1e6;
    }
    if (t <= 0.79)
    {
        return  3.016E-06 * 1e6;
    }
}//first branchstd::placeholders::_1,

Real aortaFlux6 (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00)
    {
        return  7.817E-07 * 1e6;
    }
    if (t <= 0.01)
    {
        return  7.879E-07 * 1e6;
    }
    if (t <= 0.02)
    {
        return  7.977E-07 * 1e6;
    }
    if (t <= 0.03)
    {
        return  8.077E-07 * 1e6;
    }
    if (t <= 0.04)
    {
        return  8.144E-07 * 1e6;
    }
    if (t <= 0.05)
    {
        return  8.164E-07 * 1e6;
    }
    if (t <= 0.06)
    {
        return  8.136E-07 * 1e6;
    }
    if (t <= 0.07)
    {
        return  8.129E-07 * 1e6;
    }
    if (t <= 0.08)
    {
        return  8.462E-07 * 1e6;
    }
    if (t <= 0.09)
    {
        return  9.858E-07 * 1e6;
    }
    if (t <= 0.10)
    {
        return  1.317E-06 * 1e6;
    }
    if (t <= 0.11)
    {
        return  1.854E-06 * 1e6;
    }
    if (t <= 0.12)
    {
        return  2.476E-06 * 1e6;
    }
    if (t <= 0.13)
    {
        return  3.005E-06 * 1e6;
    }
    if (t <= 0.14)
    {
        return  3.366E-06 * 1e6;
    }
    if (t <= 0.15)
    {
        return  3.626E-06 * 1e6;
    }
    if (t <= 0.16)
    {
        return  3.860E-06 * 1e6;
    }
    if (t <= 0.17)
    {
        return  4.044E-06 * 1e6;
    }
    if (t <= 0.18)
    {
        return  4.065E-06 * 1e6;
    }
    if (t <= 0.19)
    {
        return  3.824E-06 * 1e6;
    }
    if (t <= 0.20)
    {
        return  3.320E-06 * 1e6;
    }
    if (t <= 0.21)
    {
        return  2.659E-06 * 1e6;
    }
    if (t <= 0.22)
    {
        return  2.006E-06 * 1e6;
    }
    if (t <= 0.23)
    {
        return  1.503E-06 * 1e6;
    }
    if (t <= 0.24)
    {
        return  1.215E-06 * 1e6;
    }
    if (t <= 0.25)
    {
        return  1.117E-06 * 1e6;
    }
    if (t <= 0.26)
    {
        return  1.136E-06 * 1e6;
    }
    if (t <= 0.27)
    {
        return  1.190E-06 * 1e6;
    }
    if (t <= 0.28)
    {
        return  1.222E-06 * 1e6;
    }
    if (t <= 0.29)
    {
        return  1.219E-06 * 1e6;
    }
    if (t <= 0.30)
    {
        return  1.197E-06 * 1e6;
    }
    if (t <= 0.31)
    {
        return  1.179E-06 * 1e6;
    }
    if (t <= 0.32)
    {
        return  1.175E-06 * 1e6;
    }
    if (t <= 0.33)
    {
        return  1.176E-06 * 1e6;
    }
    if (t <= 0.34)
    {
        return  1.164E-06 * 1e6;
    }
    if (t <= 0.35)
    {
        return  1.129E-06 * 1e6;
    }
    if (t <= 0.36)
    {
        return  1.063E-06 * 1e6;
    }
    if (t <= 0.37)
    {
        return  9.647E-07 * 1e6;
    }
    if (t <= 0.38)
    {
        return  8.310E-07 * 1e6;
    }
    if (t <= 0.39)
    {
        return  6.635E-07 * 1e6;
    }
    if (t <= 0.40)
    {
        return  4.774E-07 * 1e6;
    }
    if (t <= 0.41)
    {
        return  3.116E-07 * 1e6;
    }
    if (t <= 0.42)
    {
        return  2.251E-07 * 1e6;
    }
    if (t <= 0.43)
    {
        return  2.627E-07 * 1e6;
    }
    if (t <= 0.44)
    {
        return  4.099E-07 * 1e6;
    }
    if (t <= 0.45)
    {
        return  5.913E-07 * 1e6;
    }
    if (t <= 0.46)
    {
        return  7.359E-07 * 1e6;
    }
    if (t <= 0.47)
    {
        return  8.403E-07 * 1e6;
    }
    if (t <= 0.48)
    {
        return  9.515E-07 * 1e6;
    }
    if (t <= 0.49)
    {
        return  1.097E-06 * 1e6;
    }
    if (t <= 0.50)
    {
        return  1.245E-06 * 1e6;
    }
    if (t <= 0.51)
    {
        return  1.331E-06 * 1e6;
    }
    if (t <= 0.52)
    {
        return  1.314E-06 * 1e6;
    }
    if (t <= 0.53)
    {
        return  1.204E-06 * 1e6;
    }
    if (t <= 0.54)
    {
        return  1.048E-06 * 1e6;
    }
    if (t <= 0.55)
    {
        return  9.004E-07 * 1e6;
    }
    if (t <= 0.56)
    {
        return  7.937E-07 * 1e6;
    }
    if (t <= 0.57)
    {
        return  7.358E-07 * 1e6;
    }
    if (t <= 0.58)
    {
        return  7.163E-07 * 1e6;
    }
    if (t <= 0.59)
    {
        return  7.186E-07 * 1e6;
    }
    if (t <= 0.60)
    {
        return  7.280E-07 * 1e6;
    }
    if (t <= 0.61)
    {
        return  7.359E-07 * 1e6;
    }
    if (t <= 0.62)
    {
        return  7.400E-07 * 1e6;
    }
    if (t <= 0.63)
    {
        return  7.423E-07 * 1e6;
    }
    if (t <= 0.64)
    {
        return  7.465E-07 * 1e6;
    }
    if (t <= 0.65)
    {
        return  7.556E-07 * 1e6;
    }
    if (t <= 0.66)
    {
        return  7.700E-07 * 1e6;
    }
    if (t <= 0.67)
    {
        return  7.871E-07 * 1e6;
    }
    if (t <= 0.68)
    {
        return  8.028E-07 * 1e6;
    }
    if (t <= 0.69)
    {
        return  8.129E-07 * 1e6;
    }
    if (t <= 0.70)
    {
        return  8.151E-07 * 1e6;
    }
    if (t <= 0.71)
    {
        return  8.101E-07 * 1e6;
    }
    if (t <= 0.72)
    {
        return  8.006E-07 * 1e6;
    }
    if (t <= 0.73)
    {
        return  7.905E-07 * 1e6;
    }
    if (t <= 0.74)
    {
        return  7.830E-07 * 1e6;
    }
    if (t <= 0.75)
    {
        return  7.792E-07 * 1e6;
    }
    if (t <= 0.76)
    {
        return  7.782E-07 * 1e6;
    }
    if (t <= 0.77)
    {
        return  7.782E-07 * 1e6;
    }
    if (t <= 0.78)
    {
        return  7.775E-07 * 1e6;
    }
    if (t <= 0.79)
    {
        return  7.754E-07 * 1e6;
    }
}//branch 1_2 smallest

Real aortaFlux7 (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00)
    {
        return  1.930E-06 * 1e6;
    }
    if (t <= 0.01)
    {
        return  1.710E-06 * 1e6;
    }
    if (t <= 0.02)
    {
        return  1.495E-06 * 1e6;
    }
    if (t <= 0.03)
    {
        return  1.289E-06 * 1e6;
    }
    if (t <= 0.04)
    {
        return  1.096E-06 * 1e6;
    }
    if (t <= 0.05)
    {
        return  9.184E-07 * 1e6;
    }
    if (t <= 0.06)
    {
        return  7.568E-07 * 1e6;
    }
    if (t <= 0.07)
    {
        return  6.240E-07 * 1e6;
    }
    if (t <= 0.08)
    {
        return  5.854E-07 * 1e6;
    }
    if (t <= 0.09)
    {
        return  8.113E-07 * 1e6;
    }
    if (t <= 0.10)
    {
        return  1.559E-06 * 1e6;
    }
    if (t <= 0.11)
    {
        return  3.032E-06 * 1e6;
    }
    if (t <= 0.12)
    {
        return  5.221E-06 * 1e6;
    }
    if (t <= 0.13)
    {
        return  7.903E-06 * 1e6;
    }
    if (t <= 0.14)
    {
        return  1.081E-05 * 1e6;
    }
    if (t <= 0.15)
    {
        return  1.376E-05 * 1e6;
    }
    if (t <= 0.16)
    {
        return  1.661E-05 * 1e6;
    }
    if (t <= 0.17)
    {
        return  1.919E-05 * 1e6;
    }
    if (t <= 0.18)
    {
        return  2.134E-05 * 1e6;
    }
    if (t <= 0.19)
    {
        return  2.292E-05 * 1e6;
    }
    if (t <= 0.20)
    {
        return  2.389E-05 * 1e6;
    }
    if (t <= 0.21)
    {
        return  2.433E-05 * 1e6;
    }
    if (t <= 0.22)
    {
        return  2.435E-05 * 1e6;
    }
    if (t <= 0.23)
    {
        return  2.406E-05 * 1e6;
    }
    if (t <= 0.24)
    {
        return  2.353E-05 * 1e6;
    }
    if (t <= 0.25)
    {
        return  2.280E-05 * 1e6;
    }
    if (t <= 0.26)
    {
        return  2.188E-05 * 1e6;
    }
    if (t <= 0.27)
    {
        return  2.078E-05 * 1e6;
    }
    if (t <= 0.28)
    {
        return  1.947E-05 * 1e6;
    }
    if (t <= 0.29)
    {
        return  1.797E-05 * 1e6;
    }
    if (t <= 0.30)
    {
        return  1.628E-05 * 1e6;
    }
    if (t <= 0.31)
    {
        return  1.445E-05 * 1e6;
    }
    if (t <= 0.32)
    {
        return  1.254E-05 * 1e6;
    }
    if (t <= 0.33)
    {
        return  1.060E-05 * 1e6;
    }
    if (t <= 0.34)
    {
        return  8.684E-06 * 1e6;
    }
    if (t <= 0.35)
    {
        return  6.838E-06 * 1e6;
    }
    if (t <= 0.36)
    {
        return  5.084E-06 * 1e6;
    }
    if (t <= 0.37)
    {
        return  3.412E-06 * 1e6;
    }
    if (t <= 0.38)
    {
        return  1.784E-06 * 1e6;
    }
    if (t <= 0.39)
    {
        return  1.534E-07 * 1e6;
    }
    if (t <= 0.40)
    {
        return  -1.494E-06 * 1e6;
    }
    if (t <= 0.41)
    {
        return  -3.093E-06 * 1e6;
    }
    if (t <= 0.42)
    {
        return  -4.495E-06 * 1e6;
    }
    if (t <= 0.43)
    {
        return  -5.521E-06 * 1e6;
    }
    if (t <= 0.44)
    {
        return  -6.086E-06 * 1e6;
    }
    if (t <= 0.45)
    {
        return  -6.252E-06 * 1e6;
    }
    if (t <= 0.46)
    {
        return  -6.160E-06 * 1e6;
    }
    if (t <= 0.47)
    {
        return  -5.908E-06 * 1e6;
    }
    if (t <= 0.48)
    {
        return  -5.508E-06 * 1e6;
    }
    if (t <= 0.49)
    {
        return  -4.944E-06 * 1e6;
    }
    if (t <= 0.50)
    {
        return  -4.234E-06 * 1e6;
    }
    if (t <= 0.51)
    {
        return  -3.441E-06 * 1e6;
    }
    if (t <= 0.52)
    {
        return  -2.639E-06 * 1e6;
    }
    if (t <= 0.53)
    {
        return  -1.873E-06 * 1e6;
    }
    if (t <= 0.54)
    {
        return  -1.160E-06 * 1e6;
    }
    if (t <= 0.55)
    {
        return  -5.019E-07 * 1e6;
    }
    if (t <= 0.56)
    {
        return  1.024E-07 * 1e6;
    }
    if (t <= 0.57)
    {
        return  6.543E-07 * 1e6;
    }
    if (t <= 0.58)
    {
        return  1.156E-06 * 1e6;
    }
    if (t <= 0.59)
    {
        return  1.614E-06 * 1e6;
    }
    if (t <= 0.60)
    {
        return  2.034E-06 * 1e6;
    }
    if (t <= 0.61)
    {
        return  2.419E-06 * 1e6;
    }
    if (t <= 0.62)
    {
        return  2.765E-06 * 1e6;
    }
    if (t <= 0.63)
    {
        return  3.065E-06 * 1e6;
    }
    if (t <= 0.64)
    {
        return  3.307E-06 * 1e6;
    }
    if (t <= 0.65)
    {
        return  3.486E-06 * 1e6;
    }
    if (t <= 0.66)
    {
        return  3.600E-06 * 1e6;
    }
    if (t <= 0.67)
    {
        return  3.654E-06 * 1e6;
    }
    if (t <= 0.68)
    {
        return  3.657E-06 * 1e6;
    }
    if (t <= 0.69)
    {
        return  3.621E-06 * 1e6;
    }
    if (t <= 0.70)
    {
        return  3.554E-06 * 1e6;
    }
    if (t <= 0.71)
    {
        return  3.465E-06 * 1e6;
    }
    if (t <= 0.72)
    {
        return  3.357E-06 * 1e6;
    }
    if (t <= 0.73)
    {
        return  3.233E-06 * 1e6;
    }
    if (t <= 0.74)
    {
        return  3.094E-06 * 1e6;
    }
    if (t <= 0.75)
    {
        return  2.941E-06 * 1e6;
    }
    if (t <= 0.76)
    {
        return  2.774E-06 * 1e6;
    }
    if (t <= 0.77)
    {
        return  2.591E-06 * 1e6;
    }
    if (t <= 0.78)
    {
        return  2.395E-06 * 1e6;
    }
    if (t <= 0.79)
    {
        return  2.185E-06 * 1e6;
    }
}//R. Brachia, branch 1_3

Real aortaFlux8 (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00)
    {
        return  2.445E-06 * 1e6;
    }
    if (t <= 0.01)
    {
        return  2.470E-06 * 1e6;
    }
    if (t <= 0.02)
    {
        return  2.519E-06 * 1e6;
    }
    if (t <= 0.03)
    {
        return  2.578E-06 * 1e6;
    }
    if (t <= 0.04)
    {
        return  2.631E-06 * 1e6;
    }
    if (t <= 0.05)
    {
        return  2.665E-06 * 1e6;
    }
    if (t <= 0.06)
    {
        return  2.674E-06 * 1e6;
    }
    if (t <= 0.07)
    {
        return  2.674E-06 * 1e6;
    }
    if (t <= 0.08)
    {
        return  2.777E-06 * 1e6;
    }
    if (t <= 0.09)
    {
        return  3.260E-06 * 1e6;
    }
    if (t <= 0.10)
    {
        return  4.539E-06 * 1e6;
    }
    if (t <= 0.11)
    {
        return  6.956E-06 * 1e6;
    }
    if (t <= 0.12)
    {
        return  1.045E-05 * 1e6;
    }
    if (t <= 0.13)
    {
        return  1.438E-05 * 1e6;
    }
    if (t <= 0.14)
    {
        return  1.769E-05 * 1e6;
    }
    if (t <= 0.15)
    {
        return  1.958E-05 * 1e6;
    }
    if (t <= 0.16)
    {
        return  1.991E-05 * 1e6;
    }
    if (t <= 0.17)
    {
        return  1.910E-05 * 1e6;
    }
    if (t <= 0.18)
    {
        return  1.775E-05 * 1e6;
    }
    if (t <= 0.19)
    {
        return  1.628E-05 * 1e6;
    }
    if (t <= 0.20)
    {
        return  1.479E-05 * 1e6;
    }
    if (t <= 0.21)
    {
        return  1.326E-05 * 1e6;
    }
    if (t <= 0.22)
    {
        return  1.165E-05 * 1e6;
    }
    if (t <= 0.23)
    {
        return  1.001E-05 * 1e6;
    }
    if (t <= 0.24)
    {
        return  8.448E-06 * 1e6;
    }
    if (t <= 0.25)
    {
        return  7.065E-06 * 1e6;
    }
    if (t <= 0.26)
    {
        return  5.930E-06 * 1e6;
    }
    if (t <= 0.27)
    {
        return  5.089E-06 * 1e6;
    }
    if (t <= 0.28)
    {
        return  4.566E-06 * 1e6;
    }
    if (t <= 0.29)
    {
        return  4.350E-06 * 1e6;
    }
    if (t <= 0.30)
    {
        return  4.382E-06 * 1e6;
    }
    if (t <= 0.31)
    {
        return  4.552E-06 * 1e6;
    }
    if (t <= 0.32)
    {
        return  4.729E-06 * 1e6;
    }
    if (t <= 0.33)
    {
        return  4.798E-06 * 1e6;
    }
    if (t <= 0.34)
    {
        return  4.698E-06 * 1e6;
    }
    if (t <= 0.35)
    {
        return  4.428E-06 * 1e6;
    }
    if (t <= 0.36)
    {
        return  4.026E-06 * 1e6;
    }
    if (t <= 0.37)
    {
        return  3.513E-06 * 1e6;
    }
    if (t <= 0.38)
    {
        return  2.875E-06 * 1e6;
    }
    if (t <= 0.39)
    {
        return  2.076E-06 * 1e6;
    }
    if (t <= 0.40)
    {
        return  1.131E-06 * 1e6;
    }
    if (t <= 0.41)
    {
        return  1.861E-07 * 1e6;
    }
    if (t <= 0.42)
    {
        return  -4.577E-0 * 1e6;
    }
    if (t <= 0.43)
    {
        return  -4.657E-0 * 1e6;
    }
    if (t <= 0.44)
    {
        return  3.124E-07 * 1e6;
    }
    if (t <= 0.45)
    {
        return  1.684E-06 * 1e6;
    }
    if (t <= 0.46)
    {
        return  3.174E-06 * 1e6;
    }
    if (t <= 0.47)
    {
        return  4.306E-06 * 1e6;
    }
    if (t <= 0.48)
    {
        return  4.873E-06 * 1e6;
    }
    if (t <= 0.49)
    {
        return  4.980E-06 * 1e6;
    }
    if (t <= 0.50)
    {
        return  4.876E-06 * 1e6;
    }
    if (t <= 0.51)
    {
        return  4.757E-06 * 1e6;
    }
    if (t <= 0.52)
    {
        return  4.681E-06 * 1e6;
    }
    if (t <= 0.53)
    {
        return  4.603E-06 * 1e6;
    }
    if (t <= 0.54)
    {
        return  4.458E-06 * 1e6;
    }
    if (t <= 0.55)
    {
        return  4.224E-06 * 1e6;
    }
    if (t <= 0.56)
    {
        return  3.928E-06 * 1e6;
    }
    if (t <= 0.57)
    {
        return  3.618E-06 * 1e6;
    }
    if (t <= 0.58)
    {
        return  3.335E-06 * 1e6;
    }
    if (t <= 0.59)
    {
        return  3.098E-06 * 1e6;
    }
    if (t <= 0.60)
    {
        return  2.913E-06 * 1e6;
    }
    if (t <= 0.61)
    {
        return  2.778E-06 * 1e6;
    }
    if (t <= 0.62)
    {
        return  2.697E-06 * 1e6;
    }
    if (t <= 0.63)
    {
        return  2.670E-06 * 1e6;
    }
    if (t <= 0.64)
    {
        return  2.692E-06 * 1e6;
    }
    if (t <= 0.65)
    {
        return  2.746E-06 * 1e6;
    }
    if (t <= 0.66)
    {
        return  2.811E-06 * 1e6;
    }
    if (t <= 0.67)
    {
        return  2.868E-06 * 1e6;
    }
    if (t <= 0.68)
    {
        return  2.904E-06 * 1e6;
    }
    if (t <= 0.69)
    {
        return  2.915E-06 * 1e6;
    }
    if (t <= 0.70)
    {
        return  2.903E-06 * 1e6;
    }
    if (t <= 0.71)
    {
        return  2.872E-06 * 1e6;
    }
    if (t <= 0.72)
    {
        return  2.826E-06 * 1e6;
    }
    if (t <= 0.73)
    {
        return  2.769E-06 * 1e6;
    }
    if (t <= 0.74)
    {
        return  2.704E-06 * 1e6;
    }
    if (t <= 0.75)
    {
        return  2.636E-06 * 1e6;
    }
    if (t <= 0.76)
    {
        return  2.569E-06 * 1e6;
    }
    if (t <= 0.77)
    {
        return  2.511E-06 * 1e6;
    }
    if (t <= 0.78)
    {
        return  2.465E-06 * 1e6;
    }
    if (t <= 0.79)
    {
        return  2.433E-06 * 1e6;
    }
}// 15, LCCA, branch 2

Real aortaFlux9 (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00)
    {
        return  7.854E-07 * 1e6;
    }
    if (t <= 0.01)
    {
        return  7.900E-07 * 1e6;
    }
    if (t <= 0.02)
    {
        return  7.988E-07 * 1e6;
    }
    if (t <= 0.03)
    {
        return  8.077E-07 * 1e6;
    }
    if (t <= 0.04)
    {
        return  8.131E-07 * 1e6;
    }
    if (t <= 0.05)
    {
        return  8.141E-07 * 1e6;
    }
    if (t <= 0.06)
    {
        return  8.120E-07 * 1e6;
    }
    if (t <= 0.07)
    {
        return  8.085E-07 * 1e6;
    }
    if (t <= 0.08)
    {
        return  8.082E-07 * 1e6;
    }
    if (t <= 0.09)
    {
        return  8.255E-07 * 1e6;
    }
    if (t <= 0.10)
    {
        return  8.927E-07 * 1e6;
    }
    if (t <= 0.11)
    {
        return  1.055E-06 * 1e6;
    }
    if (t <= 0.12)
    {
        return  1.343E-06 * 1e6;
    }
    if (t <= 0.13)
    {
        return  1.728E-06 * 1e6;
    }
    if (t <= 0.14)
    {
        return  2.118E-06 * 1e6;
    }
    if (t <= 0.15)
    {
        return  2.419E-06 * 1e6;
    }
    if (t <= 0.16)
    {
        return  2.598E-06 * 1e6;
    }
    if (t <= 0.17)
    {
        return  2.673E-06 * 1e6;
    }
    if (t <= 0.18)
    {
        return  2.642E-06 * 1e6;
    }
    if (t <= 0.19)
    {
        return  2.465E-06 * 1e6;
    }
    if (t <= 0.20)
    {
        return  2.118E-06 * 1e6;
    }
    if (t <= 0.21)
    {
        return  1.660E-06 * 1e6;
    }
    if (t <= 0.22)
    {
        return  1.210E-06 * 1e6;
    }
    if (t <= 0.23)
    {
        return  8.723E-07 * 1e6;
    }
    if (t <= 0.24)
    {
        return  6.812E-07 * 1e6;
    }
    if (t <= 0.25)
    {
        return  6.155E-07 * 1e6;
    }
    if (t <= 0.26)
    {
        return  6.323E-07 * 1e6;
    }
    if (t <= 0.27)
    {
        return  6.927E-07 * 1e6;
    }
    if (t <= 0.28)
    {
        return  7.708E-07 * 1e6;
    }
    if (t <= 0.29)
    {
        return  8.538E-07 * 1e6;
    }
    if (t <= 0.30)
    {
        return  9.397E-07 * 1e6;
    }
    if (t <= 0.31)
    {
        return  1.030E-06 * 1e6;
    }
    if (t <= 0.32)
    {
        return  1.120E-06 * 1e6;
    }
    if (t <= 0.33)
    {
        return  1.203E-06 * 1e6;
    }
    if (t <= 0.34)
    {
        return  1.263E-06 * 1e6;
    }
    if (t <= 0.35)
    {
        return  1.292E-06 * 1e6;
    }
    if (t <= 0.36)
    {
        return  1.288E-06 * 1e6;
    }
    if (t <= 0.37)
    {
        return  1.258E-06 * 1e6;
    }
    if (t <= 0.38)
    {
        return  1.207E-06 * 1e6;
    }
    if (t <= 0.39)
    {
        return  1.142E-06 * 1e6;
    }
    if (t <= 0.40)
    {
        return  1.064E-06 * 1e6;
    }
    if (t <= 0.41)
    {
        return  9.816E-07 * 1e6;
    }
    if (t <= 0.42)
    {
        return  9.133E-07 * 1e6;
    }
    if (t <= 0.43)
    {
        return  8.870E-07 * 1e6;
    }
    if (t <= 0.44)
    {
        return  9.268E-07 * 1e6;
    }
    if (t <= 0.45)
    {
        return  1.035E-06 * 1e6;
    }
    if (t <= 0.46)
    {
        return  1.183E-06 * 1e6;
    }
    if (t <= 0.47)
    {
        return  1.329E-06 * 1e6;
    }
    if (t <= 0.48)
    {
        return  1.445E-06 * 1e6;
    }
    if (t <= 0.49)
    {
        return  1.529E-06 * 1e6;
    }
    if (t <= 0.50)
    {
        return  1.587E-06 * 1e6;
    }
    if (t <= 0.51)
    {
        return  1.612E-06 * 1e6;
    }
    if (t <= 0.52)
    {
        return  1.587E-06 * 1e6;
    }
    if (t <= 0.53)
    {
        return  1.501E-06 * 1e6;
    }
    if (t <= 0.54)
    {
        return  1.370E-06 * 1e6;
    }
    if (t <= 0.55)
    {
        return  1.230E-06 * 1e6;
    }
    if (t <= 0.56)
    {
        return  1.112E-06 * 1e6;
    }
    if (t <= 0.57)
    {
        return  1.029E-06 * 1e6;
    }
    if (t <= 0.58)
    {
        return  9.776E-07 * 1e6;
    }
    if (t <= 0.59)
    {
        return  9.460E-07 * 1e6;
    }
    if (t <= 0.60)
    {
        return  9.239E-07 * 1e6;
    }
    if (t <= 0.61)
    {
        return  9.061E-07 * 1e6;
    }
    if (t <= 0.62)
    {
        return  8.917E-07 * 1e6;
    }
    if (t <= 0.63)
    {
        return  8.821E-07 * 1e6;
    }
    if (t <= 0.64)
    {
        return  8.786E-07 * 1e6;
    }
    if (t <= 0.65)
    {
        return  8.809E-07 * 1e6;
    }
    if (t <= 0.66)
    {
        return  8.871E-07 * 1e6;
    }
    if (t <= 0.67)
    {
        return  8.940E-07 * 1e6;
    }
    if (t <= 0.68)
    {
        return  8.978E-07 * 1e6;
    }
    if (t <= 0.69)
    {
        return  8.959E-07 * 1e6;
    }
    if (t <= 0.70)
    {
        return  8.874E-07 * 1e6;
    }
    if (t <= 0.71)
    {
        return  8.734E-07 * 1e6;
    }
    if (t <= 0.72)
    {
        return  8.561E-07 * 1e6;
    }
    if (t <= 0.73)
    {
        return  8.383E-07 * 1e6;
    }
    if (t <= 0.74)
    {
        return  8.221E-07 * 1e6;
    }
    if (t <= 0.75)
    {
        return  8.088E-07 * 1e6;
    }
    if (t <= 0.76)
    {
        return  7.986E-07 * 1e6;
    }
    if (t <= 0.77)
    {
        return  7.909E-07 * 1e6;
    }
    if (t <= 0.78)
    {
        return  7.852E-07 * 1e6;
    }
    if (t <= 0.79)
    {
        return  7.807E-07 * 1e6;
    }
}// 20 LVA branch 3_1

Real aortaFlux4 (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00)
    {
        return  2.581E-06 * 1e6;
    }
    if (t <= 0.01)
    {
        return  2.499E-06 * 1e6;
    }
    if (t <= 0.02)
    {
        return  2.399E-06 * 1e6;
    }
    if (t <= 0.03)
    {
        return  2.281E-06 * 1e6;
    }
    if (t <= 0.04)
    {
        return  2.147E-06 * 1e6;
    }
    if (t <= 0.05)
    {
        return  2.002E-06 * 1e6;
    }
    if (t <= 0.06)
    {
        return  1.849E-06 * 1e6;
    }
    if (t <= 0.07)
    {
        return  1.695E-06 * 1e6;
    }
    if (t <= 0.08)
    {
        return  1.565E-06 * 1e6;
    }
    if (t <= 0.09)
    {
        return  1.537E-06 * 1e6;
    }
    if (t <= 0.10)
    {
        return  1.766E-06 * 1e6;
    }
    if (t <= 0.11)
    {
        return  2.460E-06 * 1e6;
    }
    if (t <= 0.12)
    {
        return  3.773E-06 * 1e6;
    }
    if (t <= 0.13)
    {
        return  5.709E-06 * 1e6;
    }
    if (t <= 0.14)
    {
        return  8.131E-06 * 1e6;
    }
    if (t <= 0.15)
    {
        return  1.087E-05 * 1e6;
    }
    if (t <= 0.16)
    {
        return  1.379E-05 * 1e6;
    }
    if (t <= 0.17)
    {
        return  1.675E-05 * 1e6;
    }
    if (t <= 0.18)
    {
        return  1.957E-05 * 1e6;
    }
    if (t <= 0.19)
    {
        return  2.205E-05 * 1e6;
    }
    if (t <= 0.20)
    {
        return  2.403E-05 * 1e6;
    }
    if (t <= 0.21)
    {
        return  2.538E-05 * 1e6;
    }
    if (t <= 0.22)
    {
        return  2.606E-05 * 1e6;
    }
    if (t <= 0.23)
    {
        return  2.611E-05 * 1e6;
    }
    if (t <= 0.24)
    {
        return  2.562E-05 * 1e6;
    }
    if (t <= 0.25)
    {
        return  2.474E-05 * 1e6;
    }
    if (t <= 0.26)
    {
        return  2.362E-05 * 1e6;
    }
    if (t <= 0.27)
    {
        return  2.238E-05 * 1e6;
    }
    if (t <= 0.28)
    {
        return  2.111E-05 * 1e6;
    }
    if (t <= 0.29)
    {
        return  1.981E-05 * 1e6;
    }
    if (t <= 0.30)
    {
        return  1.850E-05 * 1e6;
    }
    if (t <= 0.31)
    {
        return  1.715E-05 * 1e6;
    }
    if (t <= 0.32)
    {
        return  1.576E-05 * 1e6;
    }
    if (t <= 0.33)
    {
        return  1.432E-05 * 1e6;
    }
    if (t <= 0.34)
    {
        return  1.284E-05 * 1e6;
    }
    if (t <= 0.35)
    {
        return  1.132E-05 * 1e6;
    }
    if (t <= 0.36)
    {
        return  9.768E-06 * 1e6;
    }
    if (t <= 0.37)
    {
        return  8.180E-06 * 1e6;
    }
    if (t <= 0.38)
    {
        return  6.543E-06 * 1e6;
    }
    if (t <= 0.39)
    {
        return  4.831E-06 * 1e6;
    }
    if (t <= 0.40)
    {
        return  3.030E-06 * 1e6;
    }
    if (t <= 0.41)
    {
        return  1.163E-06 * 1e6;
    }
    if (t <= 0.42)
    {
        return  -6.817E-07 * 1e6;
    }
    if (t <= 0.43)
    {
        return  -2.362E-06 * 1e6;
    }
    if (t <= 0.44)
    {
        return  -3.738E-06 * 1e6;
    }
    if (t <= 0.45)
    {
        return  -4.742E-06 * 1e6;
    }
    if (t <= 0.46)
    {
        return  -5.400E-06 * 1e6;
    }
    if (t <= 0.47)
    {
        return  -5.785E-06 * 1e6;
    }
    if (t <= 0.48)
    {
        return  -5.956E-06 * 1e6;
    }
    if (t <= 0.49)
    {
        return  -5.932E-06 * 1e6;
    }
    if (t <= 0.50)
    {
        return  -5.723E-06 * 1e6;
    }
    if (t <= 0.51)
    {
        return  -5.358E-06 * 1e6;
    }
    if (t <= 0.52)
    {
        return  -4.889E-06 * 1e6;
    }
    if (t <= 0.53)
    {
        return  -4.370E-06 * 1e6;
    }
    if (t <= 0.54)
    {
        return  -3.846E-06 * 1e6;
    }
    if (t <= 0.55)
    {
        return  -3.341E-06 * 1e6;
    }
    if (t <= 0.56)
    {
        return  -2.866E-06 * 1e6;
    }
    if (t <= 0.57)
    {
        return  -2.422E-06 * 1e6;
    }
    if (t <= 0.58)
    {
        return  -2.004E-06 * 1e6;
    }
    if (t <= 0.59)
    {
        return  -1.601E-06 * 1e6;
    }
    if (t <= 0.60)
    {
        return  -1.206E-06 * 1e6;
    }
    if (t <= 0.61)
    {
        return  -8.123E-07 * 1e6;
    }
    if (t <= 0.62)
    {
        return  -4.195E-07 * 1e6;
    }
    if (t <= 0.63)
    {
        return  -3.198E-08 * 1e6;
    }
    if (t <= 0.64)
    {
        return  3.429E-07 * 1e6;
    }
    if (t <= 0.65)
    {
        return  6.970E-07 * 1e6;
    }
    if (t <= 0.66)
    {
        return  1.024E-06 * 1e6;
    }
    if (t <= 0.67)
    {
        return  1.320E-06 * 1e6;
    }
    if (t <= 0.68)
    {
        return  1.583E-06 * 1e6;
    }
    if (t <= 0.69)
    {
        return  1.813E-06 * 1e6;
    }
    if (t <= 0.70)
    {
        return  2.011E-06 * 1e6;
    }
    if (t <= 0.71)
    {
        return  2.180E-06 * 1e6;
    }
    if (t <= 0.72)
    {
        return  2.323E-06 * 1e6;
    }
    if (t <= 0.73)
    {
        return  2.442E-06 * 1e6;
    }
    if (t <= 0.74)
    {
        return  2.539E-06 * 1e6;
    }
    if (t <= 0.75)
    {
        return  2.613E-06 * 1e6;
    }
    if (t <= 0.76)
    {
        return  2.665E-06 * 1e6;
    }
    if (t <= 0.77)
    {
        return  2.691E-06 * 1e6;
    }
    if (t <= 0.78)
    {
        return  2.690E-06 * 1e6;
    }
    if (t <= 0.79)
    {
        return  2.661E-06 * 1e6;
    }
}//21, L. Brachia, bhanch 3_2







Real f (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.;
}

Real u1 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

// Initial velocity
Real u0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real p0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}


Real E (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 1000;
    //   if(z<5)
    //       return  /*-*/10*1e3;//-29;//*e-1*5; // about [(110-60)*(133.332*10)]/[10*(2.08-1.8)] ( 10 because of mm, 133... because of mmHg)
    //   //    (see paper by Liu, Dang, etc). in their plot x and y are inverted.
    //   if(z>5 && z<6)
    //       return 10+20*(z-5)*1e3;
    //   if(z>=6)
    //       return 30*1e3;
    // if ( z<0 && z>=-3 )
    //   return 10-3*(z)*1e3;
    // if(z<-3)
    //   return 19*1e3;
}


// Initial displacement and velocity
Real d0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
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
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            break;
    }
}

Real w0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
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
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            break;
    }
}



Real linearPress2 ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaPhisPress (tii) - aortaPhisPress (ti) ) / (0.001) * (t - (ti) ) + aortaPhisPress (ti) + 115000;
}


Real linearFlux3 ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux3 (tii) - aortaFlux3 (ti) ) / (0.001) * (t - (ti) ) + aortaFlux3 (ti);
}


Real linearFluxIn (Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;

    return (aortaFluxIn (tii) - aortaFluxIn (ti) ) / (0.001) * (t - (ti) ) + aortaFluxIn (ti);
}


Real linearFlux4 (Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux4 (tii) - aortaFlux4 (ti) ) / (0.001) * (t - (ti) ) + aortaFlux4 (ti);
}



Real linearFlux5 ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux5 (tii) - aortaFlux5 (ti) ) / (0.001) * (t - (ti) ) + aortaFlux5 (ti);
}



Real linearFlux6 ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux6 (tii) - aortaFlux6 (ti) ) / (0.001) * (t - (ti) ) + aortaFlux6 (ti);
}


Real linearFlux7 ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux7 (tii) - aortaFlux7 (ti) ) / (0.001) * (t - (ti) ) + aortaFlux7 (ti);
}


Real linearFlux8 (Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux8 (tii) - aortaFlux8 (ti) ) / (0.001) * (t - (ti) ) + aortaFlux8 (ti);
}

Real linearFlux9 ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux9 (tii) - aortaFlux9 (ti) ) / (0.001) * (t - (ti) ) + aortaFlux9 (ti);
}


Real linearFlux3_ (Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux3_ (tii) - aortaFlux3_ (ti) ) / (0.001) * (t - (ti) ) + aortaFlux3_ (ti);
}


Real linearFlux6_ ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000) ) % 800) / 1000;
    }

    Real ti = floor (t * 1000) / 1000;
    Real tii = ti + 0.001;
    return (aortaFlux6_ (tii) - aortaFlux6_ (ti) ) / (0.001) * (t - (ti) ) + aortaFlux6_ (ti);
}

Real aortaFlux3_ (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00 + 0.01)
    {
        return  33.6             ;
    }
    if (t <= 0.01 + 0.01)
    {
        return  33.83            ;
    }
    if (t <= 0.02 + 0.01)
    {
        return  34.12            ;
    }
    if (t <= 0.03 + 0.01)
    {
        return  34.44            ;
    }
    if (t <= 0.04 + 0.01)
    {
        return  34.76            ;
    }
    if (t <= 0.05 + 0.01)
    {
        return  35.05            ;
    }
    if (t <= 0.06 + 0.01)
    {
        return  35.29            ;
    }
    if (t <= 0.07 + 0.01)
    {
        return  35.49            ;
    }
    if (t <= 0.08 + 0.01)
    {
        return  35.74            ;
    }
    if (t <= 0.09 + 0.01)
    {
        return  36.27            ;
    }
    if (t <= 0.10 + 0.01)
    {
        return  37.83            ;
    }
    if (t <= 0.11 + 0.01)
    {
        return  41.86            ;
    }
    if (t <= 0.12 + 0.01)
    {
        return  50.47            ;
    }
    if (t <= 0.13 + 0.01)
    {
        return  65.55999999999999;
    }
    if (t <= 0.14 + 0.01)
    {
        return  87.10999999999999;
    }
    if (t <= 0.15 + 0.01)
    {
        return  111.8            ;
    }
    if (t <= 0.16 + 0.01)
    {
        return  134.6            ;
    }
    if (t <= 0.17 + 0.01)
    {
        return  152              ;
    }
    if (t <= 0.18 + 0.01)
    {
        return  164.3            ;
    }
    if (t <= 0.19 + 0.01)
    {
        return  173.8            ;
    }
    if (t <= 0.20 + 0.01)
    {
        return  182.2            ;
    }
    if (t <= 0.21 + 0.01)
    {
        return  189.8            ;
    }
    if (t <= 0.22 + 0.01)
    {
        return  195.9            ;
    }
    if (t <= 0.23 + 0.01)
    {
        return  200              ;
    }
    if (t <= 0.24 + 0.01)
    {
        return  202              ;
    }
    if (t <= 0.25 + 0.01)
    {
        return  202.1            ;
    }
    if (t <= 0.26 + 0.01)
    {
        return  200.5            ;
    }
    if (t <= 0.27 + 0.01)
    {
        return  197.3            ;
    }
    if (t <= 0.28 + 0.01)
    {
        return  192.7            ;
    }
    if (t <= 0.29 + 0.01)
    {
        return  186.9            ;
    }
    if (t <= 0.30 + 0.01)
    {
        return  179.9            ;
    }
    if (t <= 0.31 + 0.01)
    {
        return  172              ;
    }
    if (t <= 0.32 + 0.01)
    {
        return  163.3            ;
    }
    if (t <= 0.33 + 0.01)
    {
        return  154.1            ;
    }
    if (t <= 0.34 + 0.01)
    {
        return  144.7            ;
    }
    if (t <= 0.35 + 0.01)
    {
        return  135.1            ;
    }
    if (t <= 0.36 + 0.01)
    {
        return  125.5            ;
    }
    if (t <= 0.37 + 0.01)
    {
        return  115.9            ;
    }
    if (t <= 0.38 + 0.01)
    {
        return  106.3            ;
    }
    if (t <= 0.39 + 0.01)
    {
        return  96.60999999999999;
    }
    if (t <= 0.40 + 0.01)
    {
        return  86.44            ;
    }
    if (t <= 0.41 + 0.01)
    {
        return  75.67999999999999;
    }
    if (t <= 0.42 + 0.01)
    {
        return  64.56999999999999;
    }
    if (t <= 0.43 + 0.01)
    {
        return  53.92            ;
    }
    if (t <= 0.44 + 0.01)
    {
        return  44.94            ;
    }
    if (t <= 0.45 + 0.01)
    {
        return  38.75            ;
    }
    if (t <= 0.46 + 0.01)
    {
        return  35.66            ;
    }
    if (t <= 0.47 + 0.01)
    {
        return  34.83000000000001;
    }
    if (t <= 0.48 + 0.01)
    {
        return  34.67            ;
    }
    if (t <= 0.49 + 0.01)
    {
        return  33.83            ;
    }
    if (t <= 0.50 + 0.01)
    {
        return  32.01            ;
    }
    if (t <= 0.51 + 0.01)
    {
        return  29.85            ;
    }
    if (t <= 0.52 + 0.01)
    {
        return  28.19            ;
    }
    if (t <= 0.53 + 0.01)
    {
        return  27.4             ;
    }
    if (t <= 0.54 + 0.01)
    {
        return  27.28            ;
    }
    if (t <= 0.55 + 0.01)
    {
        return  27.42            ;
    }
    if (t <= 0.56 + 0.01)
    {
        return  27.54            ;
    }
    if (t <= 0.57 + 0.01)
    {
        return  27.63            ;
    }
    if (t <= 0.58 + 0.01)
    {
        return  27.82            ;
    }
    if (t <= 0.59 + 0.01)
    {
        return  28.21            ;
    }
    if (t <= 0.60 + 0.01)
    {
        return  28.76            ;
    }
    if (t <= 0.61 + 0.01)
    {
        return  29.39            ;
    }
    if (t <= 0.62 + 0.01)
    {
        return  30               ;
    }
    if (t <= 0.63 + 0.01)
    {
        return  30.54            ;
    }
    if (t <= 0.64 + 0.01)
    {
        return  31.01            ;
    }
    if (t <= 0.65 + 0.01)
    {
        return  31.42            ;
    }
    if (t <= 0.66 + 0.01)
    {
        return  31.78            ;
    }
    if (t <= 0.67 + 0.01)
    {
        return  32.09            ;
    }
    if (t <= 0.68 + 0.01)
    {
        return  32.34            ;
    }
    if (t <= 0.69 + 0.01)
    {
        return  32.54            ;
    }
    if (t <= 0.70 + 0.01)
    {
        return  32.69            ;
    }
    if (t <= 0.71 + 0.01)
    {
        return  32.8             ;
    }
    if (t <= 0.72 + 0.01)
    {
        return  32.87            ;
    }
    if (t <= 0.73 + 0.01)
    {
        return  32.91            ;
    }
    if (t <= 0.74 + 0.01)
    {
        return  32.92            ;
    }
    if (t <= 0.75 + 0.01)
    {
        return  32.93000000000001;
    }
    if (t <= 0.76 + 0.01)
    {
        return  32.93000000000001;
    }
    if (t <= 0.77 + 0.01)
    {
        return  32.96            ;
    }
    if (t <= 0.78 + 0.01)
    {
        return  33.01000000000001;
    }
    if (t <= 0.79 + 0.01)
    {
        return  33.1             ;
    }
    RETURN_UNDEFINED;
}

Real aortaFlux6_ (const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0.00 + 0.01)
    {
        return    0.6737071250000001;
    }
    if (t <= 0.01 + 0.01)
    {
        return    0.6799071250000001;
    }
    if (t <= 0.02 + 0.01)
    {
        return    0.6897071250000001;
    }
    if (t <= 0.03 + 0.01)
    {
        return    0.6997071250000001;
    }
    if (t <= 0.04 + 0.01)
    {
        return    0.706407125       ;
    }
    if (t <= 0.05 + 0.01)
    {
        return    0.7084071250000001;
    }
    if (t <= 0.06 + 0.01)
    {
        return    0.7056071250000001;
    }
    if (t <= 0.07 + 0.01)
    {
        return    0.7049071250000001;
    }
    if (t <= 0.08 + 0.01)
    {
        return    0.7382071250000001;
    }
    if (t <= 0.09 + 0.01)
    {
        return    0.8778071249999999;
    }
    if (t <= 0.10 + 0.01)
    {
        return    1.209007125       ;
    }
    if (t <= 0.11 + 0.01)
    {
        return    1.746007125       ;
    }
    if (t <= 0.12 + 0.01)
    {
        return    2.368007125       ;
    }
    if (t <= 0.13 + 0.01)
    {
        return    2.897007125       ;
    }
    if (t <= 0.14 + 0.01)
    {
        return    3.258007125       ;
    }
    if (t <= 0.15 + 0.01)
    {
        return    3.518007125       ;
    }
    if (t <= 0.16 + 0.01)
    {
        return    3.752007125       ;
    }
    if (t <= 0.17 + 0.01)
    {
        return    3.936007124999999 ;
    }
    if (t <= 0.18 + 0.01)
    {
        return    3.957007125       ;
    }
    if (t <= 0.19 + 0.01)
    {
        return    3.716007125       ;
    }
    if (t <= 0.20 + 0.01)
    {
        return    3.212007125       ;
    }
    if (t <= 0.21 + 0.01)
    {
        return    2.551007125       ;
    }
    if (t <= 0.22 + 0.01)
    {
        return    1.898007125       ;
    }
    if (t <= 0.23 + 0.01)
    {
        return    1.395007125       ;
    }
    if (t <= 0.24 + 0.01)
    {
        return    1.107007125       ;
    }
    if (t <= 0.25 + 0.01)
    {
        return    1.009007125       ;
    }
    if (t <= 0.26 + 0.01)
    {
        return    1.028007125       ;
    }
    if (t <= 0.27 + 0.01)
    {
        return    1.082007125       ;
    }
    if (t <= 0.28 + 0.01)
    {
        return    1.114007125       ;
    }
    if (t <= 0.29 + 0.01)
    {
        return    1.111007125       ;
    }
    if (t <= 0.30 + 0.01)
    {
        return    1.089007125       ;
    }
    if (t <= 0.31 + 0.01)
    {
        return    1.071007125       ;
    }
    if (t <= 0.32 + 0.01)
    {
        return    1.067007125       ;
    }
    if (t <= 0.33 + 0.01)
    {
        return    1.068007125       ;
    }
    if (t <= 0.34 + 0.01)
    {
        return    1.056007125       ;
    }
    if (t <= 0.35 + 0.01)
    {
        return    1.021007125       ;
    }
    if (t <= 0.36 + 0.01)
    {
        return    0.955007125       ;
    }
    if (t <= 0.37 + 0.01)
    {
        return    0.856707125       ;
    }
    if (t <= 0.38 + 0.01)
    {
        return    0.7230071250000001;
    }
    if (t <= 0.39 + 0.01)
    {
        return    0.5555071250000001;
    }
    if (t <= 0.40 + 0.01)
    {
        return    0.3694071249999999;
    }
    if (t <= 0.41 + 0.01)
    {
        return    0.203607125       ;
    }
    if (t <= 0.42 + 0.01)
    {
        return    0.117107125       ;
    }
    if (t <= 0.43 + 0.01)
    {
        return    0.154707125       ;
    }
    if (t <= 0.44 + 0.01)
    {
        return    0.3019071249999999;
    }
    if (t <= 0.45 + 0.01)
    {
        return    0.4833071250000001;
    }
    if (t <= 0.46 + 0.01)
    {
        return    0.6279071250000001;
    }
    if (t <= 0.47 + 0.01)
    {
        return    0.7323071250000001;
    }
    if (t <= 0.48 + 0.01)
    {
        return    0.8435071250000001;
    }
    if (t <= 0.49 + 0.01)
    {
        return    0.989007125       ;
    }
    if (t <= 0.50 + 0.01)
    {
        return    1.137007125       ;
    }
    if (t <= 0.51 + 0.01)
    {
        return    1.223007125       ;
    }
    if (t <= 0.52 + 0.01)
    {
        return    1.206007125       ;
    }
    if (t <= 0.53 + 0.01)
    {
        return    1.096007125       ;
    }
    if (t <= 0.54 + 0.01)
    {
        return    0.9400071249999998;
    }
    if (t <= 0.55 + 0.01)
    {
        return    0.7924071250000001;
    }
    if (t <= 0.56 + 0.01)
    {
        return    0.6857071250000001;
    }
    if (t <= 0.57 + 0.01)
    {
        return    0.6278071250000001;
    }
    if (t <= 0.58 + 0.01)
    {
        return    0.6083071250000001;
    }
    if (t <= 0.59 + 0.01)
    {
        return    0.610607125       ;
    }
    if (t <= 0.60 + 0.01)
    {
        return    0.620007125       ;
    }
    if (t <= 0.61 + 0.01)
    {
        return    0.6279071250000001;
    }
    if (t <= 0.62 + 0.01)
    {
        return    0.6320071250000001;
    }
    if (t <= 0.63 + 0.01)
    {
        return    0.6343071250000001;
    }
    if (t <= 0.64 + 0.01)
    {
        return    0.6385071250000001;
    }
    if (t <= 0.65 + 0.01)
    {
        return    0.6476071250000001;
    }
    if (t <= 0.66 + 0.01)
    {
        return    0.6620071250000001;
    }
    if (t <= 0.67 + 0.01)
    {
        return    0.679107125       ;
    }
    if (t <= 0.68 + 0.01)
    {
        return    0.6948071250000001;
    }
    if (t <= 0.69 + 0.01)
    {
        return    0.7049071250000001;
    }
    if (t <= 0.70 + 0.01)
    {
        return    0.7071071250000001;
    }
    if (t <= 0.71 + 0.01)
    {
        return    0.7021071250000001;
    }
    if (t <= 0.72 + 0.01)
    {
        return    0.6926071250000001;
    }
    if (t <= 0.73 + 0.01)
    {
        return    0.6825071250000001;
    }
    if (t <= 0.74 + 0.01)
    {
        return    0.675007125       ;
    }
    if (t <= 0.75 + 0.01)
    {
        return    0.671207125       ;
    }
    if (t <= 0.76 + 0.01)
    {
        return    0.6702071250000001;
    }
    if (t <= 0.77 + 0.01)
    {
        return    0.6702071250000001;
    }
    if (t <= 0.78 + 0.01)
    {
        return    0.669507125       ;
    }
    if (t <= 0.79 + 0.01)
    {
        return    0.6674071250000001;
    }
    RETURN_UNDEFINED;
}

Real u2 (Real  /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -1e4;
}

}
