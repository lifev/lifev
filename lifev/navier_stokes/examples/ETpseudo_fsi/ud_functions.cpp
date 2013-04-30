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

Real fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}


Real E (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 500.0;
}

Real inletCylinder (Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t <= 0   + 0.02)
    {
        return      4.572358 ;
    }
    if (t <= 0.02 + 0.02)
    {
        return        4.380298   ;
    }
    if (t <= 0.04 + 0.02)
    {
        return        3.971077 ;
    }
    if (t <= 0.06 + 0.02)
    {
        return        3.789365   ;
    }
    if (t <= 0.08 + 0.02)
    {
        return        4.394491   ;
    }
    if (t <= 0.1 + 0.02)
    {
        return      6.048268 ;
    }
    if (t <= 0.12 + 0.02)
    {
        return        8.463103   ;
    }
    if (t <= 0.14 + 0.02)
    {
        return        10.925368;
    }
    if (t <= 0.16 + 0.02)
    {
        return        12.729437;
    }
    if (t <= 0.18 + 0.02)
    {
        return        13.572193;
    }
    if (t <= 0.2 + 0.02)
    {
        return      13.62635 ;
    }
    if (t <= 0.22 + 0.02)
    {
        return        13.300649;
    }
    if (t <= 0.24 + 0.02)
    {
        return        12.900495;
    }
    if (t <= 0.26 + 0.02)
    {
        return        12.473357;
    }
    if (t <= 0.28 + 0.02)
    {
        return        11.895968;
    }
    if (t <= 0.3 + 0.02)
    {
        return      11.086719    ;
    }
    if (t <= 0.32 + 0.02)
    {
        return        10.114592;
    }
    if (t <= 0.34 + 0.02)
    {
        return        9.149018   ;
    }
    if (t <= 0.36 + 0.02)
    {
        return        8.307111   ;
    }
    if (t <= 0.38 + 0.02)
    {
        return        7.585431   ;
    }
    if (t <= 0.4 + 0.02)
    {
        return      6.924356 ;
    }
    if (t <= 0.42 + 0.02)
    {
        return        6.338303   ;
    }
    if (t <= 0.44 + 0.02)
    {
        return        5.948785   ;
    }
    if (t <= 0.46 + 0.02)
    {
        return        5.880324   ;
    }
    if (t <= 0.48 + 0.02)
    {
        return        6.106962   ;
    }
    if (t <= 0.5 + 0.02)
    {
        return      6.41436  ;
    }
    if (t <= 0.52 + 0.02)
    {
        return        6.530911   ;
    }
    if (t <= 0.54 + 0.02)
    {
        return        6.329742   ;
    }
    if (t <= 0.56 + 0.02)
    {
        return        5.911398   ;
    }
    if (t <= 0.58 + 0.02)
    {
        return        5.501446   ;
    }
    if (t <= 0.6 + 0.02)
    {
        return      5.253072 ;
    }
    if (t <= 0.62 + 0.02)
    {
        return        5.138401   ;
    }
    if (t <= 0.64 + 0.02)
    {
        return        5.013553   ;
    }
    if (t <= 0.66 + 0.02)
    {
        return        4.78726    ;
    }
    if (t <= 0.68 + 0.02)
    {
        return        4.515227   ;
    }
    if (t <= 0.7 + 0.02)
    {
        return      4.341614 ;
    }
    if (t <= 0.72 + 0.02)
    {
        return        4.341829   ;
    }
    if (t <= 0.74 + 0.02)
    {
        return        4.346171   ;
    }
    if (t <= 0.76 + 0.02)
    {
        return        4.460936   ;
    }
    if (t <= 0.78 + 0.02)
    {
        return        4.535817   ;
    }
    if (t <= 0.8 + 0.02)
    {
        return      4.485449 ;
    }
    if (t <= 0.82 + 0.02)
    {
        return        4.36465    ;
    }
    if (t <= 0.84 + 0.02)
    {
        return        4.302813   ;
    }
    if (t <= 0.86 + 0.02)
    {
        return        4.363062   ;
    }
    if (t <= 0.88 + 0.02)
    {
        return        4.476906   ;
    }
    if (t <= 0.9 + 0.02)
    {
        return      4.51795  ;
    }
    if (t <= 0.92 + 0.02)
    {
        return        4.438247   ;
    }
    if (t <= 0.94 + 0.02)
    {
        return        4.319594   ;
    }
    if (t <= 0.96 + 0.02)
    {
        return        4.286588   ;
    }
    if (t <= 0.98 + 0.02)
    {
        return        4.368411   ;
    }
    if (t <= 1   + 0.02)
    {
        return      4.461585 ;
    }
    if (t <= 1.02 + 0.02)
    {
        return        4.439973   ;
    }
    if (t <= 1.04 + 0.02)
    {
        return        4.302693   ;
    }
    if (t <= 1.06 + 0.02)
    {
        return        4.194866   ;
    }
    if (t <= 1.08 + 0.02)
    {
        return        4.260641   ;
    }
    if (t <= 1.1 + 0.02)
    {
        return      4.462993 ;
    }
}

Real linearInletCylinder ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{

    Real tNew = t;
    if (t > 1.102)
    {
        tNew = ( ( (int) floor (t * 10000) ) % 11000) / 10000.0;
    }
    //Real dt =  static_cast<int>( std::floor(tNew * 100) )% 2 / 100.0;
    Real dt = ( int(floor (tNew * 1000.0) ) % 20 )/ 1000.0;
    Real ti = tNew - dt;
    Real tii = ti + 0.02;
    return - ( (inletCylinder (tii, 0, 0, 0, 0) - inletCylinder (ti, 0, 0, 0, 0) ) / (0.02) * (tNew - (ti) ) + inletCylinder (ti, 0, 0, 0, 0) );
}


Real linearVelInletCylinder ( Real  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    //Components for Simone's mesh

    Real n1 (0.0);
    Real n2 (0.0);
    Real n3 (-1.0);

    Real flux (1.0 * linearInletCylinder (t, x, y, z, i) );
    //Center for the aneurym with one IO
    Real x0 (0);
    Real y0 (0);
    Real area (0.3 * 0.3 * 3.141592653589793); //fluidBig

    Real radiusSquared (0);
    radiusSquared = 0.5 * 0.5; //* 3.1415962;

    Real peak (0);
    peak = ( 2.0  * flux ) / ( area );

    switch (i)
    {
        case 0:
            return n1 * ( peak * std::max ( 0.0, ( (radiusSquared - ( (x - x0) * (x - x0) + (y - y0) * (y - y0) ) ) / radiusSquared) ) );
        case 1:
            return n2 * ( peak * std::max ( 0.0, ( (radiusSquared - ( (x - x0) * (x - x0) + (y - y0) * (y - y0) ) ) / radiusSquared) ) );
        case 2:
            return n3 * ( peak * std::max ( 0.0, ( (radiusSquared - ( (x - x0) * (x - x0) + (y - y0) * (y - y0) ) ) / radiusSquared) ) );
        default:
            return 0.0;
    }
}



Real pont_dist (const Real t, const Real& , const Real& , const Real& , const ID& i) //inlet flux from the bypass
{
    Real newt = ( (Real) ( ( (int) (t * 1000) ) % 792) ) / 1000;
    if ( newt <= 0.    + 0.004 )
    {
        return  1.96970000e-02   * (-10);
    }
    else if ( newt <= 0.004 + 0.004 )
    {
        return  1.93150000e-02   * (-10);
    }
    else if ( newt <= 0.008 + 0.004 )
    {
        return  2.12310000e-02   * (-10);
    }
    else if ( newt <= 0.012 + 0.004 )
    {
        return  3.04300000e-02   * (-10);
    }
    else if ( newt <= 0.016 + 0.004 )
    {
        return  4.50000000e-02   * (-10);
    }
    else if ( newt <= 0.02  + 0.004 )
    {
        return  6.18670000e-02   * (-10);
    }
    else if ( newt <= 0.024 + 0.004 )
    {
        return  8.48730000e-02   * (-10);
    }
    else if ( newt <= 0.028 + 0.004 )
    {
        return  1.07874000e-01   * (-10);
    }
    else if ( newt <= 0.032 + 0.004 )
    {
        return  1.32410000e-01   * (-10);
    }
    else if ( newt <= 0.036 + 0.004 )
    {
        return  1.57714000e-01   * (-10);
    }
    else if ( newt <= 0.04  + 0.004 )
    {
        return  1.86849000e-01   * (-10);
    }
    else if ( newt <= 0.044 + 0.004 )
    {
        return  2.19821000e-01   * (-10);
    }
    else if ( newt <= 0.048 + 0.004 )
    {
        return  2.52788000e-01   * (-10);
    }
    else if ( newt <= 0.052 + 0.004 )
    {
        return  2.85760000e-01   * (-10);
    }
    else if ( newt <= 0.056 + 0.004 )
    {
        return  3.17198000e-01   * (-10);
    }
    else if ( newt <= 0.06  + 0.004 )
    {
        return  3.44798000e-01   * (-10);
    }
    else if ( newt <= 0.064 + 0.004 )
    {
        return  3.70869000e-01   * (-10);
    }
    else if ( newt <= 0.068 + 0.004 )
    {
        return  3.96959000e-01   * (-10);
    }
    else if ( newt <= 0.072 + 0.004 )
    {
        return  4.26080000e-01   * (-10);
    }
    else if ( newt <= 0.076 + 0.004 )
    {
        return  4.53675000e-01   * (-10);
    }
    else if ( newt <= 0.08  + 0.004 )
    {
        return  4.82034000e-01   * (-10);
    }
    else if ( newt <= 0.084 + 0.004 )
    {
        return  5.09630000e-01   * (-10);
    }
    else if ( newt <= 0.088 + 0.004 )
    {
        return  5.27264000e-01   * (-10);
    }
    else if ( newt <= 0.092 + 0.004 )
    {
        return  5.38798000e-01   * (-10);
    }
    else if ( newt <= 0.096 + 0.004 )
    {
        return  5.48759000e-01   * (-10);
    }
    else if ( newt <= 0.1   + 0.004 )
    {
        return  5.52572000e-01   * (-10);
    }
    else if ( newt <= 0.104 + 0.004 )
    {
        return  5.63344000e-01   * (-10);
    }
    else if ( newt <= 0.108 + 0.004 )
    {
        return  5.80978000e-01   * (-10);
    }
    else if ( newt <= 0.112 + 0.004 )
    {
        return  6.00138000e-01   * (-10);
    }
    else if ( newt <= 0.116 + 0.004 )
    {
        return  6.19298000e-01   * (-10);
    }
    else if ( newt <= 0.12  + 0.004 )
    {
        return  6.32310000e-01   * (-10);
    }
    else if ( newt <= 0.124 + 0.004 )
    {
        return  6.40746000e-01   * (-10);
    }
    else if ( newt <= 0.128 + 0.004 )
    {
        return  6.41508000e-01   * (-10);
    }
    else if ( newt <= 0.132 + 0.004 )
    {
        return  6.43081000e-01   * (-10);
    }
    else if ( newt <= 0.136 + 0.004 )
    {
        return  6.49944000e-01   * (-10);
    }
    else if ( newt <= 0.14  + 0.004 )
    {
        return  6.59190000e-01   * (-10);
    }
    else if ( newt <= 0.144 + 0.004 )
    {
        return  6.64528000e-01   * (-10);
    }
    else if ( newt <= 0.148 + 0.004 )
    {
        return  6.60716000e-01   * (-10);
    }
    else if ( newt <= 0.152 + 0.004 )
    {
        return  6.56092000e-01   * (-10);
    }
    else if ( newt <= 0.156 + 0.004 )
    {
        return  6.59190000e-01   * (-10);
    }
    else if ( newt <= 0.16  + 0.004 )
    {
        return  6.62241000e-01   * (-10);
    }
    else if ( newt <= 0.164 + 0.004 )
    {
        return  6.62241000e-01   * (-10);
    }
    else if ( newt <= 0.168 + 0.004 )
    {
        return  6.55330000e-01   * (-10);
    }
    else if ( newt <= 0.172 + 0.004 )
    {
        return  6.46131000e-01   * (-10);
    }
    else if ( newt <= 0.176 + 0.004 )
    {
        return  6.35408000e-01   * (-10);
    }
    else if ( newt <= 0.18  + 0.004 )
    {
        return  6.23111000e-01   * (-10);
    }
    else if ( newt <= 0.184 + 0.004 )
    {
        return  6.17010000e-01   * (-10);
    }
    else if ( newt <= 0.188 + 0.004 )
    {
        return  6.10862000e-01   * (-10);
    }
    else if ( newt <= 0.192 + 0.004 )
    {
        return  5.99376000e-01   * (-10);
    }
    else if ( newt <= 0.196 + 0.004 )
    {
        return  5.86316000e-01   * (-10);
    }
    else if ( newt <= 0.2   + 0.004 )
    {
        return  5.70970000e-01   * (-10);
    }
    else if ( newt <= 0.204 + 0.004 )
    {
        return  5.54145000e-01   * (-10);
    }
    else if ( newt <= 0.208 + 0.004 )
    {
        return  5.43374000e-01   * (-10);
    }
    else if ( newt <= 0.212 + 0.004 )
    {
        return  5.35700000e-01   * (-10);
    }
    else if ( newt <= 0.216 + 0.004 )
    {
        return  5.26502000e-01   * (-10);
    }
    else if ( newt <= 0.22  + 0.004 )
    {
        return  5.11965000e-01   * (-10);
    }
    else if ( newt <= 0.224 + 0.004 )
    {
        return  4.91232000e-01   * (-10);
    }
    else if ( newt <= 0.228 + 0.004 )
    {
        return  4.69022000e-01   * (-10);
    }
    else if ( newt <= 0.232 + 0.004 )
    {
        return  4.46002000e-01   * (-10);
    }
    else if ( newt <= 0.236 + 0.004 )
    {
        return  4.27605000e-01   * (-10);
    }
    else if ( newt <= 0.24  + 0.004 )
    {
        return  4.10733000e-01   * (-10);
    }
    else if ( newt <= 0.244 + 0.004 )
    {
        return  3.86950000e-01   * (-10);
    }
    else if ( newt <= 0.248 + 0.004 )
    {
        return  3.65502000e-01   * (-10);
    }
    else if ( newt <= 0.252 + 0.004 )
    {
        return  3.42501000e-01   * (-10);
    }
    else if ( newt <= 0.256 + 0.004 )
    {
        return  3.24862000e-01   * (-10);
    }
    else if ( newt <= 0.26  + 0.004 )
    {
        return  3.07994000e-01   * (-10);
    }
    else if ( newt <= 0.264 + 0.004 )
    {
        return  2.88058000e-01   * (-10);
    }
    else if ( newt <= 0.268 + 0.004 )
    {
        return  2.69656000e-01   * (-10);
    }
    else if ( newt <= 0.272 + 0.004 )
    {
        return  2.47422000e-01   * (-10);
    }
    else if ( newt <= 0.276 + 0.004 )
    {
        return  2.23653000e-01   * (-10);
    }
    else if ( newt <= 0.28  + 0.004 )
    {
        return  2.04484000e-01   * (-10);
    }
    else if ( newt <= 0.284 + 0.004 )
    {
        return  1.84547000e-01   * (-10);
    }
    else if ( newt <= 0.288 + 0.004 )
    {
        return  1.67680000e-01   * (-10);
    }
    else if ( newt <= 0.292 + 0.004 )
    {
        return  1.54644000e-01   * (-10);
    }
    else if ( newt <= 0.296 + 0.004 )
    {
        return  1.42376000e-01   * (-10);
    }
    else if ( newt <= 0.3   + 0.004 )
    {
        return  1.32410000e-01   * (-10);
    }
    else if ( newt <= 0.304 + 0.004 )
    {
        return  1.20142000e-01   * (-10);
    }
    else if ( newt <= 0.308 + 0.004 )
    {
        return  1.06340000e-01   * (-10);
    }
    else if ( newt <= 0.312 + 0.004 )
    {
        return  9.25370000e-02   * (-10);
    }
    else if ( newt <= 0.316 + 0.004 )
    {
        return  7.64370000e-02   * (-10);
    }
    else if ( newt <= 0.32  + 0.004 )
    {
        return  5.88030000e-02   * (-10);
    }
    else if ( newt <= 0.324 + 0.004 )
    {
        return  9.73060000e-03   * (-10);
    }
    else if ( newt <= 0.328 + 0.004 )
    {
        return  -3.09101000e-02  * (-10);
    }
    else if ( newt <= 0.332 + 0.004 )
    {
        return  -6.84790000e-02  * (-10);
    }
    else if ( newt <= 0.336 + 0.004 )
    {
        return  -1.58956600e-01  * (-10);
    }
    else if ( newt <= 0.34  + 0.004 )
    {
        return  -2.08029000e-01  * (-10);
    }
    else if ( newt <= 0.344 + 0.004 )
    {
        return  -2.37931000e-01  * (-10);
    }
    else if ( newt <= 0.348 + 0.004 )
    {
        return  -2.68601000e-01  * (-10);
    }
    else if ( newt <= 0.352 + 0.004 )
    {
        return  -2.71666000e-01  * (-10);
    }
    else if ( newt <= 0.356 + 0.004 )
    {
        return  -2.74735000e-01  * (-10);
    }
    else if ( newt <= 0.36  + 0.004 )
    {
        return  -2.70136000e-01  * (-10);
    }
    else if ( newt <= 0.364 + 0.004 )
    {
        return  -2.65532000e-01  * (-10);
    }
    else if ( newt <= 0.368 + 0.004 )
    {
        return  -2.59398000e-01  * (-10);
    }
    else if ( newt <= 0.372 + 0.004 )
    {
        return  -2.59398000e-01  * (-10);
    }
    else if ( newt <= 0.376 + 0.004 )
    {
        return  -2.60165000e-01  * (-10);
    }
    else if ( newt <= 0.38  + 0.004 )
    {
        return  -2.60932000e-01  * (-10);
    }
    else if ( newt <= 0.384 + 0.004 )
    {
        return  -2.60165000e-01  * (-10);
    }
    else if ( newt <= 0.388 + 0.004 )
    {
        return  -2.57101000e-01  * (-10);
    }
    else if ( newt <= 0.392 + 0.004 )
    {
        return  -2.58635000e-01  * (-10);
    }
    else if ( newt <= 0.396 + 0.004 )
    {
        return  -2.56333000e-01  * (-10);
    }
    else if ( newt <= 0.4   + 0.004 )
    {
        return  -2.49432000e-01  * (-10);
    }
    else if ( newt <= 0.404 + 0.004 )
    {
        return  -2.43298000e-01  * (-10);
    }
    else if ( newt <= 0.408 + 0.004 )
    {
        return  -2.23361000e-01  * (-10);
    }
    else if ( newt <= 0.412 + 0.004 )
    {
        return  -2.09563000e-01  * (-10);
    }
    else if ( newt <= 0.416 + 0.004 )
    {
        return  -1.71992000e-01  * (-10);
    }
    else if ( newt <= 0.42  + 0.004 )
    {
        return  -1.22151700e-01  * (-10);
    }
    else if ( newt <= 0.424 + 0.004 )
    {
        return  -9.14815000e-02  * (-10);
    }
    else if ( newt <= 0.428 + 0.004 )
    {
        return  -6.08113000e-02  * (-10);
    }
    else if ( newt <= 0.432 + 0.004 )
    {
        return  -5.69779000e-02  * (-10);
    }
    else if ( newt <= 0.436 + 0.004 )
    {
        return  -5.16103000e-02  * (-10);
    }
    else if ( newt <= 0.44  + 0.004 )
    {
        return  -3.93414000e-02  * (-10);
    }
    else if ( newt <= 0.444 + 0.004 )
    {
        return  -1.78748000e-02  * (-10);
    }
    else if ( newt <= 0.448 + 0.004 )
    {
        return  2.42960000e-02   * (-10);
    }
    else if ( newt <= 0.452 + 0.004 )
    {
        return  3.80990000e-02   * (-10);
    }
    else if ( newt <= 0.456 + 0.004 )
    {
        return  4.26980000e-02   * (-10);
    }
    else if ( newt <= 0.46  + 0.004 )
    {
        return  4.26980000e-02   * (-10);
    }
    else if ( newt <= 0.464 + 0.004 )
    {
        return  4.50000000e-02   * (-10);
    }
    else if ( newt <= 0.468 + 0.004 )
    {
        return  4.57670000e-02   * (-10);
    }
    else if ( newt <= 0.472 + 0.004 )
    {
        return  4.73020000e-02   * (-10);
    }
    else if ( newt <= 0.476 + 0.004 )
    {
        return  4.88320000e-02   * (-10);
    }
    else if ( newt <= 0.48  + 0.004 )
    {
        return  5.26690000e-02   * (-10);
    }
    else if ( newt <= 0.484 + 0.004 )
    {
        return  5.49660000e-02   * (-10);
    }
    else if ( newt <= 0.488 + 0.004 )
    {
        return  5.49660000e-02   * (-10);
    }
    else if ( newt <= 0.492 + 0.004 )
    {
        return  5.72680000e-02   * (-10);
    }
    else if ( newt <= 0.496 + 0.004 )
    {
        return  5.95700000e-02   * (-10);
    }
    else if ( newt <= 0.5   + 0.004 )
    {
        return  6.11000000e-02   * (-10);
    }
    else if ( newt <= 0.504 + 0.004 )
    {
        return  6.18670000e-02   * (-10);
    }
    else if ( newt <= 0.508 + 0.004 )
    {
        return  6.49370000e-02   * (-10);
    }
    else if ( newt <= 0.512 + 0.004 )
    {
        return  6.87680000e-02   * (-10);
    }
    else if ( newt <= 0.516 + 0.004 )
    {
        return  6.64710000e-02   * (-10);
    }
    else if ( newt <= 0.52  + 0.004 )
    {
        return  6.41690000e-02   * (-10);
    }
    else if ( newt <= 0.524 + 0.004 )
    {
        return  6.57040000e-02   * (-10);
    }
    else if ( newt <= 0.528 + 0.004 )
    {
        return  6.87680000e-02   * (-10);
    }
    else if ( newt <= 0.532 + 0.004 )
    {
        return  7.18380000e-02   * (-10);
    }
    else if ( newt <= 0.536 + 0.004 )
    {
        return  7.64370000e-02   * (-10);
    }
    else if ( newt <= 0.54  + 0.004 )
    {
        return  7.64370000e-02   * (-10);
    }
    else if ( newt <= 0.544 + 0.004 )
    {
        return  7.41350000e-02   * (-10);
    }
    else if ( newt <= 0.548 + 0.004 )
    {
        return  7.49020000e-02   * (-10);
    }
    else if ( newt <= 0.552 + 0.004 )
    {
        return  7.56700000e-02   * (-10);
    }
    else if ( newt <= 0.556 + 0.004 )
    {
        return  7.33680000e-02   * (-10);
    }
    else if ( newt <= 0.56  + 0.004 )
    {
        return  6.87680000e-02   * (-10);
    }
    else if ( newt <= 0.564 + 0.004 )
    {
        return  7.03030000e-02   * (-10);
    }
    else if ( newt <= 0.568 + 0.004 )
    {
        return  6.72340000e-02   * (-10);
    }
    else if ( newt <= 0.572 + 0.004 )
    {
        return  6.41690000e-02   * (-10);
    }
    else if ( newt <= 0.576 + 0.004 )
    {
        return  6.11000000e-02   * (-10);
    }
    else if ( newt <= 0.58  + 0.004 )
    {
        return  6.11000000e-02   * (-10);
    }
    else if ( newt <= 0.584 + 0.004 )
    {
        return  6.41690000e-02   * (-10);
    }
    else if ( newt <= 0.588 + 0.004 )
    {
        return  7.18380000e-02   * (-10);
    }
    else if ( newt <= 0.592 + 0.004 )
    {
        return  7.95020000e-02   * (-10);
    }
    else if ( newt <= 0.596 + 0.004 )
    {
        return  8.56360000e-02   * (-10);
    }
    else if ( newt <= 0.6   + 0.004 )
    {
        return  8.25710000e-02   * (-10);
    }
    else if ( newt <= 0.604 + 0.004 )
    {
        return  7.64370000e-02   * (-10);
    }
    else if ( newt <= 0.608 + 0.004 )
    {
        return  7.10710000e-02   * (-10);
    }
    else if ( newt <= 0.612 + 0.004 )
    {
        return  7.33680000e-02   * (-10);
    }
    else if ( newt <= 0.616 + 0.004 )
    {
        return  7.18380000e-02   * (-10);
    }
    else if ( newt <= 0.62  + 0.004 )
    {
        return  7.10710000e-02   * (-10);
    }
    else if ( newt <= 0.624 + 0.004 )
    {
        return  6.80010000e-02   * (-10);
    }
    else if ( newt <= 0.628 + 0.004 )
    {
        return  6.95360000e-02   * (-10);
    }
    else if ( newt <= 0.632 + 0.004 )
    {
        return  7.03030000e-02   * (-10);
    }
    else if ( newt <= 0.636 + 0.004 )
    {
        return  6.72340000e-02   * (-10);
    }
    else if ( newt <= 0.64  + 0.004 )
    {
        return  6.18670000e-02   * (-10);
    }
    else if ( newt <= 0.644 + 0.004 )
    {
        return  5.65010000e-02   * (-10);
    }
    else if ( newt <= 0.648 + 0.004 )
    {
        return  5.65010000e-02   * (-10);
    }
    else if ( newt <= 0.652 + 0.004 )
    {
        return  5.49660000e-02   * (-10);
    }
    else if ( newt <= 0.656 + 0.004 )
    {
        return  5.65010000e-02   * (-10);
    }
    else if ( newt <= 0.66  + 0.004 )
    {
        return  5.88030000e-02   * (-10);
    }
    else if ( newt <= 0.664 + 0.004 )
    {
        return  5.95700000e-02   * (-10);
    }
    else if ( newt <= 0.668 + 0.004 )
    {
        return  6.03370000e-02   * (-10);
    }
    else if ( newt <= 0.672 + 0.004 )
    {
        return  6.03370000e-02   * (-10);
    }
    else if ( newt <= 0.676 + 0.004 )
    {
        return  5.95700000e-02   * (-10);
    }
    else if ( newt <= 0.68  + 0.004 )
    {
        return  6.11000000e-02   * (-10);
    }
    else if ( newt <= 0.684 + 0.004 )
    {
        return  6.11000000e-02   * (-10);
    }
    else if ( newt <= 0.688 + 0.004 )
    {
        return  5.95700000e-02   * (-10);
    }
    else if ( newt <= 0.692 + 0.004 )
    {
        return  5.88030000e-02   * (-10);
    }
    else if ( newt <= 0.696 + 0.004 )
    {
        return  5.80350000e-02   * (-10);
    }
    else if ( newt <= 0.7   + 0.004 )
    {
        return  6.03370000e-02   * (-10);
    }
    else if ( newt <= 0.704 + 0.004 )
    {
        return  5.95700000e-02   * (-10);
    }
    else if ( newt <= 0.708 + 0.004 )
    {
        return  5.80350000e-02   * (-10);
    }
    else if ( newt <= 0.712 + 0.004 )
    {
        return  5.42030000e-02   * (-10);
    }
    else if ( newt <= 0.716 + 0.004 )
    {
        return  5.19010000e-02   * (-10);
    }
    else if ( newt <= 0.72  + 0.004 )
    {
        return  4.65350000e-02   * (-10);
    }
    else if ( newt <= 0.724 + 0.004 )
    {
        return  4.50000000e-02   * (-10);
    }
    else if ( newt <= 0.728 + 0.004 )
    {
        return  4.19350000e-02   * (-10);
    }
    else if ( newt <= 0.732 + 0.004 )
    {
        return  4.26980000e-02   * (-10);
    }
    else if ( newt <= 0.736 + 0.004 )
    {
        return  4.42330000e-02   * (-10);
    }
    else if ( newt <= 0.74  + 0.004 )
    {
        return  4.80690000e-02   * (-10);
    }
    else if ( newt <= 0.744 + 0.004 )
    {
        return  4.42330000e-02   * (-10);
    }
    else if ( newt <= 0.748 + 0.004 )
    {
        return  4.11680000e-02   * (-10);
    }
    else if ( newt <= 0.752 + 0.004 )
    {
        return  3.65640000e-02   * (-10);
    }
    else if ( newt <= 0.756 + 0.004 )
    {
        return  3.58010000e-02   * (-10);
    }
    else if ( newt <= 0.76  + 0.004 )
    {
        return  3.42670000e-02   * (-10);
    }
    else if ( newt <= 0.764 + 0.004 )
    {
        return  3.27320000e-02   * (-10);
    }
    else if ( newt <= 0.768 + 0.004 )
    {
        return  3.11970000e-02   * (-10);
    }
    else if ( newt <= 0.772 + 0.004 )
    {
        return  3.04300000e-02   * (-10);
    }
    else if ( newt <= 0.776 + 0.004 )
    {
        return  2.96670000e-02   * (-10);
    }
    else if ( newt <= 0.78  + 0.004 )
    {
        return  3.04300000e-02   * (-10);
    }
    else if ( newt <= 0.784 + 0.004 )
    {
        return  2.96670000e-02   * (-10);
    }
    else if ( newt <= 0.788 + 0.004 )
    {
        return  2.50630000e-02   * (-10);
    }
    else if ( newt <= 0.792 + 0.004 )
    {
        return  2.27660000e-02   * (-10);
    }
    else if ( newt <= 0.796 + 0.004 )
    {
        return  2.12310000e-02   * (-10);
    }
    else if ( newt <= 0.8   + 0.004 )
    {
        return  2.23800000e-02   * (-10);
    }
}


Real linearPontdist ( Real  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if (t > 0.8)
    {
        t = ( ( (int) floor (t * 1000.0) ) % 800) / 1000.0;
    }

    Real dt = ( (int) floor (t * 1000.0) ) % 4 / 1000.0;
    Real ti = floor (t * 1000.0) / 1000.0 - dt;
    Real tii = ti + 0.004;
    return (pont_dist (tii) - pont_dist (ti) ) / (0.004) * (t - (ti) ) + pont_dist (ti);
}


}
