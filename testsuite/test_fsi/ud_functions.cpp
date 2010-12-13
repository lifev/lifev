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
    @file
    @brief

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @date 00-00-0000
 */

#include "ud_functions.hpp"


namespace LifeV
{
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





Real u2(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
    case 1:
        return 0.0;
        break;
    case 3:
        if ( t <= 0.003 )
            return 1.3332e4;
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

Real PhysFlux(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    double coef = 0.001;
    double plusmoins = -1.;

    if (t < 0.01) return -1000.*t;
    if (t <= 0.02) return -10.;
    if (t <= 0.03) return 1000.*t - 30.;
    return 0;

    int    numData = 100;
    double flux[101] = { 0.,
                         0.55312181720914,
                         0.55299302643153,
                         0.55302818124406,
                         0.55317321131557,
                         0.55353364652152,
                         0.55374878962962,
                         0.55406829977313,
                         0.55585881887584,
                         0.55879983633299,
                         0.56387572718194,
                         0.57079488161935,
                         0.58817359652018,
                         0.61511673048210,
                         0.65025652077432,
                         0.79143597227331,
                         1.06959564837144,
                         1.33301653745648,
                         1.55916094914244,
                         1.69807981838757,
                         1.73941326221337,
                         1.69691789994768,
                         1.61546505715344,
                         1.51298277554169,
                         1.35636910481872,
                         1.18958468647846,
                         1.02381943290960,
                         0.89472539111700,
                         0.79994401665900,
                         0.74338513206540,
                         0.72984443940672,
                         0.74311502021666,
                         0.77155202899612,
                         0.80920592343822,
                         0.84528951107944,
                         0.87014986946118,
                         0.89016029509068,
                         0.89999452080415,
                         0.89511807514404,
                         0.87823357219620,
                         0.85663326250089,
                         0.82858260857792,
                         0.79671836139948,
                         0.75291106131325,
                         0.70711033067577,
                         0.65740190152584,
                         0.61066125251620,
                         0.56488426267136,
                         0.51713402331108,
                         0.46453623816504,
                         0.41513731950517,
                         0.47836113912116,
                         0.56452765299777,
                         0.62096051336166,
                         0.66202024502726,
                         0.69173157064612,
                         0.71835003021294,
                         0.74183126309604,
                         0.75295645424862,
                         0.75292455314576,
                         0.74514787317446,
                         0.72467414023271,
                         0.70473486985061,
                         0.68057326827129,
                         0.66194232245132,
                         0.64681425465222,
                         0.63714376881254,
                         0.62991615896879,
                         0.62662699778909,
                         0.62724985397200,
                         0.62674770176751,
                         0.62666043242736,
                         0.62617509524360,
                         0.62556258658310,
                         0.62581913341632,
                         0.62604032520998,
                         0.62582937168093,
                         0.62404471163034,
                         0.61923663804136,
                         0.61378537728592,
                         0.60976137345625,
                         0.60596975158344,
                         0.60144172708524,
                         0.59702451106965,
                         0.59319136754468,
                         0.58982107329344,
                         0.58718879911670,
                         0.58474181066352,
                         0.58208280034445,
                         0.57913818429409,
                         0.57588144074776,
                         0.57289019558638,
                         0.57076133371909,
                         0.56912637026578,
                         0.56776096894206,
                         0.56622327633393,
                         0.56376396210446,
                         0.56132345661888,
                         0.55914786876437,
                         0.55714215894326,
                         0.55457447187998
                       };


    double timescale = 1./numData;
    int     ipos = t/timescale;

    double t2 =timescale*(ipos + 1);

    double a = (flux[ipos + 1] - flux[ipos])/timescale;
    double b = flux[ipos + 1] - a*t2;

    double slope = ipos -  t;

    std::cout << "t = " << t << " f = " << t*a + b << " " << ipos << " " << a << " " << b << std::endl;


    return plusmoins*(t*a + b);
}

Real aortaPhysPress(const Real&  t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    double coef = .1;
    switch (i)
    {
    case 1:
        return 0.0;
        break;
    case 3:
        if (t<=0.00) return 110170*coef;
        if (t<=0.01) return 109540*coef;
        if (t<=0.02) return 108930*coef;
        if (t<=0.03) return 108320*coef;
        if (t<=0.04) return 107710*coef;
        if (t<=0.05) return 107120*coef;
        if (t<=0.06) return 106530*coef;
        if (t<=0.07) return 111130*coef;
        if (t<=0.08) return 115440*coef;
        if (t<=0.09) return 118690*coef;
        if (t<=0.10) return 121460*coef;
        if (t<=0.11) return 123940*coef;
        if (t<=0.12) return 126350*coef;
        if (t<=0.13) return 128890*coef;
        if (t<=0.14) return 131510*coef;
        if (t<=0.15) return 133980*coef;
        if (t<=0.16) return 136200*coef;
        if (t<=0.17) return 138330*coef;
        if (t<=0.18) return 140350*coef;
        if (t<=0.19) return 142290*coef;
        if (t<=0.20) return 144360*coef;
        if (t<=0.21) return 146130*coef;
        if (t<=0.22) return 147530*coef;
        if (t<=0.23) return 148780*coef;
        if (t<=0.24) return 149740*coef;
        if (t<=0.25) return 150320*coef;
        if (t<=0.26) return 150470*coef;
        if (t<=0.27) return 150250*coef;
        if (t<=0.28) return 149750*coef;
        if (t<=0.29) return 148990*coef;
        if (t<=0.30) return 148220*coef;
        if (t<=0.31) return 147210*coef;
        if (t<=0.32) return 145940*coef;
        if (t<=0.33) return 144960*coef;
        if (t<=0.34) return 143750*coef;
        if (t<=0.35) return 141980*coef;
        if (t<=0.36) return 139900*coef;
        if (t<=0.37) return 137260*coef;
        if (t<=0.38) return 133970*coef;
        if (t<=0.39) return 131670*coef;
        if (t<=0.40) return 131320*coef;
        if (t<=0.41) return 133150*coef;
        if (t<=0.42) return 132710*coef;
        if (t<=0.43) return 131570*coef;
        if (t<=0.44) return 130280*coef;
        if (t<=0.45) return 129750*coef;
        if (t<=0.46) return 129330*coef;
        if (t<=0.47) return 128910*coef;
        if (t<=0.48) return 128360*coef;
        if (t<=0.49) return 127680*coef;
        if (t<=0.50) return 127000*coef;
        if (t<=0.51) return 126410*coef;
        if (t<=0.52) return 125920*coef;
        if (t<=0.53) return 125480*coef;
        if (t<=0.54) return 125040*coef;
        if (t<=0.55) return 124560*coef;
        if (t<=0.56) return 124050*coef;
        if (t<=0.57) return 123530*coef;
        if (t<=0.58) return 123000*coef;
        if (t<=0.59) return 122440*coef;
        if (t<=0.60) return 121840*coef;
        if (t<=0.61) return 121220*coef;
        if (t<=0.62) return 120580*coef;
        if (t<=0.63) return 119950*coef;
        if (t<=0.64) return 119330*coef;
        if (t<=0.65) return 118710*coef;
        if (t<=0.66) return 118100*coef;
        if (t<=0.67) return 117470*coef;
        if (t<=0.68) return 116840*coef;
        if (t<=0.69) return 116200*coef;
        if (t<=0.70) return 115560*coef;
        if (t<=0.71) return 114920*coef;
        if (t<=0.72) return 114280*coef;
        if (t<=0.73) return 113650*coef;
        if (t<=0.74) return 113020*coef;
        if (t<=0.75) return 112400*coef;
        if (t<=0.76) return 111790*coef;
        if (t<=0.77) return 111200*coef;
        if (t<=0.78) return 110620*coef;
        if (t<=0.79) return 110060*coef;
        break;
    case 2:
        return 0.0;
        break;
    }

    return 0.;
}


}

