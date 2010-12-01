/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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

// ===================================================
//! User functions
// ===================================================

#ifdef TWODIM
class SourceFct
{
public:
    inline Real operator()(Real /*x*/,Real /*y*/,Real /*z*/,Real /*t*/,int /*ic*/) const
    {
        return -4.;
    }
};

class AnalyticalSol
{
    // ic stands for a component index (unuseful in the scalar case)
public:
    static Real u(Real t, Real x,Real y,Real /*z*/, UInt /*ic*/)
    {
        return t*t*(x*x+y*y);
    }

    static Real der_t(Real t, Real x,Real y,Real /*z*/, UInt /*ic*/)
    {
        return 2*t*(x*x+y*y);
    }

    static Real grad( UInt icoor, Real t, Real x,Real y,Real /*z*/, UInt /*ic*/)
    {
        switch (icoor)
        {
        case 1:
            return 2*x*t*t;
        case 2:
            return 2*y*t*t;
        default:
            return 0;
        }
    }
};

#elif defined THREEDIM
class SourceFct
{
public:
    inline Real operator()(Real /*x*/,Real /*y*/,Real /*z*/,Real /*t*/,int /*ic*/) const
    {
        return -6.;
    }
};


class AnalyticalSol
{
    // ic stands for a component index (unuseful in the scalar case)
public:
    static Real u(Real t, Real x,Real y,Real z, UInt /*ic*/)
    {
        return t*t*(x*x+y*y+z*z);
    }
    static Real grad(UInt icoor, Real t, Real x,Real y,Real z, UInt /*ic*/)
    {
        switch (icoor)
        {
        case 1: //der_x
            return 2*x*t*t;
        case 2: //der_y
            return 2*y*t*t;
        case 3: //der_z
            return 2*z*t*t;
        default:
            return 0;
        }
    }
};

#endif

Real nu(const Real& t)
{
    return 1./(t*t);
}

Real sigma(const Real& t)
{
    return -2./t;
}

}
