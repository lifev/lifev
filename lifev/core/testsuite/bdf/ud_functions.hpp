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

    @author
    @date 00-00-0000
 */

namespace LifeV
{

// ===================================================
//! User functions
// ===================================================

class SourceFct_2d
{
public:
    inline Real operator() (Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/, int /*ic*/) const
    {
        return -4.;
    }
};

class AnalyticalSol_2d
{
    // ic stands for a component index (unuseful in the scalar case)
public:
    static Real u (Real t, Real x, Real y, Real /*z*/, UInt /*ic*/)
    {
        return t * t * (x * x + y * y);
    }

    static Real der_t (Real t, Real x, Real y, Real /*z*/, UInt /*ic*/)
    {
        return 2 * t * (x * x + y * y);
    }

    static Real grad ( UInt icoor, Real t, Real x, Real y, Real /*z*/, UInt /*ic*/)
    {
        switch (icoor)
        {
            case 1:
                return 2 * x * t * t;
            case 2:
                return 2 * y * t * t;
            default:
                return 0;
        }
    }
};

class SourceFct
{
public:
    inline Real operator() (Real /*x*/, Real /*y*/, Real /*z*/, Real /*t*/, int /*ic*/) const
    {
        return -6.;
    }
};


class AnalyticalSol
{
    // ic stands for a component index (unuseful in the scalar case)
public:
    static Real u (Real t, Real x, Real y, Real z, UInt /*ic*/)
    {
        return t * t * (x * x + y * y + z * z);
    }
    static Real grad (UInt icoor, Real t, Real x, Real y, Real z, UInt /*ic*/)
    {
        switch (icoor)
        {
            case 1: //der_x
                return 2 * x * t * t;
            case 2: //der_y
                return 2 * y * t * t;
            case 3: //der_z
                return 2 * z * t * t;
            default:
                return 0;
        }
    }
};


Real nu (const Real& t)
{
    return 1. / (t * t);
}

Real sigma (const Real& t)
{
    return -2. / t;
}

}
