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


#ifndef _FP_HPP
#define _FP_HPP

#include "operFS.hpp"

namespace LifeV
{

class fixedPoint : public operFS
{
public:

    typedef operFS super;
    typedef super::fluid_type fluid_type;
    typedef super::solid_type solid_type;
    typedef super::bchandler_type bchandler_type;
    // default constructor
//     fixedPoint()
//         :
//         super()
//         {}

    // constructors
    fixedPoint();

    fixedPoint( fluid_type& fluid,
                solid_type& solid,
                GetPot    &_dataFile,
                bchandler_type &BCh_u,
                bchandler_type &BCh_d,
                bchandler_type &BCh_mesh);

    // destructor
    ~fixedPoint();

    // member functions

    void evalResidual(Vector       &_res,
                      const Vector &_disp,
                      const int     _iter);
    void solveJac    (Vector       &_muk,
                      const Vector &_res,
                      const double  _linearRelTol);

    void setUpBC     ();

    Real   defOmega() {return M_defOmega;}

    void setDataFromGetPot( GetPot const& data );

    void setup();

private:

//     void eval        (const Vector &_disp,
//                       const int     _status);
    Real                    M_defOmega;
    generalizedAitken<Vector, Real> M_aitkFS;

    void eval(Vector& dispNew, Vector& velo, const Vector& disp, int status);

    void transferOnInterface(const Vector      &_vec1,
                             const BCHandler   &_BC,
                             const std::string &_BCName,
                             Vector            &_vec2);
};

}

#endif
