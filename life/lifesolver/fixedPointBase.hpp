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

    // constructors
    fixedPoint(GetPot    &_dataFile,
               BCHandler &BCh_u,
               BCHandler &BCh_d,
               BCHandler &BCh_mesh);

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

private:

//     void eval        (const Vector &_disp,
//                       const int     _status);

    void eval(Vector& dispNew, Vector& velo, const Vector& disp, int status);

    Real                    M_defOmega;

};

}

#endif
