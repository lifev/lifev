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


#ifndef _EJ_HPP
#define _EJ_HPP

#include "operFS.hpp"

namespace LifeV
{

class exactJacobian;

class dataJacobianEJ {

 public:

  dataJacobianEJ(exactJacobian* oper):
    M_pFS(oper){}

    exactJacobian* M_pFS;

};


class exactJacobian : public operFS
{
public:

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID
                                   & )> function_type;

    // constructors
    exactJacobian(GetPot &_dataFile);

    // destructor
    ~exactJacobian();

    // member functions

    void evalResidual(const Vector &_disp,
                      const int     _iter,
                      Vector       &_res);
    void solveJac    (const Vector &_res,
                      const double  _linearRelTol,
                      Vector       &_muk);
    void solveLinearFluid();
    void solveLinearSolid();

    void setUpBC(function_type _bcf,
                 function_type _vel);

    Vector dz() {return M_dz;}

private:

    Vector            M_dz;
    Vector            M_rhs_dz;

    void eval        (const Vector &_disp,
                      const int     _status);

    dataJacobianEJ            M_dataJacobian;
};


void my_matvecJacobianEJ(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);

}

#endif
