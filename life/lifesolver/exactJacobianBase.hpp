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




class exactJacobian : public operFS
{
public:

    typedef operFS super;
    typedef super::fluid_type fluid_type;
    typedef super::solid_type solid_type;
    typedef super::bchandler_type bchandler_type;

    // default constructor
    exactJacobian()
        :
        super(),
        M_BCh_du ( new BCHandler( 0, BCHandler::HINT_BC_ONLY_ESSENTIAL) ),
        M_BCh_dz ( new BCHandler ),
        M_dz     (),
        M_rhs_dz (),
        M_dataJacobian(this)
        {}

    // constructors
    exactJacobian( fluid_type& fluid,
                   solid_type& solid,
                   GetPot &_dataFile,
                   bchandler_type &BCh_u,
                   bchandler_type &BCh_d,
                   bchandler_type &BCh_mesh);

    // destructor
    ~exactJacobian();

    // member functions

    void evalResidual(Vector       &_res,
                      const Vector &_disp,
                      const int     _iter);
    void solveJac    (Vector       &_muk,
                      const Vector &_res,
                      const double  _linearRelTol);

    void solveLinearFluid();
    void solveLinearSolid();

    void setUpBC();

    void setup();

    Vector dz() {return M_dz;}

    struct dataJacobian
    {

        dataJacobian()
            :
            M_pFS( 0 )
            {}

        dataJacobian(exactJacobian* oper)
            :
            M_pFS(oper){}

        exactJacobian* M_pFS;

    };

private:

    bchandler_type         M_BCh_du;
    bchandler_type         M_BCh_dz;

    Vector            M_dz;
    Vector            M_rhs_dz;

    void eval        (const Vector &_disp,
                      const int     _status);

    dataJacobian            M_dataJacobian;
};


void my_matvecJacobianEJ(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);

Real fzeroEJ(const Real& t,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID& i);


}

#endif
