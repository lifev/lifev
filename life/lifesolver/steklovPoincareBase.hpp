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


#ifndef _STEKLOV_HPP
#define _STEKLOV_HPP

#include "operFS.hpp"

namespace LifeV
{

class steklovPoincare : public operFS
{
public:

    typedef operFS super;
    typedef super::fluid_type fluid_type;
    typedef super::solid_type solid_type;
    typedef super::bchandler_type bchandler_type;

    // default constructor
    steklovPoincare();

    // constructors
    steklovPoincare(fluid_type& fluid,
                    solid_type& solid,
                    GetPot    &_dataFile,
                    bchandler_type &BCh_u,
                    bchandler_type &BCh_d,
                    bchandler_type &BCh_mesh);

    // destructor
    ~steklovPoincare();

    // member functions

    void evalResidual(Vector       &_res,
                      const Vector &_disp,
                      const int     _iter);

    void solveJac    (Vector        &_muk,
                      const Vector &_res,
                      const double  _linearRelTol);

    void solveLinearFluid();
    void solveLinearSolid();

    void solveInvLinearFluid();
    void solveInvLinearSolid();

    void setUpBC();

    //setters and getters

    Vector dz()       {return M_dz;}
    Real   defOmega() {return M_defOmega;}

    void setResidualS  ( Vector const& _res ){M_residualS = _res;}
    void setResidualF  ( Vector const& _res ){M_residualF = _res;}

//    void setResidualF  ( Vector const& _res );
    void setResidualFSI( double const* _res );
    void setResidualFSI( Vector const& _res );

    void setDataFromGetPot( GetPot const& data );
    Vector getResidualFSIOnSolid();
    void getResidualFSIOnSolid(Vector& _vec);
    Vector getSolidInterfaceOnSolid(double const* _vec);

    void setup();

    struct DataJacobian
    {
        DataJacobian()
            :
            M_pFS(0)
            {}
        DataJacobian(steklovPoincare* oper)
            :
            M_pFS(oper)
            {}

        steklovPoincare* M_pFS;
    };

    bchandler_type   BCh_du_inv(){return M_BCh_du_inv;}
    bchandler_type   BCh_dz_inv(){return M_BCh_dz_inv;}

    bchandler_type   BCh_du(){return M_BCh_du;}
    bchandler_type   BCh_dz(){return M_BCh_dz;}

private:

    bchandler_type          M_BCh_du;
    bchandler_type          M_BCh_dz;

    bchandler_type          M_BCh_du_inv;
    bchandler_type          M_BCh_dz_inv;

    Vector                  M_dz;
    Vector                  M_rhs_dz;

    PhysVectUnknown<Vector> M_residualS;
    PhysVectUnknown<Vector> M_residualF;
    PhysVectUnknown<Vector> M_residualFSI;

    Real                    M_defOmega;
    Real                    M_defOmegaS;
    Real                    M_defOmegaF;

    generalizedAitken<Vector, Real> M_aitkFS;


    Real                    M_linearRelTol;


    void eval           (const  Vector &disp,
                         int    status,
                         Vector &dispNew,
                         Vector &veloStruct);

    void  invSfPrime    (const Vector &res,
                         double       linear_rel_tol,
                         Vector       &step);

    void  invSsPrime    (const Vector &res,
                         double       linear_rel_tol,
                         Vector       &step);


    void  invSfSsPrime  (const Vector &res,
                         double       linear_rel_tol,
                         Vector       &step);

    void computeResidualFSI();


    Vector getSolidInterfaceOnFluid(Vector const& _vec);
    Vector getFluidInterfaceOnSolid(Vector const& _vec);

//     void transferOnInterface(const Vector      &_vec1,
//                              const BCHandler   &_BC,
//                              const std::string &_BCName,
//                              Vector            &_vec2);

    DataJacobian            M_dataJacobian;
};

Real fzeroSP(const Real& t,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID& i);

void my_matvecSfSsPrime(double *z,
                        double *Jz,
                        AZ_MATRIX* J,
                        int proc_config[]);

void my_matvecSfSsPrimePrec(double *z,
                            double *Jz,
                            AZ_MATRIX* J,
                            int proc_config[]);

}

#endif
