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

#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "vectorNorms.hpp"
#include "regionMesh3D_ALE.hpp"
#include "SolverAztec.hpp"

#ifndef _OPERFS
#define _OPERFS

namespace LifeV
{
    class operFS;
//    class GetPot;

    class DataJacobian
    {
    public:
        
        DataJacobian(operFS* oper):
            M_pFS(oper){}
        
        operFS* M_pFS;
    };
    

//
// Fluid-Structure operator Class
//
    class operFS {
        
    public:

        // constructors
        
        operFS(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
               VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid,
               BC_Handler& BCh_du, BC_Handler& BCh_dz,
               GetPot &data_file);


        // member functions
        
        void eval         (Vector &dispNew,
                           Vector &veloStruct,
                           const  Vector &disp,
                           int    status);
        
        void evalResidual (Vector &res,
                           const  Vector &sol,
                           int    iter);

        void updatePrec   (Vector& sol,
                           int     iter);

        void solvePrec    (Vector &);

        void solveLinearFluid();

        void solveLinearSolid();

        // mutators and setters

        UInt   const & nbEval() const {return M_nbEval;};

        Vector const & dz()     const {return M_dz;};
        
        NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >
        &fluid() {return M_fluid;};
        
        VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >
        &solid() {return M_solid;};

        void setTime(const Real &time) {M_time = time;};
        
    private:

        NavierStokesAleSolverPC
        < RegionMesh3D_ALE<LinearTetra> > &M_fluid;

        VenantKirchhofSolver
        < RegionMesh3D_ALE<LinearTetra> > &M_solid;

        Vector       M_dispStruct;
        
        Real         M_time;

        SolverAztec  M_solverAztec;
        
        Vector       M_velo;
        Vector       M_dz;
        Vector       M_rhs_dz;

        UInt         M_nbEval;

        BC_Handler&  M_BCh_du;
        BC_Handler&  M_BCh_dz;

        DataJacobian M_dataJacobian;
    };

    void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);
}
#endif
