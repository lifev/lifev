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
#include "generalizedAitken.hpp"
#include "bcHandler.hpp"
#include "dof.hpp"


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
               BCHandler& BCh_du, BCHandler& BCh_dz);
//               GetPot &data_file);

        // destructor

        ~operFS();

        // member functions

        void eval           (const  Vector &disp,
                             int    status,
                             Vector &dispNew,
                             Vector &veloStruct);

        void evalResidual (Vector &res,
                           Vector &sol,
                           int    iter);

        void updateJac      (Vector& sol,
                             int     iter);

        void solveJac      (Vector &,
                            const Vector &,
                            double);

        void solveLinearFluid();

        void solveLinearSolid();

        void computeResidualFSI();
//        void computeResidualFSI(const PhysVectUnknown<Vector> &_res);

        // mutators and setters

        UInt   const & nbEval()      const
            {return M_nbEval;}

        Vector const & dz()          const
            {return M_dz;}

        Vector const & residualFSI() const
            {return M_residualFSI;}

        void setResidualFSI(double *_res);
        void setResidualFSI(const Vector _res);

        Vector getResidualFSIOnSolid();
        Vector getFluidInterfaceOnSolid(Vector &_vec);
        Vector getFluidInterfaceOnSolid(Vector &_vec,
                                        BCHandler &_BCh);

//         PhysVectUnknown<Vector> const & residualS() const
//             {return M_residualS;}
        PhysVectUnknown<Vector> & residualS()             {return M_residualS;}
        PhysVectUnknown<Vector> const & residualF() const {return M_residualF;}
        PhysVectUnknown<Vector> & residualFSI()           {return M_residualFSI;}

        NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >
        &fluid() {return M_fluid;}

        VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >
        &solid() {return M_solid;}

        void setTime(const Real &time) {M_time = time;};

    private:

        NavierStokesAleSolverPC
        < RegionMesh3D_ALE<LinearTetra> > &M_fluid;

        VenantKirchhofSolver
        < RegionMesh3D_ALE<LinearTetra> > &M_solid;

        Vector                  M_dispStruct;

        Real                    M_time;

        SolverAztec             M_solverAztec;

        Vector                  M_velo;
        Vector                  M_dz;
        Vector                  M_rhs_dz;

        PhysVectUnknown<Vector> M_residualS;
        PhysVectUnknown<Vector> M_residualF;
        PhysVectUnknown<Vector> M_residualFSI;

        UInt                    M_nbEval;

        BCHandler&             M_BCh_du;
        BCHandler&             M_BCh_dz;

        DataJacobian            M_dataJacobian;

        Vector setDispOnInterface(const Vector &_disp);

        void  invSfPrime  (const Vector &res,
                           double       linear_rel_tol,
                           Vector       &step);

        void  invSsPrime  (const Vector &res,
                           double       linear_rel_tol,
                           Vector       &step);

        void  invSfSsPrime(const Vector &res,
                           double       linear_rel_tol,
                           Vector       &step);
    };




void my_matvecSfSsPrime(double *z,
                        double *Jz,
                        AZ_MATRIX* J,
                        int proc_config[]);
}
#endif
