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

#include "dofInterface3Dto3D.hpp"
#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "regionMesh3D_ALE.hpp"
#include "SolverAztec.hpp"
#include "generalizedAitken.hpp"
#include "bcHandler.hpp"
#include "bcFunction.hpp"
#include "dof.hpp"

#ifndef _OPERFS
#define _OPERFS

namespace LifeV
{
//
// Fluid-Structure operator Class
//
    class operFS {

    public:
        typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID
                                       & )> function_type;
        // constructors

        operFS(GetPot &data_file);

        // destructor

        virtual ~operFS();

        // virtual memeber functions

        virtual void evalResidual(Vector &res,
                                  const Vector &_disp,
                                  const int     _iter) = 0;

        virtual void solveJac (Vector &_muk,
                               const Vector &_res,
                                const double  _linearRelTol) = 0;

        // member functions

        virtual void setUpBC(function_type _bcf,
                             function_type _vel) = 0;
        void updateJac (Vector& sol,
                        int     iter);

        void solveLinearFluid();

        void solveLinearSolid();

        // mutators and setters

        UInt   const & nbEval()      const
            {return M_nbEval;}

        NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >
        &fluid() {return M_fluid;}

        VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >
        &solid() {return M_solid;}

        void setTime(const Real &time) {M_time = time;};
        Real time() {return M_time;};

    protected:

        BCHandler               M_BCh_u;
        BCHandler               M_BCh_d;
        BCHandler               M_BCh_mesh;

        BCHandler               M_BCh_du;
        BCHandler               M_BCh_dz;

        NavierStokesAleSolverPC
        < RegionMesh3D_ALE<LinearTetra> > M_fluid;

        VenantKirchhofSolver
        < RegionMesh3D_ALE<LinearTetra> > M_solid;

        DofInterface3Dto3D      M_dofFluidToStructure;
        DofInterface3Dto3D      M_dofStructureToSolid;
        DofInterface3Dto3D      M_dofStructureToFluidMesh;
        DofInterface3Dto3D      M_dofMeshToFluid;

        Vector                  M_dispStruct;
        Vector                  M_velo;

        SolverAztec             M_solverAztec;

        Real                    M_time;

        UInt                    M_nbEval;
    private:


        UInt                    M_method;
        UInt                    M_precond;
    };

}
#endif
