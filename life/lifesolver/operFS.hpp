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

#include <factory.hpp>
#include <singleton.hpp>

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
typedef enum OperFSPreconditioner
{
    NEUMANN_DIRICHLET = 0,
    DIRICHLET_NEUMANN,
    NEUMANN_NEUMANN
};
class operFS {

public:

    typedef boost::shared_ptr<NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > > fluid_type;
    typedef boost::shared_ptr<VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > > solid_type;

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID
                                   & )> function_type;
    // constructors
    operFS()
        :
        M_BCh_u(),
        M_BCh_d(),
        M_BCh_mesh(),
        M_fluid(),
        M_solid()
        {}

    operFS(fluid_type & fluid,
           solid_type &solid,
           GetPot    &data_file,
           BCHandler &BCh_u,
           BCHandler &BCh_d,
           BCHandler &BCh_mesh);

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

    virtual void setUpBC() = 0;

    void updateJac (Vector& sol,
                    int     iter);

    void solveLinearFluid();

    void solveLinearSolid();

    // mutators and setters

    UInt   const & nbEval()      const
        {return M_nbEval;}

    fluid_type::value_type& fluid() {return *M_fluid;}

    solid_type::value_type& solid() {return *M_solid;}

    void setPreconditioner( OperFSPreconditioner __p ) { M_precond = __p; }
    OperFSPreconditioner preconditioner() const { return M_precond; }

    void setTime(const Real &time) {M_time = time;};
    Real time() {return M_time;};

    void displacementOnInterface();

    void setFluid( fluid_type const& fluid ){ M_fluid = fluid;}
    void setSolid( solid_type const& solid ){ M_solid = solid;}

    virtual void setDataFromGetPot( GetPot const& data );

    void setBC( BCHandler const& bc_u,  BCHandler const& bc_d, BCHandler const& bc_m )
        {
            M_BCh_u = bc_u;
            M_BCh_d = bc_d;
            M_BCh_mesh = bc_m;
        }

    virtual void setup();

protected:

    BCHandler               M_BCh_u;
    BCHandler               M_BCh_d;
    BCHandler               M_BCh_mesh;

    fluid_type              M_fluid;
    solid_type              M_solid;

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
    OperFSPreconditioner    M_precond;
};
typedef boost::shared_ptr<operFS> oper_fsi_ptr;
typedef singleton<factory<operFS,  std::string> > FSIFactory;
}
#endif
