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
    NO_PRECONDITIONER=-1,
    NEUMANN_DIRICHLET,
    DIRICHLET_NEUMANN,
    NEUMANN_NEUMANN,
    NEWTON
};
class operFS {

public:

    typedef boost::shared_ptr<NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > > fluid_type;
    typedef boost::shared_ptr<VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > > solid_type;

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID
                                   & )> function_type;

    typedef boost::shared_ptr<DofInterface3Dto3D> dof_interface_type;
    typedef boost::shared_ptr<BCHandler> bchandler_type;

    // constructors
    operFS():
        M_BCh_u(),
        M_BCh_d(),
        M_BCh_mesh(),
        M_fluid(),
        M_solid(),
        M_dofFluidToStructure( new DofInterface3Dto3D ),
        M_dofStructureToSolid( new DofInterface3Dto3D ),
        M_dofStructureToFluidMesh( new DofInterface3Dto3D ),
        M_dofMeshToFluid( new DofInterface3Dto3D ),
        M_dispStruct(),
        M_dispStructOld(),
        M_velo(),
        M_nbEval( 0 ),
        M_method(),
        M_precond( NO_PRECONDITIONER )
        {}

    operFS(fluid_type & fluid,
           solid_type &solid,
           GetPot    &data_file,
           bchandler_type& BCh_u,
           bchandler_type& BCh_d,
           bchandler_type& BCh_mesh);

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

    void setPreconditioner( OperFSPreconditioner _p ) { M_precond = _p; }
    OperFSPreconditioner preconditioner() const { return M_precond; }

    void setTime(const Real &time) {M_time = time;};
    Real time() {return M_time;};

    Vector displacementOnInterface();

    void setFluid( fluid_type const& fluid ){ M_fluid = fluid;}
    void setSolid( solid_type const& solid ){ M_solid = solid;}

    virtual void setDataFromGetPot( GetPot const& data );

    void setBC( bchandler_type& bc_u,  bchandler_type& bc_d, bchandler_type& bc_m )
        {
            M_BCh_u = bc_u;
            M_BCh_d = bc_d;
            M_BCh_mesh = bc_m;
        }

    virtual void setup();

    dof_interface_type& dofMeshToFluid() { return M_dofMeshToFluid; }
    dof_interface_type const& dofMeshToFluid() const { return M_dofMeshToFluid; }

protected:

    void transferOnInterface(const Vector      &_vec1,
                             const BCHandler   &_BC,
                             const std::string &_BCName,
                             Vector            &_vec2);

    bchandler_type          M_BCh_u;
    bchandler_type          M_BCh_d;
    bchandler_type          M_BCh_mesh;

    fluid_type              M_fluid;
    solid_type              M_solid;

    dof_interface_type      M_dofFluidToStructure;
    dof_interface_type      M_dofStructureToSolid;
    dof_interface_type      M_dofStructureToFluidMesh;
    dof_interface_type      M_dofMeshToFluid;

    Vector                  M_dispStruct;
    Vector                  M_dispStructOld;
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

/*!
   \def  FOR_EACH_INTERFACE_DOF( Expr )

   \c FOR_EACH_INTERFACE_DOF is an helper macro to ease the
   the computation of quantities like the residual on the interface
   of the Fluid and Structure.

 */
#define FOR_EACH_INTERFACE_DOF( Expr )                              \
{                                                                   \
    UInt iBCf = M_fluid->BC_fluid().getBCbyName("Interface");       \
                                                                    \
    BCBase const &BC_fluidInterface = M_fluid->BC_fluid()[iBCf];    \
                                                                    \
    UInt nDofInterface = BC_fluidInterface.list_size();             \
                                                                    \
    UInt nDimF = BC_fluidInterface.numberOfComponents();            \
                                                                    \
    UInt totalDofFluid = M_fluid->uDof().numTotalDof();             \
    UInt totalDofSolid = M_solid->dDof().numTotalDof();             \
                                                                    \
    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)                 \
    {                                                               \
        ID IDfluid = BC_fluidInterface(iBC)->id();                  \
                                                                    \
        BCVectorInterface const *BCVInterface =                     \
            dynamic_cast <BCVectorInterface const *>                \
            (BC_fluidInterface.pointerToBCVector());                \
                                                                    \
        assert( BCVInterface != 0 );                                \
                                                                    \
        ID IDsolid = BCVInterface->                                 \
            dofInterface().getInterfaceDof(IDfluid);                \
                                                                    \
        for (UInt jDim = 0; jDim < nDimF; ++jDim)                   \
        {                                                           \
            ( Expr );                                               \
        }                                                           \
    }                                                               \
}

}
#endif
