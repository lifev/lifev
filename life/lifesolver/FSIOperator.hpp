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

//! previously called operFS in the CVS repository

#ifndef _OPERFS
#define _OPERFS

#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifefem/dofInterface3Dto3D.hpp>
#include <life/lifesolver/NavierStokesAleSolverPC.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
//#include <life/lifesolver/reducedLinFluid.hpp>
#include <life/lifefem/regionMesh3D_ALE.hpp>
#include <life/lifealg/SolverAztec.hpp>
#include <life/lifealg/generalizedAitken.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bcFunction.hpp>
#include <life/lifefem/dof.hpp>


namespace LifeV
{
//
// Fluid-Structure operator Class
//
typedef enum Preconditioner
{
    NO_PRECONDITIONER = -1,
    NEUMANN_DIRICHLET,
    DIRICHLET_NEUMANN,
    NEUMANN_NEUMANN,
    NEWTON
};

typedef enum DDNPreconditioner
{
    DDN_NO_PRECONDITIONER = -1,
    DDN_NEUMANN_DIRICHLET,
    DDN_DIRICHLET_NEUMANN,
    DDN_NEUMANN_NEUMANN
};

class reducedLinFluid;

class FSIOperator {

public:

    typedef boost::shared_ptr<NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > > fluid_type;
    typedef boost::shared_ptr<VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> > > solid_type;

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID
                                   & )> function_type;

    typedef boost::shared_ptr<DofInterface3Dto3D> dof_interface_type;
    typedef boost::shared_ptr<BCVectorInterface>  bc_vector_interface;
    typedef boost::shared_ptr<BCHandler>          bchandler_type;

//    typedef boost::shared_ptr<reducedLinFluid> quasi_newton_type;
    typedef boost::shared_ptr<reducedLinFluid>    quasi_newton_type;

    // constructors
    FSIOperator():
        M_BCh_u(new BCHandler),
        M_BCh_d(new BCHandler),
        M_BCh_mesh(new BCHandler),
        M_BCh_du(new BCHandler),
        M_BCh_du_inv(new BCHandler),
        M_BCh_dz(new BCHandler),
        M_BCh_dz_inv(new BCHandler),
        M_BCh_dp(new BCHandler),
        M_BCh_dp_inv(new BCHandler),
        M_fluid(),
        M_solid(),
        M_dofFluidToStructure                ( new DofInterface3Dto3D ),
//         M_dofFluidToSolid                    ( new DofInterface3Dto3D ),
//         M_dofSolidTofluid                    ( new DofInterface3Dto3D ),
        M_dofStructureToSolid                ( new DofInterface3Dto3D ),
        M_dofStructureToHarmonicExtension    ( new DofInterface3Dto3D ),
        M_dofHarmonicExtensionToFluid        ( new DofInterface3Dto3D ),
        M_dofStructureToReducedFluid         ( new DofInterface3Dto3D ),
        M_dofReducedFluidToStructure         ( new DofInterface3Dto3D ),
        // boundary vector interfaces
        M_bcvFluidLoadToStructure            ( new  BCVectorInterface ),
        M_bcvStructureDispToSolid            ( new  BCVectorInterface ),
        M_bcvStructureDispToHarmonicExtension( new  BCVectorInterface ),
        M_bcvHarmonicExtensionVelToFluid     ( new  BCVectorInterface ),
        M_bcvStructureToFluid                ( new  BCVectorInterface ),
        M_bcvStructureToReducedFluid         ( new  BCVectorInterface ),
        M_bcvReducedFluidToStructure         ( new  BCVectorInterface ),
        M_bcvDerHarmonicExtensionVelToFluid  ( new  BCVectorInterface ),
        M_bcvDerFluidLoadToStructure         ( new  BCVectorInterface ),
        M_bcvDerFluidLoadToFluid             ( new  BCVectorInterface ),
        M_bcvDerStructureDispToSolid         ( new  BCVectorInterface ),
        M_bcvDerReducedFluidLoadToStructure  ( new  BCVectorInterface ),
        M_bcvDerStructureAccToReducedFluid   ( new  BCVectorInterface ),
        M_dispStruct(),
        M_dispStructOld(),
        M_velo(),
        M_nbEval( 0 ),
        M_method(),
        M_precond( NO_PRECONDITIONER )
        {}

//     FSIOperator(fluid_type     &fluid,
//            solid_type     &solid,
//            GetPot         &data_file,
//            bchandler_type &BCh_u,
//            bchandler_type &BCh_d,
//            bchandler_type &BCh_mesh);

    // destructor

    virtual ~FSIOperator();

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

    void setPreconditioner   ( Preconditioner    _p ) { M_precond = _p; }
    void setDDNPreconditioner( DDNPreconditioner _p ) { M_DDNprecond = _p; }

    Preconditioner       preconditioner() const { return M_precond; }
    DDNPreconditioner    DDNpreconditioner() const { return M_DDNprecond; }

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

    // BC Vector Interface setters and getters

    void setHarmonicExtensionVelToFluid(PhysVectUnknown<Vector> &vel, UInt type = 0);
    bc_vector_interface bcvHarmonicExtensionVelToFluid()
        {return M_bcvHarmonicExtensionVelToFluid;}

    void setDerHarmonicExtensionVelToFluid(PhysVectUnknown<Vector> &dvel, UInt type = 0);
    bc_vector_interface bcvDerHarmonicExtensionVelToFluid()
        {return M_bcvDerHarmonicExtensionVelToFluid;}

    void setStructureDispToHarmonicExtension(PhysVectUnknown<Vector> &disp, UInt type = 0);
    bc_vector_interface bcvStructureDispToHarmonicExtension()
        {return M_bcvStructureDispToHarmonicExtension;}

    void setStructureDispToSolid(PhysVectUnknown<Vector> &disp, UInt type = 0);
    bc_vector_interface bcvStructureDispToSolid()
        {return M_bcvStructureDispToSolid;}

    void setDerStructureDispToSolid(PhysVectUnknown<Vector> &ddisp, UInt type = 0);
    bc_vector_interface bcvDerStructureDispToSolid()
        {return M_bcvDerStructureDispToSolid;}

    void setFluidLoadToStructure(Vector &load, UInt type = 0);
    bc_vector_interface bcvFluidLoadToStructure()
        {return M_bcvFluidLoadToStructure;}

    void setDerFluidLoadToStructure(Vector &dload, UInt type = 0);
    bc_vector_interface bcvDerFluidLoadToStructure()
        {return M_bcvDerFluidLoadToStructure;}

    void setDerFluidLoadToFluid(Vector &dload, UInt type = 0);
    bc_vector_interface bcvDerFluidLoadToFluid()
        {return M_bcvDerFluidLoadToFluid;}

    void setDerReducedFluidLoadToStructure(Vector &dload, UInt type = 0);
    bc_vector_interface bcvDerReducedFluidLoadToStructure()
        {return M_bcvDerReducedFluidLoadToStructure;}

    void setDerStructureAccToReducedFluid(Vector &acc, UInt type = 0);
    bc_vector_interface bcvDerStructureAccToReducedFluid()
        {return M_bcvDerStructureAccToReducedFluid;}

    //


    quasi_newton_type getReducedLinFluid(){return M_reducedLinFluid;}

    UInt reducedFluid(){return M_reducedFluid;}

    bchandler_type const& BCh_fluid(){return M_BCh_u;}

    bchandler_type const& BCh_solid(){return M_BCh_d;}

    bchandler_type const& BCh_harmonicExtension(){return M_BCh_mesh;}

    bchandler_type const& BCh_du(){return M_BCh_du;}
    bchandler_type const& BCh_du_inv(){return M_BCh_du_inv;}

    bchandler_type const& BCh_dz(){return M_BCh_dz;}
    bchandler_type const& BCh_dz_inv(){return M_BCh_dz_inv;}


protected:

    void transferOnInterface(const Vector      &_vec1,
                             const BCHandler   &_BC,
                             const std::string &_BCName,
                             Vector            &_vec2);

    bchandler_type          M_BCh_u;
    bchandler_type          M_BCh_d;
    bchandler_type          M_BCh_mesh;

    // interface operators BCs
    bchandler_type          M_BCh_du;
    bchandler_type          M_BCh_du_inv;

    bchandler_type          M_BCh_dz;
    bchandler_type          M_BCh_dz_inv;

    bchandler_type          M_BCh_dp;
    bchandler_type          M_BCh_dp_inv;

    fluid_type              M_fluid;
    solid_type              M_solid;

    quasi_newton_type       M_reducedLinFluid;

    dof_interface_type      M_dofFluidToStructure;
    dof_interface_type      M_dofStructureToSolid;
    dof_interface_type      M_dofStructureToHarmonicExtension;
    dof_interface_type      M_dofHarmonicExtensionToFluid;
    dof_interface_type      M_dofStructureToReducedFluid;
    dof_interface_type      M_dofReducedFluidToStructure;

    bc_vector_interface     M_bcvFluidLoadToStructure;
    bc_vector_interface     M_bcvStructureDispToSolid;
    bc_vector_interface     M_bcvStructureDispToHarmonicExtension;
    bc_vector_interface     M_bcvHarmonicExtensionVelToFluid;
    bc_vector_interface     M_bcvStructureToFluid;
    bc_vector_interface     M_bcvStructureToReducedFluid;
    bc_vector_interface     M_bcvReducedFluidToStructure;

    bc_vector_interface     M_bcvDerHarmonicExtensionVelToFluid;
    bc_vector_interface     M_bcvDerFluidLoadToStructure;
    bc_vector_interface     M_bcvDerFluidLoadToFluid;
    bc_vector_interface     M_bcvDerStructureDispToSolid;
    bc_vector_interface     M_bcvDerReducedFluidLoadToStructure;
    bc_vector_interface     M_bcvDerStructureAccToReducedFluid;

    Vector                  M_dispStruct;
    Vector                  M_dispStructOld;
    Vector                  M_velo;

    SolverAztec             M_solverAztec;

    Real                    M_time;

    UInt                    M_nbEval;

private:

    UInt                    M_reducedFluid;
    UInt                    M_method;
    Preconditioner          M_precond;
    DDNPreconditioner       M_DDNprecond;
};

typedef boost::shared_ptr<FSIOperator> oper_fsi_ptr;
typedef singleton<factory<FSIOperator,  std::string> > FSIFactory;

/*!
   \def  FOR_EACH_INTERFACE_DOF( Expr )

   \c FOR_EACH_INTERFACE_DOF is an helper macro to ease the
   the computation of quantities like the residual on the interface
   of the Fluid and Structure.

 */
//    UInt iBCf = M_fluid->BCh_fluid().getBCbyName("Wall");

#define FOR_EACH_INTERFACE_DOF( Expr )                              \
{\
                                                                    \
    UInt iBCf = M_fluid->BCh_HarmonicExtension().getBCbyName("Interface"); \
                                                                    \
    BCBase const &BC_fluidInterface = M_fluid->BCh_HarmonicExtension()[iBCf];   \
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
        for (UInt jDim = 0; jDim < nDimF; ++jDim)                   \
        {                                                           \
            ( Expr );                                               \
        }                                                           \
    }                                                               \
}

}
#endif
