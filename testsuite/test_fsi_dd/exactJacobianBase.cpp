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


#include "exactJacobianBase.hpp"


namespace LifeV
{
exactJacobian::exactJacobian(GetPot &_dataFile):
    operFS(_dataFile),
    M_dz     (3*M_solid.dDof().numTotalDof()),
    M_rhs_dz (3*M_solid.dDof().numTotalDof()),
    M_dataJacobian(this)
{
    setUpBC();
}

exactJacobian::~exactJacobian()
{}


//
// Residual computation
//

void exactJacobian::eval(const Vector &_disp,
                         const int     _status)
{
    this->M_solid.d() = _disp;

    for (int ii = 0; ii < M_solid().d().size(); ++ii)
        std::cout << M_solid().d()[ii] << std::endl;

    this->M_fluid.updateMesh(time());
    this->M_fluid.iterate   (time());

    this->M_fluid.postProcess();


    this->M_solid.setRecur(0);
    this->M_solid.iterate();

    this->M_solid.postProcess();
}


void exactJacobian::evalResidual(const Vector &_disp,
                                 const int     _iter,
                                 Vector       &_res)
{
    int status = 0;

    if(_iter == 0) status = 1;

    std::cout << "*** Fixed Point: Residual computation g(x_" << _iter <<")";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;

    eval(M_dispStruct, status);

    M_dispStruct = this->M_solid.d();
    M_velo       = this->M_solid.w();

    std::cout << "                ::: norm(disp     ) = " << maxnorm(_disp)  << std::endl;
    std::cout << "                ::: norm(dispNew  ) = " << maxnorm(M_dispStruct) << std::endl;
    std::cout << "                ::: norm(velo     ) = " << maxnorm(M_velo) << std::endl;

    std::cout << "Max ResidualF   = " << maxnorm(M_fluid.residual())
              << std::endl;
    std::cout << "Max ResidualS   = " << maxnorm(M_solid.residual())
              << std::endl;

    _res = _disp - M_dispStruct;
}


//
// Boundary conditions setup
//


void exactJacobian::setUpBC()
{
    std::cout << "Boundary Conditions setup ... ";

    UInt dim_solid = this->M_solid.dDof().numTotalDof();
    UInt dim_fluid = this->M_fluid.uDof().numTotalDof();

    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================
    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    BCVectorInterface g_wall(this->M_fluid.residual(),
                             dim_fluid,
                             M_dofFluidToStructure);
    //
    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    BCVectorInterface displ(this->M_solid.d(),
                            dim_solid,
                            M_dofStructureToFluidMesh);
    //
    // Passing data from structure to the fluid: solid velocity at the interface velocity
    //
    BCVectorInterface u_wall(this->M_fluid.wInterpolated(),
                             dim_fluid,
                             M_dofMeshToFluid);
    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================

    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    BCFunctionBase bcf(fZero);
    M_BCh_mesh.addBC("Interface", 1, Essential, Full, displ, 3);
    M_BCh_mesh.addBC("Top",       3, Essential, Full, bcf,   3);
    M_BCh_mesh.addBC("Base",      2, Essential, Full, bcf,   3);
    M_BCh_mesh.addBC("Edges",    20, Essential, Full, bcf,   3);


    // Boundary conditions for the fluid velocity
    BCFunctionBase in_flow(u2);
    M_BCh_u.addBC("Wall",   1,  Essential, Full, u_wall,  3);
    M_BCh_u.addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    M_BCh_u.addBC("Edges",  20, Essential, Full, bcf,     3);

    // Boundary conditions for the solid displacement
    M_BCh_d.addBC("Interface", 1, Natural,   Full, g_wall, 3);
    M_BCh_d.addBC("Top",       3, Essential, Full, bcf,    3);
    M_BCh_d.addBC("Base",      2, Essential, Full, bcf,    3);

    BCVectorInterface du_wall(M_fluid.dwInterpolated(), dim_fluid, M_dofMeshToFluid);

    // Passing data from fluid to the structure: du -> dz
    //
    BCVectorInterface dg_wall(M_fluid.residual(), dim_fluid, M_dofFluidToStructure);

    // Boundary conditions for du
    M_BCh_du.addBC("Wall",   1,  Essential, Full, du_wall,  3);
    M_BCh_du.addBC("Edges",  20, Essential, Full, bcf,      3);


    // Boundary conditions for dz
    M_BCh_dz.addBC("Interface", 1, Natural,   Full, dg_wall, 3);
    M_BCh_dz.addBC("Top",       3, Essential, Full, bcf,  3);
    M_BCh_dz.addBC("Base",      2, Essential, Full, bcf,  3);


}


//
// new step computation resolution
//


void  exactJacobian::solveJac(const Vector  &_res,
                               const double   _linearRelTol,
                               Vector        &_muk)
{
 // AZTEC specifications for the second system
  int    data_org[AZ_COMM_SIZE];   // data organisation for J
  int    proc_config[AZ_PROC_SIZE];  // Processor information:
  int    options[AZ_OPTIONS_SIZE];   // Array used to select solver options.
  double params[AZ_PARAMS_SIZE];     // User selected solver paramters.
  double status[AZ_STATUS_SIZE];     // Information returned from AZ_solve()

  AZ_set_proc_config(proc_config, AZ_NOT_MPI);

  // data_org assigned "by hands": no parallel computation is performed
  UInt dim_res = _res.size();
  data_org[AZ_N_internal]= dim_res;
  data_org[AZ_N_border]= 0;
  data_org[AZ_N_external]= 0;
  data_org[AZ_N_neigh]= 0;

 // Recovering AZTEC defaults options and params
  AZ_defaults(options,params);

  // Fixed Aztec options for this linear system
  options[AZ_solver]   = AZ_gmres;
  options[AZ_output]   = 1;
  options[AZ_poly_ord] = 5;
  options[AZ_kspace]   = 40;
  options[AZ_conv]     = AZ_rhs;
  params[AZ_tol]       = _linearRelTol;

  //AZTEC matrix for the jacobian
  AZ_MATRIX *J;
  J = AZ_matrix_create(dim_res);

  // data containing the matrices C, D, trD and H as pointers
  // are passed through A_ii and pILU_ii:
  AZ_set_MATFREE(J, &M_dataJacobian, my_matvecJacobianEJ);

  std::cout << "  o-  Solving Jacobian system... ";
  Chrono chrono;

  for (UInt i=0;i<dim_res; ++i)
    _muk[i]=0.0;

  chrono.start();
  AZ_iterate(&_muk[0], &_res[0], options, params, status, proc_config, J, NULL, NULL);
  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;

  AZ_matrix_destroy(&J);
}


void  exactJacobian::solveLinearFluid()
{
    this->M_fluid.iterateLin(time(), M_BCh_du);
}


//


void  exactJacobian::solveLinearSolid()
{
    M_rhs_dz = 0.;
    M_dz     = 0.;

    if ( !M_BCh_dz.bdUpdateDone() )
        M_BCh_dz.bdUpdate(this->M_solid.mesh(),
                          this->M_solid.feBd(),
                          this->M_solid.dof());

    bcManageVector(M_rhs_dz,
                   this->M_solid.mesh(),
                   this->M_solid.dof(),
                   M_BCh_dz,
                   this->M_solid.feBd(),
                   1., 1.);

    Real tol       = 1.e-10;

    std::cout << "rhs_dz norm = " << maxnorm(M_rhs_dz) << std::endl;
    this->M_solid.setRecur(1);
    this->M_solid.solveJac(M_rhs_dz, tol, M_dz);
    std::cout << "dz norm     = " << maxnorm(M_dz) << std::endl;
}


void my_matvecJacobianEJ(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]) {


  // Extraction of data from J
  dataJacobianEJ* my_data = static_cast< dataJacobianEJ* >(AZ_get_matvec_data(J));

  UInt dim = my_data->M_pFS->dz().size();

  double xnorm =  AZ_gvector_norm(dim,-1,z,proc_config);
  std::cout << " ***** norm (z)= " << xnorm << std::endl << std::endl;

  if ( xnorm == 0.0 ) {
    for (int i=0; i <(int)dim; ++i)
      Jz[i] =  0.0;
  }
  else {
    for (int i=0; i <(int)dim; ++i) {
      my_data->M_pFS->solid().d()[i] =  z[i];
    }
    my_data->M_pFS->fluid().updateDispVelo();
    my_data->M_pFS->solveLinearFluid();
    my_data->M_pFS->solveLinearSolid();
    for (int i=0; i <(int)dim; ++i)
      Jz[i] =  z[i] - my_data->M_pFS->dz()[i];
  }
  std::cout << " ***** norm (Jz)= " << AZ_gvector_norm(dim,-1,Jz,proc_config)
            << std::endl << std::endl;
}

}
