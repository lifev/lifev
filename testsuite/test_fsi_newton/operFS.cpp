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
#include "operFS.hpp"

namespace LifeV
{

using namespace std;

operFS::operFS(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
	       VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid,
	       BCHandler& BCh_du, BCHandler& BCh_dz):
  _fluid(fluid),
     _solid(solid),
     _dispStruct( 3*_solid.dDof().numTotalDof() ),
     _velo( 3*_solid.dDof().numTotalDof() ),
     _dz( 3*_solid.dDof().numTotalDof() ),
     _rhs_dz( 3*_solid.dDof().numTotalDof() ),
     _nbEval(0),
     _BCh_du(BCh_du),
     _BCh_dz(BCh_dz),
     _dataJacobian(this) {}


void operFS::eval(Vector& dispNew, Vector& velo, const Vector& disp, int status) {


  if(status) _nbEval = 0; // new time step
  _nbEval++ ;

  _solid.d() = disp;

  _fluid.updateMesh(_time);
  _fluid.iterate(_time);

  _solid.setRecur(0);
  _solid.iterate();

  dispNew = _solid.d();
  velo    = _solid.w();

  std::cout << "                ::: norm(disp     ) = " << maxnorm(disp) << std::endl;
  std::cout << "                ::: norm(dispNew  ) = " << maxnorm(dispNew) << std::endl;
  std::cout << "                ::: norm(velo     ) = " << maxnorm(velo) << std::endl;

}


// Residual evaluation
//
void operFS::evalResidual(const Vector& disp, int iter, Vector &res) {

  int status = 0;
  if(iter == 0) status = 1;
  std::cout << "*** Residual computation g(x_" << iter <<" )";
  if (status) std::cout << " [NEW TIME STEP] ";
  std::cout << std::endl;
  eval(_dispStruct,_velo,disp,status);
  res = disp - _dispStruct;
}


//
void  operFS::updateJac(Vector& sol,int iter) {
}


//
void  operFS::solveJac(const Vector& res, double& linear_rel_tol, Vector &step) {

 // AZTEC specifications for the second system
  int    data_org[AZ_COMM_SIZE];   // data organisation for J
  int    proc_config[AZ_PROC_SIZE];  // Processor information:
  int    options[AZ_OPTIONS_SIZE];   // Array used to select solver options.
  double params[AZ_PARAMS_SIZE];     // User selected solver paramters.
  double status[AZ_STATUS_SIZE];     // Information returned from AZ_solve()

  AZ_set_proc_config(proc_config, AZ_NOT_MPI);

  // data_org assigned "by hands": no parallel computation is performed
  UInt dim_res = res.size();
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
  params[AZ_tol]       = linear_rel_tol;

  //AZTEC matrix for the jacobian
  AZ_MATRIX *J;
  J = AZ_matrix_create(dim_res);

  // data containing the matrices C, D, trD and H as pointers
  // are passed through A_ii and pILU_ii:
  AZ_set_MATFREE(J, &_dataJacobian, my_matvecJacobian);

  std::cout << "  o-  Solving Jacobian system... ";
  Chrono chrono;

  for (UInt i=0;i<dim_res; ++i)
    step[i]=0.0;

  chrono.start();
  AZ_iterate(&step[0], &res[0], options, params, status, proc_config, J, NULL, NULL);
  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;


  AZ_matrix_destroy(&J);
}



//
void  operFS::solveLinearFluid() {

  _fluid.iterateLin(_time, _BCh_du);

}

//
void  operFS::solveLinearSolid() {

  _rhs_dz = 0.0;

  if ( !_BCh_dz.bdUpdateDone() )
    _BCh_dz.bdUpdate(_solid._mesh,_solid._feBd,_solid._dof);
  bcManageVector(_rhs_dz,_solid._mesh,_solid._dof,_BCh_dz,_solid._feBd, 1.0, 1.0);

  Real tol=1.e-10;

  _solid.setRecur(1);

  _solid.solveJac(_dz, tol, _rhs_dz);

//   for (UInt ii = 0; ii < _rhs_dz.size(); ++ii)
//       std::std::cout << _rhs_dz[ii] << " " << _dz[ii] << std::std::endl;

}


UInt operFS::nbEval() {
  return _nbEval;
}

void operFS::setTime(const Real& time) {
    _time = time;
}

void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]) {


  // Extraction of data from J
  DataJacobian* my_data = static_cast< DataJacobian* >(AZ_get_matvec_data(J));

  UInt dim = my_data->_pFS->_dz.size();

  double xnorm =  AZ_gvector_norm(dim,-1,z,proc_config);
  std::cout << " ***** norm (z)= " << xnorm << std::endl<< std::endl;

  if ( xnorm == 0.0 ) {
    for (int i=0; i <(int)dim; ++i)
      Jz[i] =  0.0;
  }
  else {
    for (int i=0; i <(int)dim; ++i) {
      my_data->_pFS->_solid.d()[i] =  z[i];
    }
    my_data->_pFS->_fluid.updateDispVelo();
    my_data->_pFS->solveLinearFluid();
    my_data->_pFS->solveLinearSolid();
    for (int i=0; i <(int)dim; ++i)
      Jz[i] =  z[i]-my_data->_pFS->_dz[i];
  }
  std::cout << " ***** norm (Jz)= " << AZ_gvector_norm(dim,-1,Jz,proc_config)<< std::endl<< std::endl;
}
}
