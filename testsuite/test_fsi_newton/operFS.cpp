#include "operFS.hpp"


operFS::operFS(NavierStokesAleSolverPC< RegionMesh3D<LinearTetra> >& fluid, 
	       VenantKirchhofSolver< RegionMesh3D<LinearTetra> >& solid, 
	       BC_Handler& BCh_du, BC_Handler& BCh_dz):
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
  
  _solid.d().vec() = disp;

  _fluid.updateMesh(); 
  _fluid.iterate();  
  
  _solid._recur=0;
  _solid.iterate(); 
 
  dispNew = _solid.d().vec(); 
  velo    = _solid.w().vec(); 
 
  cout << "                ::: norm(disp     ) = " << maxnorm(disp) << endl;
  cout << "                ::: norm(dispNew  ) = " << maxnorm(dispNew) << endl;
  cout << "                ::: norm(velo     ) = " << maxnorm(velo) << endl;

}


// Residual evaluation
// 
void operFS::evalResidual(Vector& res, const Vector& disp, int iter) {

  int status = 0;
  if(iter == 0) status = 1;
  cout << "*** Residual computation g(x_" << iter <<" )";
  if (status) cout << " [NEW TIME STEP] ";
  cout << endl;
  eval(_dispStruct,_velo,disp,status);
  res = disp - _dispStruct;

}


//  
void  operFS::updateJac(Vector& sol,int iter) {
}


//
void  operFS::solveJac(Vector& step, const Vector& res, double& linear_rel_tol) {

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

  cout << "  o-  Solving Jacobian system... ";
  Chrono chrono;

  for (UInt i=0;i<dim_res; ++i)
    step[i]=0.0;

  chrono.start();
  AZ_iterate(&step[0], &res[0], options, params, status, proc_config, J, NULL, NULL);
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;


  AZ_matrix_destroy(&J);
}



//
void  operFS::solveLinearFluid() {
 
  _fluid.iterateLin(_BCh_du);

}

//
void  operFS::solveLinearSolid() {

  _rhs_dz = 0.0;
  
  if ( !_BCh_dz.bdUpdateDone() )  
    _BCh_dz.bdUpdate(_solid._mesh,_solid._feBd,_solid._dof);
  bc_manage_vector(_rhs_dz,_solid._mesh,_solid._dof,_BCh_dz,_solid._feBd, 1.0, 1.0);
  
  Real tol=1.e-10;

  _solid._recur = 1;

  _solid.solveJac(_dz, _rhs_dz, tol);

}


UInt operFS::nbEval() {
  return _nbEval;
}


void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]) {


  // Extraction of data from J
  DataJacobian* my_data = static_cast< DataJacobian* >(AZ_get_matvec_data(J));

  UInt dim = my_data->_pFS->_dz.size();
 
  double xnorm =  AZ_gvector_norm(dim,-1,z,proc_config);
  cout << " ***** norm (z)= " << xnorm << endl<< endl;

  if ( xnorm == 0.0 ) {
    for (int i=0; i <(int)dim; ++i) 
      Jz[i] =  0.0; 
  }
  else {
    for (int i=0; i <(int)dim; ++i) {
      my_data->_pFS->_solid.d().vec()[i] =  z[i];
    } 
    my_data->_pFS->_fluid.updateDispVelo();
    my_data->_pFS->solveLinearFluid();
    my_data->_pFS->solveLinearSolid();
    for (int i=0; i <(int)dim; ++i)
      Jz[i] =  z[i]-my_data->_pFS->_dz[i];
  }
  cout << " ***** norm (Jz)= " << AZ_gvector_norm(dim,-1,Jz,proc_config)<< endl<< endl;
}
