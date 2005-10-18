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

#ifndef _OPERFS
#define _OPERFS
#include <life/lifecore/GetPot.hpp>

namespace LifeV
{

  template <typename F, typename S> 
  class operFS;

  
  template <typename F, typename S> 
  class DataJacobian {
    
  public:
    
    DataJacobian(operFS<F,S>* oper):
      _pFS(oper){}
    
    operFS<F,S>* _pFS;
    
  };


  //
  // Fluid-Structure operator Class
  //
  template <typename F, typename S> 
  class operFS {

  public:
    
    operFS(F& fluid, S& solid, BCHandler& bcHdu,  BCHandler& bcHdp, BCHandler& sum_dp, BCHandler& bcHdz, const GetPot& dataFile);

    //
    void eval(Vector& dispNew, Vector& veloStruct, const Vector& disp,int status);

    //
    void evalResidual(Vector &res, const Vector& sol, int iter);
    
    //
    void updateJacobian(Vector& sol,int iter);
    
    //
    void solveJac(Vector &step, const Vector& res, double& linear_rel_tol);
    
    //
    void solveLinearFluid();

    //
    void solveLinearSolid();
    
    //
    UInt nbEval();

    //
    Vector& da() {return _da;}
  
    F& _fluid;
  
    S& _solid;

    Vector _dispStruct;
    Vector _velo;
    Vector _dz, _da;
    Vector _rhs_dz;
    
    UInt _nbEval;
    
    BCHandler& _bcHdu;
    BCHandler& _bcHdp;
    BCHandler& _sum_dp;
    BCHandler& _bcHdz;
    
    UInt _FSIalgo;

    DataJacobian<F,S> _dataJacobian;
    
    void setTime(const Real& time);

  private:

    Real _time;
    
  };
  
  template <typename F, typename S> 
  void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);



  template <typename F, typename S> 
  operFS<F,S>::operFS(F& fluid, S& solid, BCHandler& bcHdu, BCHandler& bcHdp, BCHandler& sum_dp, BCHandler& bcHdz, const GetPot& dataFile):
    _fluid(fluid),
    _solid(solid),
    _dispStruct( 3*_solid.dDof().numTotalDof() ),
    _velo( 3*_solid.dDof().numTotalDof() ),
    _dz( 3*_solid.dDof().numTotalDof() ), 
    _da( 3*_solid.dDof().numTotalDof() ), 
    _rhs_dz( 3*_solid.dDof().numTotalDof() ),
    _nbEval(0),
    _bcHdu(bcHdu), 
    _bcHdp(bcHdp), 
    _sum_dp(sum_dp),
    _bcHdz(bcHdz),
    _FSIalgo( dataFile( "FSIalgo", 0 ) ),
    _dataJacobian(this) {}

  template <typename F, typename S> 
  void  operFS<F,S>::eval(Vector& dispNew, Vector& velo, const Vector& disp, int status) {

    if(status) _nbEval = 0; // new time step
    _nbEval++ ;
    
    _solid.disp() = disp;
    
    _fluid.updateMesh(_time);
    _fluid.iterate(_time);
    
    
    _solid.setRecur(0);
    _solid.iterate();
 
    dispNew = _solid.disp();
    velo    = _solid.w();
    
    std::cout << "                ::: norm(disp     ) = " << norm_inf(disp) << std::endl;
    std::cout << "                ::: norm(dispNew  ) = " << norm_inf(dispNew) << std::endl;
    std::cout << "                ::: norm(velo     ) = " << norm_inf(velo) << std::endl;
    
  }


  // Residual evaluation
  //
  template <typename F, typename S> 
  void  operFS<F,S>::evalResidual(Vector &res, const Vector& disp, int iter) {
    
    int status = 0;
    if(iter == 0) status = 1;
    std::cout << "*** Residual computation g(x_" << iter <<" )";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;
    eval(_dispStruct,_velo,disp,status);
    res = disp - _dispStruct;
  }
  
  
  //
  template <typename F, typename S> 
  void  operFS<F,S>::updateJacobian(Vector& sol,int iter) 
  {}
  
  
  //
  template <typename F, typename S> 
  void  operFS<F,S>::solveJac(Vector &step, const Vector& res, double& linear_rel_tol) {

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
    AZ_set_MATFREE(J, &_dataJacobian, my_matvecJacobian<F,S>);
    
    std::cout << "  o-  Solving Jacobian system... ";
    Chrono chrono;

    for (UInt i=0;i<dim_res; ++i)
      step[i]=0.0;

    chrono.start();
    AZ_iterate(&step[0], const_cast<double*>( &res[0] ), options, params, status, proc_config, J, NULL, NULL);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;


    AZ_matrix_destroy(&J);
  }



  //
  template <typename F, typename S> 
  void  operFS<F,S>::solveLinearFluid() {
    
    switch ( _FSIalgo ) 
      {
      case 1: // Newton: with shape derivative terms
	_fluid.iterateLin(_time, _bcHdu);
	break;
      case 2: // Quasi-Newton: without shape derivative terms (Linearized NS)
	_fluid.iterateLin(_time, _bcHdu, 0);
	break;
      case 3: // Quasi-Newton: simplified fluid 
	_fluid.iterateReducedLin(_time, _bcHdp, _sum_dp );
	break;
      default:
	ERROR_MSG("This FSI algorithm is not yet implemented");
      }
  }

  //
  template <typename F, typename S> 
  void  operFS<F,S>::solveLinearSolid() {

    _rhs_dz = ZeroVector( _rhs_dz.size() );
    
    if ( !_bcHdz.bdUpdateDone() )
      _bcHdz.bdUpdate(_solid.mesh(),_solid.feBd(),_solid.dof());
    bcManageVector(_rhs_dz,_solid.mesh(),_solid.dof(),_bcHdz,_solid.feBd(), 1.0, 1.0);
    
    Real tol=1.e-10;
    
    _solid.setRecur(1);
    
    _solid.solveJac(_dz, _rhs_dz, tol);
    
  }



  template <typename F, typename S> 
  UInt  operFS<F,S>::nbEval() {
    return _nbEval;
  }

  template <typename F, typename S> 
  void  operFS<F,S>::setTime(const Real& time) {
    _time = time;
  }
  
  template <typename F, typename S> 
  void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]) {


    // Extraction of data from J
    DataJacobian<F,S>* my_data = static_cast< DataJacobian<F,S>* >(AZ_get_matvec_data(J));

    UInt dim = my_data->_pFS->_dz.size();
    
    double xnorm =  AZ_gvector_norm(dim,-1,z,proc_config);
    std::cout << " ***** norm (z)= " << xnorm << std::endl<< std::endl;

    double dt = my_data->_pFS->_fluid.timestep();
    double dti2 = 1.0/( dt * dt) ;
    
    if ( xnorm == 0.0 ) {
      for (int i=0; i <(int)dim; ++i)
	Jz[i] =  0.0;
    }
    else {
      
      switch ( my_data->_pFS->_FSIalgo ) 
	{
	case 1: 
	case 2:
	  for (int i=0; i <(int)dim; ++i) 
	    my_data->_pFS->_solid.disp()[i] =  z[i];
	  my_data->_pFS->_fluid.updateDispVelo();
	  break;
	case 3:
	  for (int i=0; i <(int)dim; ++i) 
	    my_data->_pFS->_da[i] =  - my_data->_pFS->_fluid.density() * z[i] * dti2;
	  break;
	default:
	  ERROR_MSG("This FSI algorithm is not yet implemented");
	}
      
      my_data->_pFS->solveLinearFluid();
      my_data->_pFS->solveLinearSolid();
      for (int i=0; i <(int)dim; ++i)
	Jz[i] =  z[i]-my_data->_pFS->_dz[i];
    }
    std::cout << " ***** norm (Jz)= " << AZ_gvector_norm(dim,-1,Jz,proc_config)<< std::endl<< std::endl;
  }


}
#endif
