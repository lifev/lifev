/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file NavierStokesAleSolverPC.h
  \author M.A. Fernandez and A. Gauthier
  \date 05/2003 
  \version 1.0

  \brief This file contains a NavierStokes ALE solver  class which implements  a semi-implicit scheme with 
         an exact factorization. Preconditioning of the  Schur Complement is done by an algebraic 
         Chorin-Temam pressure-corrected preconditioner 

*/

#ifndef _NAVIERSTOKESALESOLVERPC_HH
#define _NAVIERSTOKESALESOLVERPC_HH

#include "NavierStokesAleHandler.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include "values.hpp"
#include "pattern.hpp"
#include "assemb.hpp"
#include "bc_manage.hpp"
#include "algebraic_facto.hpp"
#include "chrono.hpp"
#include "dataAztec.hpp"


/*! 
  \class NavierStokesSolverPC

   This class implements an NavierStokes ALE solver via exact factorization. Preconditioning of the 
   Schur Complement is done by an algebraic Chorin-Temam pressure-corrected preconditioner 

*/
template<typename Mesh>
class NavierStokesAleSolverPC:
public NavierStokesAleHandler<Mesh> {
 
 public:

  typedef typename NavierStokesHandler<Mesh>::Function Function;
  //! Constructor
  /*!
    \param data_file GetPot data file
    \param refFE_u reference FE for the velocity
    \param refFE_p reference FE for the pressure
    \param Qr_u volumic quadrature rule for the velocity
    \param bdQr_u surface quadrature rule for the velocity 
    \param Qr_p volumic quadrature rule for the pressure
    \param bdQr_p surface quadrature rule for the pressure
    \param BCh_u boundary conditions for the velocity
    \param BCh_mesh boundary conditions for the motion harmonic extension
  */
  NavierStokesAleSolverPC(const GetPot& data_file, const RefFE& refFE_u, const RefFE& refFE_p, const QuadRule& Qr_u,
			  const QuadRule& bdQr_u, const QuadRule& Qr_p, const QuadRule& bdQr_p, BC_Handler& BCh_u, 
			  BC_Handler& BCh_mesh);

  //! Update the right  hand side  for time advancing 
  /*! 
    \param source volumic source  
    \param time present time
  */
  void timeAdvance(const Function source, const Real& time);

  //! Update convective term, bc treatment and solve the linearized ns system
  void iterate(); 

  void iterateTransp();

  void iterateLin(BC_Handler& BCh_du);

  Vector& residual();

 private:


  //! Block pattern of M_u
  MSRPatt _pattM_u_block;
  
  //! Pattern for M
  MixedPattern<3,3,MSRPatt> _pattM_u;

  //! Block pattern of C
  MSRPatt _pattC;

  //! Block pattern of D: Vdiv operator
  CSRPatt _pattD_block;

  //! Pattern for D
  MixedPattern<1,3,CSRPatt> _pattD;

  //! Block  pattern of trD: transpose Vdiv operator trVdiv
  CSRPatt _pattDtr_block;
  
  //! Pattern  of trD
  MixedPattern<3,1,CSRPatt> _pattDtr;

  //! Matrix D: Vdiv operator
  MixedMatr<1,3,CSRPatt,double> _D;
 
  //! Matrix trD: transpose Vdiv operator trVdiv
  MixedMatr<3,1,CSRPatt,double> _trD;

  //! Matrix trD: transpose Vdiv operator trVdiv
  MixedMatr<3,1,CSRPatt,double> _trDAux;

  //! Matrix HinvDtr:  H^{-1}D^T
  MixedMatr<3,1,CSRPatt,double> _HinvDtr;

  //! Matrix M_u: Vmass
  MixedMatr<3,3,MSRPatt,double> _M_u;
 
  //! Matrix HinvC: H^{-1}C
  MSRMatr<double> _HinvC;
 
  //! Matrix CStokes
  MSRMatr<double> _CStokes;
  
  //! Matrix C: 1/dt*Vmass + mu*Vstiff operator + Convective_term
  MSRMatr<double> _C;
  MSRMatr<double> _CAux;
 

  //! H diag matrix: H= diag( _M_u )/sum( diag( _M_u ) ) where _M_u = mass * rho / dt
  vector<double>  _H;

  //! Elementary matrices and vectors
  ElemMat _elmatC; //velocity stiffnes 
  ElemMat _elmatM_u; //velocity mass
  ElemMat _elmatDtr; // vel_i * pres_j
  ElemVec _elvec; // Elementary right hand side
  ElemVec _elvec_du; // Elementary right hand side for the linearized velocity
  ElemVec _elvec_dp; // Elementary right hand side for the linearized pressure
  ElemVec _w_loc; // Elementary mesh velocity  
  ElemVec _uk_loc; // Elementary velocity 
  ElemVec _pk_loc; // Elementary pressure
  ElemVec _convect; // Elementary convection velocity
  ElemVec _d_loc; // Elementary displacement for right hand side
  ElemVec _dw_loc; // Elementary mesh velocity for right hand side

  //! The pressure
  ScalUnknown<Vector> _dp;

  //! The velocity
  PhysVectUnknown<Vector> _un;

  //! The velocity
  PhysVectUnknown<Vector> _du;

  //! Right  hand  side for the velocity
  PhysVectUnknown<Vector> _f_u;

  //! Right  hand  side for the velocity
  PhysVectUnknown<Vector> _f_duWithOutBC;

  //! Right  hand  side for pressure
  ScalUnknown<Vector> _f_p;

  //! Right  hand  side for the velocity
  PhysVectUnknown<Vector> _f_uWithOutBC;

  //! The residual of the momentum equations
  PhysVectUnknown<Vector> _residual_u;
   

  //!  This vector contains the product C^{-1}*trD*P where P is the pressure, solution
  //!  of the system (ii).
  PhysVectUnknown<Vector> _invCtrDP;

  DataAztec _dataAztec_i;
  DataAztec _dataAztec_ii;
  DataAztec _dataAztec_s;

  //! DataFactorisation: data passed to matrix-vector product are stored in the class
  DataFactorisation< MSRMatr<double>, MixedMatr<1,3,CSRPatt,double>,
    MixedMatr<3,1,CSRPatt,double>, vector<double>, MSRMatr<double>, Vector> _factor_data;
  
  //! DataFactorisation: data passed to matrix-vector product are stored in the class
  DataFactorisation< MSRMatr<double>, MixedMatr<1,3,CSRPatt,double>,
    MixedMatr<3,1,CSRPatt,double>, vector<double>, MSRMatr<double>, Vector> _factor_data_jacobian;

};


//
//                                         IMPLEMENTATION
//
template<typename Mesh> NavierStokesAleSolverPC<Mesh>::
NavierStokesAleSolverPC(const GetPot& data_file, const RefFE& refFE_u, const RefFE& refFE_p, const QuadRule& Qr_u,
			const QuadRule& bdQr_u, const QuadRule& Qr_p, const QuadRule& bdQr_p, BC_Handler& BCh_u, 
			BC_Handler& BCh_mesh):
     NavierStokesAleHandler<Mesh>(data_file,refFE_u,refFE_p,Qr_u,bdQr_u,Qr_p,bdQr_p,BCh_u,BCh_mesh),  
     _pattM_u_block(_dof_u),
     _pattM_u(_pattM_u_block,"diag"),
     _pattC(_dof_u,3),
     _pattD_block(_dof_p,_dof_u),
     _pattD(_pattD_block),
     _pattDtr_block(_dof_u,_dof_p),
     _pattDtr(_pattDtr_block),
     _D(_pattD),
     _trD(_pattDtr), 
     _trDAux(_pattDtr),
     _HinvDtr(_pattDtr),
     _M_u(_pattM_u),
     _HinvC(_pattC),
     _CStokes(_pattC),
     _C(_pattC),
     _CAux(_pattC),
     _H(_pattC.nRows()), 
     _elmatC(_fe_u.nbNode,nDimensions,nDimensions), 
     _elmatM_u(_fe_u.nbNode,nDimensions,nDimensions),
     _elmatDtr(_fe_u.nbNode,nDimensions,0,_fe_p.nbNode,0,1), 
     _elvec(_fe_u.nbNode,nDimensions),
     _elvec_du(_fe_u.nbNode,nDimensions),
     _elvec_dp(_fe_p.nbNode,1),
     _w_loc(_fe_u.nbNode,nDimensions),   
     _uk_loc(_fe_u.nbNode,nDimensions), 
     _pk_loc(_fe_p.nbNode,1),
     _convect(_fe_u.nbNode,nDimensions),
     _d_loc(_fe_u.nbNode,nDimensions),
     _dw_loc(_fe_u.nbNode,nDimensions), 
     _dp(_dim_p),
     _un(_dim_u), 
     _du(_dim_u),
     _f_u(_dim_u),   
     _f_duWithOutBC(_dim_u),
     _f_p(_dim_p),  
     _f_uWithOutBC(_dim_u),
     _residual_u(_dim_u),
     _invCtrDP(_dim_u),
     _dataAztec_i(data_file,"fluid/aztec_i"),
     _dataAztec_ii(data_file,"fluid/aztec_ii"),
     _dataAztec_s(data_file,"fluid/aztec_s"),
     _factor_data(_C,_D,_trD,_H,_HinvC,_HinvDtr,_invCtrDP.vec(),_dataAztec_i,_dataAztec_s,_BCh_u.fullEssential(),1), 
     _factor_data_jacobian(_C,_D,_trD,_H,_HinvC,_HinvDtr,_invCtrDP.vec(),_dataAztec_i,_dataAztec_s,_BCh_u.fullEssential(),2) {

  
  cout << endl;
  cout << "O-  Pressure unknowns: " << _dim_p     << endl; 
  cout << "O-  Velocity unknowns: " << _dim_u     << endl<<endl;
  cout << "O-  Computing mass matrix... ";  
    
  Chrono chrono;
  chrono.start();

  // Number of velocity components  
  UInt nc_u=_u.nbcomp();

  // Initializing mass matrix
  _M_u.zeros();
 
  Real dti=1.0/_dt;
  
  // loop on volumes: assembling mass term
  for(UInt i=1; i<=_mesh.numVolumes(); ++i){
    _fe_u.updateFirstDerivQuadPt(_mesh.volumeList(i));
    _elmatM_u.zero();      
    mass(_rho*dti,_elmatM_u,_fe_u,0,0,nc_u);  
    for (UInt ic=0; ic<nc_u; ++ic){ 
      assemb_mat(_M_u,_elmatM_u,_fe_u,_dof_u,ic,ic);
    }
  }
    
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
}



template<typename Mesh>  
void NavierStokesAleSolverPC<Mesh>::
timeAdvance(const Function source, const Real& time) {

  _time = time;

  cout << endl;
  cout << "O== FLUID: Now we are at time "<< _time << " s." << endl;

  // Number of velocity components  
  UInt nc_u=_u.nbcomp();

  cout << "  o-  Updating mass term on right hand side... ";
  
  Chrono chrono;
  chrono.start();

   // Right hand side for the velocity at time
  _f_uWithOutBC.vec()=0.;
  
  // loop on volumes: assembling source term
  for(UInt i=1; i<=_mesh.numVolumes(); ++i){
    _elvec.zero();
    _fe_u.updateFirstDerivQuadPt(_mesh.volumeList(i));
    for (UInt ic=0; ic<nc_u; ++ic){ 
      compute_vec(source,_elvec,_fe_u,_time,ic); // compute local vector
      assemb_vec(_f_uWithOutBC.vec(),_elvec,_fe_u,_dof_u,ic); // assemble local vector into global one       
    }
  }
  _f_uWithOutBC.vec() += _M_u * _u.vec();

  // Save last mesh displacement and fluid velocity
  _dispOld.vec() = _disp.vec();
  _un.vec() = _u.vec();

  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
}


template<typename Mesh>    
void NavierStokesAleSolverPC<Mesh>::
iterate() {

  Chrono chrono;

  // Number of velocity components  
  UInt nc_u=_u.nbcomp();

  cout << "  o-  Updating matrices... ";
 
  chrono.start();
  
  //initialize matrices
  _D.zeros();
  _trD.zeros();
  _M_u.zeros();
  _C.zeros();

  Real dti=1.0/_dt;

  // Loop on elements
  for(UInt i = 1; i <= _mesh.numVolumes(); i++){
    
    _fe_p.update( _mesh.volumeList(i) ); // just to provide the id number in the
                                        // assem_mat_mixed
    _fe_u.updateFirstDerivQuadPt(_mesh.volumeList(i));
     
    // initialization of elementary matrices
    _elmatM_u.zero();
    _elmatDtr.zero(); 
    _elmatC.zero();
           
    // mass
    mass(_rho*dti,_elmatM_u,_fe_u,0,0,nc_u);
   
    // stiffness strain
    stiff_strain(2.0*_mu,_elmatC,_fe_u);
    _elmatC.mat() += _elmatM_u.mat();

    // Non linear term, Semi-implicit approach      
    // u_loc contains the velocity values in the nodes
    for (UInt k=0 ; k<(UInt)_fe_u.nbNode ; k++){
      UInt  iloc = _fe_u.patternFirst(k);
      for (UInt ic=0; ic<nc_u; ++ic){     
	UInt ig=_dof_u.localToGlobal(i,iloc+1)-1+ic*_dim_u;       
	_elvec.vec()[iloc+ic*_fe_u.nbNode] = _rho*(_un.vec()(ig)-_wInterp.vec()(ig)); 
      }
    }

    // ALE term: 0.5 div (u^n-w) u v
    mass_divw(0.5,_elvec,_elmatC,_fe_u,0,0,nc_u);

    // loop on velocity components
    for(UInt ic=0;ic<nc_u;ic++){
      
      for (UInt jc=0;jc<nc_u;jc++) {
	grad(jc,_elvec,_elmatC,_fe_u,_fe_u,ic,ic);
        assemb_mat(_C,_elmatC,_fe_u,_dof_u,ic,jc);
      }
      
     // mass
      assemb_mat(_M_u,_elmatM_u,_fe_u,_dof_u,ic,ic);
      
      // computing  - (p, \div v) term: the minus sign is in the inner computation
      grad(ic,1.0,_elmatDtr,_fe_u,_fe_p,ic,0);
      
      // assembling p div v term and transposed
      assemb_mat_mixed(_trD,_elmatDtr,_fe_u,_fe_p,_dof_u,_dof_p,ic,0);
      assemb_tr_mat_mixed(1.0,_D,_elmatDtr,_fe_p,_fe_u,_dof_p,_dof_u,0,ic); 
    }
  }

  _trDAux = _trD;
  _CAux = _C;

  // H diag matrix: H= diag( _M_u )/sum( diag( _M_u ) ) where _M_u = mass * rho / dt
  _H=_M_u.giveDiag();
  Real sum=accumulate(_H.begin(),_H.end(),0.0);
  // Matrix equilibration
  for (UInt i=0; i<_H.size();i++) 
    _H[i]=_H[i]*dti/sum;
  
  chrono.stop();
  cout << "done in "<< chrono.diff() << "s." << endl;


  // for BC treatment (done at each time-step)
  Real tgv=1.e02; 

  cout << "  o-  Applying boundary conditions... ";
  chrono.start();  
  _f_u.vec()=_f_uWithOutBC.vec(); 
  _BCh_u.bdUpdate(_mesh, _feBd_u, _dof_u);
  bc_manage(_C, _trD, _f_u.vec(), _mesh, _dof_u, _BCh_u, _feBd_u, tgv, _time);
  chrono.stop();
  cout << "done in " << chrono.diff() << "s." << endl;
 
   //matrices HinvDtr:
  MultInvDiag(_H, _trD, _HinvDtr);
  // ---------------
  // (i) C * V = F_V
  // ---------------
  // AZTEC specifications for each system
  int    data_org_i[AZ_COMM_SIZE];   // data organisation for C
  int    proc_config_i[AZ_PROC_SIZE];// Processor information:  
  int    options_i[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_i[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_i[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.
  AZ_set_proc_config(proc_config_i, AZ_NOT_MPI);

  //AZTEC matrix and preconditioner
  AZ_MATRIX *C;
  AZ_PRECOND *prec_C;

  int N_eq_i= 3*_dim_u; // number of DOF for each component
 // data_org assigned "by hands" while no parallel computation is performed
  data_org_i[AZ_N_internal]= N_eq_i;
  data_org_i[AZ_N_border]= 0;
  data_org_i[AZ_N_external]= 0;
  data_org_i[AZ_N_neigh]= 0;
  data_org_i[AZ_name]= DATA_NAME_AZTEC;

  // create matrix and preconditionner 
  C= AZ_matrix_create(N_eq_i);
  prec_C= AZ_precond_create(C, AZ_precondition, NULL);

  AZ_set_MSR(C, (int*) _pattC.giveRaw_bindx(), (double*) _C.giveRaw_value(), data_org_i, 0, NULL, AZ_LOCAL);

  _dataAztec_i.aztecOptionsFromDataFile(options_i,params_i); 

  //keep C factorisation and preconditioner reused in my_matvec
  options_i[AZ_keep_info]= 1;       
   
  // ---------------
  // (i) C * V = F_V
  // ---------------

  // intermediate velocity computation   
  cout << "  o-  Solving system (i)... "; 
  chrono.start();
  AZ_iterate(_u.giveVec(), _f_u.giveVec(), options_i,params_i, status_i,
  	     proc_config_i, C, prec_C, NULL);
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;


  // --------------------------------------------
  // (ii) (D*C^(-1)*trD) * P = D*C^{-1}*F_V = D*V
  // --------------------------------------------
  // AZTEC specifications for the second system
  int    proc_config_ii[AZ_PROC_SIZE];// Processor information:
  int    options_ii[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_ii[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_ii[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.

  AZ_set_proc_config(proc_config_ii, AZ_NOT_MPI);

  //AZTEC matrix for A_ii=(D*C^{-1}*trD)
  AZ_MATRIX *A_ii;
  AZ_PRECOND *pILU_ii;

  int N_eq_ii= _p.size();

  A_ii= AZ_matrix_create(N_eq_ii);
  // data containing the matrices C, D, trD and H as pointers
  // are passed through A_ii and pILU_ii:
  AZ_set_MATFREE(A_ii, &_factor_data,
		 my_matvec< MixedMatr<1,3,CSRPatt,double>, MixedMatr<3,1,CSRPatt,double>,
		 vector<double>, MSRMatr<double>, Vector >);

  pILU_ii= AZ_precond_create(A_ii, my_precSchur_PC<
			     MSRMatr<double>,
			     MixedMatr<1,3,CSRPatt,double>,
			     MixedMatr<3,1,CSRPatt,double>,
			     vector<double>,
			     MSRMatr<double>,
			     Vector>, &_factor_data);

  _dataAztec_ii.aztecOptionsFromDataFile(options_ii,params_ii);
   
  // user preconditioning:
  options_ii[AZ_precond]  = AZ_user_precond;

  // RHS of the linear system (ii)
  Vector vec_DV(_p.size());

 
  //matrices HinvC (depends on time):
  MultInvDiag(_H, _C, _HinvC);
 
  
  // RHS of the linear system (ii)
  vec_DV = _D*_u.vec();

   // case of pure Dirichlet BCs:
  if (_BCh_u.fullEssential()) {
    vec_DV[_dim_p-1]   = 1.0; // correction of the right hand side.
    _p.vec()[_dim_p-1] = 1.0; // pressure value at the last node.
  }

  cout << "  o-  Solving pressure system... ";
  chrono.start();
  AZ_iterate(_p.giveVec(), &vec_DV[0],  options_ii,params_ii, status_ii,
       proc_config_ii, A_ii, pILU_ii, NULL);

  
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
  
  // ----------------------------
  // (iii) V = V-(C^(-1)*trD) * P
  // ----------------------------

  // everything is done...
  _u.vec() = _u.vec() - _invCtrDP.vec();
  cout << "  o-  Velocity updated" << endl; 

  AZ_matrix_destroy(&A_ii);
  AZ_precond_destroy(&pILU_ii);
  AZ_matrix_destroy(&C);
  AZ_precond_destroy(&prec_C);

  _residual_u.vec() = _f_uWithOutBC.vec() - _CAux * _u.vec() - _trDAux   * _p.vec(); 

 
}



template<typename Mesh>    
void NavierStokesAleSolverPC<Mesh>::
iterateTransp() {

  Chrono chrono;

  _C = _CAux;

  // for BC treatment (done at each time-step)
  Real tgv=1.e02; 

  cout << "  o-  Applying boundary conditions... ";
  chrono.start();  
  _f_u.vec()=_f_uWithOutBC.vec(); 
  _BCh_u.bdUpdate(_mesh, _feBd_u, _dof_u);
  bc_manage(_C, _trD, _f_u.vec(), _mesh, _dof_u, _BCh_u, _feBd_u, tgv, _time);
  chrono.stop();
  cout << "done in " << chrono.diff() << "s." << endl;
 

  // ---------------
  // (i) C * V = F_V
  // ---------------
  // AZTEC specifications for each system
  int    data_org_i[AZ_COMM_SIZE];   // data organisation for C
  int    proc_config_i[AZ_PROC_SIZE];// Processor information:  
  int    options_i[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_i[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_i[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.
  AZ_set_proc_config(proc_config_i, AZ_NOT_MPI);

  //AZTEC matrix and preconditioner
  AZ_MATRIX *C;
  AZ_PRECOND *prec_C;

  int N_eq_i= 3*_dim_u; // number of DOF for each component
 // data_org assigned "by hands" while no parallel computation is performed
  data_org_i[AZ_N_internal]= N_eq_i;
  data_org_i[AZ_N_border]= 0;
  data_org_i[AZ_N_external]= 0;
  data_org_i[AZ_N_neigh]= 0;
  data_org_i[AZ_name]= DATA_NAME_AZTEC;

  // create matrix and preconditionner 
  C= AZ_matrix_create(N_eq_i);
  prec_C= AZ_precond_create(C, AZ_precondition, NULL);

  AZ_set_MSR(C, (int*) _pattC.giveRaw_bindx(), (double*) _C.giveRaw_value(), data_org_i, 0, NULL, AZ_LOCAL);

  _dataAztec_i.aztecOptionsFromDataFile(options_i,params_i); 

  //keep C factorisation and preconditioner reused in my_matvec
  options_i[AZ_keep_info]= 1;       
   
  // ---------------
  // (i) C * V = F_V
  // ---------------

  // intermediate velocity computation   
  cout << "  o-  Solving system (i)... "; 
  chrono.start();
  AZ_iterate(_u.giveVec(), _f_u.giveVec(), options_i,params_i, status_i,
  	     proc_config_i, C, prec_C, NULL);
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;


  // --------------------------------------------
  // (ii) (D*C^(-1)*trD) * P = D*C^{-1}*F_V = D*V
  // --------------------------------------------
  // AZTEC specifications for the second system
  int    proc_config_ii[AZ_PROC_SIZE];// Processor information:
  int    options_ii[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_ii[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_ii[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.

  AZ_set_proc_config(proc_config_ii, AZ_NOT_MPI);

  //AZTEC matrix for A_ii=(D*C^{-1}*trD)
  AZ_MATRIX *A_ii;
  AZ_PRECOND *pILU_ii;

  int N_eq_ii= _p.size();

  A_ii= AZ_matrix_create(N_eq_ii);
  // data containing the matrices C, D, trD and H as pointers
  // are passed through A_ii and pILU_ii:
  AZ_set_MATFREE(A_ii, &_factor_data,
		 my_matvec< MixedMatr<1,3,CSRPatt,double>, MixedMatr<3,1,CSRPatt,double>,
		 vector<double>, MSRMatr<double>, Vector >);

  pILU_ii= AZ_precond_create(A_ii, my_precSchur_PC<
			     MSRMatr<double>,
			     MixedMatr<1,3,CSRPatt,double>,
			     MixedMatr<3,1,CSRPatt,double>,
			     vector<double>,
			     MSRMatr<double>,
			     Vector>, &_factor_data);

  _dataAztec_ii.aztecOptionsFromDataFile(options_ii,params_ii);
   
  // user preconditioning:
  options_ii[AZ_precond]  = AZ_user_precond;

  // RHS of the linear system (ii)
  Vector vec_DV(_p.size());

   
  // RHS of the linear system (ii)
  vec_DV = _D*_u.vec();

   // case of pure Dirichlet BCs:
  if (_BCh_u.fullEssential()) {
    vec_DV[_dim_p-1]   = 1.0; // correction of the right hand side.
    _p.vec()[_dim_p-1] = 1.0; // pressure value at the last node.
  }

  cout << "  o-  Solving pressure system... ";
  chrono.start();
  AZ_iterate(_p.giveVec(), &vec_DV[0],  options_ii,params_ii, status_ii,
       proc_config_ii, A_ii, pILU_ii, NULL);

  
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
  
  // ----------------------------
  // (iii) V = V-(C^(-1)*trD) * P
  // ----------------------------

  // everything is done...
  _u.vec() = _u.vec() - _invCtrDP.vec();
  cout << "  o-  Velocity updated" << endl; 

  AZ_matrix_destroy(&A_ii);
  AZ_precond_destroy(&pILU_ii);
  AZ_matrix_destroy(&C);
  AZ_precond_destroy(&prec_C);




}


//
// Linearized iterate for Newton FSI
//
template<typename Mesh>    
void NavierStokesAleSolverPC<Mesh>::
iterateLin(BC_Handler& BCh_du) {
  
  Chrono chrono;
  
  // Number of velocity components  
  UInt nc_u=_u.nbcomp(), iloc, ig;

  cout << "  OOO-  LINEARIZED FLUID SYSTEM\n\n";

  cout << "    o-  Updating right hand side... ";
 
  //
  // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
  //
  chrono.start();
    
  //initialize right hand side
  _f_duWithOutBC.vec()=0.0; 
  _f_p.vec()=0.0; 

  // Loop on elements
  for(UInt i = 1; i <= _mesh.numVolumes(); i++){
    
    _fe_p.update( _mesh.volumeList(i) ); 
    _fe_u.updateFirstDerivQuadPt(_mesh.volumeList(i));
    
    // initialization of elementary vectors
    _elvec_du.zero();
    _elvec_dp.zero();
    
    for (UInt k=0 ; k<(UInt)_fe_u.nbNode ; k++){
      iloc = _fe_u.patternFirst(k);
      for (UInt ic=0; ic<nc_u; ++ic){     
	ig=_dof_u.localToGlobal(i,iloc+1)-1+ic*_dim_u;       
	_convect.vec()[iloc + ic*_fe_u.nbNode] = _un.vec()(ig)-_wInterp.vec()(ig);  // u^n - w^k local
	_w_loc.vec(  )[iloc + ic*_fe_u.nbNode] = _wInterp.vec()(ig);                // w^k local
	_uk_loc.vec( )[iloc + ic*_fe_u.nbNode] = _u.vec()(ig);                      // u^k local
	_d_loc.vec(  )[iloc + ic*_fe_u.nbNode] = _dInterp.vec()(ig);                // d local
	_dw_loc.vec( )[iloc + ic*_fe_u.nbNode] = _dwInterp.vec()(ig);               // dw local
      }
    }
    
    for (UInt k=0 ; k<(UInt)_fe_p.nbNode ; k++){
      iloc = _fe_p.patternFirst(k);
      ig   = _dof_p.localToGlobal(i,iloc+1)-1;
      _pk_loc.vec()[iloc] = _p.vec()(ig);  // p^k local
    }
    
    //
    // Elementary vectors
    // 
    
    //  - \rho ( \grad( u^n-w^k ):[I\div d - (\grad d)^T] u^k + ( u^n-w^k )^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  )
    source_mass1(-_rho, _uk_loc, _convect, _d_loc, _elvec_du, _fe_u);
    
    //  + \rho * ( \grad u^k dw, v  )
    source_mass2( _rho, _uk_loc, _dw_loc, _elvec_du, _fe_u);
    
    //  - ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  )
    source_stress(-1.0, _mu, _uk_loc, _pk_loc, _d_loc, _elvec_du, _fe_u, _fe_p);
    
    // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )  
    source_stress2( _mu, _uk_loc, _d_loc, _elvec_du, _fe_u);

    //  + ( (\grad u^k):[I\div d - (\grad d)^T] , q  ) 
    source_press(1.0, _uk_loc, _d_loc, _elvec_dp, _fe_u, _fe_p);
    

   
  
    //
    // Assembling
    // 
    
    // assembling presssure right hand side
    assemb_vec(_f_p.vec(),_elvec_dp,_fe_p,_dof_p,0);
  
    // loop on velocity components
    for(UInt ic=0;ic<nc_u;ic++)
      // assembling velocity right hand side
      assemb_vec(_f_duWithOutBC.vec(),_elvec_du,_fe_u,_dof_u,ic);
  }
  
  chrono.stop();
  cout << "done in "<< chrono.diff() << "s." << endl;

  cout << "  maxnorm (_f_duWithOutBC.vec()) = " <<  maxnorm(_f_duWithOutBC.vec())   << endl; 
  
  // for BC treatment (done at each time-step)
  Real tgv=1.e02; 

  cout << "    o-  Applying boundary conditions... ";
  chrono.start();  
  _C   = _CAux;
  _trD = _trDAux;

  _f_u.vec() = _f_duWithOutBC.vec();



  BCh_du.bdUpdate(_mesh, _feBd_u, _dof_u);
  bc_manage(_C, _trD, _f_u.vec(), _mesh, _dof_u, BCh_du, _feBd_u, tgv, _time);


  chrono.stop();
  cout << "done in " << chrono.diff() << "s." << endl;
  cout << "  maxnorm (_f_du.vec()) after BC= " <<  maxnorm(_f_u.vec())   << endl; 
  cout << "  maxnorm ( difference ) after BC= " <<  maxnorm( _f_duWithOutBC.vec() - _f_u.vec())   << endl; 
  
  //matrices HinvDtr:
  MultInvDiag(_H, _trD, _HinvDtr);
  // ---------------
  // (i) C * V = F_V
  // ---------------
  // AZTEC specifications for each system
  int    data_org_i[AZ_COMM_SIZE];   // data organisation for C
  int    proc_config_i[AZ_PROC_SIZE];// Processor information:  
  int    options_i[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_i[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_i[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.

  AZ_set_proc_config(proc_config_i, AZ_NOT_MPI);

  //AZTEC matrix and preconditioner
  AZ_MATRIX *C;
  AZ_PRECOND *prec_C;

  int N_eq_i= 3*_dim_u; // number of DOF for each component
 // data_org assigned "by hands" while no parallel computation is performed
  data_org_i[AZ_N_internal]= N_eq_i;
  data_org_i[AZ_N_border]= 0;
  data_org_i[AZ_N_external]= 0;
  data_org_i[AZ_N_neigh]= 0;
  data_org_i[AZ_name]= DATA_NAME_AZTEC;

  // create matrix and preconditionner 
  C= AZ_matrix_create(N_eq_i);
  prec_C= AZ_precond_create(C, AZ_precondition, NULL);

  AZ_set_MSR(C, (int*) _pattC.giveRaw_bindx(), (double*) _C.giveRaw_value(), data_org_i, 0, NULL, AZ_LOCAL);

  _dataAztec_i.aztecOptionsFromDataFile(options_i,params_i); 

  //keep C factorisation and preconditioner reused in my_matvec
  options_i[AZ_keep_info]= 1;       
  

  // ---------------
  // (i) C * V = F_V
  // ---------------
  options_i[AZ_recursion_level]=1;

  _du.vec()=0.0;

  // intermediate velocity computation   
  cout << "  o-  Solving system (i)... "; 
  chrono.start();
  AZ_iterate(_du.giveVec(), _f_u.giveVec(), options_i,params_i, status_i,
  	     proc_config_i, C, prec_C, NULL);
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
  
  //options_i[AZ_recursion_level]=0;

  // --------------------------------------------
  // (ii) (D*C^(-1)*trD) * P = D*C^{-1}*F_V = D*V
  // --------------------------------------------
  // AZTEC specifications for the second system
  int    proc_config_ii[AZ_PROC_SIZE];// Processor information:
  int    options_ii[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_ii[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_ii[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.

  AZ_set_proc_config(proc_config_ii, AZ_NOT_MPI);

  //AZTEC matrix for A_ii=(D*C^{-1}*trD)
  AZ_MATRIX *A_ii;
  AZ_PRECOND *pILU_ii;

  int N_eq_ii= _dp.size();

  A_ii= AZ_matrix_create(N_eq_ii);
  // data containing the matrices C, D, trD and H as pointers
  // are passed through A_ii and pILU_ii:
  AZ_set_MATFREE(A_ii, &_factor_data_jacobian,
		 my_matvec< MixedMatr<1,3,CSRPatt,double>, MixedMatr<3,1,CSRPatt,double>,
		 vector<double>, MSRMatr<double>, Vector >);

  pILU_ii= AZ_precond_create(A_ii, my_precSchur_PC<
			     MSRMatr<double>,
			     MixedMatr<1,3,CSRPatt,double>,
			     MixedMatr<3,1,CSRPatt,double>,
			     vector<double>,
			     MSRMatr<double>,
			     Vector>, &_factor_data_jacobian);

  _dataAztec_ii.aztecOptionsFromDataFile(options_ii,params_ii);
   
  // user preconditioning:
  options_ii[AZ_precond]  = AZ_user_precond;

  // RHS of the linear system (ii)
  Vector vec_DV(_dp.size());

 
  //matrices HinvC (depends on time):
  MultInvDiag(_H, _C, _HinvC);
 
  
  // RHS of the linear system (ii)
  vec_DV = _D*_du.vec()-_f_p.vec();

   // case of pure Dirichlet BCs:
  if (BCh_du.fullEssential()) {
    vec_DV[_dim_p-1]    = 1.0; // correction of the right hand side.
  }
  
  _dp.vec()=0.0;

 
  cout << "  o-  Solving pressure system... \n";
  cout << "  maxnorm (vec_DV) = " <<  maxnorm(vec_DV) << endl; 
  cout << "  maxnorm (_f_p.vec()) = " <<  maxnorm(_f_p.vec())   << endl; 
  cout << "  maxnorm (_D*_du.vec() ) = " <<  maxnorm(_D*_du.vec()) << endl; 
  
  chrono.start(); 
  options_ii[AZ_recursion_level]=1;

  AZ_iterate(_dp.giveVec(), &vec_DV[0],  options_ii,params_ii, status_ii,
       proc_config_ii, A_ii, pILU_ii, NULL);

  chrono.stop();

  cout << "done in " << chrono.diff() << " s." << endl;
  //options_ii[AZ_recursion_level]=0;

  // ----------------------------
  // (iii) V = V-(C^(-1)*trD) * P
  // ----------------------------
  _du.vec() = _du.vec() - _invCtrDP.vec();
  cout << "  o-  Velocity updated" << endl; 

  AZ_matrix_destroy(&A_ii);
  AZ_precond_destroy(&pILU_ii);
  AZ_matrix_destroy(&C);
  AZ_precond_destroy(&prec_C);

  _residual_u.vec() = _f_duWithOutBC.vec() - _CAux * _du.vec() - _trDAux   * _dp.vec(); 

  cout << "  maxnorm (_residual_du ) = " <<  maxnorm( _residual_u.vec() ) << endl; 
}



template<typename Mesh>    
Vector& NavierStokesAleSolverPC<Mesh>::
residual() {
  return _residual_u.vec(); 
}


#endif
