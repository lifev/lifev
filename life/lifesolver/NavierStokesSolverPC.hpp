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
  \file NavierStokesSolverPC.h
  \author M.A. Fernandez and A. Gauthier
  \date 11/2002
  \version 1.0
  \author A. Veneziani
  \date 05/2003
  \version 1.1


  \brief This file contains a NavierStokes solver class which implements a semi-implicit scheme with an exact factorization.
         Preconditioning of the Schur Complement is done by an algebraic Chorin-Temam pressure-corrected preconditioner
         Added a class for handling high order discretization schemes (A. Veneziani)
*/

#ifndef _NAVIERSTOKESSOLVERPC_H_
#define _NAVIERSTOKESSOLVERPC_H_

#include "NavierStokesHandler.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include "values.hpp"
#include "pattern.hpp"
#include "assemb.hpp"
#include "bc_manage.hpp"
#include "algebraic_facto.hpp"
#include "bcCond.hpp"
#include "chrono.hpp"
#include "dataAztec.hpp"
#include "bdfNS.hpp"
#include "openDX_wrtrs.hpp"
#include <string>

namespace LifeV
{
/*!
  \class NavierStokesSolverPC

  This class implements an NavierStokes solver via exact factorization. Preconditioning of the
  Schur Complement is done by an algebraic Chorin-Temam pressure-corrected preconditioner

*/
template<typename Mesh>
class NavierStokesSolverPC:
        public NavierStokesHandler<Mesh> {

public:

    typedef  typename  NavierStokesHandler<Mesh>::Function Function;

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
      \param ord_bdf order of the Bdf time advancing scheme + incremental for the pressure
    */
    NavierStokesSolverPC(const GetPot& data_file, const RefFE& refFE_u, const RefFE& refFE_p, const QuadRule& Qr_u,
                         const QuadRule& bdQr_u, const QuadRule& Qr_p, const QuadRule& bdQr_p, BC_Handler& BCh_u);

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    void timeAdvance(const Function source, const Real& time);

    //! Update convective term, bc treatment and solve the linearized ns system
    void iterate(const Real& time);

    // ! Residual Computation
    PhysVectUnknown<Vector> residual();

    // ! Shear stress computation ***** Prova Agosto 2003
    void ShearStressCompute(std::string filename_sstress, std::string fe_type);


private:

    //! Block pattern of C: rho/dt*Vmass + mu*Vstiff operator
    MSRPatt _pattC_block;

    //! Pattern for C
    MixedPattern<3,3,MSRPatt> _pattC;

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

    //! Matrix HinvDtr:  H^{-1}D^T
    MixedMatr<3,1,CSRPatt,double> _HinvDtr;

    //! Matrix M_u: Vmass
    MixedMatr<3,3,MSRPatt,double> _M_u;

    //! Matrix HinvC: H^{-1}C
    MixedMatr<3,3,MSRPatt,double> _HinvC;

    //! Matrix C: rho/dt*Vmass + mu*Vstiff operator
    MixedMatr<3,3,MSRPatt,double> _CStokes;

    //! Matrix C: rho/dt*Vmass + mu*Vstiff operator + Convective_term
    MixedMatr<3,3,MSRPatt,double> _C;


    //! H diag matrix: H= diag( _M_u )/sum( diag( _M_u ) ) where _M_u = mass * rho / dt
    std::vector<double>  _H;

    //! Elementary matrices and vectors
    ElemMat _elmatC; //velocity stiffnes
    ElemMat _elmatM_u; //velocity mass
    ElemMat _elmatDtr; // vel_i * pres_j
    ElemVec _elvec; // Elementary right hand side

    //! Right  hand  side for the velocity
    PhysVectUnknown<Vector> _f_u;

    //!  This vector contains the product C^{-1}*trD*P where P is the pressure, solution
    //!  of the system (ii).
    PhysVectUnknown<Vector> _invCtrDP;

    // A Veneziani August 2003 ****************
    //! Matrix Cnobc: rho/dt*Vmass + mu*Vstiff operator + Convective_term WITHOUT BC (for the residual computation)
    MixedMatr<3,3,MSRPatt,double> _CnoBc;

    //! Matrix trD: transpose Vdiv operator trVdiv WITHOUT BC (for the residual computation)
    MixedMatr<3,1,CSRPatt,double> _trDnoBc;

    //! Right  hand  side for the velocity WITHOUT BC
    PhysVectUnknown<Vector> _f_u_noBc;
    // REM This solution is just to start: Miguel suggested a different way,
    // dealing only with the submatrices involved in the Dirichlet elimination.
    // This is a good idea....still to be done


    DataAztec _dataAztec_i;
    DataAztec _dataAztec_ii;
    DataAztec _dataAztec_s;

    //! DataFactorisation: data passed to matrix-vector product are stored in the class
    DataFactorisation<
        MixedMatr<3,3,MSRPatt,double>,
        MixedMatr<1,3,CSRPatt,double>,
        MixedMatr<3,1,CSRPatt,double>,
        std::vector<double>,
        MSRMatr<double>,
        Vector> _factor_data;
};


//
//                                         IMPLEMENTATION
//
template<typename Mesh> NavierStokesSolverPC<Mesh>::
NavierStokesSolverPC(const GetPot& data_file, const RefFE& refFE_u, const RefFE& refFE_p, const QuadRule& Qr_u,
                     const QuadRule& bdQr_u, const QuadRule& Qr_p, const QuadRule& bdQr_p, BC_Handler& BCh_u):
    NavierStokesHandler<Mesh>(data_file,refFE_u,refFE_p,Qr_u,bdQr_u,Qr_p,bdQr_p,BCh_u),
     _pattC_block(this->_dof_u),
     _pattC(_pattC_block,"diag"),
     _pattD_block(this->_dof_p,this->_dof_u),
     _pattD(_pattD_block),
     _pattDtr_block(this->_dof_u,this->_dof_p),
     _pattDtr(_pattDtr_block),
     _D(_pattD),
     _trD(_pattDtr),
     _HinvDtr(_pattDtr),
     _M_u(_pattC),
     _HinvC(_pattC),
     _CStokes(_pattC),
     _C(_pattC),
     _H(_pattC.nRows()),
     _elmatC(this->_fe_u.nbNode,nDimensions,nDimensions),
     _elmatM_u(this->_fe_u.nbNode,nDimensions,nDimensions),
     _elmatDtr(this->_fe_u.nbNode,nDimensions,0,this->_fe_p.nbNode,0,1),
     _elvec(this->_fe_u.nbNode,nDimensions),
     _f_u(this->_dim_u),
     _invCtrDP(this->_dim_u),
      _CnoBc(_pattC),_trDnoBc(_pattDtr),_f_u_noBc(this->_dim_u),
     _dataAztec_i(data_file,"fluid/aztec_i"),
     _dataAztec_ii(data_file,"fluid/aztec_ii"),
     _dataAztec_s(data_file,"fluid/aztec_s"),
     _factor_data(_C,_D,_trD,_H,_HinvC,_HinvDtr,_invCtrDP,_dataAztec_i,_dataAztec_s,this->_BCh_u.fullEssential()) {

  std::cout << std::endl;
  std::cout << "O-  Pressure unknowns: " << this->_dim_p     << std::endl;
  std::cout << "O-  Velocity unknowns: " << this->_dim_u     << std::endl<<std::endl;
  std::cout << "O-  Computing mass and Stokes matrices... ";

  Chrono chrono;
  chrono.start();

  // Matrices initialization
  _D.zeros();
  _trD.zeros();
  _M_u.zeros();
  _CStokes.zeros();
  _C.zeros();

  _CnoBc.zeros();
  _trDnoBc.zeros();


  // Number of velocity components
  UInt nc_u=this->_u.nbcomp();

  //inverse of dt:
  Real dti=1./this->_dt;


  // *******************************************************
  // Coefficient of the mass term at time t^{n+1}
  Real first_coeff = this->_bdf.bdf_u().coeff_der(0);
  std::cout << std::endl;
  std::cout << "Bdf NS first coeff " << first_coeff << std::endl;

  this->_bdf.bdf_u().showMe();
  this->_bdf.bdf_p().showMe();

  // Auxiliary matrix
  ElemMat elmatM_u_St(this->_fe_u.nbNode,nDimensions,nDimensions);

  // Elementary computation and matrix assembling
  // Loop on elements
  for(UInt i = 1; i <= this->_mesh.numVolumes(); i++){

    this->_fe_p.update( this->_mesh.volumeList(i) ); // just to provide the id number in the assem_mat_mixed
    this->_fe_u.updateFirstDerivQuadPt(this->_mesh.volumeList(i));

    _elmatC.zero();
    _elmatM_u.zero();
     elmatM_u_St.zero();
    _elmatDtr.zero();

    stiff(this->_mu,_elmatC,this->_fe_u,0,0,nDimensions);
  // *******************************************************
    mass(first_coeff*this->_rho*dti,elmatM_u_St,this->_fe_u,0,0,nDimensions);
    mass(this->_rho*dti,_elmatM_u,this->_fe_u,0,0,nDimensions);
    // stiffness + mass

    _elmatC.mat() += elmatM_u_St.mat();

    for(UInt ic=0;ic<nc_u;ic++){

      // stiffness
      assemb_mat(_CStokes,_elmatC,this->_fe_u,this->_dof_u,ic,ic);

      // mass
      assemb_mat(_M_u,_elmatM_u,this->_fe_u,this->_dof_u,ic,ic);

      // div
      grad(ic,1.0,_elmatDtr,this->_fe_u,this->_fe_p,ic,0);

      // assembling
      assemb_mat_mixed(_trD,_elmatDtr,this->_fe_u,this->_fe_p,this->_dof_u,this->_dof_p,ic,0);
      assemb_tr_mat_mixed(1.0,_D,_elmatDtr,this->_fe_p,this->_fe_u,this->_dof_p,this->_dof_u,0,ic);
    }
  }

  // H diag matrix: H= diag( _M_u )/sum( diag( _M_u ) ) where _M_u = mass * rho / dt * first_coeff_bdf_scheme
  _H=_M_u.giveDiag();
  double sum=accumulate(_H.begin(),_H.end(),0.0);
  // *******************************************************
  // Matrix equilibration:
  for (UInt i=0; i<_H.size();i++)
    _H[i]=_H[i]*first_coeff*dti/sum; // H = lumping of first_coeff/dt*Mass

  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;
}

template<typename Mesh>
void NavierStokesSolverPC<Mesh>::
timeAdvance(const Function source, const Real& time) {

  std::cout << std::endl;
  std::cout << "O== Now we are at time "<< time << " s." << std::endl;

  // Number of velocity components
  UInt nc_u=this->_u.nbcomp();

  std::cout << "  o-  Updating mass term on right hand side... ";


  Chrono chrono;
  chrono.start();

  // Right hand side for the velocity at time
  _f_u=0.;

  // loop on volumes: assembling source term
  for(UInt i=1; i<=this->_mesh.numVolumes(); ++i){
    _elvec.zero();
    this->_fe_u.update(this->_mesh.volumeList(i));

    for (UInt ic=0; ic<nc_u; ++ic){
      compute_vec(source,_elvec,this->_fe_u,time,ic); // compute local vector
      assemb_vec(_f_u,_elvec,this->_fe_u,this->_dof_u,ic); // assemble local vector into global one
    }
  }

  // *******************************************************
  _f_u += _M_u*this->_bdf.bdf_u().time_der(); //_M_u is the mass matrix divided by the time step
  //  _f_u += _M_u * _u;
  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;
}


template<typename Mesh>
void NavierStokesSolverPC<Mesh>::
iterate(const Real& time) {

  // Number of velocity components
  UInt nc_u=this->_u.nbcomp();
  Vector u_extrap = this->_bdf.bdf_u().extrap();

  Chrono  chrono;

  // C = CStokes + convective term
  chrono.start();
  _C=_CStokes;
  chrono.stop();


  std::cout << "  o-  Stokes matrix was copied in " << chrono.diff() << "s." << std::endl;
  std::cout << "  o-  Updating convective term... ";

  chrono.start();

  // loop on volumes
  for(UInt i=1; i<=this->_mesh.numVolumes(); ++i){

    this->_fe_u.updateFirstDeriv(this->_mesh.volumeList(i)); // as updateFirstDer

    _elmatC.zero();

    UInt eleID = this->_fe_u.currentId();
    // Non linear term, Semi-implicit approach
    // ULoc contains the velocity values in the nodes
    for (UInt k=0 ; k<(UInt)this->_fe_u.nbNode ; k++){
	UInt  iloc = this->_fe_u.patternFirst(k);
	for (UInt ic=0; ic<nc_u; ++ic){
	  UInt ig=this->_dof_u.localToGlobal(eleID,iloc+1)-1+ic*this->_dim_u;
	  _elvec[iloc+ic*this->_fe_u.nbNode] = u_extrap(ig);
	}
    }

    // loop on components
    for (UInt ic=0; ic<nc_u; ++ic){
      // compute local convective term and assembling
      grad(0,_elvec,_elmatC,this->_fe_u,this->_fe_u,ic,ic);
      grad(1,_elvec,_elmatC,this->_fe_u,this->_fe_u,ic,ic);
      grad(2,_elvec,_elmatC,this->_fe_u,this->_fe_u,ic,ic);
      assemb_mat(_C,_elmatC,this->_fe_u,this->_dof_u,ic,ic);
    }
  }
  chrono.stop();
  std::cout << "done in " << chrono.diff() << "s." << std::endl;

  // QUI VENGONO APPLICATE LE BC
  _CnoBc=_C;
  _trDnoBc=_trD;
  _f_u_noBc=_f_u;
  //for (UInt myindex=0;myindex<_dim_u;myindex++) _f_u_noBc[myindex] = _f_u[myindex];

  _f_u -= _trD*this->_bdf.bdf_p().extrap(); // INCREMENTAL correction for the pressure: IT MUST BE AFTER THE INITIALIZATION of _f_u_noBC

  // for BC treatment (done at each time-step)
  Real tgv=1.e02;

  std::cout << "  o-  Applying boundary conditions... ";
  chrono.start();
  // BC manage for the velocity
  if ( !this->_BCh_u.bdUpdateDone() )
      this->_BCh_u.bdUpdate(this->_mesh, this->_feBd_u, this->_dof_u);
  bc_manage(_C, _trD, _f_u, this->_mesh, this->_dof_u, this->_BCh_u, this->_feBd_u, tgv, time);
  chrono.stop();
  std::cout << "done in " << chrono.diff() << "s." << std::endl;

  //matrices HinvDtr:
  MultInvDiag(_H, _trD, _HinvDtr);

  // AZTEC specifications for each system
  int    data_org_i1[AZ_COMM_SIZE];   // data organisation for C1
  int    data_org_i2[AZ_COMM_SIZE];   // data organisation for C2
  int    data_org_i3[AZ_COMM_SIZE];   // data organisation for C3

  // identical for each system
  int    proc_config_i[AZ_PROC_SIZE];// Processor information:
  int    options_i[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_i[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_i[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.
  AZ_set_proc_config(proc_config_i, AZ_NOT_MPI);

  //AZTEC matrix and preconditioner
  AZ_MATRIX *C1, *C2, *C3;
  AZ_PRECOND *prec_C1, *prec_C2, *prec_C3;

  int N_eq_i= this->_dim_u; // number of DOF for each component

  // set each block
  C1= AZ_matrix_create(N_eq_i);
  C2= AZ_matrix_create(N_eq_i);
  C3= AZ_matrix_create(N_eq_i);

  // create preconditioner for each block
  prec_C1= AZ_precond_create(C1, AZ_precondition, NULL);
  prec_C2= AZ_precond_create(C2, AZ_precondition, NULL);
  prec_C3= AZ_precond_create(C3, AZ_precondition, NULL);

  // data_org assigned "by hands" while no parallel computation is performed
  data_org_i1[AZ_N_internal]= N_eq_i;
  data_org_i1[AZ_N_border]= 0;
  data_org_i1[AZ_N_external]= 0;
  data_org_i1[AZ_N_neigh]= 0;
  data_org_i1[AZ_name]= DATA_NAME_AZTEC1;

  data_org_i2[AZ_N_internal]= N_eq_i;
  data_org_i2[AZ_N_border]= 0;
  data_org_i2[AZ_N_external]= 0;
  data_org_i2[AZ_N_neigh]= 0;
  data_org_i2[AZ_name]= DATA_NAME_AZTEC2;

  data_org_i3[AZ_N_internal]= N_eq_i;
  data_org_i3[AZ_N_border]= 0;
  data_org_i3[AZ_N_external]= 0;
  data_org_i3[AZ_N_neigh]= 0;
  data_org_i3[AZ_name]= DATA_NAME_AZTEC3;

  // --------------------------------------
  // (ii) (D*C^{-1}*trD) * \delta P =  D*V
  // --------------------------------------

  // AZTEC specifications for the second system
  int    proc_config_ii[AZ_PROC_SIZE];// Processor information:
  //  proc_config[AZ_node] = node name
  //  proc_config[AZ_N_procs] = # of nodes
  int    options_ii[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_ii[AZ_PARAMS_SIZE];   // User selected solver paramters.
  double status_ii[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                     // indicating success or failure.
  //

  AZ_set_proc_config(proc_config_ii, AZ_NOT_MPI);

  //AZTEC matrix for A_ii=(D*C^{-1}*trD)
  AZ_MATRIX *A_ii;
  AZ_PRECOND *pILU_ii;

  int N_eq_ii= this->_p.size();

  A_ii= AZ_matrix_create(N_eq_ii);
  // data containing the matrices C, D, trD and H as pointers
  // are passed through A_ii and pILU_ii:
  AZ_set_MATFREE(A_ii, &_factor_data,
		 my_matvec_block<
		 MixedMatr<3,3,MSRPatt,double>,
		 MixedMatr<1,3,CSRPatt,double>,
		 MixedMatr<3,1,CSRPatt,double>,
		 std::vector<double>,
		 MSRMatr<double>,
		 Vector>);

  pILU_ii= AZ_precond_create(A_ii, my_precSchur_PC<
			     MixedMatr<3,3,MSRPatt,double>,
			     MixedMatr<1,3,CSRPatt,double>,
			     MixedMatr<3,1,CSRPatt,double>,
			     std::vector<double>,
			     MSRMatr<double>,
			     Vector>, &_factor_data);


  _dataAztec_ii.aztecOptionsFromDataFile(options_ii,params_ii);

  // user preconditioning:
  options_ii[AZ_precond]  = AZ_user_precond;

  // RHS of the linear system (ii)
  Vector vec_DV( this->_p.size() );

  //matrices HinvC (depends on time):
  MultInvDiag(_H, _C, _HinvC);

  // ---------------
  // (i) C * V = F_V
  // ---------------
  std::cout << "  o-  Solving first system... ";

  AZ_set_MSR(C1, (int*) _pattC_block.giveRaw_bindx(),
	     (double*) _C.giveRaw_value(0,0),
	     data_org_i1, 0, NULL, AZ_LOCAL);
  AZ_set_MSR(C2, (int*) _pattC_block.giveRaw_bindx(),
	     (double*) _C.giveRaw_value(1,1),
	     data_org_i2, 0, NULL, AZ_LOCAL);
  AZ_set_MSR(C3, (int*) _pattC_block.giveRaw_bindx(),
	     (double*) _C.giveRaw_value(2,2),
	     data_org_i3, 0, NULL, AZ_LOCAL);

  _dataAztec_i.aztecOptionsFromDataFile(options_i,params_i);

  //keep C factorisation and preconditioner reused in my_matvec
  options_i[AZ_keep_info]= 1;               // keep information

  // ---------------
  // (i) C * V = F_V (= forcing term + bc - D^T*P^n)
  // ---------------
  // intermediate velocity computation

  chrono.start();

  //for each block
  AZ_iterate( this->_u.giveVec(), _f_u.giveVec(), options_i,params_i, status_i,
             proc_config_i, C1, prec_C1, NULL);
  AZ_iterate( this->_u.giveVec()+this->_dim_u, _f_u.giveVec()+this->_dim_u,  options_i,params_i,
              status_i, proc_config_i, C2, prec_C2, NULL);
  AZ_iterate( this->_u.giveVec()+2*this->_dim_u, _f_u.giveVec()+2*this->_dim_u,  options_i,params_i,
              status_ii, proc_config_i, C3, prec_C3, NULL);

  //
  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;

  // ---------------------------------------------------
  // (ii) (D*C^(-1)*trD) * \delta P = D*C^{-1}*F_V = D*V
  // ---------------------------------------------------

  // RHS of the linear system (ii)
  vec_DV = _D*this->_u;
  this->_p = 0.0; // AT this point, this vector stands for the "pressure increment"

  // case of pure Dirichlet BCs:
  if ( this->_BCh_u.fullEssential() )
  {
      vec_DV[this->_dim_p-1]   = 0.0; // correction of the right hand side.
      this->_p[this->_dim_p-1] = 0.0; // pressure value at the last node.
  }

  std::cout << "  o-  Solving second system... ";

  chrono.start();
  AZ_iterate( this->_p.giveVec(), &vec_DV[0],  options_ii,params_ii, status_ii,
             proc_config_ii, A_ii, pILU_ii, NULL);

  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;

  // ------------------------------------
  // (iii) V = V-(C^(-1)*trD) * \delta P
  // ------------------------------------


  // everything is done...
  this->_u = this->_u - _invCtrDP;
  std::cout << "  o-  Velocity updated" << std::endl;

  // *******************************************************
  // This is the REAL pressure (not the increment)
  this->_p += this->_bdf.bdf_p().extrap();
  std::cout << "  o-  Pressure updated" << std::endl;

  // *******************************************************
  // update the array of the previous solutions
  this->_bdf.bdf_u().shift_right(this->_u);
  this->_bdf.bdf_p().shift_right(this->_p);


  // destroy Ci and prec_Ci
  AZ_matrix_destroy(&C1);
  AZ_precond_destroy(&prec_C1);
  AZ_matrix_destroy(&C2);
  AZ_precond_destroy(&prec_C2);
  AZ_matrix_destroy(&C3);
  AZ_precond_destroy(&prec_C3);
  // destroy A_ii and pILU_ii
  AZ_matrix_destroy(&A_ii);
  AZ_precond_destroy(&pILU_ii);
}


////////////////

template<typename Mesh>
PhysVectUnknown<Vector> NavierStokesSolverPC<Mesh>::residual()
 {
  Chrono chrono;
  PhysVectUnknown<Vector> r(this->_dim_u);
  std::cout << "  o- Computing the residual...";
  chrono.start();
  r = _f_u_noBc-_CnoBc*this->_u - _trDnoBc*this->_p;
  chrono.stop();
  std::cout << "done in " << chrono.diff() << "s" << std::endl;
  return r;
  }

template<typename Mesh>
void NavierStokesSolverPC<Mesh>::ShearStressCompute(std::string filename_stress, std::string fe_type)
{
  Vector residual(this->residual());
  UInt ss = residual.size()/NDIM;
  Vector sstress(this->_ns_post_proc.compute_sstress(residual,(UInt)NDIM));
  // just a stupid way for writing the shear stress in OpenDx or Medit formats,
  // exploting the existent subroutines
  UInt s=sstress.size()/NDIM;
  UInt where;
  for (UInt i=0;i<s;i++){
    where=this->_ns_post_proc.fBdToIn()[i];
    //    std::cout << residual.size() << " " << s << " " << where << " " << i << std::endl;
   for (UInt j=0;j<NDIM;j++) residual[where-1+j*ss]=sstress[i+j*s];
  }
  wr_opendx_header(filename_stress, this->_mesh, this->_dof_u,this->_fe_u,fe_type);
  wr_opendx_vector(filename_stress,"ShearStress",residual,NDIM);
}

}

#endif
