/*!
  \file convDiffReactSolverPC.h
  \author M. Prosi
  \date 03/2004
  \version 1.0

  \brief This file contains a solver class for the Convection-Diffusion-Reaction equation 
*/

#ifndef _CONVDIFFREACTSOLVERPC_H_
#define _CONVDIFFREACTSOLVERPC_H_

#include "convDiffReactHandler.hpp"
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
#include "bdf.hpp"
#include "openDX_wrtrs.hpp"
#include <string>

/*! 
  \class convDiffReactSolverPC

   This class implements a solver for the Convection-Diffusion-Reaction equation

*/
template<typename Mesh>
class ConvDiffReactSolverPC:
public ConvDiffReactHandler<Mesh> {
 
 public:
  
  typedef  typename  ConvDiffReactHandler<Mesh>::Function Function; 

  //! Constructor 
  /*!
    \param data_file GetPot data file
    \param refFE_c reference FE for the concentration
    \param Qr_c volumic quadrature rule for the concentration
    \param bdQr_c surface quadrature rule for the concentration
    \param BCh_c boundary conditions for the concentration
  */
  ConvDiffReactSolverPC(const GetPot& data_file, const RefFE& refFE_c, const QuadRule& Qr_c,
	    const QuadRule& bdQr_c, BC_Handler& BCh_c);

  //! Update the right  hand side  for time advancing 
  /*! 
    \param source volumic source  
    \param time present time
  */
  void timeAdvance(const Function source, const Real& time);

  //! Update convective term, bc treatment and solve the linearized ns system
  void iterate(const Real& time);

  //! Projection of the velocity on grid of concentration discretization
  template <typename RegionMesh3D>
  void getvel(RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BC_Handler& BCh_u, const Real& time);

 private:

  //! Pattern of M
  MSRPatt _pattM;
  
  //! Matrix C:  1/dt*Cmass + D*Cstiff operator (+ r*Cmass ... reaction Term TODO !!!)
  MSRMatr<double> _DR;
 
  //! Matrix C:  1/dt*Cmass + D*Cstiff operator + Convective_transport term (+ r*Cmass ... reaction Term TODO !!!)
  MSRMatr<double> _CDR;

  //! Matrix C_u: Cmass
  MSRMatr<double> _M_c;
 
  //! Elementary matrices and vectors
  ElemMat _elmatC; //Concentration stiffnes 
  ElemMat _elmatM_c; //Concentration mass
  ElemVec _elvec; // Elementary right hand side
  ElemVec _elvec_u; // Elementary velocity for convection term

  //! Right  hand  side for the concentration
  ScalUnknown<Vector> _f_c;

  //! velocity vector on the concentration nodes
  PhysVectUnknown<Vector> _u_c;

  DataAztec _dataAztec_o;
};

//
//                                         IMPLEMENTATION
//
template<typename Mesh> ConvDiffReactSolverPC<Mesh>::
ConvDiffReactSolverPC(const GetPot& data_file, const RefFE& refFE_c, const QuadRule& Qr_c,
		      const QuadRule& bdQr_c, BC_Handler& BCh_c):  
  ConvDiffReactHandler<Mesh>(data_file,refFE_c,Qr_c, bdQr_c, BCh_c),
     _pattM(_dof_c),
     _DR(_pattM),
     _CDR(_pattM),
     _M_c(_pattM),
     _elmatC(_fe_c.nbNode,1,1), 
     _elmatM_c(_fe_c.nbNode,1,1),
     _elvec(_fe_c.nbNode,1), 
     _elvec_u(_fe_c.nbNode,nDimensions),
     _f_c(_dim_c),
     _u_c(_dim_c),
     _dataAztec_o(data_file,"masstransport/aztec_o"){
  
  cout << endl;
  cout << "O-  Concentration unknowns: " << _dim_c     << endl; 
  cout << "O-  Computing mass and stiffness matrices... ";  
    
  Chrono chrono;
  chrono.start();

  // Matrices initialization 
  _DR.zeros();
  _CDR.zeros();
  _M_c.zeros();
  
  //inverse of dt:
  Real dti=1./_dt;

  // *******************************************************
  // Coefficient of the mass term at time t^{n+1}
  Real first_coeff = _bdf.coeff_der(0);
  cout << endl;
  cout << "Bdf CDR first coeff " << first_coeff << endl; 

  _bdf.showMe();

  // Elementary computation and matrix assembling  

  for(UInt i = 1; i <= _mesh.numVolumes(); i++){          // Loop on elements

    _fe_c.updateFirstDerivQuadPt(_mesh.volumeList(i));
    
    _elmatC.zero();
    _elmatM_c.zero();
 
    stiff(_diffusivity,_elmatC,_fe_c);
//    _elmatC.showMe();
    mass(first_coeff*dti,_elmatM_c,_fe_c);
//    _elmatM_c.showMe();

    // stiffness + mass
    _elmatC.mat() += _elmatM_c.mat();
    
    
    // stiffness
    assemb_mat(_DR,_elmatC,_fe_c,_dof_c);
      
    // mass
    assemb_mat(_M_c,_elmatM_c,_fe_c,_dof_c);
     
  }
   
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;

}

template<typename Mesh>  
void ConvDiffReactSolverPC<Mesh>::
timeAdvance(const Function source, const Real& time) {

  cout << "  o-  Updating mass term on right hand side (concentration)... ";

  Chrono chrono;
  chrono.start();

  // Right hand side for the velocity at time
  _f_c=0.;

  // loop on volumes: assembling source term
  for(UInt i=1; i<=_mesh.numVolumes(); ++i){
     _elvec.zero();
     _fe_c.update(_mesh.volumeList(i));

      compute_vec(source,_elvec,_fe_c,time,0); // compute local vector
      assemb_vec(_f_c,_elvec,_fe_c,_dof_c,0); // assemble local vector into global one       
  }

  // ******************************************************* 
  _f_c += _M_c*_bdf.time_der(); //_M_u is the mass matrix divided by the time step
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
}


template<typename Mesh>  
void ConvDiffReactSolverPC<Mesh>::
iterate(const Real& time) {

  UInt nc_u=_u_c.nbcomp();

  Chrono  chrono;

  // CDR = DR + convective term (C)
  chrono.start();
  _CDR=_DR;
  chrono.stop();


  cout << "  o-  Diffusion-Reaction matrix was copied in " << chrono.diff() << "s." << endl;
  cout << "  o-  Updating convective transport... ";
 
  chrono.start();

  // loop on volumes
  for(UInt i=1; i<=_mesh.numVolumes(); ++i){
  
    _fe_c.updateFirstDeriv(_mesh.volumeList(i)); // as updateFirstDer

    _elmatC.zero();

    UInt eleID = _fe_c.currentId();

// ********** copy global velocity vector to local velocity vector *************
// ********** assuming velocity is given on concentration mesh *****************

    for (UInt k=0 ; k<(UInt)_fe_c.nbNode ; k++){
       UInt  iloc = _fe_c.patternFirst(k);
       for (UInt ic=0; ic<nc_u; ++ic){
	  UInt ig=_dof_c.localToGlobal(eleID,iloc+1)-1+ic*_dim_c;     
	  _elvec_u[iloc+ic*_fe_c.nbNode] = _u_c(ig);
       }
     }

    grad(0,_elvec_u,_elmatC,_fe_c,_fe_c,_fe_c);
    grad(1,_elvec_u,_elmatC,_fe_c,_fe_c,_fe_c);
    grad(2,_elvec_u,_elmatC,_fe_c,_fe_c,_fe_c);


// *************************************** Upwind ******************************

      Real VLoc_infty=0.;
      Real VLoc_mean=0.;
      Real VLoc_c=0.;
      for (UInt ih_c=0 ; ih_c<(UInt)_fe_c.nbNode ; ih_c++){
         UInt  iloc = _fe_c.patternFirst(ih_c);
	 for (UInt ic=0; ic<nc_u;++ic){
	   UInt ig=_dof_c.localToGlobal(eleID,iloc+1)-1+ic*_dim_c;
           _elvec_u[iloc+ic*_fe_c.nbNode] = _u_c(ig);
	   VLoc_c+=_u_c(ig)*_u_c(ig);}
	 VLoc_c=sqrt(VLoc_c);
	 VLoc_mean += VLoc_c;
	if (VLoc_c>VLoc_infty) VLoc_infty=VLoc_c;
      }
      VLoc_mean=VLoc_mean/_fe_c.nbNode;

      Real coef_stab, Pe_loc;
//      coef_stab=_fe_c.diameter()*VLoc_infty; // Alessandro - method

      Pe_loc=VLoc_infty*_fe_c.diameter()/(2.0*_diffusivity);

//      coef_stab=(1.0/tanh(Pe_loc))-(1.0/Pe_loc); // classical approach

      if(Pe_loc < -3.0)
	 coef_stab= -1.0;
      else {
	 if(Pe_loc > 3.0)
	    coef_stab=1.0;
	 else	
	    coef_stab=Pe_loc/3.0;}

// ******************************* STREAMLINEUPWIND ****************************
     stiff_sd(coef_stab/(VLoc_mean*VLoc_mean),_elvec_u,_elmatC,_fe_c,_fe_c);

// ************************* Assembling ****************************************

     assemb_mat(_CDR,_elmatC,_fe_c,_dof_c);

  }
  chrono.stop();
  cout << "done in " << chrono.diff() << "s." << endl;

  // for BC treatment (done at each time-step)
  Real tgv=1.e02; 

  cout << "  o-  Applying boundary conditions... ";
  chrono.start(); 
  // BC manage for the concentration
  if ( !_BCh_c.bdUpdateDone() )  
    _BCh_c.bdUpdate(_mesh, _feBd_c, _dof_c);
  bc_manage(_CDR, _f_c, _mesh, _dof_c, _BCh_c, _feBd_c, tgv, time);
  chrono.stop();

  cout << "done in " << chrono.diff() << "s." << endl;


  int    proc_config_o[AZ_PROC_SIZE];// Processor information:                 
  //  proc_config[AZ_node] = node name      
  //  proc_config[AZ_N_procs] = # of nodes  
  int    options_o[AZ_OPTIONS_SIZE]; // Array used to select solver options.     
  double params_o[AZ_PARAMS_SIZE];   // User selected solver paramters.          
  int    *data_org_o;                // Array to specify data layout   
  double status_o[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                   // indicating success or failure.           
  // altre dichiarazioni per AZTEC  
  int    *update_o,                  // vector elements updated on this node. 
         *external_o;                // vector elements needed by this node.    
  int    *update_index_o;            // ordering of update[] and external[]     
  int    *extern_index_o;            // locally on this processor.              
  //  int    *bindx;                 // Sparse matrix to be solved is stored    
  //  double *val;                   // in these MSR arrays.                    
  int    N_update_o;                 // # of unknowns updated on this node      

  AZ_set_proc_config(proc_config_o, AZ_NOT_MPI );

    AZ_read_update(&N_update_o, &update_o, proc_config_o, _dim_c, 1, AZ_linear);
    AZ_defaults(options_o,params_o);
    _dataAztec_o.aztecOptionsFromDataFile(options_o,params_o);
    AZ_transform(proc_config_o, &external_o, 
	       (int *)_pattM.giveRaw_bindx(), _CDR.giveRaw_value(), 
	       update_o, &update_index_o,
	       &extern_index_o, &data_org_o, N_update_o, NULL, NULL, NULL, NULL,
	       AZ_MSR_MATRIX);
  
    chrono.start();
//    init_options_c(options_o,params_o);
  
    AZ_solve(_c.giveVec(),_f_c.giveVec(), options_o, params_o, NULL, 
	   (int *)_pattM.giveRaw_bindx(), NULL, NULL, NULL, 
	   _CDR.giveRaw_value(), data_org_o,status_o, proc_config_o);
  //
    chrono.stop();
    cout << "*** Solution (Concentration) computed in " << chrono.diff() << "s." << endl;
  _bdf.shift_right(_c);
  
}

template<typename Mesh> template<typename RegionMesh3D>  
void ConvDiffReactSolverPC<Mesh>::
getvel(RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BC_Handler& BCh_u, const Real& time){

   for (UInt j=0; j<3; j++){
       for(UInt i=0; i< _dim_c; i++){
	  _u_c(i+j*_dim_c)=u(i+j*u.size()/3);
       }}


}


#endif
