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
/*!
  \file oneDModelSolver.cpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file implements a solver for 1D model.
*/

#include "oneDModelSolver.hpp"



OneDModelSolver::OneDModelSolver(const GetPot& data_file):  
  OneDModelHandler(data_file),
  _M_elmatMass (_M_fe.nbNode,1,1), 
  _M_elmatStiff(_M_fe.nbNode,1,1), 
  _M_elmatGrad (_M_fe.nbNode,1,1), 
  _M_elmatDiv  (_M_fe.nbNode,1,1), 
  _M_elvec     (_M_fe.nbNode,1), 
  _M_rhs(_M_dimDof),
  _M_massMatrix(_M_dimDof)
{
  
  cout << endl;
  cout << "O-  Nb of unknowns: " << _M_dimDof     << endl; 
  cout << "O-  Computing mass matrix... \n";  
    
  Chrono chrono;
  chrono.start();

  //! Matrices initialization 
  _M_massMatrix.zero();
  _M_massMatrix.showMe(std::cout, _M_verbose);

  //inverse of the time step:
  double dti=1./_M_time_step;

  _M_coeffMass  = 1.;
  _M_coeffStiff = 1.;
  _M_coeffGrad  = 1.;
  _M_coeffDiv   = 1.;

  //! Elementary computation and matrix assembling  
  //! Loop on elements
  for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){          

    //! update _M_elmat*
    _updateElemMatrices( iedge );


    // stiffness + mass
    //    _elmatC.mat() += _elmatM_c.mat();
    
    
    // stiffness
    // assemb_mat(_DR,_elmatC,_M_fe,_dof_c);
      
    // mass
    //assemb_mat(_M_c,_elmatM_c,_M_fe,_dof_c);
     
  }
   
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;

}

//! Update the element matrices with the current element (from 1)
void OneDModelSolver::_updateElemMatrices( const UInt& iedge )
{
  //! set the elementary matrices to 0.
  _M_elmatMass.zero();
  _M_elmatStiff.zero();
  _M_elmatGrad.zero();
  _M_elmatDiv.zero();

  //! update the current element 
  _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge)); 
  std::cout << _M_fe.currentId() << std::endl;

  //! update the mass matrix
  mass( _M_coeffMass, _M_elmatMass, _M_fe,0, 0 );
  std::cout << "Elem Mass matrix :" << std::endl;
  _M_elmatMass.showMe( std::cout );

  //! update the stiffness matrix
  stiff( _M_coeffStiff, _M_elmatStiff, _M_fe,0 ,0 );
  std::cout << "Elem Stiff matrix :" << std::endl;
  _M_elmatStiff.showMe( std::cout );

  /*! update the gradient matrix 
      gradient operator: 
      grad_{ij} = \int_{fe} coeff \phi_j \frac{d \phi_i}{d x}
      
      BEWARE :
      \param 0: the first argument "0" corresponds to the first
      and only coordinate (1D!), and HERE it starts from 0... (Damm'!)

      \param - _M_coeffGrad: the sign "-" in the second argument 
      is added to correspond to the described operator. 
      (There is a minus in the elemOper implementation).
  */
  grad( 0 , - _M_coeffGrad, _M_elmatGrad, _M_fe, _M_fe, 0, 0 );
  std::cout << "Elem Grad matrix :" << std::endl;
  _M_elmatGrad.showMe( std::cout );

  /*! update the divergence matrix 
      divergence operator: (transpose of the gradient) 
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this _M_elmatDiv is not necessary
      as it is the transpose of the _M_elmatGrad.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
  */
  div( 0 , - _M_coeffDiv, _M_elmatDiv, _M_fe, _M_fe, 0, 0 );
  std::cout << "Elem Div matrix :" << std::endl;
  _M_elmatDiv.showMe( std::cout );
}


/*
template<typename Mesh>  
void OneDModelSolver<Mesh>::
timeAdvance(const Function source, const double& time) {

  cout << "  o-  Updating mass term on right hand side (concentration)... ";

  Chrono chrono;
  chrono.start();

  // Right hand side for the velocity at time
  _f_c=0.;

  // loop on volumes: assembling source term
  for(UInt i=1; i<=_mesh.numVolumes(); ++i){
     _elvec.zero();
     _M_fe.update(_mesh.volumeList(i));

      compute_vec(source,_elvec,_M_fe,time,0); // compute local vector
      assemb_vec(_f_c,_elvec,_M_fe,_dof_c,0); // assemble local vector into global one       
  }

  // ******************************************************* 
  _f_c += _M_c*_bdf.time_der(); //_M_u is the mass matrix divided by the time step
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
}


template<typename Mesh>  
void OneDModelSolver<Mesh>::
iterate(const double& time, PhysVectUnknown<Vector> & u) {

  UInt nc_u=u.nbcomp();

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
  
    _M_fe.updateFirstDeriv(_mesh.volumeList(i)); // as updateFirstDer

    _elmatC.zero();

    UInt eleID = _M_fe.currentId();

// ********** copy global velocity vector to local velocity vector *************
// ********** assuming velocity is given on concentration mesh *****************

    for (UInt k=0 ; k<(UInt)_M_fe.nbNode ; k++){
       UInt  iloc = _M_fe.patternFirst(k);
       for (UInt ic=0; ic<nc_u; ++ic){
	  UInt ig=_dof_c.localToGlobal(eleID,iloc+1)-1+ic*_dim_c;     
	  _elvec_u[iloc+ic*_M_fe.nbNode] = u(ig);
       }
     }

    grad(0,_elvec_u,_elmatC,_M_fe,_M_fe,_M_fe);
    grad(1,_elvec_u,_elmatC,_M_fe,_M_fe,_M_fe);
    grad(2,_elvec_u,_elmatC,_M_fe,_M_fe,_M_fe);


// *************************************** Upwind ******************************

      double VLoc_infty=0.;
      double VLoc_mean=0.;
      double VLoc_c=0.;
      for (UInt ih_c=0 ; ih_c<(UInt)_M_fe.nbNode ; ih_c++){
         UInt  iloc = _M_fe.patternFirst(ih_c);
	 for (UInt ic=0; ic<nc_u;++ic){
	   UInt ig=_dof_c.localToGlobal(eleID,iloc+1)-1+ic*_dim_c;
           _elvec_u[iloc+ic*_M_fe.nbNode] = u(ig);
	   VLoc_c+=u(ig)*u(ig);}
	 VLoc_c=sqrt(VLoc_c);
	 VLoc_mean += VLoc_c;
	if (VLoc_c>VLoc_infty) VLoc_infty=VLoc_c;
      }
      VLoc_mean=VLoc_mean/_M_fe.nbNode;

      double coef_stab, Pe_loc;
//      coef_stab=_M_fe.diameter()*VLoc_infty; // Alessandro - method

      Pe_loc=VLoc_infty*_M_fe.diameter()/(2.0*_diffusivity);

//      coef_stab=(1.0/tanh(Pe_loc))-(1.0/Pe_loc); // classical approach

      if(Pe_loc < -3.0)
	 coef_stab= -1.0;
      else {
	 if(Pe_loc > 3.0)
	    coef_stab=1.0;
	 else	
	    coef_stab=Pe_loc/3.0;}

// ******************************* STREAMLINEUPWIND ****************************
     stiff_sd(coef_stab/(VLoc_mean*VLoc_mean),_elvec_u,_elmatC,_M_fe,_M_fe);

// ************************* Assembling ****************************************

     assemb_mat(_CDR,_elmatC,_M_fe,_dof_c);

  }
  chrono.stop();
  cout << "done in " << chrono.diff() << "s." << endl;

  // for BC treatment (done at each time-step)
  double tgv=1.e02; 

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

*/
