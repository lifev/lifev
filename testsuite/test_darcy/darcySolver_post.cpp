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
#include "darcySolver.hpp"
#include "medit_wrtrs.hpp"

void DarcySolver::postProcessPressureQ0()
{
  if(verbose) cout << "Postprocessing of pressure (constant by element)\n";
  if(post_proc_format == "medit"){
    wr_medit_ascii_scalar(post_dir + "/presQ0.bb",globalP.giveVec(),globalP.size(),1);
  } else {
    cerr
      <<"Warning: Solution constant by element is possible only with medit for the moment\n";
  }
}
  
void DarcySolver::postProcessPressureQ1()
{
  if(verbose)
    cout << "Postprocessing of pressure (L2 projection on the nodes)\n";
  // Q1 elements
  const RefFE& refFE    = feHexaQ1;
  CurrentFE fe_q1(refFE,geoMap,qr);
  Dof dof_q1(refFE); 
  dof_q1.update(mesh);
  UInt dim_q1 = dof_q1.numTotalDof();
  ScalUnknown<Vector> p_q1(dim_q1), f_q1(dim_q1);
  p_q1.vec()=0.0;
  f_q1.vec()=0.0;
  MSRPatt pattA_q1(dof_q1);
  MSRMatr<double> A_q1(pattA_q1);
  ElemMat elmat(fe_q1.nbNode,1,1);
  ElemVec elvec(fe_q1.nbNode,1);
  for(UInt i = 1; i<=mesh.numVolumes(); i++){
    fe_q1.updateJac(mesh.volumeList(i));
    elmat.zero();
    elvec.zero();
    mass(1.,elmat,fe_q1);
    source(globalP.vec()(i-1),elvec,fe_q1,0);
    assemb_mat(A_q1,elmat,fe_q1,dof_q1,0,0);
    assemb_vec(f_q1.vec(),elvec,fe_q1,dof_q1,0);
  }
  int    options[AZ_OPTIONS_SIZE]; 
  double params[AZ_PARAMS_SIZE];
  // we first initialize Aztec with its defaults and user's parameters
  aztecOptionsFromDataFile(options,params);
  // next, we overload some of them
  params[AZ_tol] = 1e-9;
  options[AZ_solver] = AZ_cg;
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_icc;
  aztecSolveLinearSyst(A_q1,p_q1.giveVec(),f_q1.giveVec(),p_q1.size(),
		       pattA_q1,options,params);  
  //
  if(post_proc_format == "medit"){
    wr_medit_ascii_scalar(post_dir + "/presQ1.bb",p_q1.giveVec(),p_q1.size());
  } else if (post_proc_format == "vtk"){
    wr_vtk_ascii_header(post_dir + "/pres.vtk","Title",mesh, dof_q1, fe_q1);
    wr_vtk_ascii_scalar(post_dir + "/pres.vtk","scal",p_q1.giveVec(), p_q1.size());
  }
}

void DarcySolver::postProcessVelocityQ1()
{
  if(verbose)
    cout << "Postprocessing of velocity (L2 projection on the nodes)\n";
  // Q1 elements
  const RefFE& refFE    = feHexaQ1;
  CurrentFE fe_q1(refFE,geoMap,qr);
  Dof dof_q1(refFE);
  dof_q1.update(mesh);
  UInt dim_q1 = dof_q1.numTotalDof();
  PhysVectUnknown<Vector> u_q1(dim_q1), f_q1(dim_q1);
  u_q1.vec()=0.0;
  f_q1.vec()=0.0;
  MSRPatt pattA_q1(dof_q1,nbCoor);
  MSRMatr<double> A_q1(pattA_q1);
  ElemMat elmat_hdiv(fe_q1.nbNode,nbCoor,0,
		     vfe.nbNode,0,1);  
  ElemMat elmat(fe_q1.nbNode,nbCoor,nbCoor);
  ElemVec elvec(fe_q1.nbNode,nbCoor);
  ElemVec elvec_hdiv(vfe.nbNode,1);
  Tab1dView elvec_hdiv_vec = elvec_hdiv.block(0);
  
  for(UInt i = 1; i<=mesh.numVolumes(); i++){
    fe_q1.updateJac(mesh.volumeList(i));
    vfe.updatePiola(mesh.volumeList(i));
    elmat.zero();
    elmat_hdiv.zero();
    mass(1.,elmat,fe_q1,0,0,nbCoor);
    mass_Mixed_Hdiv(1.,elmat_hdiv,fe_q1,vfe,0,0);
    extract_vec(globalFlux.vec(),elvec_hdiv,refVFE,vdof,mesh.volumeList(i).id(),0);
    //
    for(int j=0;j<(int) mesh.volumeList(i).numLocalFaces;j++){
      elvec_hdiv_vec[j] *= signLocalFace( (int)mesh.volumeList(i).id() - 1, j);
    }
    //
    elvec.vec() = elmat_hdiv.mat() * elvec_hdiv.vec();
    for(UInt icoor = 0; icoor < nbCoor; icoor++){
      assemb_mat(A_q1,elmat,fe_q1,dof_q1,icoor,icoor);
      assemb_vec(f_q1.vec(),elvec,fe_q1,dof_q1,icoor);
    }
  }
  int    options[AZ_OPTIONS_SIZE];
  double params[AZ_PARAMS_SIZE];
  // we first initialize Aztec with its defaults and user's parameters
  aztecOptionsFromDataFile(options,params);
  // next, we overload some of them
  params[AZ_tol] = 1e-9;
  options[AZ_solver] = AZ_cg;
  options[AZ_precond] = AZ_dom_decomp;
  options[AZ_subdomain_solve] = AZ_icc;
  aztecSolveLinearSyst(A_q1,u_q1.giveVec(),f_q1.giveVec(),u_q1.size(),
		       pattA_q1,options,params);  
  //
  if(post_proc_format == "medit"){
    wr_medit_ascii_vector(post_dir + "/velQ1.bb",u_q1.giveVec(),u_q1.size());
  } else if (post_proc_format == "vtk"){
    wr_vtk_ascii_header(post_dir + "/vel.vtk","Title",mesh, dof_q1, fe_q1);
    wr_vtk_ascii_vector(post_dir + "/vel.vtk","Velocity",u_q1.giveVec(), u_q1.size());
  }
}
