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
#include "sobolevNorms.hpp"

#define ANALYTICAL_SOL 0

namespace LifeV
{
void DarcySolver::postProcessTraceOfPressureRT0()
{
  if(verbose) std::cout << "Postprocessing of TP (RT0 per element)\n";
  if(post_proc_format == "medit"){
    wr_medit_ascii_scalar(post_dir + "/presTP0.bb",globalTP.giveVec(),globalTP.size(),1);
  } else {
    std::cerr
      <<"Warning: Solution constant by element is possible only with medit for the moment\n";
  }
}

void DarcySolver::postProcessVelocityRT0()
{
  if(verbose) std::cout << "Postprocessing of velocity (RT0 per element)\n";
  if(post_proc_format == "medit"){
    wr_medit_ascii_scalar(post_dir + "/velocRT0.bb",globalFlux.giveVec(),globalFlux.size(),1);
  } else {
    std::cerr
      <<"Warning: Solution constant by element is possible only with medit for the moment\n";
  }
}

void DarcySolver::postProcessPressureQ0()
{
  if(verbose) std::cout << "Postprocessing of pressure (constant by element)\n";
  if(post_proc_format == "medit"){
    wr_medit_ascii_scalar(post_dir + "/presQ0.bb",globalP.giveVec(),globalP.size(),1);
  } else {
    std::cerr
      <<"Warning: Solution constant by element is possible only with medit for the moment\n";
  }


#if ANALYTICAL_SOL
  std::cout <<"Compute L2 pressure error:\n";
  AnalyticalSolPres analyticSol;

  double normL2=0., normL2diff=0., normL2sol=0.;
  double normL2sq=0., normL2diffsq=0., normL2solsq=0.;

  for(UInt i=1; i<=mesh.numVolumes(); ++i){

    pfe.updateFirstDeriv(mesh.volumeList(i));

    normL2sq     += elem_L2_2(globalP,pfe,pdof);
    normL2solsq  += elem_L2_2(analyticSol,pfe);
    normL2diffsq += elem_L2_diff_2(globalP,analyticSol,pfe,pdof);

  }

  normL2     = sqrt(normL2sq);
  normL2sol  = sqrt(normL2solsq);
  normL2diff = sqrt(normL2diffsq);

  std::string errname = post_dir + "/errQ0Pres.txt";
  std::ofstream ofile(errname.c_str());

  ASSERT(ofile,"Error: Output file cannot be opened.");
  ofile << "PRESSION ERROR (Q0)" << std::endl;
  ofile << "|| p       ||_{L^2}                   = " << normL2 << std::endl;
  ofile << "|| p_ex     ||_{L^2}                   = " << normL2sol << std::endl;
  ofile << "|| p - p_ex ||_{L^2}                   = " << normL2diff<< std::endl;
  ofile << "|| P - p_ex ||_{L^2} / || p_ex ||_{L^2} = "
	<< normL2diff / normL2sol << "\n" << std::endl;
  ofile << "SQUARE of PRESSION ERROR (Q0)" << std::endl;
  ofile << "|| p       ||^2_{L^2}                   = " << normL2sq << std::endl;
  ofile << "|| p_ex     ||^2_{L^2}                   = " << normL2solsq << std::endl;
  ofile << "|| p - p_ex ||^2_{L^2}                   = " << normL2diffsq << std::endl;
  ofile << "|| P - p_ex ||^2_{L^2} / || p_ex ||^2_{L^2} = "
	<< normL2diffsq / normL2solsq << std::endl;
#endif

}

void DarcySolver::postProcessPressureQ1()
{
  if(verbose)
    std::cout << "Postprocessing of pressure (L2 projection on the nodes)\n";
  // Q1 elements
  const RefFE& refFE    = feHexaQ1;
  CurrentFE fe_q1(refFE,geoMap,qr);
  Dof dof_q1(refFE);
  dof_q1.update(mesh);
  UInt dim_q1 = dof_q1.numTotalDof();
  ScalUnknown<Vector> p_q1(dim_q1), f_q1(dim_q1);
  p_q1=0.0;
  f_q1=0.0;
  MSRPatt pattA_q1(dof_q1);
  MSRMatr<double> A_q1(pattA_q1);
  ElemMat elmat(fe_q1.nbNode,1,1);
  ElemVec elvec(fe_q1.nbNode,1);
  for(UInt i = 1; i<=mesh.numVolumes(); i++){
    fe_q1.updateJac(mesh.volumeList(i));
    elmat.zero();
    elvec.zero();
    mass(1.,elmat,fe_q1);
    source(globalP(i-1),elvec,fe_q1,0);
    assemb_mat(A_q1,elmat,fe_q1,dof_q1,0,0);
    assemb_vec(f_q1,elvec,fe_q1,dof_q1,0);
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
  std::string vtkname,bbname;
  /*
  char str_iter[10],str_time[10];
  static int iter_post=0;
  sprintf(str_time,"t=%f",time);
  sprintf(str_iter,".%03d",iter_post);
  */
  //
  if(post_proc_format == "medit"){
    //bbname = post_dir + "/presQ1" + (string) str_iter + ".bb";
    bbname = post_dir + "/presQ1.bb";
    wr_medit_ascii_scalar(bbname,p_q1.giveVec(),p_q1.size());
  } else if (post_proc_format == "vtk"){
    //vtkname = post_dir + "/presQ1" + (string) str_iter + ".vtk";
    vtkname = post_dir + "/presQ1.vtk";
    wr_vtk_ascii_header(vtkname,"Pressure",mesh, dof_q1, fe_q1);
    wr_vtk_ascii_scalar(vtkname,"P",p_q1.giveVec(), p_q1.size());
  }
  //iter_post ++;
  //---------------------------------------------


#if ANALYTICAL_SOL
  std::cout <<"Compute pressure error:\n";
  AnalyticalSolPres analyticSol;

  double normL2=0., normL2diff=0., normL2sol=0.;
  double normH1=0., normH1diff=0., normH1sol=0.;

  for(UInt i=1; i<=mesh.numVolumes(); ++i){

    fe_q1.updateFirstDeriv(mesh.volumeList(i));

    normL2     += elem_L2_2(p_q1,fe_q1,dof_q1);
    normL2sol  += elem_L2_2(analyticSol,fe_q1);
    normL2diff += elem_L2_diff_2(p_q1,analyticSol,fe_q1,dof_q1);

    normH1     += elem_H1_2(p_q1,fe_q1,dof_q1);
    normH1sol  += elem_H1_2(analyticSol,fe_q1);
    normH1diff += elem_H1_diff_2(p_q1,analyticSol,fe_q1,dof_q1);
  }

  normL2     = sqrt(normL2);
  normL2sol  = sqrt(normL2sol);
  normL2diff = sqrt(normL2diff);

  normH1     = sqrt(normH1);
  normH1sol  = sqrt(normH1sol);
  normH1diff = sqrt(normH1diff);

  std::string errname = post_dir + "/errQ1Pres.txt";
  std::ofstream ofile(errname.c_str());

  ASSERT(ofile,"Error: Output file cannot be opened.");
  ofile << "PRESSION ERROR (Q1)" << std::endl;
  ofile << "|| p       ||_{L^2}                   = " << normL2 << std::endl;
  ofile << "|| p_ex     ||_{L^2}                   = " << normL2sol << std::endl;
  ofile << "|| p - p_ex ||_{L^2}                   = " << normL2diff<< std::endl;
  ofile << "|| P - p_ex ||_{L^2} / || p_ex ||_{L^2} = " << normL2diff/normL2sol
       << std::endl;

  ofile << "|| U       ||_{H^1}                   = " << normH1 << std::endl;
  ofile << "|| sol     ||_{H^1}                   = " << normH1sol << std::endl;
  ofile << "|| U - sol ||_{H^1}                   = " << normH1diff<< std::endl;
  ofile << "|| U - sol ||_{H^1} / || sol ||_{H^1} = " << normH1diff/normH1sol
       << std::endl;
#endif
}

void DarcySolver::postProcessVelocityQ1()
{
  if(verbose)
    std::cout << "Postprocessing of velocity (L2 projection on the nodes)\n";
  // Q1 elements
  const RefFE& refFE    = feHexaQ1;
  CurrentFE fe_q1(refFE,geoMap,qr);
  Dof dof_q1(refFE);
  dof_q1.update(mesh);
  UInt dim_q1 = dof_q1.numTotalDof();
  PhysVectUnknown<Vector> u_q1(dim_q1), f_q1(dim_q1);
  u_q1=0.0;
  f_q1=0.0;
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
    extract_vec(globalFlux,elvec_hdiv,refVFE,vdof,mesh.volumeList(i).id(),0);
    //
    for(int j=0;j<(int) mesh.volumeList(i).numLocalFaces;j++){
      elvec_hdiv_vec[j] *= signLocalFace( (int)mesh.volumeList(i).id() - 1, j);
    }
    //
    elvec.vec() = elmat_hdiv.mat() * elvec_hdiv.vec();
    for(UInt icoor = 0; icoor < nbCoor; icoor++){
      assemb_mat(A_q1,elmat,fe_q1,dof_q1,icoor,icoor);
      assemb_vec(f_q1,elvec,fe_q1,dof_q1,icoor);
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
  std::string vtkname,bbname;
  /*
  char str_iter[10],str_time[10];
  static int iter_post=0;
  sprintf(str_time,"t=%f",time);
  sprintf(str_iter,".%03d",iter_post);
  */
  //
  if(post_proc_format == "medit"){
    //bbname = post_dir + "/vel" + (string) str_iter + ".bb";
    bbname = post_dir + "/velQ1.bb";
    wr_medit_ascii_vector(bbname,u_q1.giveVec(),u_q1.size());
  } else if (post_proc_format == "vtk"){
    //vtkname = post_dir + "/vel" + (string) str_iter + ".vtk";
    vtkname = post_dir + "/velQ1.vtk";
    wr_vtk_ascii_header(vtkname,"Velocity",mesh, dof_q1, fe_q1);
    wr_vtk_ascii_vector(vtkname,"U",u_q1.giveVec(), u_q1.size());
  }
  //iter_post ++;

  //---------------------------------------------

#if ANALYTICAL_SOL
  std::cout <<"Velocity error:\n";
  AnalyticalSolFlux analyticSol;

  /*
  analyticSol.init(diffusion_tensor(0,0)  / diffusion_scalar,
		   diffusion_tensor(2,0) / diffusion_scalar,
		   diffusion_tensor(2,2) / diffusion_scalar,
		   diffusion_tensor(1,1) / diffusion_scalar);
  */

  double normL2=0., normL2diff=0., normL2sol=0.;
  double normH1=0., normH1diff=0., normH1sol=0.;

  for(UInt i=1; i<=mesh.numVolumes(); ++i){

    fe_q1.updateFirstDeriv(mesh.volumeList(i));

    normL2     += elem_L2_2(u_q1,fe_q1,dof_q1,3);
    normL2sol  = -1.; //! elem_L2_2 is not reckognized
    //! (confusion with another templated function)
    // normL2sol  += elem_L2_2<AnalyticalSolFlux>(analyticSol,fe_q1,0.0,3);
    normL2diff += elem_L2_diff_2(u_q1,analyticSol,fe_q1,dof_q1,0.,3);
    /*
    normH1     += elem_H1_2(u_q1,fe_q1,dof_q1,0,3);
    normH1sol  += elem_H1_2(analyticSol,fe_q1,0,3);
    normH1diff += elem_H1_diff_2(u_q1,analyticSol,fe_q1,dof_q1,0,3);
    */
  }

  normL2     = sqrt(normL2);
  normL2sol  = sqrt(normL2sol);
  normL2diff = sqrt(normL2diff);

  normH1     = sqrt(normH1);
  normH1sol  = sqrt(normH1sol);
  normH1diff = sqrt(normH1diff);

  std::string errname = post_dir + "/errQ1Vel.txt";
  std::ofstream ofile(errname.c_str());

  ASSERT(ofile,"Error: Output file cannot be opened.");
  ofile << "VELOCITY ERROR (Q1)" << std::endl;

  ofile << "|| U         ||_{L^2}                   = " << normL2 << std::endl;
  ofile << "|| exact     ||_{L^2}                   = " << normL2sol << std::endl;
  ofile << "|| U - exact ||_{L^2}                   = " << normL2diff<< std::endl;
  ofile << "|| U - exact ||_{L^2}/|| exact ||_{L^2} = " << normL2diff/normL2sol
       << std::endl;

  //  std::cerr << "|| U       ||_{H^1}                   = " << normH1 << std::endl;
  //  std::cerr << "|| exact     ||_{H^1}                   = " << normH1sol << std::endl;
  //  std::cerr << "|| U - exact ||_{H^1}                   = " << normH1diff<< std::endl;
  // std::cerr << "|| U - exact ||_{H^1} / || exact ||_{H^1} = " << normH1diff/normH1sol
  //   << std::endl;

#endif
}
}
