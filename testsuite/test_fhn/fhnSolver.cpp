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
  \brief A Fitzhugh-Nagumo solver (see the details in fhnHandler.hpp)
  \file fhnSolver.cpp
  \author J.-F. Gerbeau
  \date 09/2004
*/

#include "fhnSolver.hpp"
#include <life/lifefem/bcManage.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <life/lifefem/sobolevNorms.hpp>

namespace LifeV{
//====================================================

  FhNSolver::FhNSolver(const GetPot& data_file):
    FhNHandler(data_file),
    msrPattern(dof,2),
    mat(msrPattern),
    rhs(2*dimdof,2),
    uv(2*dimdof,2),
    uv_1(2*dimdof,2),
    elvec(refFE.nbDof,2),
    elvec_uv(refFE.nbDof,2),
    elmat(refFE.nbDof,2,2),
    vol_source(*this),
    uv_ref(nb_ref_time,dimdof),
    uv_stored(nb_ref_time,dimdof),
    ref_time(nb_ref_time)
  {
    switch(test_case){
    case 1:
      mesh.list_of_points(select_points1,pts_proc);
      break;
    }
    string refname = post_dir + "/" + ref_file;
    ifstream reffile(refname.c_str());
    std::cout << "Read a reference solution on the file " << refname << std::endl;
    if(reffile){
      flag_ref_solution = true;
      double nbt,nbd;
      reffile >> nbt;
      reffile >> nbd;
      if ( (nbt != nb_ref_time) || (nbd != dimdof) ){
    flag_ref_solution = false;
    std::cout << "WARNING : " << std::endl;
    std::cout << "   On your reference file : nb_ref_time = " << nbt
         << " and dim = " << nbd << std::endl;
    std::cout << "   It should be : nb_ref_time = " << nb_ref_time
         << " and dim = " << dimdof << std::endl;
      } else {
    double dummy;
    reffile >> dummy; // t=0.
    for(int ip=0;ip<(int)dimdof;ip++){
      reffile >> dummy; //  0.
    }
    for(int it=0;it<nb_ref_time;it++){
      reffile >> ref_time(it);
      for(int ip=0;ip<(int)dimdof;ip++){
        reffile >> uv_ref(it,ip);
      }
    }
    std::cout << "Reference solution ok." << std::endl;
      }
    } else {
      std::cout << "WARNING : cannot read the reference file " << refname << std::endl;
      flag_ref_solution = false;
    }    
  }

  void FhNSolver::compute_mat()
  {
    std::cout <<"*** Matrix computation\n";
    for(UInt i = 1; i<=mesh.numVolumes(); i++){
      fe.updateFirstDerivQuadPt(mesh.volumeList(i));
      elmat.zero();
      // first equation:
      // \dt u - div (fhn_diff \grad u) + v = fhn_f0 u(1-u)(u-fhn_alpha) 
      mass(1./time_step,elmat,fe,0,0);
      mass(1.,elmat,fe,0,1);
      switch(fhn_diff_fct){
      case 0:
    stiff(fhn_diff,elmat,fe,0,0);
    break;
      case 1:
    stiff(fhn_diff,diff_fct_1,elmat,fe,0,0);
    break;
      default:
    std::cout << "FhNSolver::compute_mat(): Unknown diffusion function "<< std::endl;
      }
      // second equation:
      // \dt v + fhn_eps fhn_gamma v - fhn_eps fhn_beta u = 0
      mass(-fhn_beta*fhn_eps,elmat,fe,1,0);
      mass((1./time_step + fhn_eps*fhn_gamma),elmat,fe,1,1);
      //
      assemb_mat(mat,elmat,fe,dof,0,0);
      assemb_mat(mat,elmat,fe,dof,1,0);
      assemb_mat(mat,elmat,fe,dof,0,1);
      assemb_mat(mat,elmat,fe,dof,1,1);
    }
    bcManageMatrix(mat,mesh,dof,bc,feBd,1.,0.);
  }
  
  void FhNSolver::compute_rhs()
  {
    std::cout <<"*** RHS computation..." << std::flush;
    rhs = ZeroVector(2*dimdof);
    for(UInt i=1; i<=mesh.numVolumes(); ++i){
      elvec.zero();
      fe.updateFirstDerivQuadPt(mesh.volumeList(i));
      extract_vec(uv_1,elvec_uv,refFE,dof,i,0);
      source(1./time_step,elvec_uv,elvec,fe,0,0); 
      source_fhn(fhn_f0,fhn_alpha,elvec_uv,elvec,fe,0,0);
      switch(test_case){
      case 2:case 3:
    source(vol_source,elvec,fe,time,0);
    break;
      }
      assemb_vec(rhs,elvec,fe,dof,0);

      extract_vec(uv_1,elvec_uv,refFE,dof,i,1);
      source(1./time_step,elvec_uv,elvec,fe,1,1); 
      assemb_vec(rhs,elvec,fe,dof,1);
    }
    std::cout <<"done." << std::endl;
  }

  void FhNSolver::initial_data()
  {
    std::cout << "*** Initial data " << std::endl;
    // be careful : we impose the value only on the points of the mesh (ok for P1)
    switch(init_data){
    case 0: // initial solution computed 
      for(UInt i=1;i<=mesh.numPoints();i++){
    Geo0D& pt = mesh.pointList(i);
    switch(test_case){
    case 1:
      uv[i-1] = u_init1(pt.x(),pt.y(),pt.z());
      break;
    case 2:case 3:
      uv[i-1] = u_init2(pt.x(),pt.y(),pt.z());
      break;
    default:
      std::cout << "initial_data() : test case = " << test_case << " ?? " << std::endl;
      exit(1);
    }
      }
      break;
    case 1:
      // initial solution read on a file
      {
    UInt dummy;
    string init_file_name_u = post_dir +"/"+ init_file_name;
    string init_file_name_v =  post_dir +"/" +init_file_name;
    init_file_name_u.replace(init_file_name_u.find(".bb"),3,"_u.bb");
    init_file_name_v.replace(init_file_name_v.find(".bb"),3,"_v.bb");
    std::cout << "*** The initial solutions are read on the files " ;
    std::cout << init_file_name_u << " and " << init_file_name_v << std::endl;
    //--------
    ifstream init_file_u(init_file_name_u.c_str());
    init_file_u >> dummy ; // dimension
    init_file_u >> dummy ; // 1 
    init_file_u >> dummy ; // size
    if(dummy != dimdof){
      std::cout << "FhNSolver::initial_data() Critical error : "
           << dimdof << " != " << dummy << std::endl;
      exit(1);
    }
    init_file_u >> dummy ; // type
    for(UInt i=0;i<mesh.numPoints();i++){
      init_file_u >> uv[i];
    }
    //--------
    ifstream init_file_v(init_file_name_v.c_str());
    init_file_v >> dummy ; // dimension
    init_file_v >> dummy ; // 1 
    init_file_v >> dummy ; // size
    if(dummy != dimdof){
      std::cout << "FhNSolver::initial_data() Critical error : "
           << dimdof << " != " << dummy << std::endl;
      exit(1);
    }
    init_file_v >> dummy ; // type
    std::cout << "dimdof = " << dimdof << std::endl;
    std::cout << "meshnumpts = " << mesh.numPoints() << std::endl;
    
    for(UInt i=0;i<mesh.numPoints();i++){
      init_file_v >> uv[dimdof + i];
    }
    //--------
    break;
      }
    }
  }
  
  void FhNSolver::time_advance()
  {
    //
    static bool flag = true;
    static int az_name;
    compute_rhs();
    bcManageVector(rhs,mesh,dof,bc,feBd,time,1.);
    //
    int options[AZ_OPTIONS_SIZE]; 
    double params[AZ_PARAMS_SIZE];
    // we first initialize Aztec with its defaults and user's parameters
    aztecOptionsFromDataFile(options,params);
    if(flag){
      std::cout << "*** Aztec: factorization of the preconditioner" << std::endl;
      options[AZ_keep_info] = 1;
    } else {
      std::cout << "*** Aztec: old preconditioner" << std::endl;
      options[AZ_keep_info] = 1;
      options[AZ_pre_calc] = AZ_reuse;
    }                         
    // next, we may overload some of them
    /*
      params[AZ_tol] = 1e-9;
      options[AZ_solver] = AZ_gmres;
      options[AZ_precond] = AZ_dom_decomp;
      options[AZ_subdomain_solve] = AZ_ilu;
    */
    aztecSolveLinearSyst(mat,uv.giveVec(),rhs.giveVec(),uv.size(),
             msrPattern,options,params,az_name,flag);  
    flag = false; // to avoid the re-computation of the preconditionner
    //
  }
  
  void FhNSolver::solve()
  {
    initial_data();
    time = init_time;
    compute_mat();
    post_process();
    write_adapt();
    for(int iter=1;iter<=max_time_iter && time<max_time;iter++){
      time += time_step;
      timeBanner(iter,time);
      uv_1 = uv;
      time_advance();
      if  (!(iter % post_proc_period)) post_process();
      if  (!(iter % adapt_period)) write_adapt();
      if( max_time_iter / nb_ref_time ){
    if  (!(iter % ( max_time_iter / nb_ref_time))) store_solution();
      }
    }
    
    save_stored_solution();
    if(flag_ref_solution) compare_stored_with_ref_sol();
    //
    
  }

  void FhNSolver::post_process()
  {
    static int iter_post=0;
    int ipts;
    double x,y,z;
    string mtvname;
    ofstream fmtv;
    
    char str_iter[10],str_time[10];
    string vtkname,bbname;
    sprintf(str_time,"t=%f",time);
    sprintf(str_iter,".%03d",iter_post);
    //
    if(post_proc_format == "medit"){
      bbname = post_dir + "/u" + (string) str_iter + ".bb";
      if(verbose) std::cout << "*** Postprocessing on " << bbname << std::endl;
      wr_medit_ascii_scalar(bbname,uv.giveVec(),dimdof);
    } else if (post_proc_format == "vtk"){
      vtkname = post_dir + "/u" + (string) str_iter + ".vtk";
      if(verbose) std::cout << "*** Postprocessing on " << vtkname << std::endl;
      wr_vtk_ascii_header(vtkname,"Action potential",mesh,dof,fe);
      wr_vtk_ascii_scalar(vtkname,"u",uv.giveVec(),dimdof);
    }
    //
    switch(test_case){
    case 1:
      mtvname = post_dir + "/u.mtv";
      if(iter_post==0){
    fmtv.open(mtvname.c_str());
      } else {
    fmtv.open(mtvname.c_str(),ios::app);
      }
      fmtv << "$ DATA = CURVE2D\n\n";
      fmtv << "%toplabel = 't=" << time << "'\n\n";
      fmtv << "% boundary = True \n";
      fmtv << "% ymax = 1.2 ymin = -0.5 \n";
      for(UInt i=0;i<pts_proc.size();i++){
    ipts = pts_proc[i];
    Geo0D& pt = mesh.pointList(ipts);
    fmtv << pt.z() << "\t" << uv[ipts-1] << std::endl;
      }
      //
      fmtv.close();

      mtvname = post_dir + "/v.mtv";
      if(iter_post==0){
    fmtv.open(mtvname.c_str());
      } else {
    fmtv.open(mtvname.c_str(),ios::app);
      }
      fmtv << "$ DATA = CURVE2D\n\n";
      fmtv << "%toplabel = 't=" << time << "'\n\n";
      fmtv << "% boundary = True \n";
      fmtv << "% ymax = 1.2 ymin = -0.5 \n";
      for(UInt i=0;i<pts_proc.size();i++){
    ipts = pts_proc[i];
    Geo0D& pt = mesh.pointList(ipts);
    fmtv << pt.z() << "\t" << uv[dimdof+ipts-1] << std::endl;
      }
      //
      fmtv.close();
      break;
    case 2:case 3:
      //-------
      ipts = 719;
      x = mesh.pointList(ipts).x();
      y= mesh.pointList(ipts).y();
      z= mesh.pointList(ipts).z();
      mtvname = post_dir + "/u_appex.mtv";
      if(iter_post==0){
    fmtv.open(mtvname.c_str());
    fmtv << "$ DATA = CURVE2D\n\n";
    fmtv << "% toplabel = 'points (" << x <<","<< y <<"," << z << ")'\n";
    fmtv << "% boundary = True \n";
    fmtv << "% ymax = 1.1 ymin = -0.3 \n";
      } else {
    fmtv.open(mtvname.c_str(),ios::app);
      }
      fmtv << time << "\t" << uv[ipts-1] << std::endl;
      fmtv.close();
      //-------
      ipts = 761;
      x = mesh.pointList(ipts).x();
      y= mesh.pointList(ipts).y();
      z= mesh.pointList(ipts).z();
      mtvname = post_dir + "/u_mid.mtv";
      if(iter_post==0){
    fmtv.open(mtvname.c_str());
    fmtv << "$ DATA = CURVE2D\n\n";
    fmtv << "% toplabel = 'points (" << x <<","<< y <<"," << z << ")'\n";
    fmtv << "% boundary = True \n";
    fmtv << "% ymax = 1.1 ymin = -0.3 \n";
      } else {
    fmtv.open(mtvname.c_str(),ios::app);
      }
      fmtv << time << "\t" << uv[ipts-1] << std::endl;
      fmtv.close();      
      //-------
      ipts = 625;
      x = mesh.pointList(ipts).x();
      y= mesh.pointList(ipts).y();
      z= mesh.pointList(ipts).z();
      mtvname = post_dir + "/u_base.mtv";
      if(iter_post==0){
    fmtv.open(mtvname.c_str());
    fmtv << "$ DATA = CURVE2D\n\n";
    fmtv << "% toplabel = 'points (" << x <<","<< y <<"," << z << ")'\n";
    fmtv << "% boundary = True \n";
    fmtv << "% ymax = 1.1 ymin = -0.3 \n";
      } else {
    fmtv.open(mtvname.c_str(),ios::app);
      }
      fmtv << time << "\t" << uv[ipts-1] << std::endl;
      fmtv.close();
      break;
    }
    iter_post ++;
  }

  void FhNSolver::write_adapt()
  {
    static int iter_adapt=0;
    char str_iter[10],str_time[10];
    string bbname;
    sprintf(str_time,"t=%f",time);
    sprintf(str_iter,".%03d",iter_adapt);
    //
    bbname = post_dir + "/u_adap" + (string) str_iter + ".bb";
    if(verbose) std::cout << "*** Save for adaption on " << bbname << std::endl;
    wr_medit_ascii_scalar(bbname,uv.giveVec(),dimdof);
    bbname = post_dir + "/v_adap" + (string) str_iter + ".bb";
    if(verbose) std::cout << "*** Save for adaption on " << bbname << std::endl;
    wr_medit_ascii_scalar(bbname,uv.giveVec()+dimdof,dimdof);
    iter_adapt ++;
  }

  void FhNSolver::store_solution()
  {
    std::cout << "*** store the solution" << std::endl;
    static int iter_save=0;
    for(UInt ipts=0;ipts<mesh.numPoints();ipts++){
      if(iter_save == nb_ref_time){
    std::cout << "pb de dim" << std::endl;
    exit(1);
      }
      if(flag_ref_solution && fabs(time-1000*ref_time(iter_save))>1.e-6){// ref time in milliseconds
    std::cout << "t        = " << time << std::endl;
    std::cout << "ref_time = " << ref_time(iter_save)*1000 << std::endl;
    std::cout << "WARNING : the current time is different from the corresponding reference time. The current solution will not be compared to the reference one." << std::endl;
    flag_ref_solution = false;
      }
      if(!flag_ref_solution) ref_time(iter_save) = 0.001*time; // convert in milliseconds
      uv_stored(iter_save,(int)ipts) = uv[ipts];
    }
    iter_save++;
  }
  
  void FhNSolver::save_stored_solution()
  {
    ofstream fref;
    string refname = post_dir + "/" + store_file;
    std::cout << "*** save the stored solution on " << refname << std::endl;
    fref.open(refname.c_str());
    fref << nb_ref_time << " " << dimdof << std::endl;
    // we write 0. at time t=0.
    fref << 0. ;
    for(int ipts=0;ipts<(int)mesh.numPoints();ipts++){
      fref.precision(9);
      fref << " " << 0.;
    }
    fref << std::endl;
    for(int i=0;i<nb_ref_time; i++){
      fref << ref_time(i);
      for(int ipts=0;ipts<(int)mesh.numPoints();ipts++){
    fref.precision(9);
    fref << " " << uv_stored(i,ipts);
      }
      fref << std::endl;
    }
  }

  void FhNSolver::compare_stored_with_ref_sol()
  {
    ofstream fcomp;
    string compname = post_dir + "/" + "diff_" + store_file;
    double L2x_ref,H1x_ref,L2x,H1x;
    double L2tL2x=0.,L2tH1x=0.,L2x_1=0.,H1x_1=0.;
    double L2tL2x_ref=0.,L2tH1x_ref=0.,L2x_ref_1=0.,H1x_ref_1=0.;
    KN<double> tmp(mesh.numPoints());
    fcomp.open(compname.c_str());
    fcomp << "# time \tL2 norm \tRel diff L2 \tRel diff semi H1" << std::endl;
    for(int it=0;it<nb_ref_time;it++){
      fcomp << ref_time(it) << "\t";
      for(int ip=0;ip<(int)mesh.numPoints();ip++){
    tmp(ip) = (uv_stored(it,ip) - uv_ref(it,ip));
      }
      L2_and_H1_norms(tmp,L2x,H1x);
      for(int ip=0;ip<(int)mesh.numPoints();ip++){
    tmp(ip) = uv_ref(it,ip);
      }
      L2_and_H1_norms(tmp,L2x_ref,H1x_ref);
      if(it){
    L2tL2x += 0.5*(L2x+L2x_1)*(ref_time(it) - ref_time(it-1));
    L2tL2x_ref += 0.5*(L2x_ref+L2x_ref_1)*(ref_time(it) - ref_time(it-1));
    L2tH1x += 0.5*(H1x+H1x_1)*(ref_time(it) - ref_time(it-1));
    L2tH1x_ref += 0.5*(H1x_ref+H1x_ref_1)*(ref_time(it) - ref_time(it-1));
      }
      L2x_1 = L2x;
      H1x_1 = H1x;
      L2x_ref_1 = L2x_ref;
      H1x_ref_1 = H1x_ref;
      fcomp << sqrt(L2x)
        <<"\t" << sqrt(L2x/L2x_ref) << "\t" << sqrt(H1x/H1x_ref) << std::endl;
    }
    fcomp << "# Relative difference in L^2_t(0,T;L^2_x) " <<  sqrt(L2tL2x/L2tL2x_ref)
      << std::endl;
    fcomp << "# Relative difference in L^2_t(0,T;semi H^1_x) "
      <<  sqrt(L2tH1x/L2tH1x_ref) << std::endl;
    fcomp.close();
  }

  void FhNSolver::L2_and_H1_norms(const KN<double>& vec,double& L2_2,double& H1_2)
  {
    //return the *square* of the L2 and H1 norms 
    L2_2 = 0.;
    H1_2 = 0.;
    for(UInt i = 1; i<=mesh.numVolumes(); i++){
      fe.updateFirstDeriv(mesh.volumeList(i));
      L2_2 += elem_L2_2(vec,fe,dof);
      H1_2 += elem_H1_2(vec,fe,dof);
    }
  }
}

