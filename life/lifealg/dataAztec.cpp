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
#include <iostream>
#include <vector>
#include <fstream>

#include "dataAztec.hpp"

namespace LifeV
{

DataAztec::DataAztec( const GetPot& dfile, const std::string& section ) :
        aztec_solver_list( section + "/solver" ),
        aztec_precond_list( section + "/precond" ),
        aztec_scaling_list( section + "/scaling" ),
        aztec_conv_list( section + "/conv" ),
        aztec_output_list( section + "/output" ),
        aztec_subdomain_solve_list( section + "/subdomain_solve" )
{
    //
    // solver
    //
    aztec_solver_list.add( "cg", AZ_cg, "preconditioned conjugate gradient method" );
    aztec_solver_list.add( "gmres", AZ_gmres,
                           "preconditioned gmres method (default)" );
    aztec_solver_list.add( "cgs", AZ_cgs, "preconditioned cg squared method" );
    aztec_solver_list.add( "tfqmr", AZ_tfqmr,
                           "preconditioned transpose-free qmr method" );
    aztec_solver_list.add( "bicgstab", AZ_bicgstab,
                           "preconditioned stabilized bi-cg method" );
    aztec_solver_list.add( "slu", AZ_slu, "super LU direct method" );
    aztec_solver_list.add( "symmlq", AZ_symmlq, "indefinite symmetric like symmlq" );
    aztec_solver_list.add( "fixed_pt", AZ_fixed_pt, "fixed point iteration" );
    aztec_solver_list.add( "lu", AZ_lu, "sparse LU direct method" );
    //
    aztec_solver_str = dfile( ( section + "/solver" ).data(), "gmres" );
    aztec_solver = aztec_solver_list.value( aztec_solver_str );
    //
    // scaling
    //
    aztec_scaling_list.add( "none", AZ_none, "no scaling (default)" );
    aztec_scaling_list.add( "jacobi", AZ_Jacobi, "Jacobi scaling" );
    aztec_scaling_list.add( "bjacobi", AZ_BJacobi, "block Jacobi scaling" );
    aztec_scaling_list.add( "row_sum", AZ_row_sum, "point row-sum scaling" );
    aztec_scaling_list.add( "sym_diag", AZ_sym_diag, "symmetric diagonal scaling" );
    aztec_scaling_list.add( "sym_row_sum", AZ_sym_row_sum,
                            "symmetric row-sum scaling" );
    aztec_scaling_list.add( "sym_bjacobi", AZ_sym_BJacobi,
                            "symmetric block Jacobi scaling" );
    //
    aztec_scaling_str = dfile( ( section + "/scaling" ).data(), "none" );
    aztec_scaling = aztec_scaling_list.value( aztec_scaling_str );
    //
    // precond
    //
    aztec_precond_list.add( "none", AZ_none, "no preconditioning (default)" );
    aztec_precond_list.add( "jacobi", AZ_Jacobi, "Jacobi preconditioning" );
    aztec_precond_list.add( "sym_gs", AZ_sym_GS,
                            "symmetric Gauss-Siedel preconditioning" );
    aztec_precond_list.add( "neumann", AZ_Neumann,
                            "Neumann series polynomial preconditioning" );
    aztec_precond_list.add( "ls", AZ_ls,
                            "least-squares polynomial preconditioning" );
    aztec_precond_list.add( "smoother", AZ_smoother,
                            "Recursive call to AZ_iterate()" );
    aztec_precond_list.add( "dom_decomp", AZ_dom_decomp,
                            "Domain decomposition using subdomain solver given by options[AZ_subdomain_solve] = ilu, ilut, icc, rilu, bilu or lu" );
    aztec_precond_list.add( "user_precond", AZ_user_precond,
                            "user's preconditioning" );
    //
    aztec_precond_str = dfile( ( section + "/precond" ).data(), "none" );
    aztec_precond = aztec_precond_list.value( aztec_precond_str );
    //
    // Convergence type
    //
    aztec_conv_list.add( "r0", AZ_r0, " ||r||_2 / ||r^{(0)}||_2 (default)" );
    aztec_conv_list.add( "rhs", AZ_rhs, "||r||_2 / ||b||_2 " );
    aztec_conv_list.add( "anorm", AZ_Anorm, "||r||_2 / ||A||_infty " );
    aztec_conv_list.add( "sol", AZ_sol,
                         "||r||_infty/(||A||_infty ||x||_1+||b||_infty)" );
    aztec_conv_list.add( "weighted", AZ_weighted, "||r||_WRMS" );
    aztec_conv_list.add( "expected_values", AZ_expected_values,
                         "||r||_WRMS with weights taken as |A||x0|" );
    aztec_conv_list.add( "noscaled", AZ_noscaled, "||r||_2 " );
    //
    aztec_conv_str = dfile( ( section + "/conv" ).data(), "r0" );
    aztec_conv = aztec_conv_list.value( aztec_conv_str );
    //
    // output
    //
    aztec_output_list.add( "all", AZ_all, "Print out everything including matrix" );
    aztec_output_list.add( "none", AZ_none,
                           "Print out no results (not even warnings)" );
    aztec_output_list.add( "last", AZ_last,
                           "Print out final residual and warnings" );
    aztec_output_list.add( "warnings", AZ_warnings,
                           "Print out only warning messages" );
    aztec_output_list.add( "all_res", 1,
                           "Print out all residuals and warnings (default)" );

    aztec_output_str = dfile( ( section + "/output" ).data(), "all_res" );
    aztec_output = aztec_output_list.value( aztec_output_str );
    aztec_pre_calc = dfile( ( section + "/pre_calc" ).data(), AZ_calc );
    aztec_max_iter = dfile( ( section + "/max_iter" ).data(), 500 );
    aztec_poly_ord = dfile( ( section + "/poly_ord" ).data(), 3 );
    aztec_overlap = dfile( ( section + "/overlap" ).data(), AZ_none );
    aztec_type_overlap = dfile( ( section + "/type_overlap" ).data(), AZ_standard );
    aztec_kspace = dfile( ( section + "/kspace" ).data(), 30 );
    aztec_orthog = dfile( ( section + "/orthog" ).data(), AZ_classic );
    aztec_aux_vec = dfile( ( section + "/aux_vec" ).data(), AZ_resid );
    aztec_reorder = dfile( ( section + "/reorder" ).data(), 1 );
    aztec_keep_info = dfile( ( section + "/keep_info" ).data(), 0 );
    //
    // subdomain solver
    //
    aztec_subdomain_solve_list.add( "bilu", AZ_bilu,
                                    "domain decomp with block ilu in subdomains" );
    aztec_subdomain_solve_list.add( "ilu", AZ_ilu,
                                    "domain decomp with ilu in subdomains" );
    aztec_subdomain_solve_list.add( "lu", AZ_lu,
                                    "domain decomp with lu in subdomains" );
    aztec_subdomain_solve_list.add( "icc", AZ_icc,
                                    "domain decomp with incomp Choleski in domains" );
    aztec_subdomain_solve_list.add( "ilut", AZ_ilut,
                                    "domain decomp with ilut in subdomains" );
    aztec_subdomain_solve_list.add( "rilu", AZ_rilu,
                                    "domain decomp with rilu in subdomains" );
    //
    aztec_subdomain_solve_str = dfile( ( section + "/subdomain_solve" ).data(),
                                       "ilut" );
    aztec_subdomain_solve = aztec_subdomain_solve_list.value( aztec_subdomain_solve_str );
    //
    aztec_graph_fill = dfile( ( section + "/graph_fill" ).data(), 0 );
    aztec_init_guess = dfile( ( section + "/init_guess" ).data(), AZ_NOT_ZERO );
    aztec_keep_kvecs = dfile( ( section + "/keep_kvecs" ).data(), 0 );
    aztec_apply_kvecs = dfile( ( section + "/apply_kvecs" ).data(), AZ_FALSE );
    aztec_orth_kvecs = dfile( ( section + "/orth_kvecs" ).data(), AZ_FALSE );
    aztec_ignore_scaling = dfile( ( section + "/ignore_scaling" ).data(), AZ_FALSE );
    aztec_check_update_size = dfile( ( section + "/check_update_size" ).data(),
                                     AZ_FALSE );
    //
    aztec_tol = dfile( ( section + "/tol" ).data(), 1.0e-06 );
    aztec_drop = dfile( ( section + "/drop" ).data(), 0.0 );
    aztec_ilut_fill = dfile( ( section + "/ilut_fill" ).data(), 1. );
    aztec_omega = dfile( ( section + "/omega" ).data(), 1. );
    aztec_update_reduction = dfile( ( section + "/update_reduction" ).data(), 10e10 );
}


void DataAztec::aztecOptionsFromDataFile( int* options, double* params )
{
    /*
      We first initialize by the defaults parameters.  Note that this step is
      not strictly necessary since the GetPot default parameters are
      precisely the parameters set by AZ_defaults. Nevertheless, we prefer to
      call this function, by security (if, in a new version of aztec, this
      function initialize new variables that are not consider at the moment,
      there will be no bad surprise, even if we do not update accordingly the
      GetPot default parameters).
    */
    AZ_defaults( options, params );

    // Next, we overload with the users parameters from the data file:

    options[ AZ_solver ] = aztec_solver;
    options[ AZ_scaling ] = aztec_scaling;
    options[ AZ_precond ] = aztec_precond;
    options[ AZ_conv ] = aztec_conv;
    options[ AZ_output ] = aztec_output;
    options[ AZ_pre_calc ] = aztec_pre_calc;
    options[ AZ_max_iter ] = aztec_max_iter;
    options[ AZ_poly_ord ] = aztec_poly_ord;
    options[ AZ_overlap ] = aztec_overlap;
    options[ AZ_type_overlap ] = aztec_type_overlap;
    options[ AZ_kspace ] = aztec_kspace;
    options[ AZ_orthog ] = aztec_orthog;
    options[ AZ_aux_vec ] = aztec_aux_vec;
    options[ AZ_reorder ] = aztec_reorder;
    options[ AZ_keep_info ] = aztec_keep_info;
    options[ AZ_subdomain_solve ] = aztec_subdomain_solve;
    options[ AZ_graph_fill ] = aztec_graph_fill;
    options[ AZ_init_guess ] = aztec_init_guess;
    options[ AZ_keep_kvecs ] = aztec_keep_kvecs;
    options[ AZ_apply_kvecs ] = aztec_apply_kvecs;
    options[ AZ_orth_kvecs ] = aztec_orth_kvecs;
    options[ AZ_ignore_scaling ] = aztec_ignore_scaling;
    options[ AZ_check_update_size ] = aztec_check_update_size;
    //
    params[ AZ_tol ] = aztec_tol;
    params[ AZ_drop ] = aztec_drop;
    params[ AZ_ilut_fill ] = aztec_ilut_fill;
    params[ AZ_omega ] = aztec_omega;
    params[ AZ_update_reduction ] = aztec_update_reduction;
}

void DataAztec::aztecSolveLinearSyst( MSRMatr<double>& mat,
                                      Real* unknown, Real* rhs,
                                      int unknown_size, MSRPatt& pattern )
{
    int proc_config[ AZ_PROC_SIZE ]; // Processor information:
    int options[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    int *data_org;                // Array to specify data layout
    double status[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    int *update;                  // vector elements updated on this node.
    int *external;                // vector elements needed by this node.
    int *update_index;            // ordering of update[] and external[]
    int *extern_index;            // locally on this processor.
    int N_update;                 // # of unknowns updated on this node
    AZ_set_proc_config( proc_config, AZ_NOT_MPI );
    AZ_read_update( &N_update, &update, proc_config, unknown_size, 1, AZ_linear );
    AZ_transform( proc_config, &external,
                  ( int * ) pattern.giveRaw_bindx(), mat.giveRaw_value(),
                  update, &update_index,
                  &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
                  AZ_MSR_MATRIX );
    // We initialize Aztec options and parameters from the data file:
    aztecOptionsFromDataFile( options, params );
    AZ_solve( unknown, rhs, options, params, NULL,
              ( int * ) pattern.giveRaw_bindx(), NULL, NULL, NULL,
              mat.giveRaw_value(), data_org, status, proc_config );
}

void DataAztec::aztecSolveLinearSyst( MSRMatr<double>& mat,
                                      Real* unknown, Real* rhs,
                                      int unknown_size, MSRPatt& pattern,
                                      int* options, double* params )
{
  int    proc_config[AZ_PROC_SIZE];// Processor information:
  int    *data_org;                // Array to specify data layout
  double status[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
  int    *update;                  // vector elements updated on this node.
  int    *external;                // vector elements needed by this node.
  int    *update_index;            // ordering of update[] and external[]
  int    *extern_index;            // locally on this processor.
  int    N_update;                 // # of unknowns updated on this node
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
  AZ_read_update(&N_update, &update, proc_config, unknown_size,1,AZ_linear);

  //  cout << "1) data_org[AZ_name] = " << data_org[AZ_name] << endl;
  AZ_transform(proc_config, &external,
	       (int *)pattern.giveRaw_bindx(), mat.giveRaw_value(),
	       update, &update_index,
	       &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
	       AZ_MSR_MATRIX);
  //  cout << "2) data_org[AZ_name] = " << data_org[AZ_name] << endl;
  // data_org[AZ_name] = 1;
  //Here, Aztec options and parameters are provided by the calling function
  AZ_solve(unknown, rhs, options, params, NULL,
	   (int *)pattern.giveRaw_bindx(), NULL, NULL, NULL,
	   mat.giveRaw_value(), data_org,status, proc_config);
  //cout << "3) data_org[AZ_name] = " << data_org[AZ_name] << endl;
}


void DataAztec::aztecSolveLinearSyst(MSRMatr<double>& mat,
				     Real* unknown,Real* rhs,
				     int unknown_size,MSRPatt& pattern,
				     int* options,double* params,
				     int& az_name,bool flag)
{
  /*
    if flag=true  : return az_name
    if flag=false : use the given az_name (to reuse previous information like
    preconditioner factorization)

    Typical use in the calling function:
    
    ...
    static bool flag = true;
    static int az_name;
    ...
    if(flag){
      options[AZ_keep_info] = 1;
    } else {
      options[AZ_keep_info] = 1;
      options[AZ_pre_calc] = AZ_reuse;
    }					     
    aztecSolveLinearSyst(mat,u.giveVec(),rhs.giveVec(),u.size(),
			 msrPattern,options,params,az_name,flag);  
    flag = false;
    
  */
  int    proc_config[AZ_PROC_SIZE];// Processor information:
  int    *data_org;                // Array to specify data layout
  double status[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
  int    *update;                  // vector elements updated on this node.
  int    *external;                // vector elements needed by this node.
  int    *update_index;            // ordering of update[] and external[]
  int    *extern_index;            // locally on this processor.
  int    N_update;                 // # of unknowns updated on this node
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
  AZ_read_update(&N_update, &update, proc_config, unknown_size,1,AZ_linear);

  AZ_transform(proc_config, &external,
	       (int *)pattern.giveRaw_bindx(), mat.giveRaw_value(),
	       update, &update_index,
	       &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
	       AZ_MSR_MATRIX);
  if(flag) az_name = data_org[AZ_name];
  else data_org[AZ_name]=az_name;
  //Here, Aztec options and parameters are provided by the calling function
  AZ_solve(unknown, rhs, options, params, NULL,
	   (int *)pattern.giveRaw_bindx(), NULL, NULL, NULL,
	   mat.giveRaw_value(), data_org,status, proc_config);
}

void DataAztec::dataAztecHelp( std::ostream& c )
{
    c << "\n*** Help for data [aztec]\n\n";
    c << std::endl;
    aztec_solver_list.showMe();
    c << std::endl;
    aztec_precond_list.showMe();
    c << std::endl;
    aztec_subdomain_solve_list.showMe();
    c << std::endl;
    c << "to be completed...\n";
}

void DataAztec::dataAztecShowMe( std::ostream& c )
{
    c << "\n*** Values for data [aztec]\n\n";
    c << "aztec_solver            = " << aztec_solver << " ("
    << aztec_solver_str << ")" << std::endl;
    c << "aztec_scaling           = " << aztec_scaling << " ("
    << aztec_scaling_str << ")" << std::endl;
    c << "aztec_precond           = " << aztec_precond << " ("
    << aztec_precond_str << ")" << std::endl;
    c << "aztec_tol               = " << aztec_tol << std::endl;
    c << "aztec_conv              = " << aztec_conv << " ("
    << aztec_conv_str << ")" << std::endl;
    c << "aztec_output            = " << aztec_output << " ("
    << aztec_output_str << ")" << std::endl;
    c << "aztec_pre_calc          = " << aztec_pre_calc << std::endl;
    c << "aztec_max_iter          = " << aztec_max_iter << std::endl;
    c << "aztec_poly_ord          = " << aztec_poly_ord << std::endl;
    c << "aztec_overlap           = " << aztec_overlap << std::endl;
    c << "aztec_type_overlap      = " << aztec_type_overlap << std::endl;
    c << "aztec_kspace            = " << aztec_kspace << std::endl;
    c << "aztec_orthog            = " << aztec_orthog << std::endl;
    c << "aztec_aux_vec           = " << aztec_aux_vec << std::endl;
    c << "aztec_reorder           = " << aztec_reorder << std::endl;
    c << "aztec_keep_info         = " << aztec_keep_info << std::endl;
    c << "aztec_subdomain_solve   = " << aztec_subdomain_solve
    << " (" << aztec_subdomain_solve_str << ")" << std::endl;
    c << "aztec_graph_fill        = " << aztec_graph_fill << std::endl;
    c << "aztec_init_guess        = " << aztec_init_guess << std::endl;
    c << "aztec_keep_kvecs        = " << aztec_keep_kvecs << std::endl;
    c << "aztec_apply_kvecs       = " << aztec_apply_kvecs << std::endl;
    c << "aztec_orth_kvecs        = " << aztec_orth_kvecs << std::endl;
    c << "aztec_ignore_scaling    = " << aztec_ignore_scaling << std::endl;
    c << "aztec_check_update_size = " << aztec_check_update_size << std::endl;
    c << "aztec_drop              = " << aztec_drop << std::endl;
    c << "aztec_ilut_fill         = " << aztec_ilut_fill << std::endl;
    c << "aztec_omega             = " << aztec_omega << std::endl;
    c << "aztec_update_reduction  = " << aztec_update_reduction << std::endl;
}
}
