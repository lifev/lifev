#ifndef _DATAAZTEC_H_
#define _DATAAZTEC_H_
#include <string>
#include <iostream>

//! Aztec include
#include "az_aztec.h"


#include "GetPot.hpp"
#include "values.hpp"
#include "dataString.hpp"

/*!
  \file dataAztec.h
  \author J.-F. Gerbeau
  \date 11/2002
        01/2003 M.A. Fernandez, GetPot sections as parameter
  
  \brief To use Aztec with GetPot

  All the parameters and options can be selected from a GetPot
  file. Some of them (solvers, preconditionner,...) can be
  given with a string (e.g.: gmres, ilut, ...)

  To see the values (including the default ones): dataAztecShowMe()

  To see all the possible choices:  dataAztecHelp() 
  
  \todo dataAztecHelp() is not complete, and other items could
  be named with a string rather than with an integer
*/
using namespace std;

class DataAztec
{
public:
  // Aztec Options
  int aztec_solver;
  string aztec_solver_str;
  int aztec_scaling;
  string aztec_scaling_str;
  int aztec_precond;
  string aztec_precond_str;
  int aztec_conv;
  string aztec_conv_str;
  int aztec_output;
  string aztec_output_str;
  int aztec_pre_calc;
  int aztec_max_iter;
  int aztec_poly_ord;
  int aztec_overlap;
  int aztec_type_overlap;
  int aztec_kspace;
  int aztec_orthog;
  int aztec_aux_vec;
  int aztec_reorder;
  int aztec_keep_info;
  int aztec_subdomain_solve;
  string aztec_subdomain_solve_str;
  int aztec_graph_fill;
  int aztec_init_guess;
  int aztec_keep_kvecs;
  int aztec_apply_kvecs;
  int aztec_orth_kvecs;
  int aztec_ignore_scaling;
  int aztec_check_update_size;
  // Aztec Parameters
  double aztec_tol;
  double aztec_drop;
  double aztec_ilut_fill;
  double aztec_omega;
  double aztec_update_reduction;
  //--------------------------
  DataStringList aztec_solver_list;
  DataStringList aztec_precond_list;
  DataStringList aztec_scaling_list;
  DataStringList aztec_conv_list;
  DataStringList aztec_output_list;
  DataStringList aztec_subdomain_solve_list;
  // 
  DataAztec(const GetPot& dfile, const string& section="aztec");
 
  void aztecOptionsFromDataFile(int* option,double* param);
  /*! solve a linear system with the parameters given in the data file */
  void aztecSolveLinearSyst(MSRMatr<double>& mat,Real* unknown,Real* rhs,
			    int unknown_size,MSRPatt& pattern); 
  /*! solve a linear system with the parameters given by the users and the data file */
  void aztecSolveLinearSyst(MSRMatr<double>& mat,
			    Real* unknown,Real* rhs,int unknown_size,
			    MSRPatt& pattern,int* options,double* params);
  /*! to see a little help */
  void dataAztecHelp(ostream& c=cout);
  /*! to see the values of the items (including the default ones) */
  void dataAztecShowMe(ostream& c=cout);
};
#endif
