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
#ifndef _FHN_SOLVER_H
#define _FHN_SOLVER_H
#include "fhnHandler.hpp"
#include "user_fct.hpp"
#include "user_functor.hpp"

namespace LifeV
{
/*!
  \brief A simple Fitzhugh-Nagumo solver
  \file fhnSolver.hpp
  \author J.-F. Gerbeau
  \date 09/2004
*/

class FhNSolver:
  public FhNHandler
{
public:
  FhNSolver(const GetPot& data_file);
  MSRPatt msrPattern;
  MSRMatr<double> mat;
  GenericVecHdl<Vector> rhs;
  GenericVecHdl<Vector> uv;
  GenericVecHdl<Vector> uv_1;
  ElemVec elvec;//
  ElemVec elvec_uv;//
  ElemMat elmat;// 

  vector<UInt> pts_proc;
  Vol_source vol_source;
  KNM<double> uv_ref;   // reference solution (read on a file)
  KNM<double> uv_stored; // current solution stored for comparison with reference 
  KN<double> ref_time;  // times at which the solution is stored 
  
  void compute_mat();
  void compute_rhs();
  void time_advance();
  void initial_data();
  void solve();
  void post_process();
  void store_solution();
  void save_stored_solution();
  bool flag_ref_solution;
  void compare_stored_with_ref_sol();
  void L2_and_H1_norms(const KN<double>& vec,double& L2_2,double& H1_2);
  void write_adapt();
  
};
}
#endif
