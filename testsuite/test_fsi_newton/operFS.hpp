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
#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "norm.hpp"

#ifndef _OPERFS
#define _OPERFS


class operFS;


class DataJacobian {

 public:

  DataJacobian(operFS* oper):
    _pFS(oper){}

  operFS* _pFS;

};


//
// Fluid-Structure operator Class
//
class operFS {

 public:

  operFS(NavierStokesAleSolverPC< RegionMesh3D<LinearTetra> >& fluid, 
	 VenantKirchhofSolver< RegionMesh3D<LinearTetra> >& solid,
	 BC_Handler& BCh_du, BC_Handler& BCh_dz);
  
  //
  void eval(Vector& dispNew, Vector& veloStruct, const Vector& disp,int status);

  // 
  void evalResidual(Vector&res, const Vector& sol, int iter);

  //  
  void updateJac(Vector& sol,int iter);

  //
  void solveJac(Vector& step, const Vector& res, double& linear_rel_tol);

  //
  void solveLinearFluid();
  
  //
  void solveLinearSolid();

  //
  UInt nbEval(); 


  NavierStokesAleSolverPC< RegionMesh3D<LinearTetra> >& _fluid;
  
  VenantKirchhofSolver< RegionMesh3D<LinearTetra> >& _solid;
  
  Vector _dispStruct;
  Vector _velo; 
  Vector _dz;
  Vector _rhs_dz;

  UInt _nbEval;
  
  BC_Handler& _BCh_du;
  BC_Handler& _BCh_dz;
  
  DataJacobian _dataJacobian;

};

void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);

#endif
