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
#include "regionMesh3D_ALE.hpp"
#ifndef _OPERFS
#define _OPERFS

namespace LifeV
{
  class operFS;
  
  
  class DataJacobian {
    
  public:
    
    DataJacobian(operFS* oper):
      _pFS(oper){}
    
    operFS* _pFS;
    
  };
  
  void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);
  
  
  // Fluid-Structure operator Class
  class operFS {
    
  public:
    
    operFS(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
	   VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid,
	   BCHandler& BCh_dp, BCHandler& BCh_dz, const GetPot& data_file);
    
    //
    void eval(Vector& dispNew, Vector& veloStruct, const Vector& disp,int status);
    
    //
    void evalResidual(Vector &res, const Vector& sol, int iter);
    
    //
    void updateJac(Vector& sol,int iter);
    
    //
    void solveJac(Vector &step, const Vector& res, double& linear_rel_tol);
    
    //
    void solveReducedLinearFluid();
    
    //
    void solveLinearSolid();
    
    //
    UInt nbEval();
    
    void setTime(const Real& time);
    
    Vector& da() {return _da;}
    
    Vector& minusdp() {return _minusdp;}
    
    friend void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);

  private:

    //! fluid and solid solvers
    NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& _fluid;
    VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& _solid;
    
    //! aux vectors
    Vector _dispStruct, _velo;
    
    //! displacement variation
    Vector _dz;
    
    //! righthand side of the solid tangent problem
    Vector _rhs_dz;
    
    //! number of residual evaluations
    UInt _nbEval;
    
    //! linearized displaceemnt BC's 
    BCHandler& _BCh_dz;
    
    //! data for the coupled tangent problem
    DataJacobian _dataJacobian;

    //! present time
    Real _time;

    //! linearized acceleration
    Vector _da;

    //
    // laplacian staff
    //
    
    //! reference FE
    const RefFE& _refFE;

    //! volumic quadrature rule
    const QuadRule& _Qr;

    //! boundary quadrature ru;e
    const QuadRule& _bdQr;

    //! BC for pressure
    BCHandler& _BCh_dp;

    //! dof 
    Dof _dof;
    
    //! number of total dof 
    UInt _dim;
    
    //! pattern of C
    MSRPatt _pattC;
    
    //! global matrix
    MSRMatr<double> _C, _CAux;
    
    //! Current FE
    CurrentFE _fe;
    CurrentBdFE _feBd;

    //! elementary matrix and vector
    ElemMat _elmatC; 

    //! unknown: presure variation
    Vector _dp, _minusdp;

    //! right handside 
    Vector _f;

    //! linear solver
    SolverAztec _linearSolver;

    //! computed matrix C?
    bool _computedC;

};

void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);
}
#endif
