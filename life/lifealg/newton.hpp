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
#ifndef _NEWTON_H
#define _NEWTON_H

#include <algorithm> // for min and max
#include "linesearch_parabolic.hpp"
#include "linesearch_cubic.hpp"

namespace LifeV
{
template<class Fct,class Vector,class Real, class Norm>
int newton(Vector& sol,Fct& f,Norm& norm,Real abstol,Real reltol,int& maxit,
	   Real eta_max,int linesearch=0, ofstream& out_res, const Real& time)
{
  /*
    sol            :  the solution
    maxit          :  input: maximum iterations, output: nb of iterations
    abstol, reltol :  the stoping criteria is abstol+reltol*norm(residual_0),
    eta_max        :  Maximum error tolerance for residual in linear solver.
                      The linear solver terminates when the relative
		      linear residual is smaller than eta*| f(sol) |.
		      The value linear_rel_tol send for the relative tolerance
		      to the linear solver is therefore eta. eta is determined
		      by the modified Eisenstat-Walker formula if etamax > 0.
		      If eta_max < 0, then eta = |etamax| for the entire
		      iteration (e.g. etamax = -1e-6 ensures that the linear
		      tolerance would be always 1e-6). Default value = 0.9
  */

  /*
    max_increase_res: maximum number of successive increases in residual
                      before failure is reported
  */
  const int max_increase_res=5;
  /*
     Parameters for the linear solver, gamma: Default value = 0.9
  */
  const Real gamma   = 0.9;
  //----------------------------------------------------------------------
  Real linres;

  int iter=0,increase_res=0;
  Vector residual=sol,step=sol;
  step=0.;
  Real  normResOld=1,lambda,slope;
  f.evalResidual(residual,sol,iter);
  Real normRes = norm(residual),normStep=0,
    stop_tol = abstol + reltol*normRes,
    ratio;
  Real eta_old,eta_new,linear_rel_tol=fabs(eta_max);
  //
  cout << "------------------------------------------------------------------"
       << endl;
  cout << "    Newton 0: residual=" << normRes << ", stoping tolerance = "
       << stop_tol << endl;
  cout << "------------------------------------------------------------------"
       << endl;

  out_res << time << "    " << iter << "   " << normRes << endl;

  while( normRes > stop_tol && iter < maxit){
    iter++;
    ratio = normRes / normResOld;
    normResOld = normRes;
    normRes = norm(residual);
    f.updateJac(sol,iter);
    linres = linear_rel_tol;
    f.solveJac(step,-1.*residual,linres); // residual = f(sol)
    /*
      linres contains the relative linear tolerance achieved by the
      linear solver, i.e linear_rel_tol = | -f(sol) - J step | / |-f(sol)|
    */
    slope = normRes*normRes*(linres*linres - 1);
    cout << "### slope = " << slope << endl;
    /*
      slope denotes the quantity f^T J step, which is generally used by
      line search algorithms. This formula comes form Brown & Saad (1990),
      formula (3.4). BE CAREFUL: it assumes that the linear solver is GMRES,
      with zero as initial guess (in particular it does not work with restart
      gmres, see formula (3.7) and (3.9) of Brown & Saad (1990))
    */
    lambda = 1.;
    //
    // -- line search
    //
    switch(linesearch){
    case 0:// no linesearch
      sol += step;
      f.evalResidual(residual,sol,iter);
      normRes = norm(residual);
      break;
    case 1:
      lineSearch_parab(f,norm,residual,sol,step,normRes,lambda,iter);
      break;
    case 2: // recommended
      lineSearch_cubic(f,norm,residual,sol,step,normRes,lambda,slope,iter);
      break;
    default:
      cout << "Unknown linesearch \n";
      exit(1);
    }
    //
    //-- end of line search
    //
    normStep = lambda*norm(step);
    ratio = normRes/normResOld;
    if(ratio > 1){
      increase_res ++;
      cout << "!!! Newton warning: increase in residual \n";
      if(increase_res == max_increase_res){
	cout << "!!! Newton:" << max_increase_res
	     << " consecutive increases in residual" << endl;
	maxit = iter;
	return 1;
      }
    } else {
      increase_res=0;
    }
    cout << "------------------------------------------------------------------"
	 << endl;
    cout << "    Newton " << iter << ": residual=" << normRes << ",  step="
	 << normStep << endl;
    cout << "------------------------------------------------------------------"
	 << endl;

    out_res << time << "    " << iter << "   " << normRes << endl;

    //
    //-- forcing term computation (Eisenstat-Walker)
    //
    if (eta_max>0){
      eta_old=linear_rel_tol;
      eta_new=gamma*ratio*ratio;
      if(gamma*eta_old*eta_old>.1) eta_new=max(eta_new,gamma*eta_old*eta_old);
      linear_rel_tol=min(eta_new,eta_max);
      linear_rel_tol=min(eta_max,max(linear_rel_tol,.5*stop_tol/normRes));
      cout <<"    Newton: forcing term eta = " << linear_rel_tol << endl;
    }
    //
    //-- end of forcing term computation
    //
  }
  if(normRes > stop_tol){
    cout << "!!! Newton: convergence fails" << endl;
    maxit = iter;
    return 1;
  }
  cout << "--- Newton: convergence in " << iter << " iterations\n\n";
  maxit = iter;
  return 0;
}
}
#endif


