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
#ifndef _NONLINRICHARDSON_HPP
#define _NONLINRICHARDSON_HPP

#include <algorithm> // for min and max
#include "linesearch_parabolic.hpp"
#include "linesearch_cubic.hpp"
#include "generalizedAitken.hpp"

namespace LifeV
{
    template<class Fct,class Vector,class Real, class Norm>
    int nonLinRichardson(Vector&     sol,
                         Fct&        f,
                         Norm&       norm,
                         Real        abstol,
                         Real        reltol,
                         int&        maxit,
                         Real        eta_max,
                         int         linesearch=0,
                         std::ofstream&   out_res,
                         const Real& time,
                         const Real  omega)
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

          linesearch     :  for now consider only the case linesearch=0
          (coded but not theoretically analysed)

          omega          :  default relaxation parameter to be passed to Aitken.
          if omega is negative, then its absolute value is
          taken as constant relaxation parameter
        */

        /*
          max_increase_res: maximum number of successive increases in residual
          before failure is reported
        */
//        const int max_increase_res = 5;
        /*
          Parameters for the linear solver, gamma: Default value = 0.9
        */
//        const Real gamma   = 0.9;

//----------------------------------------------------------------------

        Real linres = 1.e-3;

        int    iter         = 0;

        int    nDofFS       = sol.size();

        Vector residual     = sol;
        Vector step         = sol;

        Vector muS(f.solid().dDof().numTotalDof());
        Vector muF(f.fluid().uDof().numTotalDof());

        muS                 = -2.;
        muF                 = -1.;
        step                = 0.;

        Real   normResOld   = 1;

        Real   omegaS       = omega;
        Real   omegaF       = omega;

        f.evalResidual(sol, iter, residual);

        Real normRes        = norm(residual);
        Real stop_tol       = abstol + reltol*normRes;
        Real linear_rel_tol = fabs(eta_max);

        generalizedAitken<Vector,Real> aitken(nDofFS, omegaS, omegaF);

//

        out_res << time << "    " << iter << "   " << normRes << std::endl;
        

        while( normRes > stop_tol && iter < maxit)
        {

            std::cout << std::endl;
            std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << "  NonLinRichardson: residual = " << normRes
                      << ", stoping tolerance = "          << stop_tol << std::endl;
            std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << std::endl;
            
            iter++;

            normResOld = normRes;
            normRes    = norm(residual);

//            muS        = residual;
            
            step       = aitken.computeDeltaLambda(sol, residual);

            std::cout << "Step norm = " << norm(step) << std::endl;

            out_res << "Step norm = " << norm(step);
            out_res << "size d = " << f.solid().dDof().numTotalDof() << std::endl;
            out_res << "size f = " << f.fluid().uDof().numTotalDof() << std::endl;
            
            muS        = f.residualS();
            muF        = f.residualF();

            sol       += step;
            
            f.evalResidual(sol, iter, residual);

            out_res << " Res norm = " << norm(residual) << std::endl;
            
            normRes    = norm(residual);

        }

        if(normRes > stop_tol)
        {
            std::cout << "!!! NonLinRichardson: convergence fails" << std::endl;
            maxit = iter;
            return 1;
        }

        std::cout << "--- NonLinRichardson: convergence in " << iter << " iterations\n\n";

        maxit = iter;

        return 0;
    }
}
#endif


