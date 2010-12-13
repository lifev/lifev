//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Preconditioned relaxed solver for non linear problems.

    @author
    @date 01-01-2004

    Given a functional this method tries to find a root near an
    initial guess sol. The two methods that the functional must
    have are evalRes and solveJac. If solveJac inverts the exact
    Jacobian of the functional, than this is equivalent to a Newton
    method. If convergence do not appens in the expected ratio,
    a linesearch algorithm is called (uadratic or parabolic).

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>


 */

#ifndef _NONLINRICHARDSON_HPP
#define _NONLINRICHARDSON_HPP

#include <algorithm> // for min and max
#include <life/lifealg/linesearch_parabolic.hpp>
#include <life/lifealg/linesearch_cubic.hpp>
#include <life/lifealg/generalizedAitken.hpp>
#include <life/lifearray/EpetraVector.hpp>

namespace LifeV
{
//! Preconditioned relaxed solver for non linear problems.
   /*!
       Add more details about the constructor.
       NOTE: short description is automatically added before this part.
       @param sol            :  the solution
       @param maxit          :  input: maximum iterations, output: nb of iterations
       @param abstol, reltol :  the stoping criteria is abstol+reltol*norm(residual_0),
       @param eta_max        :  Maximum error tolerance for residual in linear solver.

       The linear solver terminates when the relative
       linear residual is smaller than eta*| f(sol) |.

       The value linear_rel_tol send for the relative tolerance
       to the linear solver is therefore eta. eta is determined
       by the modified Eisenstat-Walker formula if etamax > 0.

       @param linesearch     :  for now consider only the case linesearch=0
         (coded but not theoretically analysed)

       @param omega          :  default relaxation parameter to be passed to Aitken.
       if omega is negative, then its absolute value is
       taken as constant relaxation parameter

    */

template < class Fct >
Int nonLinRichardson( EpetraVector& sol,
                      Fct&        functional,
                      Real        abstol,
                      Real        reltol,
                      UInt&       maxit,
                      Real        eta_max,
                      Int         linesearch,
                      std::ofstream& out_res,
                      const Real& time,
                      UInt iter = UInt(0) )
{
    /*
        */

    /*
      max_increase_res: maximum number of successive increases in residual
      before failure is reported
    */

//    const Int max_increase_res = 5;

    /*
      Parameters for the linear solver, gamma: Default value = 0.9
    */

    const Real gamma   = 0.9;

    //----------------------------------------------------------------------

    bool const verbose(sol.Comm().MyPID() == 0);

    //UInt iter = 0;

    EpetraVector residual ( sol.getMap() );
    EpetraVector step     ( sol.getMap() );

    step *= 0.;

    Real normResOld = 1;

    if (verbose)
    {
        //std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        std::cout << "      Non-Linear Richardson: starting          " << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        std::cout << std::endl;
        //std::cout << "------------------------------------------------------------------" << std::endl;
    }
    functional.evalResidual( residual, sol, iter );

    Real normRes      = residual.NormInf();
    Real stop_tol     = abstol + reltol*normRes;
    Real linearRelTol = fabs(eta_max);
    Real eta_old;
    Real eta_new;
    Real ratio;
    Real slope;
    Real linres;
    Real lambda;

    //

    Real solNormInf(sol.NormInf());
    Real stepNormInf;
    if (verbose)
    {
        out_res << std::scientific;
        out_res << "# time = ";
        out_res << time << "   " << "initial norm_res " <<  normRes
        << " stop tol = " << stop_tol
        << "initial norm_sol "
        << solNormInf << std::endl;
        out_res << "#iter      disp_norm       step_norm       residual_norm" << std::endl;
    }
    while ( normRes > stop_tol && iter < maxit )
    {
        if (verbose)
        {
            std::cout << std::endl;
            //std::cout << "------------------------------------------------------------------" << std::endl;
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            std::cout << "      Non-Linear Richardson: iteration  =      " << iter << std::endl
                      << "                             residual   =      " << normRes << std::endl
                      << "                             tolerance  =      " << stop_tol << std::endl;
            std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            std::cout << std::endl;
            //std::cout << "------------------------------------------------------------------" << std::endl;
            //std::cout << std::endl;
        }

        iter++;

        ratio      = normRes/normResOld;
        normResOld = normRes;
        normRes    = residual.NormInf();

        residual *= -1;
        functional.solveJac(step, residual, linearRelTol); // J*step = -R

        solNormInf = sol.NormInf();
        stepNormInf = step.NormInf();
        if (verbose)
        {
            out_res   << std::setw(5) << iter
            << std::setw(15) << solNormInf
            << std::setw(15) << stepNormInf;
        }
        linres = linearRelTol;

        lambda = 1.;
        slope  = normRes * normRes * ( linres * linres - 1 );

        Int status(EXIT_SUCCESS);
        switch ( linesearch )
        {
        case 0: // no linesearch
            sol += step;
            functional.evalResidual( residual, sol, iter);
//                normRes = residual.NormInf();
            break;
        case 1:
            status = lineSearch_parab( functional, residual, sol, step, normRes, lambda, iter, verbose );
            break;
        case 2:  // recommended
            status = lineSearch_cubic( functional, residual, sol, step, normRes, lambda, slope, iter, verbose );
            break;
        default:
            std::cout << "Unknown linesearch \n";
            status = EXIT_FAILURE;
        }

        if (status == EXIT_FAILURE)
            return status;



        normRes = residual.NormInf();

        if (verbose)
            out_res << std::setw(15) << normRes << std::endl;

        if ( eta_max > 0 )
        {
            eta_old = linearRelTol;
            eta_new = gamma * ratio * ratio;
            if ( gamma * eta_old * eta_old > .1 )
            {
                eta_new = std::max<Real>( eta_new, gamma * eta_old * eta_old );
            }
            linearRelTol = std::min<Real>( eta_new, eta_max );
            linearRelTol = std::min<Real>( eta_max,
                                           std::max<Real>( linearRelTol,
                                                           .5 * stop_tol / normRes ) );
            //if (verbose)
            //    std::cout << "    Newton: forcing term eta = " << linearRelTol << std::endl;
        }

    }

    if ( normRes > stop_tol )
    {
        if (verbose)
            std::cout << "!!! NonLinRichardson: convergence fails" << std::endl;
        maxit = iter;
        return EXIT_FAILURE;
    }

    if (verbose)
    {
        std::cout << std::endl;
        //std::cout << "------------------------------------------------------------------" << std::endl;
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        std::cout << "      Non-Linear Richardson: convergence =     " << normRes << std::endl
                  << "                             iterations  =     " << iter << std::endl;
        std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        std::cout << std::endl;
        //std::cout << "------------------------------------------------------------------" << std::endl;
    }
    maxit = iter;

    return EXIT_SUCCESS;
}
}
#endif
