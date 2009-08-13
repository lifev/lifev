/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
#ifndef _LINESEARCH_CUBIC_H_
#define _LINESEARCH_CUBIC_H_

#include <limits>

#undef min

namespace LifeV
{
/*
  This line search algorithm comes from Dennis & Schnabel

  (i)   lambda given (usually 1 when Newton method is used)
  (ii)  sol_test = sol + lambda step_test
  (iii) Goldstein - Price + cubic interpolation

  sigma0, sigma1: safeguarding bounds (default values 0.1 and 0.5)
  m1, m2
  max_linesearch: maximum number of steplength reductions before
  failure is reported (default 50)


*/

template <class Fct, class Vector, class Real, class Norm>
void lineSearch_cubic( Fct& f, Norm& /*norm*/, Vector& residual, Vector& sol, Vector& step,
                       Real& normRes, Real& lambda, Real slope, UInt iter )
{
    //----------------------------------------------------------------------

    std::cout << "Cubic line search ..." << std::endl;

    const Real sigma0 = 0.1;
    const Real sigma1 = 0.5;
    const Real m1 = 0.25;
    const Real m2 = 0.75;
    const int max_linesearch = 50;
    //----------------------------------------------------------------------
    static Vector sol_cur = sol; // static to avoid too many alloc/de-alloc
    int iter_linesearch;
    bool first_time = true;
    Real lambda2, lambda_old, lambda_old2, lambda_tmp,
    normRes_test, f0, ftest, fold, c, c11, c12, c21, c22, a, b, disc, g1, g2, gprev = 0;
    //
    f0 = 0.5 * normRes * normRes;
    lambda_old = lambda;
    sol_cur = sol;
    sol += lambda * step;
    f.evalResidual( residual, sol, iter );
//    f.evalResidual( sol, iter, residual );
    normRes_test = residual.NormInf();
    ftest = 0.5 * normRes_test * normRes_test;
    fold = ftest;
    iter_linesearch = 0;
    while ( ftest < f0 + m2 * slope * lambda // lambda is too small: extrapolation
            && iter_linesearch < max_linesearch )
    {
        iter_linesearch++;
        lambda *= 2;
        sol =  sol_cur;
        sol += lambda * step;
        std::cout << "--- line search (extrapolation, Goldstein rule)" << std::endl;
        f.evalResidual( residual, sol, iter );
//        f.evalResidual( sol, iter, residual );
        std::cout << "    line search iter : " << iter_linesearch << " residual test = "
        << normRes_test << ", lambda = " << lambda << std::endl;
        normRes_test = residual.NormInf();
        ftest = 0.5 * normRes_test * normRes_test;
    }
    if ( iter_linesearch == max_linesearch )
    {
        std::cout << "line search: too many extrapolations" << std::endl;
        exit( 1 );
    }
    lambda_old = lambda;
    while ( ftest > f0 + m1 * slope * lambda // Armijo's rule: lambda is too large
            && iter_linesearch < max_linesearch )
    {
        iter_linesearch++;
        //-- cubic interpolation of lambda:
        lambda2 = lambda * lambda;
        g1 = ftest - f0 - slope * lambda;
        if ( first_time )
        {
            lambda_tmp = - slope * lambda2 / ( 2. * g1 );
            first_time = false;
        }
        else
        {
            lambda_old2 = lambda_old * lambda_old;
            g2 = gprev - f0 - lambda_old * slope;
            c = 1. / ( lambda - lambda_old );
            c11 = 1. / lambda2;
            c12 = -1. / lambda_old2;
            c21 = - lambda_old / lambda2;
            c22 = lambda / lambda_old2;
            a = c * ( c11 * g1 + c12 * g2 );
            b = c * ( c21 * g1 + c22 * g2 );
            disc = b * b - 3. * a * slope;
            if ( ( fabs( a ) > std::numeric_limits<Real>::min() ) &&
                    ( disc > std::numeric_limits<Real>::min() )
               )
                lambda_tmp = ( - b + sqrt( disc ) ) / ( 3. * a );
            else
                lambda_tmp = slope * lambda2 / ( 2. * g1 ) ;
            if ( lambda_tmp >= sigma1 * lambda )
                lambda_tmp = sigma1 * lambda;
        }
        lambda_old = lambda;
        gprev = ftest;
        if ( lambda_tmp < sigma0 * lambda )
            lambda *= sigma0;
        else
            lambda = lambda_tmp;
        //--
        sol =  sol_cur;
        sol += lambda * step;
        std::cout << "--- line search (cubic interpolation, Armijo rule)" << std::endl;
//        f.evalResidual( sol, iter, residual );
        f.evalResidual( residual, sol, iter );
        normRes_test = residual.NormInf();
        std::cout << "    line search iter : " << iter_linesearch << " residual test = "
        << normRes_test << ", lambda = " << lambda << std::endl;
        ftest = 0.5 * normRes_test * normRes_test;
    }
    if ( iter_linesearch == max_linesearch )
    {
        std::cout << "line search: too many interpolations" << std::endl;
        exit( 1 );
    }
    normRes = normRes_test;
    std::cout << "ok." << std::endl;

}
}
#endif
