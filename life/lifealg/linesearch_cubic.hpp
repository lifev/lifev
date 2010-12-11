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
    @brief Line-search algorithm with cubic interpolation.

    @author
    @date

    @contributor Alessandro Melani <alessandro.melani@mail.polimi.it>
    @mantainer Alessandro Melani <alessandro.melani@mail.polimi.it>

    @see Chapter 6 of J.E. Dennis, R.B.Schnabel
    "Numerical Methods for Unconstrained Optimization and Nonlinear Equations",
    no. 16 in Classics in Applied Mathematics, SIAM, Philadelphia, 1996

*/

#ifndef _LINESEARCH_CUBIC_H_
#define _LINESEARCH_CUBIC_H_ 1

#include <limits>

#undef min

namespace LifeV
{
/*!

    @author
    @date

    This line search algorithm comes from chapter 6 of J.E. Dennis, R.B.Schnabel
    "Numerical Methods for Unconstrained Optimization and Nonlinear Equations",
    no. 16 in Classics in Applied Mathematics, SIAM, Philadelphia, 1996

    (i)   lambda given (usually 1 when Newton method is used)
    (ii)  solTest = sol + lambda stepTest
    (iii) Goldstein - Price + cubic interpolation

    sigma0, sigma1:     safeguarding bounds (default values 0.1 and 0.5)
    m1, m2                    parameters in Goldstein conditions
    maxIterations:      maximum number of steplength reductions before
    failure is reported (default 50)

    @param f          Function
    @param residual    Residual
    @param sol        Solution
    @param step        Step to update the solution
    @param normRes     Norm of the residual
    @param lambda      Length of the Step
    @param slope       Slope value in linesearch algorithm
    @param iter        Iterations
    @param verbose     Option for detailed description
*/

template <class Fct, class VectorType>
Int lineSearch_cubic( Fct& f, VectorType& residual, VectorType& sol, VectorType& step,
                      Real& normRes, Real& lambda, Real slope, UInt iter, bool const verbose = true )
{

    //----------------------------------------------------------------------
    if (verbose)
        std::cout << "Cubic line search ..." << std::endl;

    const Real sigma0 = 0.1;
    const Real sigma1 = 0.5;
    const Real m1 = 0.25;
    const Real m2 = 0.75;
    const Int maxIterations = 50;
    //----------------------------------------------------------------------
    static VectorType S_solCurrent = sol; // static to avoid too many alloc/de-alloc
    Int iterLinesearch;
    bool firstTime = true;
    Real lambda2, lambdaOld, lambdaOld2, lambdaTemporary,
    normResTest, f0, ftest, fold, c, c11, c12, c21, c22, a, b, disc, g1, g2, gprev = 0;
    //
    f0 = 0.5 * normRes * normRes;
    lambdaOld = lambda;
    S_solCurrent = sol;
    sol += lambda * step;
    f.evalResidual( residual, sol, iter );
    normResTest = residual.NormInf();
    ftest = 0.5 * normResTest * normResTest;
    fold = ftest;
    iterLinesearch = 0;
    while ( ftest < f0 + m2 * slope * lambda // lambda is too small: extrapolation
            && iterLinesearch < maxIterations )
    {
        iterLinesearch++;
        lambda *= 2;
        sol =  S_solCurrent;
        sol += lambda * step;
        if (verbose)
            std::cout << "--- line search (extrapolation, Goldstein rule)" << std::endl;
        f.evalResidual( residual, sol, iter );

        if (verbose)
            std::cout << "    line search iter : " << iterLinesearch << " residual test = "
                      << normResTest << ", lambda = " << lambda << std::endl;
        normResTest = residual.NormInf();
        ftest = 0.5 * normResTest * normResTest;
    }
    if ( iterLinesearch == maxIterations )
    {
        if (verbose)
            std::cout << "line search: too many extrapolations" << std::endl;
        return EXIT_FAILURE;
    }
    lambdaOld = lambda;
    while ( ftest > f0 + m1 * slope * lambda // Armijo's rule: lambda is too large
            && iterLinesearch < maxIterations )
    {
        iterLinesearch++;
        //-- cubic interpolation of lambda:
        lambda2 = lambda * lambda;
        g1 = ftest - f0 - slope * lambda;
        if ( firstTime )
        {
            lambdaTemporary = - slope * lambda2 / ( 2. * g1 );
            firstTime = false;
        }
        else
        {
            lambdaOld2 = lambdaOld * lambdaOld;
            g2 = gprev - f0 - lambdaOld * slope;
            c = 1. / ( lambda - lambdaOld );
            c11 = 1. / lambda2;
            c12 = -1. / lambdaOld2;
            c21 = - lambdaOld / lambda2;
            c22 = lambda / lambdaOld2;
            a = c * ( c11 * g1 + c12 * g2 );
            b = c * ( c21 * g1 + c22 * g2 );
            disc = b * b - 3. * a * slope;
            if ( ( fabs( a ) > std::numeric_limits<Real>::min() ) &&
                    ( disc > std::numeric_limits<Real>::min() )
               )
                lambdaTemporary = ( - b + sqrt( disc ) ) / ( 3. * a );
            else
                lambdaTemporary = slope * lambda2 / ( 2. * g1 ) ;
            if ( lambdaTemporary >= sigma1 * lambda )
                lambdaTemporary = sigma1 * lambda;
        }
        lambdaOld = lambda;
        gprev = ftest;
        if ( lambdaTemporary < sigma0 * lambda )
            lambda *= sigma0;
        else
            lambda = lambdaTemporary;
        //--
        sol =  S_solCurrent;
        sol += lambda * step;
        if (verbose)
            std::cout << "--- line search (cubic interpolation, Armijo rule)" << std::endl;
        f.evalResidual( residual, sol, iter );
        normResTest = residual.NormInf();
        if (verbose)
            std::cout << "    line search iter : " << iterLinesearch << " residual test = "
                      << normResTest << ", lambda = " << lambda << std::endl;
        ftest = 0.5 * normResTest * normResTest;
    }
    if ( iterLinesearch == maxIterations )
    {
        if (verbose)
            std::cout << "line search: too many interpolations" << std::endl;
        return EXIT_FAILURE;
    }
    normRes = normResTest;
    if (verbose)
        std::cout << "Parabolic line search: final residual = " << normRes << std::endl;
    return EXIT_SUCCESS;
}

} // Namespace LifeV

#endif /* _LINESEARCH_CUBIC_H_ */
