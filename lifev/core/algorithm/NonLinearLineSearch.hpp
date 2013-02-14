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
    @brief Line-search algorithm with parabolic interpolation.

    @author
    @date

    @contributor Alessandro Melani <alessandro.melani@mail.polimi.it>
    @mantainer Alessandro Melani <alessandro.melani@mail.polimi.it>

    @see Chapter 8 of C.T. Kelley,
    "Iterative methods for linear and nonlinear equations", SIAM 1995

  @see Chapter 6 of J.E. Dennis, R.B.Schnabel
    "Numerical Methods for Unconstrained Optimization and Nonlinear Equations",
    no. 16 in Classics in Applied Mathematics, SIAM, Philadelphia, 1996

*/

#ifndef _NonLinearLineSearch_H_
#define _NonLinearLineSearch_H_ 1

namespace LifeV
{
//! Implementation of Line Search method with parabolic interpolation
/*!

    @author
    @date

    This line search algorithm comes from chapter 8 of C.T. Kelley,
    "Iterative methods for linear and nonlinear equations", SIAM 1995

    (i)   lambda given (usually 1 when Newton method is used)
    (ii)  solTest = sol + lambda stepTest
    (iii) if residuTest < (1 - alpha lambda) residu
          then sol = solTest
    else choose lambda via a three point parabolic interpolation
        (that does not net the derivative) and apply a safeguarding
        step: if lambda < sigma0 lambdaCurrent then lambda = sigma0 lambdaCurrent
                     if lambda > sigma1 lambdaCurrent then lambda = sigma1 lambdaCurrent

    Constant parameters:

    sigma0, sigma1:     safeguarding bounds (default values 0.1 and 0.5)
    alpha:                     parameter to measure sufficient decrease (default 1e-4)
    maxIterations:      maximum number of steplength reductions before
                                   failure is reported (default 50)

    @param f            Function
    @param residual     Residual
    @param sol          Solution
    @param step         Step to update the solution
    @param normRes      Norm of the residual
    @param lambda       Length of the Step
    @param iter         Iterations
    @param verboseLevel Option for detailed description
*/
template <class Fct, class VectorType>
Int NonLinearLineSearchParabolic ( Fct& f, VectorType& residual, VectorType& sol, VectorType& step, Real& normRes,
                                   Real& lambda, UInt iter, UInt const verboseLevel = 1)
{

    //----------------------------------------------------------------------
    if ( verboseLevel > 0 )
    {
        std::cout << "Parabolic line search ..." << std::endl;
    }
    const Real sigma0 = 0.1;
    const Real sigma1 = 0.5;
    const Real alpha = 1.e-4;
    const Int maxIterations = 50;
    //----------------------------------------------------------------------
    static VectorType S_solCurrent = sol; // static to avoid too many alloc/de-alloc
    Int iterNonLinearLineSearch;
    Real lambdaCurrent, lambdaOld, normResTest, res2, resTest2, resTestOld2, c1, c2;
    //
    res2 = normRes * normRes;
    lambdaOld = lambda;
    lambdaCurrent = lambda;
    S_solCurrent = sol;
    sol += lambda * step;
    f.evalResidual ( residual, sol, iter );
    normResTest = residual.normInf();
    resTest2 = normResTest * normResTest;
    resTestOld2 = resTest2;
    iterNonLinearLineSearch = 0;
    while ( normResTest >= ( 1 - alpha * lambda ) * normRes )
    {
        iterNonLinearLineSearch++;
        // parabolic interpolation of lambda:
        c2 = lambdaOld * ( resTest2 - res2 ) - lambdaCurrent * ( resTestOld2 - res2 );
        if ( c2 >= 0 )
        {
            lambda = sigma1 * lambdaCurrent;
        }
        else
        {
            c1 = lambdaCurrent * lambdaCurrent * ( resTestOld2 - res2 )
                 - lambdaOld * lambdaOld * ( resTest2 - res2 );
            lambda = -c1 * .5 / c2;
            if ( lambda < sigma0 * lambdaCurrent )
            {
                lambda = sigma0 * lambdaCurrent;
            }
            if ( lambda > sigma1 * lambdaCurrent )
            {
                lambda = sigma1 * lambdaCurrent;
            }
        }
        if ( verboseLevel > 0 )
            std::cout << "--- line search " << iterNonLinearLineSearch << " : residual test = "
                      << normResTest << ", reduction = " << lambda << std::endl;
        // update solTest
        sol =  S_solCurrent;
        sol += lambda * step;
        lambdaOld = lambdaCurrent;
        lambdaCurrent = lambda;
        // eval norms
        f.evalResidual ( residual, sol, iter );
        normResTest = residual.normInf();
        resTestOld2 = resTest2;
        resTest2 = normResTest * normResTest;
        if ( iterNonLinearLineSearch > maxIterations )
        {
            if ( verboseLevel > 0 )
            {
                std::cout << "!!! Too many iterations in the line search algorithm" << std::endl;
            }
            return EXIT_FAILURE;
        }
    }
    normRes = normResTest;
    if ( verboseLevel > 0 )
    {
        std::cout << "Parabolic line search: final residual = " << normRes << std::endl;
    }

    return EXIT_SUCCESS;

}

//! Implementation of Line Search method with cubic interpolation
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

    @param f            Function
    @param residual     Residual
    @param sol          Solution
    @param step         Step to update the solution
    @param normRes      Norm of the residual
    @param lambda       Length of the Step
    @param slope        Slope value in linesearch algorithm
    @param iter         Iterations
    @param verboseLevel Option for detailed description
*/

template <class Fct, class VectorType>
Int NonLinearLineSearchCubic ( Fct& f, VectorType& residual, VectorType& sol, VectorType& step,
                               Real& normRes, Real& lambda, Real& slope, UInt iter, UInt const verboseLevel = 1 )
{

    //----------------------------------------------------------------------
    if ( verboseLevel > 0 )
    {
        std::cout << "Cubic line search ..." << std::endl;
    }

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
         normResTest, f0, ftest, c, c11, c12, c21, c22, a, b, disc, g1, g2, gprev = 0;
    //
    f0 = 0.5 * normRes * normRes;
    lambdaOld = lambda;
    S_solCurrent = sol;
    sol += lambda * step;
    f.evalResidual ( residual, sol, iter );
    normResTest = residual.normInf();
    ftest = 0.5 * normResTest * normResTest;
    iterLinesearch = 0;
    while ( ftest < f0 + m2 * slope * lambda // lambda is too small: extrapolation
            && iterLinesearch < maxIterations )
    {
        iterLinesearch++;
        lambda *= 2;
        sol =  S_solCurrent;
        sol += lambda * step;
        if ( verboseLevel > 0 )
        {
            std::cout << "--- line search (extrapolation, Goldstein rule)" << std::endl;
        }
        f.evalResidual ( residual, sol, iter );

        if ( verboseLevel > 0 )
            std::cout << "    line search iter : " << iterLinesearch << " residual test = "
                      << normResTest << ", lambda = " << lambda << std::endl;
        normResTest = residual.normInf();
        ftest = 0.5 * normResTest * normResTest;
    }
    if ( iterLinesearch == maxIterations )
    {
        if ( verboseLevel > 0 )
        {
            std::cout << "line search: too many extrapolations" << std::endl;
        }
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
            if ( ( std::fabs ( a ) > std::numeric_limits<Real>::min() ) &&
                    ( disc > std::numeric_limits<Real>::min() )
               )
            {
                lambdaTemporary = ( - b + std::sqrt ( disc ) ) / ( 3. * a );
            }
            else
            {
                lambdaTemporary = slope * lambda2 / ( 2. * g1 ) ;
            }
            if ( lambdaTemporary >= sigma1 * lambda )
            {
                lambdaTemporary = sigma1 * lambda;
            }
        }
        lambdaOld = lambda;
        gprev = ftest;
        if ( lambdaTemporary < sigma0 * lambda )
        {
            lambda *= sigma0;
        }
        else
        {
            lambda = lambdaTemporary;
        }
        //--
        sol =  S_solCurrent;
        sol += lambda * step;
        if ( verboseLevel > 0 )
        {
            std::cout << "--- line search (cubic interpolation, Armijo rule)" << std::endl;
        }
        f.evalResidual ( residual, sol, iter );
        normResTest = residual.normInf();
        if ( verboseLevel > 0 )
            std::cout << "    line search iter : " << iterLinesearch << " residual test = "
                      << normResTest << ", lambda = " << lambda << std::endl;
        ftest = 0.5 * normResTest * normResTest;
    }
    if ( iterLinesearch == maxIterations )
    {
        if ( verboseLevel > 0 )
        {
            std::cout << "line search: too many interpolations" << std::endl;
        }
        return EXIT_FAILURE;
    }
    normRes = normResTest;
    if ( verboseLevel > 0 )
    {
        std::cout << "Parabolic line search: final residual = " << normRes << std::endl;
    }
    return EXIT_SUCCESS;
}

} // Namespace LifeV

#endif /* _NonLinearLineSearch_H_ */
