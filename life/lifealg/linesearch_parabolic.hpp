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

*/

#ifndef _LINESEARCH_PARAB_H_
#define _LINESEARCH_PARAB_H_ 1

namespace LifeV
{
//! brent Implementation of Brent's method for root finding.
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

    @param f          Function
    @param residual   Residual
    @param sol        Solution
    @param step       Step to update the solution
    @param normRes    Norm of the residual
    @param lambda     Length of the Step
    @param iter       Iterations
    @param verbose    Option for detailed description
*/
template <class Fct, class VectorType>
Int LineSearchParabolic ( Fct& f, VectorType& residual, VectorType& sol, VectorType& step, Real& normRes,
                      Real& lambda, UInt iter, bool const verbose = true)
{

    //----------------------------------------------------------------------
    if (verbose)
        std::cout << "Parabolic line search ..." << std::endl;
    const Real sigma0 = 0.1;
    const Real sigma1 = 0.5;
    const Real alpha = 1.e-4;
    const Int maxIterations = 50;
    //----------------------------------------------------------------------
    static VectorType S_solCurrent = sol; // static to avoid too many alloc/de-alloc
    Int iterLinesearch;
    Real lambdaCurrent, lambdaOld, normResTest, res2, resTest2, resTestOld2, c1, c2;
    //
    res2 = normRes * normRes;
    lambdaOld = lambda;
    lambdaCurrent = lambda;
    S_solCurrent = sol;
    sol += lambda * step;
    f.evalResidual( residual, sol, iter );
    normResTest = residual.NormInf();
    resTest2 = normResTest * normResTest;
    resTestOld2 = resTest2;
    iterLinesearch = 0;
    while ( normResTest >= ( 1 - alpha * lambda ) * normRes )
    {
        iterLinesearch++;
        // parabolic interpolation of lambda:
        c2 = lambdaOld * ( resTest2 - res2 ) - lambdaCurrent * ( resTestOld2 - res2 );
        if ( c2 >= 0 )
            lambda = sigma1 * lambdaCurrent;
        else
        {
            c1 = lambdaCurrent * lambdaCurrent * ( resTestOld2 - res2 )
                 - lambdaOld * lambdaOld * ( resTest2 - res2 );
            lambda = -c1 * .5 / c2;
            if ( lambda < sigma0 * lambdaCurrent )
                lambda = sigma0 * lambdaCurrent;
            if ( lambda > sigma1 * lambdaCurrent )
                lambda = sigma1 * lambdaCurrent;
        }
        if (verbose)
            std::cout << "--- line search " << iterLinesearch << " : residual test = "
                      << normResTest << ", reduction = " << lambda << std::endl;
        // update solTest
        sol =  S_solCurrent;
        sol += lambda * step;
        lambdaOld = lambdaCurrent;
        lambdaCurrent = lambda;
        // eval norms
        f.evalResidual( residual, sol, iter );
        normResTest = residual.NormInf();
        resTestOld2 = resTest2;
        resTest2 = normResTest * normResTest;
        if ( iterLinesearch > maxIterations )
        {
            if (verbose)
                std::cout << "!!! Too many iterations in the line search algorithm" << std::endl;
            return EXIT_FAILURE;
        }
    }
    normRes = normResTest;
    if (verbose)
        std::cout << "Parabolic line search: final residual = " << normRes << std::endl;

    return EXIT_SUCCESS;

}



template <class Fct, class VectorType>
Int __attribute__ ((__deprecated__)) lineSearch_parab( Fct& f, VectorType& residual, VectorType& sol, VectorType& step, Real& normRes,
                      Real& lambda, UInt iter, bool const verbose = true)
{

    // you should replace any call to lineSearch_parab with a call to LineSearchParabolic
    return  LineSearchParabolic (  f, residual, sol, step,  normRes, lambda, iter, verbose );
}

} // Namespace LifeV

#endif /* _LINESEARCH_PARAB_H_ */
