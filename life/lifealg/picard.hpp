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
 *  @file
 *  @brief fixed point algorithms with relaxation
 *
 *  @date 00-00-0000
 *  @author Name Surname <name.surname@email.org>
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

namespace LifeV
{
/*

  Fixed point algorithms with relaxation

  input: f  : operator that must have a function
              eval(Real fx0,Real gx0,Real x0,int status)
       where fx0 = f(x0), gx0 is not used here, and status
       indicate whether it is the first iterate (status=1)
       or not (status=0).
  x0 : initial guess (modified)

  method:
       0: constant relaxation factor given by omega
              (if omega=1 : stantard fixed point without relaxation)
       1: a heuristic extension of Aitken formula (omega not used)
       2: Aitken formula due to Irons & Tuck (omega not used)
                 (in exact arithmetic 1 and 2 coincide)

   output: x1: the last iterate (i.e. |x1 - fx1| < stopping tolerance)
               (x0 is the iterate before the last)

   N.B.: the initial values of x1, fx0, fx1  are not used at
         the first iteration. x0 is the only initial guess
*/
template <typename Vector, typename Oper, typename Norm>
Int picard( Oper* f, Norm norm, Vector& fx1, Vector& fx0,
            Vector& gx1, Vector& gx0, Vector& x1, Vector& x0,
            Real abstol, Real reltol, Int& maxit, Int method, Real omega )
{
    Real xxnorm;
    Real mu = 0;
    Int iter = 1;
    Vector tmp = x0;
    f->eval( fx0, gx0, x0, 1 );
    x1 = fx0; // the first iteration is not relaxed
    Real normRes = norm( fx0 - x0 );
    Real stop_tol = abstol + reltol * normRes;
    std::string methodName;
    switch ( method )
    {
    case 0:
        methodName = "constant relaxation";
        break;
    case 1:
        methodName = "Aitken";
        break;
    case 2:
        methodName = "Irons & Tuck";
        break;
    default:
        std::cout << "Unknown Picard method\n";
        exit( 1 );
    }
    std::cout << "--------------------------=----------------------------------------"
              << std::endl;
    std::cout << "    Picard 1 : residual=" << normRes << ", stoping tolerance = "
              << stop_tol << std::endl;
    std::cout << "-------------------------------------------------------------------"
              << std::endl;
    if ( normRes <= stop_tol )
    {
        // the algorithm has converged in 1 iteration without relaxation
        std::cout << "--- Picard : convergence in 1 iteration\n\n";
        maxit = iter;
        return 0;
    }
    // One iteration was not enough : a second iteration is necessary
    iter++;
    f->eval( fx1, gx1, x1, 0 );
    normRes = norm( x1 - fx1 );
    /*
      if the above normRes is small enough, the algorithm
       has converged in 2 iterations without any relaxation. Otherwise,
       we enter in the loop with relaxation
    */
    while ( normRes > stop_tol && iter < maxit )
    {
        // computation of the relaxation parameter omega
        switch ( method )
        {
        case 1:
            tmp = x1 - fx1 - x0 + fx0;
            omega = dot( x1 - x0, tmp );
            xxnorm = dot( tmp, tmp );
            if ( xxnorm != 0 )
                omega = omega / xxnorm;
            else
                omega = 1.;
            break;
        case 2:
            tmp = x0 - fx0 - x1 + fx1;
            xxnorm = dot( tmp, tmp );
            if ( xxnorm != 0 )
                mu = mu + ( mu - 1 ) * dot( x1 - fx1 , tmp ) / xxnorm;
            omega = 1 - mu;
            break;
        }
        std::cout << "-----------------------------------------------------------------"
                  << std::endl;
        std::cout << "Picard " << iter << " ( " << methodName
                  << ") " << "omega=" << omega << ", residual=" << normRes << std::endl;
        std::cout << "-----------------------------------------------------------------"
                  << std::endl;
        iter++;
        // save the old iterate
        fx0 = fx1;
        x0 = x1;
        // compute the new iterate with relaxation
        x1 = ( 1 - omega ) * x1 + omega * fx1;
        // evaluate the function on the new iterate
        f->eval( fx1, gx1, x1, 0 );
        normRes = norm( fx1 - x1 );
    }
    if ( normRes > stop_tol )
    {
        std::cout << "!!! picard: convergence fails" << std::endl;
        maxit = iter;
        return 1;
    }
    std::cout << "--- picard: convergence in " << iter << " iterations\n\n";
    maxit = iter;
    return 0;
}
}

