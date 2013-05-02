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
    @brief Implementation of Brent's method for root finding.

    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @date

    @contributor Alessandro Melani <alessandro.melani@mail.polimi.it>
    @mantainer Alessandro Melani <alessandro.melani@mail.polimi.it>

    Brent's method  is a root-finding algorithm combining the bisection method, the secant method
    and inverse quadratic interpolation. It has the reliability of bisection but it can be as quick
    as some of the less reliable methods. The idea is to use the secant method or inverse quadratic
    interpolation if possible, because they converge faster, but to fall back to the more robust bisection
    method if necessary.

    @see  Chapter 4 of R.P. Brent, "Algorithms for Minimization without Derivatives", Prentice-Hall, Englewood Cliffs, NJ. 1973
 */
#ifndef _NonLinearBrent_H
#define _NonLinearBrent_H 1

#include <limits>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{
//! Implementation of Brent's method for root finding.
/*!
    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @date

    Brent's method  is a root-finding algorithm combining the bisection method, the secant method and inverse quadratic interpolation.
    It has the reliability of bisection but it can be as quick as some of the less reliable methods.
    The idea is to use the secant method or inverse quadratic interpolation if possible, because they converge faster, but to fall back
    to the more robust bisection method if necessary.

    See  Chapter 4 of R.P. Brent, "Algorithms for Minimization without Derivatives", Prentice-Hall, Englewood Cliffs, NJ. 1973

    @param f                     Function
    @param leftExtremeBase       Left extreme of the interval
    @param rightExtremeBase      Right extreme of the interval
    @param toll                  Tollerance
    @param maxIter               Maximum number of iterations
*/

template <class Function>
Real NonLinearBrent ( const Function& f, const Real& leftExtremeBase, const Real& rightExtremeBase, const Real& toll, const UInt& maxIter )
{

    // Trivial case
    if ( leftExtremeBase == rightExtremeBase )
    {
        return leftExtremeBase;
    }

    // Current left and right extreme of the interval
    Real leftExtreme ( leftExtremeBase ), rightExtreme ( rightExtremeBase );

    if ( leftExtreme > rightExtreme )
    {
        std::swap ( leftExtreme, rightExtreme );
    }

    // Current iteration
    UInt numIter ( static_cast<UInt> (0) );

    // Medium point of the current interval
    Real midpoint ( ( leftExtreme + rightExtreme ) / static_cast<Real> (2.) );

    // Gold
    Real gold ( static_cast<Real> ( (3. - std::sqrt (5.) ) / 2. ) );

    // Ausiliar variables
    Real p (0), q (0), r (0);
    Real x ( leftExtreme + gold * ( rightExtreme - leftExtreme ) );
    Real u (0), e (0), w (0), v (0), d (0);
    Real fx ( f (x) ), fv ( fx ), fw ( fx ), fu (0);

    // Relative tollerance
    Real tollRelative ( std::numeric_limits<Real>::epsilon() * std::fabs ( x ) + toll );


    while ( std::fabs ( x - midpoint) > ( static_cast<Real> (2.) * tollRelative - ( rightExtreme - leftExtreme ) / static_cast<Real> (2.) ) && numIter < maxIter )
    {

        // Clear some ausiliar variables
        p = q = r = static_cast<Real> (0);

        if ( std::fabs ( e ) > tollRelative )
        {
            r = ( x - w ) * ( fx - fv );
            q = ( x - v ) * ( fx - fw );
            p = ( x - v ) * q - ( x - w ) * r;
            q = static_cast<Real> (2.) * ( q - r );

            if ( q > static_cast<Real> (0.) )
            {
                p = - p;
            }
            else
            {
                q = - q;
            }

            r = e;
            e = d;
        }


        if ( std::fabs ( p ) < std::fabs ( q * r / static_cast<Real> (2.) ) && p > q * ( leftExtreme - x ) && p < q * ( rightExtreme - x ) )
        {
            d = p / q;
            u = x + d;

            if ( ( u - leftExtreme) < static_cast<Real> (2.) * tollRelative || ( rightExtreme - u ) < static_cast<Real> (2.) * tollRelative )
            {
                if ( x < midpoint )
                {
                    d = tollRelative;
                }
                else
                {
                    d = -tollRelative;
                }
            }
        }
        else
        {
            if ( x < midpoint )
            {
                e = rightExtreme - x;
            }
            else
            {
                e = leftExtreme - x;
            }

            d = gold * e;
        }

        if ( std::fabs ( d ) >= tollRelative )
        {
            u = x + d;
        }
        else
        {
            if ( d > static_cast<Real> (0.) )
            {
                u = x + tollRelative;
            }
            else
            {
                u = x - tollRelative;
            }
        }

        // Compute the value of f in the point u
        fu = f (u);

        if ( fu <= fx )
        {
            if ( u < x )
            {
                rightExtreme = x;
            }
            else
            {
                leftExtreme = x;
            }

            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if ( u < x )
            {
                leftExtreme = u;
            }
            else
            {
                rightExtreme = u;
            }

            if ( fu <= fw || x == w )
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else
            {
                if ( fu <= fv || v == x || v == w )
                {
                    v = u;
                    fv = fu;
                }
            }
        }

        // Compute the midpoint of the interval
        midpoint = ( leftExtreme + rightExtreme ) / static_cast<Real> (2.);

        // Compute the relative tollerance
        tollRelative = std::numeric_limits<Real>::epsilon() * std::fabs ( x ) + toll;

        // Increase the iterations number
        ++numIter;

    }
    std::ostringstream os;
    os << "Attention the brent scheme does not reach the convergence in "
       << numIter << ", with tollerance "
       << tollRelative << std::endl;

    // Check if the method reach the tollerance.
    ASSERT ( maxIter > numIter, os.str().c_str() );

    return x;

}

} // Namespace LifeV

#endif /* _NonLinearBrent_ */
