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
#ifndef _brentAlgorithm_H
#define _brentAlgorithm_H

#include <limits>
//#include <cmath>
//#include <iostream>


//#include <boost/function.hpp>
//#include <boost/bind.hpp>


// Policy for scalar functions
//typedef boost::function<Real ( const Real& )> fct_type;

namespace LifeV
{

template <class Function>
Real brent( const Function& f/* const fct_type& f*/, const Real& leftExtremeBase, const Real& rightExtremeBase, const Real& toll, const UInt& maxIter )
{
    // Current left and right extreme of the interval
    Real leftExtreme( leftExtremeBase ), rightExtreme( rightExtremeBase );

    if ( leftExtreme > rightExtreme )
    {
        std::swap( leftExtreme, rightExtreme );
    }

    // Current iteration
    UInt numIter( static_cast<UInt>(0) );

    // Medium point of the current interval
    Real midpoint( ( leftExtreme + rightExtreme ) / static_cast<Real>(2.) );

    // Gold (???)
    Real gold( static_cast<Real>( (3. - sqrt(5.)) / 2. ) );

    // Ausiliar variables
    Real p(0), q(0), r(0);
    Real x( leftExtreme + gold * ( rightExtreme - leftExtreme ) );
    Real u(0), e(0), w(0), v(0), d(0);
    Real fx( f(x) ), fv( fx ), fw( fx ), fu(0);

    // Relative tollerance
    Real tollRelative( std::numeric_limits<Real>::epsilon() * fabs( x ) + toll );


    while ( fabs( x - midpoint) > ( static_cast<Real>(2.) * tollRelative - ( rightExtreme - leftExtreme ) / static_cast<Real>(2.) ) && numIter < maxIter )
    {

        // Clear some ausiliar variables
        p = q = r = static_cast<Real>(0);

        if ( fabs( e ) > tollRelative )
        {
            r = ( x - w ) * ( fx - fv );
            q = ( x - v ) * ( fx - fw );
            p = ( x - v ) * q - ( x - w ) * r;
            q = static_cast<Real>(2.) * ( q - r );

            if ( q > static_cast<Real>(0.) )
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


        if ( fabs( p ) < fabs( q * r / static_cast<Real>(2.) ) && p > q * ( leftExtreme - x ) && p < q * ( rightExtreme - x ) )
        {
            d = p / q;
            u = x + d;

            if ( ( u - leftExtreme) < static_cast<Real>(2.) * tollRelative || ( rightExtreme - u ) < static_cast<Real>(2.) * tollRelative )
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

        if ( fabs( d ) >= tollRelative )
        {
            u = x + d;
        }
        else
        {
            if ( d > static_cast<Real>(0.) )
            {
                u = x + tollRelative;
            }
            else
            {
                u = x - tollRelative;
            }
        }

        // Compute the value of f in the point u
        fu = f(u);

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
        midpoint = ( leftExtreme + rightExtreme ) / static_cast<Real>(2.);

        // Compute the relative tollerance
        tollRelative = std::numeric_limits<Real>::epsilon() * fabs( x ) + toll;

        // Increase the iterations number
        ++numIter;
    }
    std::ostringstream os;
    os << "Attention the brent scheme does not reach the convergence in "
    << numIter
    << ", with tollerance "
    << tollRelative << std::endl;


    // Check if the method reach the tollerance.
    ASSERT( maxIter > numIter, os.str().c_str() );

    return x;

}

}
#endif // _brentAlgorithm_H
