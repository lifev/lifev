/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-11-04

  Copyright (C) 2004 EPFL

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
/**
   \file test_boost.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-04
 */
#include <cmath>
#include <iostream>
#include <fstream>


#include <lifeconfig.h>
#include <lifeV.hpp>

#include <boost/function.hpp>
#include <boost/timer.hpp>

double f1( double, double, double, double, LifeV::ID const& )
{
    return 1;
};

double f2( double t, double x , double y, double z, LifeV::ID const& i)
{
    switch ( i )
    {
        case 1:
            return x*x*cos( t );
            break;
        case 2:
            return y*y*cos( t );
            break;
        case 3:
            return z*z*cos( t );
            break;
    }
};
void
test_function()
{
    const ulong TGV = 10000001;
    std::ofstream __out( "bench.txt" );
    __out.precision( 8 );
    for ( ulong N = 10;N < TGV;N*=10 )
    {
        __out << N << " ";
        typedef boost::function<double ( double, double, double, double, LifeV::ID const& )> f_type;
        double t, x, y, z;
        LifeV::ID id = 1;

        std::cout << "testing dummy function with " << N << " calls\n";
        boost::timer __timer;
        for ( ulong __i = 0;__i < N;++__i )
        {
            double a = f1(t, x, y, z, id);
        }
        std::cout << "Elapsed time for pure function pointer: " << __timer.elapsed() << "\n";
        __out << __timer.elapsed() << " ";

        f_type myfunctor( f1 );

        __timer.restart();
        for ( ulong __i = 0;__i < N;++__i )
        {
            double a = myfunctor(t, x, y, z, id);
        }
        std::cout << "Elapsed time for boost::function: " << __timer.elapsed() << "\n";
        __out << __timer.elapsed() << " ";

        std::cout << "testing not so dummy function with " << N << " calls\n";
        __timer.restart();

        for ( ulong __i = 0;__i < N;++__i )
        {
            id = __i%3+1;
            double a = f2(t, x, y, z, id);
        }
        std::cout << "Elapsed time for pure function pointer: " << __timer.elapsed() << "\n";
        __out << __timer.elapsed() << " ";

        f_type myfunctor2( f2 );
        __timer.restart();

        for ( ulong __i = 0;__i < N;++__i )
        {
            id = __i%3+1;
            double a = myfunctor2(t, x, y, z, id);
        }
        std::cout << "Elapsed time for boost::function: " << __timer.elapsed() << "\n";
        __out << __timer.elapsed() << "\n";
    }

}

int main()
{
    test_function();
}
