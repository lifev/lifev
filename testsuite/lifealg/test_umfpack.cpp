/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-28

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
   \file test_umfpack.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-28
 */
#include <stdio.h>

#include <vecUnknown.hpp>
#include <lifeconfig.h>
#include <debug.hpp>

#if defined(HAVE_BOOST_TEST)
// Boost.Test
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>

using boost::unit_test_framework::test_suite;

#if defined(HAVE_UMFPACK_H)

#include <SolverUMFPACK.hpp>


void test_umfpack()
{
    using namespace LifeV;

    int    n=5;

    /**
       2  3  0 0 0
       3  0  4 0 6
    A= 0 -1 -3 2 0 .
       0  0  1 0 0
       0  4  2 0 1
    */
    uint    Ap [ ] = {0, 2, 5, 9, 10, 12} ;
    uint   Ai [ ] = { 0, 1, 0,     2, 4, 1, 2, 3,       4, 2, 1, 4} ;
    double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;

    SolverUMFPACK __solver;
    //__solver.setMatrix( n, Ap, Ai, Ax );
    Vector x( n );
    Vector b( n );
    b[0] = 8.;
    b[1] = 45.;
    b[2] = -3.;
    b[3] = 3.;
    b[4] = 19. ;

    // solution should be x=[1 2 3 4 5]^T

    __solver.solve( x, b );

    Vector sol( n );
    sol[0]=1;
    sol[1]=2;
    sol[2]=3;
    sol[3]=4;
    sol[4]=5;
    Debug( 10000 ) << "norm_2( sol - x ) = " << norm_2( sol - x ) << "\n";
    BOOST_REQUIRE( norm_2( sol - x ) < 1e-10 );
}

test_suite*
init_unit_test_suite( int argc, char** argv )
{
    test_suite* test= BOOST_TEST_SUITE( "UMFPACK Unit Test" );

    // this example will pass cause we know ahead of time number of expected failures
    //test->add( BOOST_TEST_CASE( &test_umfpack ), 0 );

    return test;
}
#else
test_suite*
init_unit_test_suite( int argc, char** argv )
{
    test_suite* test= BOOST_TEST_SUITE( "UMFPACK Unit Test" );
    return test;
}
#endif  /* HAVE_UMFPACK_H */
#else
int main()
{
    return EXIT_SUCCESS;
}
#endif

