/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-27

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
   \file test_petsc.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-27
 */
#include <cstdlib>
#include <cassert>

#include <lifeconfig.h>

#include <vecUnknown.hpp>

#if defined(HAVE_BOOST_TEST)
// Boost.Test
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>

#include <debug.hpp>

#if defined(HAVE_PETSC_H)
#include <SolverPETSC.hpp>

using boost::unit_test_framework::test_suite;

void petsc_manager()
{
    using namespace LifeV;

    // petsc should be already initialized (singleton)
    BOOST_REQUIRE( PETSC::instance().isInitialized() );

    // should not do anything : already initialized
    PETSC::instance().initialize();
    BOOST_REQUIRE( PETSC::instance().isInitialized() );

    // finalize
    PETSC::instance().finalize();
    BOOST_REQUIRE( PETSC::instance().isInitialized() == false );

}

test_suite*
init_unit_test_suite( int argc, char** argv )
{
    test_suite* test= BOOST_TEST_SUITE( "PETSC Unit Test" );

    // this example will pass cause we know ahead of time number of expected failures
    test->add( BOOST_TEST_CASE( &petsc_manager ), 0 );

    return test;
}
#else
test_suite*
init_unit_test_suite( int argc, char** argv )
{
    test_suite* test= BOOST_TEST_SUITE( "PETSC Unit Test" );
    return test;
}
#endif /* HAVE_PETSC_H */
#else
int main()
{
    return EXIT_SUCCESS;
}
#endif
