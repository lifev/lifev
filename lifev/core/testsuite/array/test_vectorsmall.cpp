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

/* ========================================================

Simple VectorSmall class test

*/


/**
   @file test_vectorsmall.cpp
   @author A. Cervone <ant.cervone@gmail.com>
   @date 2011-06-15
*/


// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorSmall.hpp>

using namespace LifeV;

// ===================================================
//! Main
// ===================================================
int main()
{
    // test for dim = 3
    Vector3D v1 ( 1., 1., 2. ), v2 ( 0., 1., 0. ), v3;

    std::cout << v1              << std::endl << std::endl;
    std::cout << v2              << std::endl << std::endl;
    std::cout << v3              << std::endl << std::endl;
    std::cout << v1[ 0 ]         << std::endl << std::endl;
    v1 [ 0 ] = 0.;
    std::cout << v1[ 0 ]         << std::endl << std::endl;
    std::cout << v1 ( 0 )         << std::endl << std::endl;
    v1 ( 0 ) = 1.;
    std::cout << v1 ( 0 )         << std::endl << std::endl;
    std::cout << v1 + v2         << std::endl << std::endl;
    std::cout << v1 - v2         << std::endl << std::endl;
    std::cout << 2. * v1         << std::endl << std::endl;
    std::cout << v1 / 2.         << std::endl << std::endl;
    std::cout << v1.dot ( v2 )    << std::endl << std::endl;
    std::cout << v1.cross ( v2 )  << std::endl << std::endl;
    std::cout << v1.normalized() << std::endl << std::endl;
    v1.normalize();
    std::cout << v1              << std::endl << std::endl;

    std::vector<Real> v4 ( 3, 1. );
    std::cout << castToVector3D ( v4 ) << std::endl << std::endl;
    KN<Real> v5 ( 3, 2. );
    std::cout << castToVector3D ( v5 ) << std::endl << std::endl;


    // test for dim = 5
    VectorSmall<5> a, b;
    for ( UInt i = 0; i < 5; i++ )
    {
        a[ i ] = i;
        b ( i ) = 4 - i;
    }

    std::cout << a              << std::endl << std::endl;
    std::cout << b              << std::endl << std::endl;
    std::cout << a + b          << std::endl << std::endl;
    std::cout << a - b          << std::endl << std::endl;
    std::cout << 0.5 * a        << std::endl << std::endl;
    std::cout << b / 2.         << std::endl << std::endl;
    std::cout << a.dot (b)       << std::endl << std::endl;
    std::cout << a.normalized() << std::endl << std::endl;
    b.normalize();
    std::cout << b              << std::endl << std::endl;

    std::vector<Real> c ( 5, 1. );
    std::cout << castToVectorSmall<5> ( c ) << std::endl << std::endl;
    KN<Real> d ( 5, 2. );
    std::cout << castToVectorSmall<5> ( d ) << std::endl << std::endl;

    VectorSmall<10> v6 = VectorSmall<10>::Constant ( 3. );
    std::cout << v6 << std::endl << std::endl;
    v6 = VectorSmall<10>::Zero();
    std::cout << v6 << std::endl << std::endl;

    return 0;
}
