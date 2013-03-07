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
 *  @brief File containing the Vector Container Test
 *
 *  @date 29-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <iomanip>
#include <string>


// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorContainer.hpp>

using namespace LifeV;

template < class VectorType >
UInt
TestFunction ( boost::shared_ptr<VectorType> A1, boost::shared_ptr<VectorType> B1,
               boost::shared_ptr<VectorType> A2, boost::shared_ptr<VectorType> B2,
               boost::shared_ptr<VectorType> A3, boost::shared_ptr<VectorType> B3,
               boost::shared_ptr<VectorType> A4, boost::shared_ptr<VectorType> B4 )
{
    //CONTAINER OF BASE VECTORS
    typedef VectorContainer< VectorType >                ContainerOfBaseVectors;
    typedef boost::shared_ptr<ContainerOfBaseVectors>    ContainerOfBaseVectors_ptr;

    ContainerOfBaseVectors_ptr V1, V2, V3, V4;
    Real scalar = 1.0;

    V1.reset ( new ContainerOfBaseVectors() );
    V1->push_back ( A1 );
    V1->push_back ( B1 );
    std::cout << "V1" << std::endl;
    V1->showMe();

    V2.reset ( new ContainerOfBaseVectors() );
    V2->push_back ( A2 );
    V2->push_back ( B2 );
    std::cout << "V2" << std::endl;
    V2->showMe();

    V3.reset ( new ContainerOfBaseVectors() );
    V3->push_back ( A3 );
    V3->push_back ( B3 );
    std::cout << "V3" << std::endl;
    V3->showMe();

    V4.reset ( new ContainerOfBaseVectors() );
    V4->push_back ( A4 );
    V4->push_back ( B4 );
    std::cout << "V4" << std::endl;
    V4->showMe();

    // Test operator= (initialize a vector as a copy of another vector)
    ContainerOfBaseVectors_ptr VV1, VV2, VV3;
    VV1.reset ( new ContainerOfBaseVectors() );
    *VV1 = *V1;
    std::cout << "VV1 = V1" << std::endl;
    VV1->showMe();

    // Test operator= (for initialize the vector with a scalar - cannot be done for an empty vector!)
    *VV1 = scalar;
    std::cout << "VV1 = 1.0" << std::endl;
    VV1->showMe();

    // Test operator+=
    *VV1 += *V1;
    std::cout << "VV1 += V1" << std::endl;
    VV1->showMe();

    // Test operator-
    *VV1 = *VV1 - *V2;
    std::cout << "VV1 = VV1 - V2" << std::endl;
    VV1->showMe();

    // Test operator*= (multiplication: scalar * vector)
    scalar = 1.23456789;
    *VV1 = *V1;
    *VV1 = scalar * *V1;
    std::cout << "VV1 = 1.23456789 * V1" << std::endl;
    VV1->showMe();

    // Test operator* (scalar product multiplication)
    scalar = VV1->dot ( *V1 );
    std::cout << "scalarProduct = VV1.Dot(V1) = " << scalar << std::endl << std::endl;

    // Concatenate two vector of vectors
    VV1->push_back ( *V1 );
    VV1->push_back ( *V2 );
    std::cout << "VV1->push_back( V1 ); VV1->push_back( V2 )" << std::endl;
    VV1->showMe();

    VV2.reset ( new ContainerOfBaseVectors() );
    VV2->push_back ( *V3 );
    VV2->push_back ( *V4 );
    VV2->push_back ( *V3 );
    std::cout << "VV2->push_back( V3 ); VV2->push_back( V4 ); VV2->push_back( V3 )" << std::endl;
    VV2->showMe();

    // Element by element Multiplication
    std::cout << "VV2 *= VV1" << std::endl;
    *VV2 *= *VV1;
    VV2->showMe();
    std::cout << "VV1 =" << std::endl;
    VV1->showMe();

    // Element by element Division
    std::cout << "VV2 /= VV1" << std::endl;
    *VV2 /= *VV1;
    VV2->showMe();
    std::cout << "VV1 =" << std::endl;
    VV1->showMe();

    // Element by element Multiplication (NaN is interpreted as a zero)
    std::cout << "VV2 && VV1" << std::endl;
    (*VV2 && *VV1).showMe();

    // Element by element Comparison ==
    std::cout << "VV1 == 0" << std::endl;
    VV3.reset ( new ContainerOfBaseVectors() );
    *VV3 = ( *VV1 == 0.0 );
    VV3->showMe();

    // Element by element Comparison >=
    std::cout << "VV1 >= 2.2" << std::endl;
    *VV3 = ( *VV1 >= 2.2 );
    VV3->showMe();

    // Abs
    std::cout << "VV2 *= -1" << std::endl;
    scalar = -1.0;
    *VV2 *= scalar;
    VV2->showMe();

    std::cout << "VV2->abs( *VV3 )" << std::endl;
    VV2->abs ( *VV3 );
    std::cout << "VV2:" << std::endl;
    VV2->showMe();
    std::cout << "VV3:" << std::endl;
    VV3->showMe();







    // Replace a vector
    UInt pos = 2;
    (*VV1) ( pos ) = B1;
    std::cout << "VV1( 2 ) = B1" << std::endl;
    VV1->showMe();

    // Compute the weight norm2 of the vector
    VV1->weightNorm2();
    std::cout << "WeightNorm2(VV1) = " << VV1->weightNorm2() << std::endl << std::endl;

    // operator!
    std::cout << "!VV1 " << std::endl;
    (! (*VV1) ).showMe();

    std::cout << "!VV3 " << std::endl;
    (! (*VV3) ).showMe();

    return EXIT_SUCCESS;
}
