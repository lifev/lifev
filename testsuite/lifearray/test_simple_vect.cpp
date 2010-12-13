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
    @brief Test SimpleVect class

    @author
    @contributor
    @maintainer

    @date 00-00-0000

	Test if the template class SimpleVect compiles and works correctly.
 */

#include<iostream>

#include <life/lifearray/SimpleVect.hpp>

int
main()
{
    using namespace LifeV;

    SimpleVect<int> a;
    SimpleVect<float> b(10);

    for (SimpleVect<float>::iterator p=b.begin(); p!= b.end(); ++p)
    {
        *p=10.0;
    }
    std::cout << b(4) << "\n";
}


