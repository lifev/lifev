/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#include <algorithm>

#include "bdfNS.hpp"

namespace LifeV
{
BdfNS::BdfNS( const UInt n )
    :
    _bdf_u( n ),
    _bdf_p( std::max( UInt( 1 ), n - 1 ) )
{}

Bdf& BdfNS::bdf_u()
{
    return _bdf_u;
}

Bdf& BdfNS::bdf_p()
{
    return _bdf_p;
}
BdfNS::~BdfNS()
{}
}
