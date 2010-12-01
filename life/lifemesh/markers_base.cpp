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
#include <climits>
#include <life/lifemesh/markers_base.hpp>

namespace LifeV
{
//  ***********************************************************************************************************
//                                           IMPLEMENTATION
//  ***********************************************************************************************************
//MarkerTraits_Base
const MarkerTraits_Base::EntityFlag MarkerTraits_Base::NULLFLAG = LONG_MIN;

//MM: if you modity these changes here recheck function readNetgenMesh
//        becouse it uses this changes

MarkerTraits_Base::EntityFlag MarkerTraits_Base::strongerFlag( EntityFlag const & a, EntityFlag const & b )
{
    return a > b ? a : b ;
}

MarkerTraits_Base::EntityFlag MarkerTraits_Base::weakerFlag( EntityFlag const & a, EntityFlag const & b )
{
    if (a==NULLFLAG)return b;
    if (b==NULLFLAG)return a;
    return a < b ? a : b ;
}

bool MarkerTraits_Base::EqualFlags(const EntityFlag& a, const EntityFlag& b)
{return a==b;}

}
