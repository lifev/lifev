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
#include "assemb.hpp"

namespace LifeV
{
//
void compute_vec( Real constant, ElemVec& elvec, const CurrentFE& fe, int iblock )
{
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    for ( i = 0;i < fe.nbNode;i++ )
    {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
            s += fe.phi( i, ig ) * fe.weightDet( ig );
        vec( i ) += constant * s;
    }
}
}
