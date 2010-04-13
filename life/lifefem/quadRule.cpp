/*-*- mode: c++ -*-
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
#include <life/lifefem/quadRule.hpp>

namespace LifeV
{
QuadRule::QuadRule( const QuadPoint* pt, int _id, std::string _name,
                    ReferenceShapes _shape, UInt _nbQuadPt, UInt _degOfExact ) :
        M_pt( pt ), M_shape( _shape ), M_id( _id ), M_name( _name ),
        M_nbQuadPt( _nbQuadPt ), M_degOfExact( _degOfExact )
{
    CONSTRUCTOR( "QuadRule" );
}

QuadRule::QuadRule( const QuadRule& qr ) :
        M_pt( qr.M_pt ), M_shape( qr.M_shape ), M_id( qr.M_id ), M_name( qr.M_name ),
        M_nbQuadPt( qr.M_nbQuadPt ), M_degOfExact( qr.M_degOfExact )
{
    CONSTRUCTOR( "QuadRule" );
}

QuadRule::~QuadRule()
{
    DESTRUCTOR( "QuadRule" );
}

std::ostream& operator << ( std::ostream& c, const QuadRule& qr )
{
    c << " name: " << qr.M_name << std::endl;
    c << " shape:" << ( int ) qr.M_shape << std::endl;
    c << " id: " << qr.M_id << std::endl;
    c << " nbQuadPt: " << qr.M_nbQuadPt << std::endl;
    c << " Points: \n";
    for ( UInt i (0);i < qr.M_nbQuadPt;++i )
        c << qr.M_pt[ i ] << std::endl;
    return c;
}


}
