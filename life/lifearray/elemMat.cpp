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
#include <life/lifearray/elemMat.hpp>

namespace LifeV
{
ElemMat::~ElemMat()
{
    DESTRUCTOR( "ElemMat" );
}

ElemMat::ElemMat( UInt nNode1, UInt nbr1, UInt nbc1 ) :
        _mat( nNode1*nbr1, nNode1*nbc1 )
{
    CONSTRUCTOR( "ElemMat 1" );
    //
    _nBlockRow = nbr1;
    _nBlockCol = nbc1;
    //
    _nRow.resize( _nBlockRow );
    _firstRow.resize( _nBlockRow );
    _nCol.resize( _nBlockCol );
    _firstCol.resize( _nBlockCol );
    //
    UInt first = 0;
    for ( UInt n = 0;n < nbr1;n++ )
    {
        _nRow[ n ] = nNode1;
        _firstRow[ n ] = first;
        first += nNode1;
    }
    //
    first = 0;
    for ( UInt n = 0;n < nbc1;n++ )
    {
        _nCol[ n ] = nNode1;
        _firstCol[ n ] = first;
        first += nNode1;
    }
}

ElemMat::ElemMat( UInt nNode1, UInt nbr1, UInt nbc1,
                  UInt nNode2, UInt nbr2, UInt nbc2 ) :
        _mat( nNode1*nbr1 + nNode2*nbr2, nNode1*nbc1 + nNode2*nbc2 )
{
    CONSTRUCTOR( "ElemMat 2" );
    //
    _nBlockRow = nbr1 + nbr2;
    _nBlockCol = nbc1 + nbc2;
    //
    _nRow.resize( _nBlockRow );
    _firstRow.resize( _nBlockRow );
    _nCol.resize( _nBlockCol );
    _firstCol.resize( _nBlockCol );
    //
    UInt first = 0, n;
    for ( n = 0;n < nbr1;n++ )
    {
        _nRow[ n ] = nNode1;
        _firstRow[ n ] = first;
        first += nNode1;
    }
    for ( n = nbr1;n < nbr1 + nbr2;n++ )
    {
        _nRow[ n ] = nNode2;
        _firstRow[ n ] = first;
        first += nNode2;
    }
    //
    first = 0;
    for ( n = 0;n < nbc1;n++ )
    {
        _nCol[ n ] = nNode1;
        _firstCol[ n ] = first;
        first += nNode1;
    }
    for ( n = nbc1;n < nbc1 + nbc2;n++ )
    {
        _nCol[ n ] = nNode2;
        _firstCol[ n ] = first;
        first += nNode2;
    }
}


ElemMat::ElemMat( UInt nNode1, UInt nbr1, UInt nbc1,
                  UInt nNode2, UInt nbr2, UInt nbc2,
                  UInt nNode3, UInt nbr3, UInt nbc3 ) :
        _mat( nNode1*nbr1 + nNode2*nbr2 + nNode3*nbr3,
              nNode1*nbc1 + nNode2*nbc2 + nNode3*nbc3 )
{
    CONSTRUCTOR( "ElemMat 3" );
    //
    _nBlockRow = nbr1 + nbr2 + nbr3;
    _nBlockCol = nbc1 + nbc2 + nbc3;
    //
    _nRow.resize( _nBlockRow );
    _firstRow.resize( _nBlockRow );
    _nCol.resize( _nBlockCol );
    _firstCol.resize( _nBlockCol );
    //
    UInt first = 0, n;
    for ( n = 0;n < nbr1;n++ )
    {
        _nRow[ n ] = nNode1;
        _firstRow[ n ] = first;
        first += nNode1;
    }
    for ( n = nbr1;n < nbr1 + nbr2;n++ )
    {
        _nRow[ n ] = nNode2;
        _firstRow[ n ] = first;
        first += nNode2;
    }
    for ( n = nbr1 + nbr2;n < nbr1 + nbr2 + nbr3;n++ )
    {
        _nRow[ n ] = nNode3;
        _firstRow[ n ] = first;
        first += nNode3;
    }
    //
    first = 0;
    for ( n = 0;n < nbc1;n++ )
    {
        _nCol[ n ] = nNode1;
        _firstCol[ n ] = first;
        first += nNode1;
    }
    for ( n = nbc1;n < nbc1 + nbc2;n++ )
    {
        _nCol[ n ] = nNode2;
        _firstCol[ n ] = first;
        first += nNode2;
    }
    for ( n = nbc1 + nbc2;n < nbc1 + nbc2 + nbc3;n++ )
    {
        _nCol[ n ] = nNode3;
        _firstCol[ n ] = first;
        first += nNode3;
    }
}

void ElemMat::showMe( std::ostream& c )
{
    UInt i, j;
    for ( i = 0;i < _nBlockRow;i++ )
    {
        for ( j = 0;j < _nBlockCol;j++ )
        {
            c << "Block (" << i << "," << j << "), ";
            c << block( i, j ) << std::endl;
        }
    }
}
}
