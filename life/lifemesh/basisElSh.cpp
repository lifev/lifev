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
//! \file basisElSh.cc

#include <life/lifemesh/basisElSh.hpp>

namespace LifeV
{

UInt getReferenceDimension(const ReferenceShapes& shape)
{
    switch(shape)
    {
    case NONE:
    case POINT:
        return 0;
    case LINE:
        return 1;
    case TRIANGLE:
    case QUAD:
        return 2;
    case HEXA:
    case PRISM:
    case TETRA:
        return 3;
    default:
        ERROR_MSG(" Unknown shape " );
    }
    // Avoid warning
    return 0;
}

/*******************************************************************
          IMPLEMENTATION
*******************************************************************/

/*
          -- LinearTriangle
                3
               /
              /  \
             /    \
            /      \
           /        \
         1 ----------2

*/

id_type
LinearTriangle::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// id_type eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 2 * numEdges ] =
        {
            1, 2, 2, 3, 3, 1
        };
    ASSERT_BD( _point > 0 && _point < 3 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 2 * _localEdge + _point - 3 ];
}

/*
          -- QuadraticTriangle
                3
               /
              /  \
             5    \
            /      4
           /        \
         1 -----3----2
*/
id_type
QuadraticTriangle::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 3 * numEdges ] =
        {
            1, 2, 4, 2, 3, 5, 3, 1, 6
        };
    ASSERT_BD( _point > 0 && _point < 4 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 3 * _localEdge + _point - 4 ];
}

/*
          -- LinearQuad

        4-------3
         !     !
         !     !
        1-------2


*/
id_type
LinearQuad::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 2 * numEdges ] =
        {
            1, 2, 2, 3, 3, 4, 4, 1
        };
    ASSERT_BD( _point > 0 && _point < 3 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 2 * _localEdge + _point - 3 ];
}

/*
          -- QuadraticQuad

        4---7---3
        !       !
        8   9   6
        !       !
        1---5---2

*/
id_type
QuadraticQuad::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 3 * numEdges ] =
        {
            1, 2, 5, 2, 3, 6, 3, 4, 7, 4, 1, 8
        };
    ASSERT_BD( _point > 0 && _point < 4 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 3 * _localEdge + _point - 4 ];
}

/*
          -- LinearTetra

                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2


*/
id_type
LinearTetra::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 2 * numEdges ] =
        {
            1, 2, 2, 3, 3, 1, 1, 4, 2, 4, 3, 4
        };
    ASSERT_BD( _point > 0 && _point < 3 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 2 * _localEdge + _point - 3 ];
}

id_type
LinearTetra::fToP( id_type const _localFace, id_type const _point )   // indexing from 1!!!
// fToP(i,j) = localId of jth point on ith local face
{
    static const id_type _fToP[ 3 * numFaces ] =
        {
            1, 3, 2, 1, 2, 4, 2, 3, 4, 1, 4, 3
        }
        ; // AV - November 2000: fixed a little bug
    //  {1,3,2, 1,2,4, 2,3,4, 3,1,4};
    ASSERT_BD( _point > 0 && _point < 4 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;
    return _fToP[ 3 * _localFace + _point - 4 ];
}
pair<id_type, bool>
LinearTetra::fToE( id_type const _localFace, id_type const _edge )   // indexing from 1!!!
// fToE(i,j) = localId of jth edge on ith local face
{
    static const id_type _fToE[ 3 * numFaces ] =
        {
            3, 2, 1, 1, 5, 4, 2, 6, 5, 4, 6, 3
        };
    static const bool _orient[ 3 * numFaces ] =
        {
            0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1
        };

    ASSERT_BD( _edge > 0 && _edge < 4 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;

    return make_pair( _fToE[ 3 * _localFace + _edge - 4 ], _orient[ 3 * _localFace + _edge - 4 ] );
}

///////////////////////

/*
          -- LinearTetraBubble

                4
               / .
              /  \.3
             /  . \\
            / . .5 \\
           /.       \!
         1 ----------2


*/
id_type
LinearTetraBubble::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 2 * numEdges ] =
        {
            1, 2, 2, 3, 3, 1, 1, 4, 2, 4, 3, 4
        };
    ASSERT_BD( _point > 0 && _point < 3 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 2 * _localEdge + _point - 3 ];
}

id_type
LinearTetraBubble::fToP( id_type const _localFace, id_type const _point )   // indexing from 1!!!
// fToP(i,j) = localId of jth point on ith local face
{
    static const id_type _fToP[ 3 * numFaces ] =
        {
            1, 3, 2, 1, 2, 4, 2, 3, 4, 1, 4, 3
        }
        ; // AV - November 2000: fixed a little bug
    //  {1,3,2, 1,2,4, 2,3,4, 3,1,4};
    ASSERT_BD( _point > 0 && _point < 4 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;
    return _fToP[ 3 * _localFace + _point - 4 ];
}
pair<id_type, bool>
LinearTetraBubble::fToE( id_type const _localFace, id_type const _edge )   // indexing from 1!!!
// fToE(i,j) = localId of jth edge on ith local face
{
    static const id_type _fToE[ 3 * numFaces ] =
        {
            3, 2, 1, 1, 5, 4, 2, 6, 5, 4, 6, 3
        };
    static const bool _orient[ 3 * numFaces ] =
        {
            0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1
        };

    ASSERT_BD( _edge > 0 && _edge < 4 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;

    return make_pair( _fToE[ 3 * _localFace + _edge - 4 ], _orient[ 3 * _localFace + _edge - 4 ] );
}

///////////////////////

/*
          -- QuadraticTetra


                4
               / .10
              /  \.3
             8  . 9\
            / 7    \6
           /.       \!
         1 -----5----2
*/
id_type
QuadraticTetra::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 3 * numEdges ] =
        {
            1, 2, 5, 2, 3, 6, 3, 1, 7, 1, 4, 8, 2, 4, 9, 3, 4, 10
        };
    ASSERT_BD( _point > 0 && _point < 4 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 3 * _localEdge + _point - 4 ];
}

id_type
QuadraticTetra::fToP( id_type const _localFace, id_type const _point )   // indexing from 1!!!
// fToP(i,j) = localId of jth point on ith local face
{
    static const id_type _fToP[ 6 * numFaces ] =
        {
            1, 3, 2, 7, 6, 5,
            1, 2, 4, 5, 9, 8,
            2, 3, 4, 6, 10, 9,
            1, 4, 3, 8, 10, 7
        };
    ASSERT_BD( _point > 0 && _point < 7 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;
    return _fToP[ 6 * _localFace + _point - 7 ];
}

pair<id_type, bool>
QuadraticTetra::fToE( id_type const _localFace, id_type const _edge )   // indexing from 1!!!
{
    return LinearTetra::fToE( _localFace, _edge );
}

/*
          -- LinearHexa
*/
/*
                      8-------7
                     /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2
*/
id_type
LinearHexa::eToP( id_type const _localEdge, id_type const _point )   // indexing from 1!!!
// eToP(i,j) = localId of jth point on ith local edge
{
    static const id_type _eToP[ 2 * numEdges ] =
        {
            1, 2, 2, 3, 3, 4, 4, 1,
            1, 5, 2, 6, 3, 7, 4, 8,
            5, 6, 6, 7, 7, 8, 8, 5
        };
    ASSERT_BD( _point > 0 && _point < 3 ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 2 * _localEdge + _point - 3 ];
}

id_type
LinearHexa::fToP( id_type const _localFace, id_type const _point )   // indexing from 1!!!
// fToP(i,j) = localId of jth point on ith local face
{
    static const id_type _fToP[ 4 * numFaces ] =
        {
            1, 4, 3, 2,
            1, 5, 8, 4,
            1, 2, 6, 5,
            2, 3, 7, 6,
            3, 4, 8, 7,
            5, 6, 7, 8
        };
    ASSERT_BD( _point > 0 && _point < 5 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;
    return _fToP[ 4 * _localFace + _point - 5 ];
}

pair<id_type, bool>
LinearHexa::fToE( id_type const _localFace, id_type const _edge )   // indexing from 1!!!
// fToE(i,j) = localId of jth edge on ith local face
{
    static const id_type _fToE[ 4 * numFaces ] =
        {
            4, 3, 2, 1, 5, 12, 8, 4, 1, 6, 9, 5,
            2, 7, 10, 6, 3, 8, 11, 7, 9, 10, 11, 12
        };

    static const bool _orient[ 4 * numFaces ] =
        {
            0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0,
            1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1
        };


    ASSERT_BD( _edge > 0 && _edge < 5 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;

    return
        make_pair( _fToE[ 4 * _localFace + _edge - 5 ], _orient[ 4 * _localFace + _edge - 5 ] );
}

/*
          -- QuadraticHexa
*/

id_type
QuadraticHexa::eToP( id_type const _localEdge, id_type const _point )
{
    static const id_type _eToP[ 3 * numEdges ] =
        {
            1, 2, 9, 2, 3, 10, 3, 4, 11, 4, 1, 12,
            1, 5, 13, 2, 6, 14, 3, 7, 15, 4, 8, 16,
            5, 6, 17, 6, 7, 18, 7, 8, 19, 8, 5, 20
        };
    ASSERT_BD( ( _point > 0 && _point < 4 ) ) ;
    ASSERT_BD( _localEdge > 0 && _localEdge <= numEdges ) ;
    return _eToP[ 3 * _localEdge + _point - 4 ];
}

id_type
QuadraticHexa::fToP( id_type const _localFace, id_type const _point )
{
    static const id_type _fToP[ 9 * numFaces ] =
        {
            1, 4, 3, 2, 12, 11, 10, 9, 21,
            1, 5, 8, 4, 13, 20, 16, 12, 22,
            1, 2, 6, 5, 9, 14, 17, 13, 23,
            2, 3, 7, 6, 10, 15, 18, 14, 24,
            3, 4, 8, 7, 11, 16, 19, 15, 25,
            5, 6, 7, 8, 17, 18, 19, 20, 26
        };
    ASSERT_BD( _point > 0 && _point < 10 ) ;
    ASSERT_BD( _localFace > 0 && _localFace <= numFaces ) ;
    return _fToP[ 9 * _localFace + _point - 10 ];
}

pair<id_type, bool>
QuadraticHexa::fToE( id_type const _localFace, id_type const _edge )   // indexing from 1!!!
// fToE(i,j) = localId of jth edge on ith local face
{
    return LinearHexa::fToE( _localFace, _edge );
}
}
