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
    @brief Contains the basic element shapes, to be used by Geometric
           and Finite Elements

    @author Luca Formaggia
    @contributor Zhen Wang <zwang26@emory.edu>
    @contributor Tiziano Passerini <tiziano@mathcs.emory.edu>
*/

#include <life/lifemesh/basisElSh.hpp>

namespace LifeV
{
UInt __attribute__ ((__deprecated__)) getReferenceDimension(const ReferenceShapes& shape)
{
	// You should substitute any call to getReferenceDimension with a call to getReferenceShapeDimension
	return getReferenceShapeDimension(shape);
}

UInt getReferenceShapeDimension(const ReferenceShapes& shape)
{
    switch (shape)
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
ID
LinearTriangle::eToP( ID const iEdge, ID const jPoint )
// ID eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 2 * numEdges ] =
    {
        1, 2, 2, 3, 3, 1
    };
    ASSERT_BD( jPoint > 0 && jPoint < 3 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 2 * iEdge + jPoint - 3 ];
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
ID
QuadraticTriangle::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 3 * numEdges ] =
    {
        1, 2, 4, 2, 3, 5, 3, 1, 6
    };
    ASSERT_BD( jPoint > 0 && jPoint < 4 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 3 * iEdge + jPoint - 4 ];
}


/*
          -- LinearQuad

        4-------3
         !     !
         !     !
        1-------2


*/
ID
LinearQuad::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 2 * numEdges ] =
    {
        1, 2, 2, 3, 3, 4, 4, 1
    };
    ASSERT_BD( jPoint > 0 && jPoint < 3 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 2 * iEdge + jPoint - 3 ];
}


/*
          -- QuadraticQuad

        4---7---3
        !       !
        8   9   6
        !       !
        1---5---2

*/
ID
QuadraticQuad::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 3 * numEdges ] =
    {
        1, 2, 5, 2, 3, 6, 3, 4, 7, 4, 1, 8
    };
    ASSERT_BD( jPoint > 0 && jPoint < 4 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 3 * iEdge + jPoint - 4 ];
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
ID
LinearTetra::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 2 * numEdges ] =
    {
        1, 2, 2, 3, 3, 1, 1, 4, 2, 4, 3, 4
    };
    ASSERT_BD( jPoint > 0 && jPoint < 3 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 2 * iEdge + jPoint - 3 ];
}

ID
LinearTetra::fToP( ID const iFace, ID const jPoint )
// fToP(i,j) = localId of jth point on ith local face
{
    static const ID _fToP[ 3 * numFaces ] =
    {
        1, 3, 2, 1, 2, 4, 2, 3, 4, 1, 4, 3
    }
    ; // AV - November 2000: fixed a little bug
    //  {1,3,2, 1,2,4, 2,3,4, 3,1,4};
    ASSERT_BD( jPoint > 0 && jPoint < 4 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;
    return _fToP[ 3 * iFace + jPoint - 4 ];
}

std::pair<ID, bool>
LinearTetra::fToE( ID const iFace, ID const jEdge )
// fToE(i,j) = localId of jth edge on ith local face
{
    static const ID _fToE[ 3 * numFaces ] =
    {
        3, 2, 1, 1, 5, 4, 2, 6, 5, 4, 6, 3
    };
    static const bool _orient[ 3 * numFaces ] =
    {
        0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1
    };

    ASSERT_BD( jEdge > 0 && jEdge < 4 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;

    return std::make_pair( _fToE[ 3 * iFace + jEdge - 4 ], _orient[ 3 * iFace + jEdge - 4 ] );
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
ID
LinearTetraBubble::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 2 * numEdges ] =
    {
        1, 2, 2, 3, 3, 1, 1, 4, 2, 4, 3, 4
    };
    ASSERT_BD( jPoint > 0 && jPoint < 3 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 2 * iEdge + jPoint - 3 ];
}

ID
LinearTetraBubble::fToP( ID const iFace, ID const jPoint )
// fToP(i,j) = localId of jth point on ith local face
{
    static const ID _fToP[ 3 * numFaces ] =
    {
        1, 3, 2, 1, 2, 4, 2, 3, 4, 1, 4, 3
    }
    ; // AV - November 2000: fixed a little bug
    //  {1,3,2, 1,2,4, 2,3,4, 3,1,4};
    ASSERT_BD( jPoint > 0 && jPoint < 4 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;
    return _fToP[ 3 * iFace + jPoint - 4 ];
}

std::pair<ID, bool>
LinearTetraBubble::fToE( ID const iFace, ID const jEdge )
// fToE(i,j) = localId of jth edge on ith local face
{
    static const ID _fToE[ 3 * numFaces ] =
    {
        3, 2, 1, 1, 5, 4, 2, 6, 5, 4, 6, 3
    };
    static const bool _orient[ 3 * numFaces ] =
    {
        0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1
    };

    ASSERT_BD( jEdge > 0 && jEdge < 4 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;

    return std::make_pair( _fToE[ 3 * iFace + jEdge - 4 ], _orient[ 3 * iFace + jEdge - 4 ] );
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
ID
QuadraticTetra::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 3 * numEdges ] =
    {
        1, 2, 5, 2, 3, 6, 3, 1, 7, 1, 4, 8, 2, 4, 9, 3, 4, 10
    };
    ASSERT_BD( jPoint > 0 && jPoint < 4 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 3 * iEdge + jPoint - 4 ];
}

ID
QuadraticTetra::fToP( ID const iFace, ID const jPoint )
// fToP(i,j) = localId of jth point on ith local face
{
    static const ID _fToP[ 6 * numFaces ] =
    {
        1, 3, 2, 7, 6, 5,
        1, 2, 4, 5, 9, 8,
        2, 3, 4, 6, 10, 9,
        1, 4, 3, 8, 10, 7
    };
    ASSERT_BD( jPoint > 0 && jPoint < 7 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;
    return _fToP[ 6 * iFace + jPoint - 7 ];
}

std::pair<ID, bool>
QuadraticTetra::fToE( ID const iFace, ID const jEdge )
// fToE(i,j) = localId of jth edge on ith local face
{
    return LinearTetra::fToE( iFace, jEdge );
}


/*
          -- LinearHexa


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
ID
LinearHexa::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 2 * numEdges ] =
    {
        1, 2, 2, 3, 3, 4, 4, 1,
        1, 5, 2, 6, 3, 7, 4, 8,
        5, 6, 6, 7, 7, 8, 8, 5
    };
    ASSERT_BD( jPoint > 0 && jPoint < 3 ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 2 * iEdge + jPoint - 3 ];
}

ID
LinearHexa::fToP( ID const iFace, ID const jPoint )
// fToP(i,j) = localId of jth point on ith local face
{
    static const ID _fToP[ 4 * numFaces ] =
    {
        1, 4, 3, 2,
        1, 5, 8, 4,
        1, 2, 6, 5,
        2, 3, 7, 6,
        3, 4, 8, 7,
        5, 6, 7, 8
    };
    ASSERT_BD( jPoint > 0 && jPoint < 5 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;
    return _fToP[ 4 * iFace + jPoint - 5 ];
}

std::pair<ID, bool>
LinearHexa::fToE( ID const iFace, ID const jEdge )
// fToE(i,j) = localId of jth edge on ith local face
{
    static const ID _fToE[ 4 * numFaces ] =
    {
        4, 3, 2, 1, 5, 12, 8, 4, 1, 6, 9, 5,
        2, 7, 10, 6, 3, 8, 11, 7, 9, 10, 11, 12
    };

    static const bool _orient[ 4 * numFaces ] =
    {
        0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0,
        1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1
    };


    ASSERT_BD( jEdge > 0 && jEdge < 5 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;

    return
        std::make_pair( _fToE[ 4 * iFace + jEdge - 5 ], _orient[ 4 * iFace + jEdge - 5 ] );
}


/*
          -- QuadraticHexa
*/
ID
QuadraticHexa::eToP( ID const iEdge, ID const jPoint )
// eToP(i,j) = localId of jth point on ith local edge
{
    static const ID _eToP[ 3 * numEdges ] =
    {
        1, 2, 9, 2, 3, 10, 3, 4, 11, 4, 1, 12,
        1, 5, 13, 2, 6, 14, 3, 7, 15, 4, 8, 16,
        5, 6, 17, 6, 7, 18, 7, 8, 19, 8, 5, 20
    };
    ASSERT_BD( ( jPoint > 0 && jPoint < 4 ) ) ;
    ASSERT_BD( iEdge > 0 && iEdge <= numEdges ) ;
    return _eToP[ 3 * iEdge + jPoint - 4 ];
}

ID
QuadraticHexa::fToP( ID const iFace, ID const jPoint )
// fToP(i,j) = localId of jth point on ith local face
{
    static const ID _fToP[ 9 * numFaces ] =
    {
        1, 4, 3, 2, 12, 11, 10, 9, 21,
        1, 5, 8, 4, 13, 20, 16, 12, 22,
        1, 2, 6, 5, 9, 14, 17, 13, 23,
        2, 3, 7, 6, 10, 15, 18, 14, 24,
        3, 4, 8, 7, 11, 16, 19, 15, 25,
        5, 6, 7, 8, 17, 18, 19, 20, 26
    };
    ASSERT_BD( jPoint > 0 && jPoint < 10 ) ;
    ASSERT_BD( iFace > 0 && iFace <= numFaces ) ;
    return _fToP[ 9 * iFace + jPoint - 10 ];
}

std::pair<ID, bool>
QuadraticHexa::fToE( ID const iFace, ID const jEdge )
// fToE(i,j) = localId of jth edge on ith local face
{
    return LinearHexa::fToE( iFace, jEdge );
}
}
