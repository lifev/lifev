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

#include <lifev/core/mesh/ElementShapes.hpp>

namespace LifeV
{

UInt shapeDimension (const ReferenceShapes& shape)
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
            ERROR_MSG (" Unknown shape " );
    }
    // Avoid warning
    return 0;
}



/*******************************************************************
          IMPLEMENTATION
*******************************************************************/

/*
         Tetra
 */

std::pair<ID, bool>
Tetra::faceToEdge ( ID const& iFace, ID const& jEdge )
// faceToEdge(i,j) = localId of jth edge on ith local face
{
    static const ID _faceToEdge[ 3 * S_numFaces ] =
    {
        2, 1, 0,
        0, 4, 3,
        1, 5, 4,
        3, 5, 2
    };

    static const bool _orient[ 3 * S_numFaces ] =
    {
        0, 0, 0,
        1, 1, 0,
        1, 1, 0,
        1, 0, 1
    };


    ASSERT_BD ( jEdge < 3 ) ;
    ASSERT_BD ( iFace < S_numFaces ) ;

    return std::make_pair ( _faceToEdge[ 3 * iFace + jEdge ], _orient[ 3 * iFace + jEdge ] );
}

/*
         Hexa
 */

std::pair<ID, bool>
Hexa::faceToEdge ( ID const& iFace, ID const& jEdge )
// faceToEdge(i,j) = localId of jth edge on ith local face
{
    static const ID _faceToEdge[ 4 * S_numFaces ] =
    {
        3, 2, 1, 0,
        4, 11, 7, 3,
        0, 5, 8, 4,
        1, 6, 9, 5,
        2, 7, 10, 6,
        8, 9, 10, 11
    };

    static const bool _orient[ 4 * S_numFaces ] =
    {
        0, 0, 0, 0,
        1, 0, 0, 1,
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 1, 1
    };

    ASSERT_BD ( jEdge < 4 ) ;
    ASSERT_BD ( iFace < S_numFaces ) ;

    return
        std::make_pair ( _faceToEdge[ 4 * iFace + jEdge ] , _orient[ 4 * iFace + jEdge ] );
}



/*
          -- LinearTriangle
                2
               /
              /  \
             /    \
            /      \
           /        \
         0 ----------1

*/
ID
LinearTriangle::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// ID edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 2 * S_numEdges ] =
    {
        1, 2,
        2, 3,
        3, 1
    };
    ASSERT_BD ( jPoint < 2 ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 2 * iEdge + jPoint ] - 1;
}


/*
          -- QuadraticTriangle
                2
               /
              /  \
             5    \
            /      4
           /        \
         0 -----3----1
*/
ID
QuadraticTriangle::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 3 * S_numEdges ] =
    {
        0, 1, 3,
        1, 2, 4,
        2, 0, 5
    };
    ASSERT_BD ( jPoint < 3 ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 3 * iEdge + jPoint ];
}


/*
          -- LinearQuad

        3-------2
         !     !
         !     !
        0-------1


*/
ID
LinearQuad::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 2 * S_numEdges ] =
    {
        0, 1,
        1, 2,
        2, 3,
        3, 0
    };
    ASSERT_BD ( jPoint < 2 ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 2 * iEdge + jPoint];
}

/*
          -- QuadraticQuad

        3---6---2
        !       !
        7   8   5
        !       !
        0---4---1

*/
ID
QuadraticQuad::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 3 * S_numEdges ] =
    {
        0, 1, 4,
        1, 2, 5,
        2, 3, 6,
        3, 0, 7
    };
    ASSERT_BD ( jPoint < 3 ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 3 * iEdge + jPoint ];
}


/*
          -- LinearTetra

                3
               / .
              /  \.2
             /  . \\
            / .    \\
           /.       \!
         0 ----------1


*/
ID
LinearTetra::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 2 * S_numEdges ] =
    {
        0, 1,
        1, 2,
        2, 0,
        0, 3,
        1, 3,
        2, 3
    };
    ASSERT_BD ( jPoint < 2 ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 2 * iEdge + jPoint ];
}

ID
LinearTetra::faceToPoint ( ID const& iFace, ID const& jPoint )
// faceToPoint(i,j) = localId of jth point on ith local face
{
    static const ID _faceToPoint[ 3 * S_numFaces ] =
    {
        0, 2, 1,
        0, 1, 3,
        1, 2, 3,
        0, 3, 2
    };
    ASSERT_BD ( jPoint < 3 ) ;
    ASSERT_BD ( iFace < S_numFaces ) ;
    return _faceToPoint[ 3 * iFace + jPoint ];
}



///////////////////////

/*
          -- LinearTetraBubble

                3
               / .
              /  \.2
             /  . \\
            / . .4 \\
           /.       \!
         0 ----------1


*/
ID
LinearTetraBubble::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    return LinearTetra::faceToPoint (iEdge, jPoint);
}


///////////////////////

/*
          -- QuadraticTetra


                3
               / .9
              /  \.2
             7  . 8\
            / 6    \5
           /.       \!
         0 -----4----1
*/
ID
QuadraticTetra::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 3 * S_numEdges ] =
    {
        0, 1, 4,
        1, 2, 5,
        2, 0, 6,
        0, 3, 7,
        1, 3, 8,
        2, 3, 9
    };
    ASSERT_BD ( jPoint < 3 ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 3 * iEdge + jPoint ];
}

ID
QuadraticTetra::faceToPoint ( ID const& iFace, ID const& jPoint )
// faceToPoint(i,j) = localId of jth point on ith local face
{
    static const ID _faceToPoint[ 6 * S_numFaces ] =
    {
        0, 2, 1, 6, 5, 4,
        0, 1, 3, 4, 8, 7,
        1, 2, 3, 5, 9, 8,
        0, 3, 2, 7, 9, 6
    };
    ASSERT_BD ( jPoint < 6 ) ;
    ASSERT_BD ( iFace < S_numFaces ) ;
    return _faceToPoint[ 6 * iFace + jPoint ];
}


/*
          -- LinearHexa


        7-------6
       /.      /|
      / .     / |
     4_______5  |
     |  .    |  |
     |  3....|..2
     | .     | /
     |.      |/
     0_______1
*/
ID
LinearHexa::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 2 * S_numEdges ] =
    {
        0, 1,
        1, 2,
        2, 3,
        3, 0,
        0, 4,
        1, 5,
        2, 6,
        3, 7,
        4, 5,
        5, 6,
        6, 7,
        7, 4
    };
    ASSERT_BD ( jPoint < 2 ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 2 * iEdge + jPoint ];
}

ID
LinearHexa::faceToPoint ( ID const& iFace, ID const& jPoint )
// faceToPoint(i,j) = localId of jth point on ith local face
{
    static const ID _faceToPoint[ 4 * S_numFaces ] =
    {
        0, 3, 2, 1,
        0, 4, 7, 3,
        0, 1, 5, 4,
        1, 2, 6, 5,
        2, 3, 7, 6,
        4, 5, 6, 7
    };
    ASSERT_BD ( jPoint < 4 ) ;
    ASSERT_BD ( iFace < S_numFaces ) ;
    return _faceToPoint[ 4 * iFace + jPoint ];
}




/*
          -- QuadraticHexa
*/
ID
QuadraticHexa::edgeToPoint ( ID const& iEdge, ID const& jPoint )
// edgeToPoint(i,j) = localId of jth point on ith local edge
{
    static const ID _edgeToPoint[ 3 * S_numEdges ] =
    {
        0, 1, 8,
        1, 2, 9,
        2, 3, 10,
        3, 0, 11,
        0, 4, 12,
        1, 5, 13,
        2, 6, 14,
        3, 7, 15,
        4, 5, 16,
        5, 6, 17,
        6, 7, 18,
        7, 4, 19
    };
    ASSERT_BD ( ( jPoint < 3 ) ) ;
    ASSERT_BD ( iEdge < S_numEdges ) ;
    return _edgeToPoint[ 3 * iEdge + jPoint ];
}

ID
QuadraticHexa::faceToPoint ( ID const& iFace, ID const& jPoint )
// faceToPoint(i,j) = localId of jth point on ith local face
{
    static const ID _faceToPoint[ 9 * S_numFaces ] =
    {
        0, 3, 2, 1, 11, 10, 9, 8, 20,
        0, 4, 7, 3, 12, 19, 15, 11, 21,
        0, 1, 5, 4, 8, 13, 16, 12, 22,
        1, 2, 6, 5, 9, 14, 17, 13, 23,
        2, 3, 7, 6, 10, 15, 18, 14, 24,
        4, 5, 6, 7, 16, 17, 18, 19, 25
    };
    ASSERT_BD ( jPoint < 9 ) ;
    ASSERT_BD ( iFace < S_numFaces ) ;
    return _faceToPoint[ 9 * iFace + jPoint ];
}

}
