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
#include <life/lifefem/quadRule.hpp>
#include <life/lifefem/refFEScalar.hpp>
#include <life/lifefem/refFEHdiv.hpp>
#include <life/lifefem/refFEHybrid.hpp>
#include <life/lifefem/geoMap.hpp>

namespace LifeV
{
//! UNKNOWN means: replace me with the degree of exactness of the quadrature rule !
const size_t UNKNOWN = size_t( -1 );




/*======================================================================
 *
 *                          Quadrature Rules on Nodes
 *
 *=======================================================================*/
//! total number of quadrature rules on segments

const size_t NB_QUAD_RULE_NODE = 3;
//! id of the quadrature rules on nodes
const int QUAD_RULE_NODE_1PT = 1;


static const QuadPoint pt_node_1pt[ 1 ] =
{
    QuadPoint( 0., 1. )
};
const QuadRule quadRuleNode1pt( pt_node_1pt,
                                QUAD_RULE_NODE_1PT,
                                "Gauss Legendre 1 point on a node", POINT, 1, 1 );


/*======================================================================
 *
 *                          Quadrature Rules on segments
 *
 *=======================================================================*/
//! total number of quadrature rules on segments
const size_t NB_QUAD_RULE_SEG = 3;
//! id of the quadrature rules on segments
const size_t QUAD_RULE_SEG_1PT = 1;
const size_t QUAD_RULE_SEG_2PT = 2;
const size_t QUAD_RULE_SEG_3PT = 3;
//----------------------------------------------------------------------

static const QuadPoint pt_seg_1pt[ 1 ] =
{
    QuadPoint( 0.5, 1. )
};
const QuadRule quadRuleSeg1pt( pt_seg_1pt,
                               QUAD_RULE_SEG_1PT,
                               "Gauss Legendre 1 point on a segment", LINE, 1, 1 );
//
const QuadRule quadRuleDummy( pt_seg_1pt, 1, "Dummy quadrature rule", LINE, 1, 1 );
//----------------------------------------------------------------------
const Real q2ptx1 = ( 1 - sqrt( 1. / 3. ) ) / 2., q2ptx2 = ( 1 + sqrt( 1. / 3. ) ) / 2.;
const Real q2ptw1 = 0.5, q2ptw2 = 0.5;

static const QuadPoint pt_seg_2pt[ 2 ] =
{
    QuadPoint( q2ptx1 , q2ptw1 ),
    QuadPoint( q2ptx2 , q2ptw2 )
};
const QuadRule quadRuleSeg2pt( pt_seg_2pt,
                               QUAD_RULE_SEG_2PT,
                               "Gauss Legendre 2 points on a segment", LINE, 2, 3 );
//----------------------------------------------------------------------
const Real q3ptx1 = 0.5, q3ptx2 = ( 1 - sqrt( 3. / 5. ) ) / 2., q3ptx3 = ( 1 + sqrt( 3. / 5. ) ) / 2.;
const Real q3ptw1 = 8. / 18., q3ptw2 = 5. / 18., q3ptw3 = 5. / 18.;

static const QuadPoint pt_seg_3pt[ 3 ] =
{
    QuadPoint( q3ptx1, q3ptw1 ),
    QuadPoint( q3ptx2, q3ptw2 ),
    QuadPoint( q3ptx3, q3ptw3 )
};

const QuadRule quadRuleSeg3pt( pt_seg_3pt,
                               QUAD_RULE_SEG_3PT,
                               "Gauss Legendre 3 points on a segment", LINE, 3, 5 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on segments
  ----------------------------------------------------------------------*/
static const QuadRule quad_rule_seg[ NB_QUAD_RULE_SEG ] =
{
    quadRuleSeg1pt,
    quadRuleSeg2pt,
    quadRuleSeg3pt
};
/*======================================================================
 *
 *                     Quadrature Rules 2D on triangles
 *
 *=======================================================================*/
//! total number of quadrature rules in 2D on triangle
#define NB_QUAD_RULE_TRIA 5
//! id of the quadrature rules on triangles
#define QUAD_RULE_TRIA_1PT     1
#define QUAD_RULE_TRIA_3PT     2
#define QUAD_RULE_TRIA_4PT     3
#define QUAD_RULE_TRIA_6PT     4
#define QUAD_RULE_TRIA_7PT     5
//----------------------------------------------------------------------

static const QuadPoint pt_tria_1pt[ 1 ] =
{
    QuadPoint( 1./3., 1./3., 1./2. )
};
const QuadRule quadRuleTria1pt( pt_tria_1pt,
                                QUAD_RULE_TRIA_1PT,
                                "Quadrature rule 1 point on a triangle", TRIANGLE, 1, 1 );
//----------------------------------------------------------------------
static const QuadPoint pt_tria_3pt[ 3 ] =
{
    QuadPoint( 0.5, 0. , 1. / 6. ),
    QuadPoint( 0. , 0.5, 1. / 6. ),
    QuadPoint( 0.5, 0.5, 1. / 6. )
};
const QuadRule quadRuleTria3pt( pt_tria_3pt,
                                QUAD_RULE_TRIA_3PT,
                                "Quadrature rule 3 points on a triangle", TRIANGLE, 3, 2 );
//----------------------------------------------------------------------
// 4 points Integration rule for triangle (Ref. e.g. Comincioli pag. 234) D of Ex = 3
const Real t4pt_xb1 = 3. / 5.,
                      t4pt_xb2 =  1. / 5.,
                                  t4pt_w1  = 25. / 96.,
                                             t4pt_w2  = -9. / 32.,
                                                        t4pt_a   =  1. / 3.;

static const QuadPoint pt_tria_4pt[ 4 ] =
{
    QuadPoint( t4pt_xb1, t4pt_xb2, t4pt_w1 ),
    QuadPoint( t4pt_xb2, t4pt_xb1, t4pt_w1 ),
    QuadPoint( t4pt_xb2, t4pt_xb2, t4pt_w1 ),
    QuadPoint( t4pt_a, t4pt_a, t4pt_w2 )
};

const QuadRule quadRuleTria4pt( pt_tria_4pt,
                                QUAD_RULE_TRIA_4PT,
                                "Quadrature rule 4 points on a triangle", TRIANGLE, 4, 3 );
//----------------------------------------------------------------------
// 6 points Integration rule for triangle, D of Ex = 4
// Ref: G.R. Cowper,  Gaussian quadrature formulas for triangles,
//      Internat. J. Numer. Methods Engrg.  7 (1973), 405--408.
const Real t6pt_x1 = 0.091576213509770743;
const Real t6pt_x2 = 0.44594849091596488;
const Real t6pt_w1 = 0.054975871827660933;
const Real t6pt_w2 = 0.11169079483900573;

static const QuadPoint pt_tria_6pt[ 6 ] =
{
    QuadPoint(     t6pt_x1,     t6pt_x1, t6pt_w1 ),
    QuadPoint(     t6pt_x1, 1-2*t6pt_x1, t6pt_w1 ),
    QuadPoint( 1-2*t6pt_x1,     t6pt_x1, t6pt_w1 ),
    QuadPoint(     t6pt_x2,     t6pt_x2, t6pt_w2 ),
    QuadPoint(     t6pt_x2, 1-2*t6pt_x2, t6pt_w2 ),
    QuadPoint( 1-2*t6pt_x2,     t6pt_x2, t6pt_w2 ),
};
const QuadRule quadRuleTria6pt( pt_tria_6pt, QUAD_RULE_TRIA_6PT,
                                "Quadrature rule 6 points on a triangle",
                                TRIANGLE, 6, 4 );
//----------------------------------------------------------------------
// 7 points Integration rule for triangle (Ref. Stroud) D of Ex = 5
const Real t7pt_x0 = 1./3.;
const Real t7pt_x1 = 0.10128650732345633;
const Real t7pt_x2 = 0.47014206410511508;
const Real t7pt_w0 = 0.1125;
const Real t7pt_w1 = 0.062969590272413576;
const Real t7pt_w2 = 0.066197076394253090;

static const QuadPoint pt_tria_7pt[ 7 ] =
{
    QuadPoint(     t7pt_x0,     t7pt_x0, t7pt_w0 ),
    QuadPoint(     t7pt_x1,     t7pt_x1, t7pt_w1 ),
    QuadPoint(     t7pt_x1, 1-2*t7pt_x1, t7pt_w1 ),
    QuadPoint( 1-2*t7pt_x1,     t7pt_x1, t7pt_w1 ),
    QuadPoint(     t7pt_x2,     t7pt_x2, t7pt_w2 ),
    QuadPoint(     t7pt_x2, 1-2*t7pt_x2, t7pt_w2 ),
    QuadPoint( 1-2*t7pt_x2,     t7pt_x2, t7pt_w2 ),
};
const QuadRule quadRuleTria7pt( pt_tria_7pt, QUAD_RULE_TRIA_7PT,
                                "Quadrature rule 7 points on a triangle",
                                TRIANGLE, 7, 5 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on triangle
  ----------------------------------------------------------------------*/
static const QuadRule quad_rule_tria[ NB_QUAD_RULE_TRIA ] =
{
    quadRuleTria1pt,
    quadRuleTria3pt,
    quadRuleTria4pt,
    quadRuleTria6pt,
    quadRuleTria7pt
};
//----------------------------------------------------------------------
/*======================================================================
 *
 *                     Quadrature Rules 2D on quadrangles
 *
 *=======================================================================*/
//! total number of quadrature rules in 2D on quadrangle
#define NB_QUAD_RULE_QUAD 3
//! id of the quadrature rules on quadrangles
#define QUAD_RULE_QUAD_1PT     1
#define QUAD_RULE_QUAD_4PT     2
#define QUAD_RULE_QUAD_9PT     3
//----------------------------------------------------------------------

static const QuadPoint pt_quad_1pt[ 1 ] =
{
    QuadPoint( .5, .5, 1. )
};
const QuadRule quadRuleQuad1pt( pt_quad_1pt,
                                QUAD_RULE_QUAD_1PT,
                                "Quadrature rule 1 point on a quadrangle", QUAD, 1, 1 );
//----------------------------------------------------------------------
static const QuadPoint pt_quad_4pt[ 4 ] =
{
    QuadPoint( q2ptx1, q2ptx1, q2ptw1 * q2ptw1 ),
    QuadPoint( q2ptx1, q2ptx2, q2ptw1 * q2ptw2 ),
    QuadPoint( q2ptx2, q2ptx1, q2ptw2 * q2ptw1 ),
    QuadPoint( q2ptx2, q2ptx2, q2ptw2 * q2ptw2 )
};
const QuadRule quadRuleQuad4pt( pt_quad_4pt,
                                QUAD_RULE_QUAD_4PT,
                                "Quadrature rule 4 points on a quadrangle", QUAD, 4, 3 );
//----------------------------------------------------------------------
// 4 points Integration rule for quadrangle

static const QuadPoint pt_quad_9pt[ 9 ] =
{
    QuadPoint( q3ptx1, q3ptx1, q3ptw1 * q3ptw1 ),
    QuadPoint( q3ptx2, q3ptx1, q3ptw2 * q3ptw1 ),
    QuadPoint( q3ptx3, q3ptx1, q3ptw3 * q3ptw1 ),
    QuadPoint( q3ptx1, q3ptx2, q3ptw1 * q3ptw2 ),
    QuadPoint( q3ptx2, q3ptx2, q3ptw2 * q3ptw2 ),
    QuadPoint( q3ptx3, q3ptx2, q3ptw3 * q3ptw2 ),
    QuadPoint( q3ptx1, q3ptx3, q3ptw1 * q3ptw3 ),
    QuadPoint( q3ptx2, q3ptx3, q3ptw2 * q3ptw3 ),
    QuadPoint( q3ptx3, q3ptx3, q3ptw3 * q3ptw3 )
};

const QuadRule quadRuleQuad9pt( pt_quad_9pt,
                                QUAD_RULE_QUAD_9PT,
                                "Quadrature rule 9 points on a quadrangle", QUAD, 9, 5 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on quadrangle
  ----------------------------------------------------------------------*/
static const QuadRule quad_rule_quad[ NB_QUAD_RULE_QUAD ] =
{
    quadRuleQuad1pt,
    quadRuleQuad4pt,
    quadRuleQuad9pt
};
//----------------------------------------------------------------------
/*======================================================================
 *
 *                     Quadrature Rules 3D on tetraedras
 *
 *=======================================================================*/
//! total number of quadrature rules in 3D on tetraedra
#define NB_QUAD_RULE_TETRA 5
//! id of the quadrature rules on tetraedra
#define QUAD_RULE_TETRA_1PT     1
#define QUAD_RULE_TETRA_4PT     2
#define QUAD_RULE_TETRA_5PT     3
#define QUAD_RULE_TETRA_15PT    4
#define QUAD_RULE_TETRA_64PT    5
//----------------------------------------------------------------------

static const QuadPoint pt_tetra_1pt[ 1 ] =
{
    QuadPoint( 1. / 4., 1. / 4., 1. / 4., 1. / 6. )
};
const QuadRule quadRuleTetra1pt( pt_tetra_1pt,
                                 QUAD_RULE_TETRA_1PT,
                                 "Quadrature rule 1 point on a tetraedra", TETRA, 1, 1 );
//----------------------------------------------------------------------
const Real tet4ptx1 = ( 5. - sqrt( 5. ) ) / 20., tet4ptx2 = ( 5. + 3*sqrt( 5. ) ) / 20.;

static const QuadPoint pt_tetra_4pt[ 4 ] =
{
    QuadPoint( tet4ptx1, tet4ptx1, tet4ptx1, 1. / 24. ),
    QuadPoint( tet4ptx1, tet4ptx1, tet4ptx2, 1. / 24. ),
    QuadPoint( tet4ptx1, tet4ptx2, tet4ptx1, 1. / 24. ),
    QuadPoint( tet4ptx2, tet4ptx1, tet4ptx1, 1. / 24. )
};
const QuadRule quadRuleTetra4pt( pt_tetra_4pt,
                                 QUAD_RULE_TETRA_4PT,
                                 "Quadrature rule 4 points on a tetraedra", TETRA, 4, 2 );
//----------------------------------------------------------------------
// 5 points Integration rule for tetraedra (Ref. e.g. Comincioli pag. 236)
const Real tet5ptx1 = 1. / 6. , tet5ptx2 = 1. / 2., tet5ptx3 = 1. / 4.;

static const QuadPoint pt_tetra_5pt[ 5 ] =
{
    QuadPoint( tet5ptx1, tet5ptx1, tet5ptx1, 9. / 120. ),
    QuadPoint( tet5ptx1, tet5ptx1, tet5ptx2, 9. / 120. ),
    QuadPoint( tet5ptx1, tet5ptx2, tet5ptx1, 9. / 120. ),
    QuadPoint( tet5ptx2, tet5ptx1, tet5ptx1, 9. / 120. ),
    QuadPoint( tet5ptx3, tet5ptx3, tet5ptx3, -16. / 120. )
};

const QuadRule quadRuleTetra5pt( pt_tetra_5pt,
                                 QUAD_RULE_TETRA_5PT,
                                 "Quadrature rule 5 points on a tetraedra", TETRA, 5, 3 );
//
//----------------------------------------------------------------------
//                     15 points integration rule for tetra.
//                   D o E = 5 (Stroud, T3:5-1 pag. 315)
// r
const Real r5 = 0.25;
// s
const Real s5[ 4 ] =
{
    0.09197107805272303, 0.3197936278296299
}
; // (7 \mp \sqrt(15))/34
// t
const Real t5[ 4 ] =
{
    0.7240867658418310, 0.04061911651111023
}
; // (13 \pm 3*sqrt(15))/34
// u
const Real u5 = 0.05635083268962915; // (10-2*sqrt(15))/40
// v
const Real v5 = 0.4436491673103708; // (10+2*sqrt(15))/40
// A
const Real A5 = 0.01975308641975309; // 16/135*1/6
// B
const Real B5[ 2 ] =
{
    0.01198951396316977, 0.01151136787104540
}
; // 1/6*(2665 \pm 14*sqrt(15))/37800
// C
const Real C5 = 0.008818342151675485; // 20/378*1/6
//
static const QuadPoint pt_tetra_15pt[ 15 ] =
{
    QuadPoint( r5, r5, r5, A5 ),
    QuadPoint( s5[ 0 ], s5[ 0 ], s5[ 0 ], B5[ 0 ] ),
    QuadPoint( t5[ 0 ], s5[ 0 ], s5[ 0 ], B5[ 0 ] ),
    QuadPoint( s5[ 0 ], t5[ 0 ], s5[ 0 ], B5[ 0 ] ),
    QuadPoint( s5[ 0 ], s5[ 0 ], t5[ 0 ], B5[ 0 ] ),
    QuadPoint( s5[ 1 ], s5[ 1 ], s5[ 1 ], B5[ 1 ] ),
    QuadPoint( t5[ 1 ], s5[ 1 ], s5[ 1 ], B5[ 1 ] ),
    QuadPoint( s5[ 1 ], t5[ 1 ], s5[ 1 ], B5[ 1 ] ),
    QuadPoint( s5[ 1 ], s5[ 1 ], t5[ 1 ], B5[ 1 ] ),
    QuadPoint( u5, u5, v5, C5 ),
    QuadPoint( u5, v5, u5, C5 ),
    QuadPoint( v5, u5, u5, C5 ),
    QuadPoint( v5, v5, u5, C5 ),
    QuadPoint( v5, u5, v5, C5 ),
    QuadPoint( u5, v5, v5, C5 )
};
//
const QuadRule quadRuleTetra15pt( pt_tetra_15pt,
                                  QUAD_RULE_TETRA_15PT,
                                  "Quadrature rule 15 points on a tetraedra",
                                  TETRA, 15, 5 );
//----------------------------------------------------------------------
//                     64 points integration rule for tetra.
//                   D o E = 7 (Stroud, T3:7-1 pag. 315)
//
// t
const Real t[ 4 ] =
{
    0.0485005494, 0.2386007376, 0.5170472951, 0.7958514179
};
// s
const Real s[ 4 ] =
{
    0.0571041961, 0.2768430136, 0.5835904324, 0.8602401357
};
// r
const Real r[ 4 ] =
{
    0.0694318422, 0.3300094782, 0.6699905218, 0.9305681558
};
// A
const Real A[ 4 ] =
{
    0.1739274226, 0.3260725774, 0.3260725774, 0.1739274226
};
// B
const Real B[ 4 ] =
{
    0.1355069134, 0.2034645680, 0.1298475476, 0.0311809709
};
// C
const Real C[ 4 ] =
{
    0.1108884156, 0.1434587898, 0.0686338872, 0.0103522407
};
//
/*
for (i=0;i<4;i++){
  for (j=0;j<4;j++){
    for (k=0;k<4;k++){
      pt_tetra_64pt[k+4*j+16*i]=
 QuadPoint(t[k], s[j]*(1-t[k]), r[i]*(1-s[j])*(1-t[k]), A[i]*B[j]*C[k]);
    }
  }
}
*/
static const QuadPoint pt_tetra_64pt[ 64 ] =
{
    QuadPoint( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 0 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 0 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 0 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 0 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 1 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 1 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 1 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 1 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 2 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 2 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 2 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 2 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 0 ] * B[ 3 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 0 ] * B[ 3 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 0 ] * B[ 3 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 0 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 0 ] * B[ 3 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 0 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 0 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 0 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 0 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 1 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 1 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 1 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 1 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 2 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 2 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 2 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 2 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 1 ] * B[ 3 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 1 ] * B[ 3 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 1 ] * B[ 3 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 1 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 1 ] * B[ 3 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 0 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 0 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 0 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 0 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 1 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 1 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 1 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 1 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 2 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 2 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 2 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 2 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 2 ] * B[ 3 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 2 ] * B[ 3 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 2 ] * B[ 3 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 2 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 2 ] * B[ 3 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 0 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 0 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 0 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 0 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 0 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 0 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 0 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 0 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 0 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 1 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 1 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 1 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 1 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 1 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 1 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 1 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 1 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 1 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 2 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 2 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 2 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 2 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 2 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 2 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 2 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 2 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 2 ] * C[ 3 ] ),
    QuadPoint( t[ 0 ], s[ 3 ] * ( 1 - t[ 0 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 0 ] ), A[ 3 ] * B[ 3 ] * C[ 0 ] ),
    QuadPoint( t[ 1 ], s[ 3 ] * ( 1 - t[ 1 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 1 ] ), A[ 3 ] * B[ 3 ] * C[ 1 ] ),
    QuadPoint( t[ 2 ], s[ 3 ] * ( 1 - t[ 2 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 2 ] ), A[ 3 ] * B[ 3 ] * C[ 2 ] ),
    QuadPoint( t[ 3 ], s[ 3 ] * ( 1 - t[ 3 ] ), r[ 3 ] * ( 1 - s[ 3 ] ) * ( 1 - t[ 3 ] ), A[ 3 ] * B[ 3 ] * C[ 3 ] )
};
//
const QuadRule quadRuleTetra64pt( pt_tetra_64pt,
                                  QUAD_RULE_TETRA_64PT,
                                  "Quadrature rule 64 points on a tetraedra",
                                  TETRA, 64, 7 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on tetraedra
  ----------------------------------------------------------------------*/
static const QuadRule quad_rule_tetra[ NB_QUAD_RULE_TETRA ] =
{
    quadRuleTetra1pt,
    quadRuleTetra4pt,
    quadRuleTetra5pt,
    quadRuleTetra15pt,
    quadRuleTetra64pt
};
//----------------------------------------------------------------------
/*======================================================================
 *
 *                     Quadrature Rules 3D on hexaedras
 *
 *=======================================================================*/
//! total number of quadrature rules in 3D on hexa
#define NB_QUAD_RULE_HEXA 2
//! id of the quadrature rules on quadrangles
#define QUAD_RULE_HEXA_1PT     1
#define QUAD_RULE_HEXA_8PT     2
//----------------------------------------------------------------------

static const QuadPoint pt_hexa_1pt[ 1 ] =
{
    QuadPoint( .5, .5, .5, 1. )
};
const QuadRule quadRuleHexa1pt( pt_hexa_1pt,
                                QUAD_RULE_HEXA_1PT,
                                "Quadrature rule 1 point on a hexa", HEXA, 1, 1 );
//----------------------------------------------------------------------
static const QuadPoint pt_hexa_8pt[ 8 ] =
{
    QuadPoint( q2ptx1, q2ptx1, q2ptx1, q2ptw1 * q2ptw1 * q2ptw1 ),
    QuadPoint( q2ptx1, q2ptx2, q2ptx1, q2ptw1 * q2ptw2 * q2ptw1 ),
    QuadPoint( q2ptx2, q2ptx1, q2ptx1, q2ptw2 * q2ptw1 * q2ptw1 ),
    QuadPoint( q2ptx2, q2ptx2, q2ptx1, q2ptw2 * q2ptw2 * q2ptw1 ),
    QuadPoint( q2ptx1, q2ptx1, q2ptx2, q2ptw1 * q2ptw1 * q2ptw2 ),
    QuadPoint( q2ptx1, q2ptx2, q2ptx2, q2ptw1 * q2ptw2 * q2ptw2 ),
    QuadPoint( q2ptx2, q2ptx1, q2ptx2, q2ptw2 * q2ptw1 * q2ptw2 ),
    QuadPoint( q2ptx2, q2ptx2, q2ptx2, q2ptw2 * q2ptw2 * q2ptw2 )
};
const QuadRule quadRuleHexa8pt( pt_hexa_8pt,
                                QUAD_RULE_HEXA_8PT,
                                "Quadrature rule 8 points on a hexa", HEXA, 8, 3 );
/*----------------------------------------------------------------------
  Set of all quadrature rules on hexa
  ----------------------------------------------------------------------*/
static const QuadRule quad_rule_hexa[ NB_QUAD_RULE_HEXA ] =
{
    quadRuleHexa1pt,
    quadRuleHexa8pt
};



//----------------------------------------------------------------------
//
//                       GeoMaps
//
//----------------------------------------------------------------------

const GeoMap geoLinearNode( "Mapping of a point", POINT,
                            1, 1,
                            fct_P0_0D, derfct_P0_0D, der2fct_P0_0D,
                            refcoor_P0_0D,
                            ( GeoMap* ) NULL );

const GeoMap geoLinearSeg( "Linear mapping on a segment", LINE,
                           2, 1,
                           fct_P1_1D, derfct_P1_1D, der2fct_P1_1D,
                           refcoor_P1_1D,
                           &geoLinearNode );

const GeoMap geoLinearTria( "Linear mapping on a triangle", TRIANGLE,
                            3, 2,
                            fct_P1_2D, derfct_P1_2D, der2fct_P1_2D,
                            refcoor_P1_2D,
                            &geoLinearSeg );

const GeoMap geoBilinearQuad( "Bilinear mapping on a quadrangle", QUAD,
                              4, 2,
                              fct_Q1_2D, derfct_Q1_2D, der2fct_Q1_2D,
                              refcoor_Q1_2D,
                              &geoLinearSeg );

const GeoMap geoLinearTetra( "Linear mapping on a tetraedra", TETRA,
                             4, 3,
                             fct_P1_3D, derfct_P1_3D, der2fct_P1_3D,
                             refcoor_P1_3D,
                             &geoLinearTria );

const GeoMap geoBilinearHexa( "Bilinear mapping on an hexaedra", HEXA,
                              8, 3,
                              fct_Q1_3D, derfct_Q1_3D, der2fct_Q1_3D,
                              refcoor_Q1_3D,
                              &geoBilinearQuad );



//======================================================================
//
//                            P0  (0D)
//
//======================================================================
/*
                           1
*/
Real fct1_P0_0D( const GeoVector& )
{
    return 1.;
}
Real derfct1_P0_0D( const GeoVector& )
{
    return 0.;
}
Real der2fct1_P0_0D( const GeoVector& )
{
    return 0.;
}
//======================================================================
//
//                            P1  (1D)
//
//======================================================================
/*
                           1-----2
*/
Real fct1_P1_1D( const GeoVector& v )
{
    return 1 - v[0];
}
Real fct2_P1_1D( const GeoVector& v )
{
    return v[0];
}

Real derfct1_1_P1_1D( const GeoVector& )
{
    return -1;
}
Real derfct2_1_P1_1D( const GeoVector& )
{
    return 1;
}

Real der2fct1_P1_1D( const GeoVector& )
{
    return 0;
}

//======================================================================
//
//                            P2  (1D)
//
//======================================================================
/*
                           1--3--2
*/
Real fct1_P2_1D( const GeoVector& v )
{
    return 2. * ( v[0] - 1. ) * ( v[0] - 0.5 );
}
Real fct3_P2_1D( const GeoVector& v )
{
    return 4. * v[0] * ( 1. - v[0] );
}
Real fct2_P2_1D( const GeoVector& v )
{
    return 2. * v[0] * ( v[0] - 0.5 );
}

Real derfct1_1_P2_1D( const GeoVector& v )
{
    return 4. * v[0] - 3.;
}
Real derfct3_1_P2_1D( const GeoVector& v )
{
    return -8. * v[0] + 4.;
}
Real derfct2_1_P2_1D( const GeoVector& v )
{
    return 4. * v[0] - 1.;
}

Real der2fct1_11_P2_1D( const GeoVector& )
{
    return 4;
}
Real der2fct3_11_P2_1D(  const GeoVector& )
{
    return -8;
}
Real der2fct2_11_P2_1D(  const GeoVector& )
{
    return 4;
}


//======================================================================
//
//                            P0  (2D)
//
//======================================================================
/*

                           |\
                           | \
                           | 1\
                            ---
*/
Real fct1_P0_2D( const GeoVector& )
{
    return 1. ;
}   //check this : 1. or 2. (\int fct1 = 0.5 or 1.)  ???
// First and Second derivatives are both equal (to 0).
Real derfct1_P0_2D( const GeoVector& )
{
    return 0. ;
}
Real der2fct1_P0_2D( const GeoVector& )
{
    return 0. ;
}

//======================================================================
//
//                            P1  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |  \
                           1---2
*/
Real fct1_P1_2D( const GeoVector& v )
{
    return ( 1. - v[0] - v[1] );
}
Real fct2_P1_2D( const GeoVector& v )
{
    return v[0] ;
}
Real fct3_P1_2D( const GeoVector& v )
{
    return v[1] ;
}

Real derfct1_1_P1_2D( const GeoVector& )
{
    return -1 ;
}
Real derfct1_2_P1_2D( const GeoVector& )
{
    return -1 ;
}
Real derfct2_1_P1_2D( const GeoVector& )
{
    return 1 ;
}
Real derfct2_2_P1_2D( const GeoVector& )
{
    return 0 ;
}
Real derfct3_1_P1_2D( const GeoVector& )
{
    return 0 ;
}
Real derfct3_2_P1_2D( const GeoVector& )
{
    return 1 ;
}

// Second derivatives
Real der2fctx_xx_P1_2D( const GeoVector& )
{
    return 0;
}



//======================================================================
//
//                            P2  (2D)
//
//======================================================================
/*
                           3
                           |\
                           6 5
                           |  \
                           1-4-2
*/
Real fct1_P2_2D( const GeoVector& v )
{
    return ( 1 -v[0] - v[1] ) * ( 1 - v[0] - v[0] - v[1] - v[1] );
}
Real fct2_P2_2D( const GeoVector& v )
{
    return -v[0] * ( 1 - v[0] - v[0] );
}
Real fct3_P2_2D( const GeoVector& v )
{
    return -v[1] * ( 1 - v[1] - v[1] );
}
Real fct4_P2_2D( const GeoVector& v )
{
    return 4 * v[0] * ( 1 - v[0] - v[1] );
}
Real fct5_P2_2D( const GeoVector& v )
{
    return 4 * v[0] * v[1];
}
Real fct6_P2_2D( const GeoVector& v )
{
    return 4 * v[1] * ( 1 - v[0] - v[1] );
}

Real derfct1_1_P2_2D( const GeoVector& v )
{
    return 4 * ( v[0] + v[1] ) - 3;
}
Real derfct1_2_P2_2D( const GeoVector& v )
{
    return 4 * ( v[0] + v[1] ) - 3;
}
Real derfct2_1_P2_2D( const GeoVector& v )
{
    return 4 * v[0] - 1;
}
Real derfct2_2_P2_2D( const GeoVector& )
{
    return 0;
}
Real derfct3_1_P2_2D( const GeoVector& )
{
    return 0;
}
Real derfct3_2_P2_2D( const GeoVector& v )
{
    return 4 * v[1] - 1;
}
Real derfct4_1_P2_2D( const GeoVector& v )
{
    return 4 * ( 1 - v[0] - v[0] - v[1] );
}
Real derfct4_2_P2_2D( const GeoVector& v )
{
    return -4 * v[0];
}
Real derfct5_1_P2_2D( const GeoVector& v )
{
    return 4 * v[1];
}
Real derfct5_2_P2_2D( const GeoVector& v )
{
    return 4 * v[0];
}
Real derfct6_1_P2_2D( const GeoVector& v )
{
    return -4 * v[1];
}
Real derfct6_2_P2_2D( const GeoVector& v )
{
    return 4 * ( 1 - v[0] - v[1] - v[1] );
}

Real der2fct1_11_P2_2D( const GeoVector& )
{
    return 4;
}
Real der2fct1_12_P2_2D( const GeoVector& )
{
    return 4;
}
Real der2fct1_21_P2_2D( const GeoVector& )
{
    return 4;
}
Real der2fct1_22_P2_2D( const GeoVector& )
{
    return 4;
}

Real der2fct2_11_P2_2D( const GeoVector& )
{
    return 4;
}
Real der2fct2_12_P2_2D( const GeoVector& )
{
    return 0;
}
Real der2fct2_21_P2_2D( const GeoVector& )
{
    return 0;
}
Real der2fct2_22_P2_2D( const GeoVector& )
{
    return 0;
}

Real der2fct3_11_P2_2D( const GeoVector& )
{
    return 0;
}
Real der2fct3_12_P2_2D( const GeoVector& )
{
    return 0;
}
Real der2fct3_21_P2_2D( const GeoVector& )
{
    return 0;
}
Real der2fct3_22_P2_2D( const GeoVector& )
{
    return 4;
}

Real der2fct4_11_P2_2D( const GeoVector& )
{
    return -8;
}
Real der2fct4_12_P2_2D( const GeoVector& )
{
    return -4;
}
Real der2fct4_21_P2_2D( const GeoVector& )
{
    return -4;
}
Real der2fct4_22_P2_2D( const GeoVector& )
{
    return 0;
}

Real der2fct5_11_P2_2D( const GeoVector& )
{
    return 0;
}
Real der2fct5_12_P2_2D( const GeoVector& )
{
    return 4;
}
Real der2fct5_21_P2_2D( const GeoVector& )
{
    return 4;
}
Real der2fct5_22_P2_2D( const GeoVector& )
{
    return 0;
}

Real der2fct6_11_P2_2D( const GeoVector& )
{
    return 0;
}
Real der2fct6_12_P2_2D( const GeoVector& )
{
    return -4;
}
Real der2fct6_21_P2_2D( const GeoVector& )
{
    return -4;
}
Real der2fct6_22_P2_2D( const GeoVector& )
{
    return -8;
}
//======================================================================
//
//                            Q0  (2D)
//
//======================================================================
/*
                            -------
                           |       |
                           |   1   |
                           |       |
                            -------
*/
Real fct1_Q0_2D( const GeoVector& )
{
    return 1. ;
}
Real derfct1_Q0_2D( const GeoVector& v )
{
    return 0. ;
}
// The second derivative is equal to the first : both are equal to 0.
Real der2fct1_Q0_2D( const GeoVector& v )
{
    return 0. ;
}

//======================================================================
//
//                            Q1  (2D)
//
//======================================================================
/*
                           4-------3
                           |       |
                           |       |
                           |       |
                           1-------2
*/
Real fct1_Q1_2D( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] );
}
Real fct2_Q1_2D( const GeoVector& v )
{
    return ( 1. - v[1] ) * v[0];
}
Real fct3_Q1_2D( const GeoVector& v )
{
    return v[0] * v[1];
}
Real fct4_Q1_2D( const GeoVector& v )
{
    return v[1] * ( 1. - v[0] );
}

Real derfct1_1_Q1_2D( const GeoVector& v )
{
    return -( 1. - v[1] );
}
Real derfct1_2_Q1_2D( const GeoVector& v )
{
    return -( 1. - v[0] );
}
Real derfct2_1_Q1_2D( const GeoVector& v )
{
    return ( 1. - v[1] );
}
Real derfct2_2_Q1_2D( const GeoVector& v )
{
    return -v[0];
}
Real derfct3_1_Q1_2D( const GeoVector& v )
{
    return v[1];
}
Real derfct3_2_Q1_2D( const GeoVector& v )
{
    return v[0];
}
Real derfct4_1_Q1_2D( const GeoVector& v )
{
    return -v[1];
}
Real derfct4_2_Q1_2D( const GeoVector& v )
{
    return ( 1. - v[0] );
}

// Second derivatives
Real der2fctx_xx_Q1_2D( const GeoVector& )
{
    return 0;
}
//======================================================================
//
//                            Q2  (2D)
//
//======================================================================
/*
                           4---7---3
                           |       |
                           8   9   6
                           |       |
                           1---5---2
*/
Real fct1_Q2_2D( const GeoVector& v )
{
    return 4. * ( 1 - v[0] ) * ( 0.5 - v[0] ) * ( 1 - v[1] ) * ( 0.5 - v[1] );
}
Real fct5_Q2_2D( const GeoVector& v )
{
    return 8. * v[0] * ( 1 - v[0] ) * ( 1 - v[1] ) * ( 0.5 - v[1] );
}
Real fct2_Q2_2D( const GeoVector& v )
{
    return 4. * v[0] * ( v[0] - 0.5 ) * ( 1 - v[1] ) * ( 0.5 - v[1] );
}
Real fct6_Q2_2D( const GeoVector& v )
{
    return 8. * v[0] * ( v[0] - 0.5 ) * v[1] * ( 1 - v[1] );
}
Real fct3_Q2_2D( const GeoVector& v )
{
    return 4. * v[0] * ( v[0] - 0.5 ) * v[1] * ( v[1] - 0.5 );
}
Real fct7_Q2_2D( const GeoVector& v )
{
    return 8. * v[0] * ( 1 - v[0] ) * v[1] * ( v[1] - 0.5 );
}
Real fct4_Q2_2D( const GeoVector& v )
{
    return 4. * ( 1 - v[0] ) * ( 0.5 - v[0] ) * v[1] * ( v[1] - 0.5 );
}
Real fct8_Q2_2D( const GeoVector& v )
{
    return 8. * ( 0.5 - v[0] ) * ( 1 - v[0] ) * v[1] * ( 1 - v[1] );
}
Real fct9_Q2_2D( const GeoVector& v )
{
    return 16. * v[0] * ( 1 - v[0] ) * v[1] * ( 1 - v[1] );
}

Real derfct1_1_Q2_2D( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * ( 4. * v[0] - 3. );
}
Real derfct1_2_Q2_2D( const GeoVector& v )
{
    return ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * ( 4. * v[1] - 3. );
}
Real derfct5_1_Q2_2D( const GeoVector& v )
{
    return -4. * ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * ( 2. * v[0] - 1. );
}
Real derfct5_2_Q2_2D( const GeoVector& v )
{
    return -4. * v[0] * ( v[0] - 1. ) * ( 4. * v[1] - 3. );
}
Real derfct2_1_Q2_2D( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * ( 4. * v[0] - 1. );
}
Real derfct2_2_Q2_2D( const GeoVector& v )
{
    return v[0] * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 3. );
}
Real derfct6_1_Q2_2D( const GeoVector& v )
{
    return -4. * v[1] * ( 4. * v[0] - 1. ) * ( v[1] - 1. );
}
Real derfct6_2_Q2_2D( const GeoVector& v )
{
    return -4. * v[0] * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct3_1_Q2_2D( const GeoVector& v )
{
    return v[1] * ( 4. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct3_2_Q2_2D( const GeoVector& v )
{
    return v[0] * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real derfct7_1_Q2_2D( const GeoVector& v )
{
    return -4. * v[1] * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct7_2_Q2_2D( const GeoVector& v )
{
    return -4. * v[0] * ( v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real derfct4_1_Q2_2D( const GeoVector& v )
{
    return v[1] * ( 4. * v[0] - 3. ) * ( 2. * v[1] - 1. );
}
Real derfct4_2_Q2_2D( const GeoVector& v )
{
    return ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real derfct8_1_Q2_2D( const GeoVector& v )
{
    return -4. * v[1] * ( 4. * v[0] - 3. ) * ( v[1] - 1. );
}
Real derfct8_2_Q2_2D( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real derfct9_1_Q2_2D( const GeoVector& v )
{
    return 16. * v[1] * ( 2. * v[0] - 1. ) * ( v[1] - 1. );
}
Real derfct9_2_Q2_2D( const GeoVector& v )
{
    return 16. * v[0] * ( v[0] - 1. ) * ( 2. * v[1] - 1. );
}

Real der2fct1_11_Q2_2D( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * 4.;
}
Real der2fct1_12_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[1] - 3. ) * ( 4. * v[0] - 3. );
}
Real der2fct1_21_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[1] - 3. ) * ( 4. * v[0] - 3. );
}
Real der2fct1_22_Q2_2D( const GeoVector& v )
{
    return ( 2. * v[0] - 1. ) * ( v[0] - 1. ) * 4.;
}

Real der2fct5_11_Q2_2D( const GeoVector& v )
{
    return -8. * ( 2. * v[1] - 1. ) * ( v[1] - 1. );
}
Real der2fct5_12_Q2_2D( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1 ) * ( 4. * v[1] - 3 );
}
Real der2fct5_21_Q2_2D( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1 ) * ( 4. * v[1] - 3 );
    ;
}
Real der2fct5_22_Q2_2D( const GeoVector& v )
{
    return -16. * v[0] * ( v[0] - 1. );
}

Real der2fct2_11_Q2_2D( const GeoVector& v )
{
    return ( 2. * v[1] - 1. ) * ( v[1] - 1. ) * 4.;
}
Real der2fct2_12_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[0] - 1 ) * ( 4. * v[1] - 3. );
}
Real der2fct2_21_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[1] - 3. ) * ( 4. * v[0] - 1. );
}
Real der2fct2_22_Q2_2D( const GeoVector& v )
{
    return v[0] * ( 2. * v[0] - 1. ) * 4.;
}

Real der2fct6_11_Q2_2D( const GeoVector& v )
{
    return -16. * v[1] * ( v[1] - 1. );
}
Real der2fct6_12_Q2_2D( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct6_21_Q2_2D( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct6_22_Q2_2D( const GeoVector& v )
{
    return -8. * v[0] * ( 2. * v[0] - 1. );
}

Real der2fct3_11_Q2_2D( const GeoVector& v )
{
    return 4. * v[1] * ( 2. * v[1] - 1. );
}
Real der2fct3_12_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct3_21_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct3_22_Q2_2D( const GeoVector& v )
{
    return 4. * v[0] * ( 2. * v[0] - 1. );
}

Real der2fct7_11_Q2_2D( const GeoVector& v )
{
    return -8. * v[1] * ( 2. * v[1] - 1. );
}
Real der2fct7_12_Q2_2D( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct7_21_Q2_2D( const GeoVector& v )
{
    return -4. * ( 2. * v[0] - 1. ) * ( 4. * v[1] - 1. );
}
Real der2fct7_22_Q2_2D( const GeoVector& v )
{
    return -16. * v[0] * ( v[0] - 1. );
}

Real der2fct4_11_Q2_2D( const GeoVector& v )
{
    return 4. * v[1] * ( 2. * v[1] - 1. );
}
Real der2fct4_12_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[0] - 3. ) * ( 4. * v[1] - 1. );
}
Real der2fct4_21_Q2_2D( const GeoVector& v )
{
    return ( 4. * v[0] - 3. ) * ( 4. * v[1] - 1. );
}
Real der2fct4_22_Q2_2D( const GeoVector& v )
{
    return 4. * ( 2. * v[0] - 1. ) * ( v[0] - 1. );
}

Real der2fct8_11_Q2_2D( const GeoVector& v )
{
    return -16. * v[1] * ( v[1] - 1. );
}
Real der2fct8_12_Q2_2D( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 3. ) * ( 2. * v[1] - 1. );
}
Real der2fct8_21_Q2_2D( const GeoVector& v )
{
    return -4. * ( 4. * v[0] - 3. ) * ( 2. * v[1] - 1. );
}
Real der2fct8_22_Q2_2D( const GeoVector& v )
{
    return -8. * ( 2. * v[0] - 1. ) * ( v[0] - 1. );
}

Real der2fct9_11_Q2_2D( const GeoVector& v )
{
    return 32. * v[1] * ( v[1] - 1. );
}
Real der2fct9_12_Q2_2D( const GeoVector& v )
{
    return 16. * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct9_21_Q2_2D( const GeoVector& v )
{
    return 16. * ( 2. * v[0] - 1. ) * ( 2. * v[1] - 1. );
}
Real der2fct9_22_Q2_2D( const GeoVector& v )
{
    return 32. * v[0] * ( v[0] - 1. );
}

//======================================================================
//
//                            P0  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
Real fct1_P0_3D( const GeoVector& )
{
    return 1.;
}

Real derfct1_P0_3D( const GeoVector& )
{
    return 0.;
}

// Second derivatives
Real der2fct1_P0_3D( const GeoVector& )
{
    return 0;
}


//======================================================================
//
//                            P1  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
Real fct1_P1_3D( const GeoVector& v )
{
    return 1 -v[0] - v[1] - v[2];
}
Real fct2_P1_3D( const GeoVector& v )
{
    return v[0];
}
Real fct3_P1_3D( const GeoVector& v )
{
    return v[1];
}
Real fct4_P1_3D( const GeoVector& v )
{
    return v[2];
}

Real derfct1_1_P1_3D( const GeoVector& )
{
    return -1;
}
Real derfct1_2_P1_3D( const GeoVector& )
{
    return -1;
}
Real derfct1_3_P1_3D( const GeoVector& )
{
    return -1;
}
Real derfct2_1_P1_3D( const GeoVector& )
{
    return 1;
}
Real derfct2_2_P1_3D( const GeoVector& )
{
    return 0;
}
Real derfct2_3_P1_3D( const GeoVector& )
{
    return 0;
}
Real derfct3_1_P1_3D( const GeoVector& )
{
    return 0;
}
Real derfct3_2_P1_3D( const GeoVector& )
{
    return 1;
}
Real derfct3_3_P1_3D( const GeoVector& )
{
    return 0;
}
Real derfct4_1_P1_3D( const GeoVector& )
{
    return 0;
}
Real derfct4_2_P1_3D( const GeoVector& )
{
    return 0;
}
Real derfct4_3_P1_3D( const GeoVector& )
{
    return 1;
}

// Second derivatives
Real der2fctx_xx_P1_3D( const GeoVector& )
{
    return 0;
}
//======================================================================
//======================================================================
//
//                            P1bubble  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / . .5 \\
           /.       \!
         1 ----------2
*/
Real fct1_P1bubble_3D( const GeoVector& v )
{
    return 1 -v[0] - v[1] - v[2];
}
Real fct2_P1bubble_3D( const GeoVector& v )
{
    return v[0];
}
Real fct3_P1bubble_3D( const GeoVector& v )
{
    return v[1];
}
Real fct4_P1bubble_3D( const GeoVector& v )
{
    return v[2];
}
Real fct5_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -v[0] - v[1] - v[2] ) * v[0] * v[1] * v[2];
}

Real derfct1_1_P1bubble_3D( const GeoVector& )
{
    return -1;
}
Real derfct1_2_P1bubble_3D( const GeoVector& )
{
    return -1;
}
Real derfct1_3_P1bubble_3D( const GeoVector& )
{
    return -1;
}
Real derfct2_1_P1bubble_3D( const GeoVector& )
{
    return 1;
}
Real derfct2_2_P1bubble_3D( const GeoVector& )
{
    return 0;
}
Real derfct2_3_P1bubble_3D( const GeoVector& )
{
    return 0;
}
Real derfct3_1_P1bubble_3D( const GeoVector& )
{
    return 0;
}
Real derfct3_2_P1bubble_3D( const GeoVector& )
{
    return 1;
}
Real derfct3_3_P1bubble_3D( const GeoVector& )
{
    return 0;
}
Real derfct4_1_P1bubble_3D( const GeoVector& )
{
    return 0;
}
Real derfct4_2_P1bubble_3D( const GeoVector& )
{
    return 0;
}
Real derfct4_3_P1bubble_3D( const GeoVector& )
{
    return 1;
}
Real derfct5_1_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -2 * v[0] - v[1] - v[2] ) * v[1] * v[2];
}
Real derfct5_2_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -v[0] - 2 * v[1] - v[2] ) * v[0] * v[2];
}
Real derfct5_3_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -v[0] - v[1] - 2 * v[2] ) * v[0] * v[1];
}

// Second derivatives
Real der2fctx_xx_P1bubble_3D( const GeoVector& )
{
    return 0;
}
Real der2fct5_11_P1bubble_3D( const GeoVector& v )
{
    return -2 * v[1] * v[2];
}
Real der2fct5_12_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -2 * v[0] - 2 * v[1] - v[2] ) * v[2];
}
Real der2fct5_13_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -2 * v[0] - v[1] - 2 * v[2] ) * v[1];
}
Real der2fct5_21_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -2 * v[0] - 2 * v[1] - v[2] ) * v[2];
}
Real der2fct5_22_P1bubble_3D( const GeoVector& v )
{
    return -2 * v[0] * v[2];
}
Real der2fct5_23_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -v[0] - 2 * v[1] - 2 * v[2] ) * v[0];
}
Real der2fct5_31_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -2 * v[0] - v[1] - 2 * v[2] ) * v[1];
}
Real der2fct5_32_P1bubble_3D( const GeoVector& v )
{
    return ( 1 -v[0] - 2 * v[1] - 2 * v[2] ) * v[0];
}
Real der2fct5_33_P1bubble_3D( const GeoVector& v )
{
    return -2 * v[0] * v[1];
}

//======================================================================
//
//                            P2  (3D)
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7    \6
           /.       \!
         1 -----5----2
*/
Real fct1_P2_3D( const GeoVector& v )
{
    return -( 1 - v[0] - v[1] - v[2] ) * ( 1 - 2 * ( 1 - v[0] - v[1] - v[2] ) );
}
Real fct2_P2_3D( const GeoVector& v )
{
    return -v[0] * ( 1 - 2 * v[0] );
}
Real fct3_P2_3D( const GeoVector& v )
{
    return -v[1] * ( 1 - 2 * v[1] );
}
Real fct4_P2_3D( const GeoVector& v )
{
    return -v[2] * ( 1 - 2 * v[2] );
}
Real fct5_P2_3D( const GeoVector& v )
{
    return 4 * v[0] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct6_P2_3D( const GeoVector& v )
{
    return 4 * v[0] * v[1];
}
Real fct7_P2_3D( const GeoVector& v )
{
    return 4 * v[1] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct8_P2_3D( const GeoVector& v )
{
    return 4 * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct9_P2_3D( const GeoVector& v )
{
    return 4 * v[0] * v[2];
}
Real fct10_P2_3D( const GeoVector& v )
{
    return 4 * v[1] * v[2];
}


Real derfct1_1_P2_3D( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2];
}
Real derfct1_2_P2_3D( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2];
}
Real derfct1_3_P2_3D( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2];
}

Real derfct2_1_P2_3D( const GeoVector& v )
{
    return -1 + 4 * v[0];
}
Real derfct2_2_P2_3D( const GeoVector& )
{
    return 0.;
}
Real derfct2_3_P2_3D( const GeoVector& )
{
    return 0.;
}

Real derfct3_1_P2_3D( const GeoVector& )
{
    return 0.;
}
Real derfct3_2_P2_3D( const GeoVector& v )
{
    return -1 + 4 * v[1];
}
Real derfct3_3_P2_3D( const GeoVector& )
{
    return 0.;
}

Real derfct4_1_P2_3D( const GeoVector& )
{
    return 0.;
}
Real derfct4_2_P2_3D( const GeoVector& )
{
    return 0.;
}
Real derfct4_3_P2_3D( const GeoVector& v )
{
    return -1 + 4 * v[2];
}

Real derfct5_1_P2_3D( const GeoVector& v )
{
    return 4 - 8 * v[0] - 4 * v[1] - 4 * v[2];
}
Real derfct5_2_P2_3D( const GeoVector& v )
{
    return -4 * v[0];
}
Real derfct5_3_P2_3D( const GeoVector& v )
{
    return -4 * v[0];
}

Real derfct6_1_P2_3D( const GeoVector& v )
{
    return 4 * v[1];
}
Real derfct6_2_P2_3D( const GeoVector& v )
{
    return 4 * v[0];
}
Real derfct6_3_P2_3D( const GeoVector& )
{
    return 0.;
}

Real derfct7_1_P2_3D( const GeoVector& v )
{
    return -4 * v[1];
}
Real derfct7_2_P2_3D( const GeoVector& v )
{
    return 4 - 4 * v[0] - 8 * v[1] - 4 * v[2];
}
Real derfct7_3_P2_3D( const GeoVector& v )
{
    return -4 * v[1];
}

Real derfct8_1_P2_3D( const GeoVector& v )
{
    return -4 * v[2];
}
Real derfct8_2_P2_3D( const GeoVector& v )
{
    return -4 * v[2];
}
Real derfct8_3_P2_3D( const GeoVector& v )
{
    return 4 - 4 * v[0] - 4 * v[1] - 8 * v[2];
}

Real derfct9_1_P2_3D( const GeoVector& v )
{
    return 4 * v[2];
}
Real derfct9_2_P2_3D( const GeoVector& )
{
    return 0.;
}
Real derfct9_3_P2_3D( const GeoVector& v )
{
    return 4 * v[0];
}

Real derfct10_1_P2_3D( const GeoVector& )
{
    return 0.;
}
Real derfct10_2_P2_3D( const GeoVector& v )
{
    return 4 * v[2];
}
Real derfct10_3_P2_3D( const GeoVector& v )
{
    return 4 * v[1];
}


Real der2fct1_11_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_12_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_13_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_21_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_22_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_23_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_31_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_32_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct1_33_P2_3D( const GeoVector& )
{
    return 4;
}

Real der2fct2_11_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct2_12_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_13_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_21_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_22_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_23_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_31_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_32_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_33_P2_3D( const GeoVector& )
{
    return 0;
}

Real der2fct3_11_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_12_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_13_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_21_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_22_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct3_23_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_31_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_32_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_33_P2_3D( const GeoVector& )
{
    return 0;
}

Real der2fct4_11_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_12_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_13_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_21_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_22_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_23_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_31_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_32_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_33_P2_3D( const GeoVector& )
{
    return 4;
}

Real der2fct5_11_P2_3D( const GeoVector& )
{
    return -8;
}
Real der2fct5_12_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct5_13_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct5_21_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct5_22_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct5_23_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct5_31_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct5_32_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct5_33_P2_3D( const GeoVector& )
{
    return 0;
}

Real der2fct6_11_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_12_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct6_13_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_21_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct6_22_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_23_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_31_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_32_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_33_P2_3D( const GeoVector& )
{
    return 0;
}

Real der2fct7_11_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct7_12_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct7_13_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct7_21_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct7_22_P2_3D( const GeoVector& )
{
    return -8;
}
Real der2fct7_23_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct7_31_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct7_32_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct7_33_P2_3D( const GeoVector& )
{
    return 0;
}

Real der2fct8_11_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct8_12_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct8_13_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct8_21_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct8_22_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct8_23_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct8_31_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct8_32_P2_3D( const GeoVector& )
{
    return -4;
}
Real der2fct8_33_P2_3D( const GeoVector& )
{
    return -8;
}

Real der2fct9_11_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct9_12_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct9_13_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct9_21_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct9_22_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct9_23_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct9_31_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct9_32_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct9_33_P2_3D( const GeoVector& )
{
    return 0;
}

Real der2fct10_11_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct10_12_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct10_13_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct10_21_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct10_22_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct10_23_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct10_31_P2_3D( const GeoVector& )
{
    return 0;
}
Real der2fct10_32_P2_3D( const GeoVector& )
{
    return 4;
}
Real der2fct10_33_P2_3D( const GeoVector& )
{
    return 0;
}
//======================================================================
//
//                            P2tilde  (3D)
// NAVIER-STOKES P2 Basis Oriented to the mass lumping
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7 .11 \6
           /.       \!
         1 -----5----2
*/
Real fct1_P2tilde_3D( const GeoVector& v )
{
    return -( 1 - v[0] - v[1] - v[2] ) * ( 1 - 2 * ( 1 - v[0] - v[1] - v[2] ) ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct2_P2tilde_3D( const GeoVector& v )
{
    return -v[0] * ( 1 - 2 * v[0] ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct3_P2tilde_3D( const GeoVector& v )
{
    return -v[1] * ( 1 - 2 * v[1] ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct4_P2tilde_3D( const GeoVector& v )
{
    return -v[2] * ( 1 - 2 * v[2] ) + 32 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}

Real fct5_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[0] * ( 1 - v[0] - v[1] - v[2] ) - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct6_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[0] * v[1] - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct7_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[1] * ( 1 - v[0] - v[1] - v[2] ) - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct8_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[2] * ( 1 - v[0] - v[1] - v[2] ) - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct9_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[0] * v[2] - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}
Real fct10_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[1] * v[2] - 64 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}

Real fct11_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[0] * v[1] * v[2] * ( 1 - v[0] - v[1] - v[2] );
}


Real derfct1_1_P2tilde_3D( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2] + 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct1_2_P2tilde_3D( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2] + 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct1_3_P2tilde_3D( const GeoVector& v )
{
    return -3 + 4 * v[0] + 4 * v[1] + 4 * v[2] + 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct2_1_P2tilde_3D( const GeoVector& v )
{
    return -1 + 4 * v[0] + 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct2_2_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct2_3_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct3_1_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct3_2_P2tilde_3D( const GeoVector& v )
{
    return -1 + 4 * v[1] + 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct3_3_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct4_1_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct4_2_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct4_3_P2tilde_3D( const GeoVector& v )
{
    return -1 + 4 * v[2] + 32 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct5_1_P2tilde_3D( const GeoVector& v )
{
    return 4 - 8 * v[0] - 4 * v[1] - 4 * v[2] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct5_2_P2tilde_3D( const GeoVector& v )
{
    return -4 * v[0] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct5_3_P2tilde_3D( const GeoVector& v )
{
    return -4 * v[0] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct6_1_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[1] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct6_2_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[0] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct6_3_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct7_1_P2tilde_3D( const GeoVector& v )
{
    return -4 * v[1] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct7_2_P2tilde_3D( const GeoVector& v )
{
    return 4 - 4 * v[0] - 8 * v[1] - 4 * v[2] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct7_3_P2tilde_3D( const GeoVector& v )
{
    return -4 * v[1] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct8_1_P2tilde_3D( const GeoVector& v )
{
    return -4 * v[2] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct8_2_P2tilde_3D( const GeoVector& v )
{
    return -4 * v[2] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct8_3_P2tilde_3D( const GeoVector& v )
{
    return 4 - 4 * v[0] - 4 * v[1] - 8 * v[2] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct9_1_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[2] - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct9_2_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct9_3_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[0] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct10_1_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct10_2_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[2] - 64 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct10_3_P2tilde_3D( const GeoVector& v )
{
    return 4 * v[1] - 64 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real derfct11_1_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[1] * v[2] * ( 1 - 2 * v[0] - v[1] - v[2] );
}
Real derfct11_2_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[0] * v[2] * ( 1 - v[0] - 2 * v[1] - v[2] );
}
Real derfct11_3_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[0] * v[1] * ( 1 - v[0] - v[1] - 2 * v[2] );
}

Real der2fct1_11_P2tilde_3D( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct1_12_P2tilde_3D( const GeoVector& v )
{
    return 4 + 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct1_13_P2tilde_3D( const GeoVector& v )
{
    return 4 + 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct1_21_P2tilde_3D( const GeoVector& v )
{
    return 4 + 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct1_22_P2tilde_3D( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct1_23_P2tilde_3D( const GeoVector& v )
{
    return 4 + 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct1_31_P2tilde_3D( const GeoVector& v )
{
    return 4 + 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct1_32_P2tilde_3D( const GeoVector& v )
{
    return 4 + 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
    ;
}
Real der2fct1_33_P2tilde_3D( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}

Real der2fct2_11_P2tilde_3D( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct2_12_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct2_13_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct2_21_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct2_22_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}
Real der2fct2_23_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct2_31_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct2_32_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct2_33_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}

Real der2fct3_11_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}
Real der2fct3_12_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct3_13_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct3_21_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct3_22_P2tilde_3D( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}
Real der2fct3_23_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct3_31_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct3_32_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct3_33_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}

Real der2fct4_11_P2tilde_3D( const GeoVector& v )
{
    return -64 * v[0] * v[1] * v[2];
}
Real der2fct4_12_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct4_13_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct4_21_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct4_22_P2tilde_3D( const GeoVector& v )
{
    return - 64 * v[0] * v[1] * v[2];
}
Real der2fct4_23_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct4_31_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct4_32_P2tilde_3D( const GeoVector& v )
{
    return 32 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct4_33_P2tilde_3D( const GeoVector& v )
{
    return 4 - 64 * v[0] * v[1] * v[2];
}

Real der2fct5_11_P2tilde_3D( const GeoVector& v )
{
    return -8 - 128 * v[0] * v[1] * v[2];
}
Real der2fct5_12_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct5_13_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct5_21_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct5_22_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct5_23_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct5_31_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct5_32_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct5_33_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct6_11_P2tilde_3D( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct6_12_P2tilde_3D( const GeoVector& v )
{
    return 4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct6_13_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct6_21_P2tilde_3D( const GeoVector& v )
{
    return 4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct6_22_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct6_23_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct6_31_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct6_32_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct6_33_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct7_11_P2tilde_3D( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct7_12_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct7_13_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct7_21_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct7_22_P2tilde_3D( const GeoVector& v )
{
    return -8 - 128 * v[0] * v[1] * v[2];
}
Real der2fct7_23_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct7_31_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct7_32_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct7_33_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct8_11_P2tilde_3D( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct8_12_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct8_13_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct8_21_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct8_22_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct8_23_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct8_31_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct8_32_P2tilde_3D( const GeoVector& v )
{
    return -4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct8_33_P2tilde_3D( const GeoVector& v )
{
    return -8 - 128 * v[0] * v[1] * v[2];
}

Real der2fct9_11_P2tilde_3D( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct9_12_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct9_13_P2tilde_3D( const GeoVector& v )
{
    return 4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct9_21_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct9_22_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct9_23_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct9_31_P2tilde_3D( const GeoVector& v )
{
    return 4 + 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct9_32_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct9_33_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct10_11_P2tilde_3D( const GeoVector& v )
{
    return -128 * v[0] * v[1] * v[2];
}
Real der2fct10_12_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct10_13_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct10_21_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct10_22_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}
Real der2fct10_23_P2tilde_3D( const GeoVector& v )
{
    return 4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct10_31_P2tilde_3D( const GeoVector& v )
{
    return 64 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct10_32_P2tilde_3D( const GeoVector& v )
{
    return 4 + 64 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct10_33_P2tilde_3D( const GeoVector& v )
{
    return - 128 * v[0] * v[1] * v[2];
}

Real der2fct11_11_P2tilde_3D( const GeoVector& v )
{
    return -512 * v[0] * v[1] * v[2];
}
Real der2fct11_12_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct11_13_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct11_21_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[2] * ( 1 - 2 * v[0] - 2 * v[1] - v[2] );
}
Real der2fct11_22_P2tilde_3D( const GeoVector& v )
{
    return -512 * v[0] * v[1] * v[2];
}
Real der2fct11_23_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct11_31_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[1] * ( 1 - 2 * v[0] - v[1] - 2 * v[2] );
}
Real der2fct11_32_P2tilde_3D( const GeoVector& v )
{
    return 256 * v[0] * ( 1 - v[0] - 2 * v[1] - 2 * v[2] );
}
Real der2fct11_33_P2tilde_3D( const GeoVector& v )
{
    return -512 * v[0] * v[1] * v[2];
}

//======================================================================
//
//                            Q0  (3D)
//
//======================================================================
/*
        ________
       /.      /|
      / .     / |
     /_______/  |
     |  .  1 |  |
     |  .....|..|
     | .     | /
     |.      |/
     |_______|

*/
Real fct1_Q0_3D( const GeoVector& )
{
    return 1.;
}
Real derfct1_Q0_3D( const GeoVector& )
{
    return 0.;
}
// The second derivative is equal to the first : both are equal to 0.
Real der2fct1_Q0_3D( const GeoVector& )
{
    return 0.;
}

//======================================================================
//
//                            Q1  (3D)
//
//======================================================================
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
Real fct1_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] ) * ( 1. - v[2] );
}
Real fct2_Q1_3D( const GeoVector& v )
{
    return v[0] * ( 1. - v[1] ) * ( 1. - v[2] );
}
Real fct3_Q1_3D( const GeoVector& v )
{
    return v[0] * v[1] * ( 1. - v[2] );
}
Real fct4_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[1] * ( 1. - v[2] );
}
Real fct5_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] ) * v[2];
}
Real fct6_Q1_3D( const GeoVector& v )
{
    return v[0] * ( 1. - v[1] ) * v[2];
}
Real fct7_Q1_3D( const GeoVector& v )
{
    return v[0] * v[1] * v[2];
}
Real fct8_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[1] * v[2];
}

Real derfct1_1_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[1] ) * ( 1. - v[2] );
}
Real derfct1_2_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] ) * ( 1. - v[2] );
}
Real derfct1_3_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] ) * ( 1. - v[1] );
}
Real derfct2_1_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[1] ) * ( 1. - v[2] );
}
Real derfct2_2_Q1_3D( const GeoVector& v )
{
    return -v[0] * ( 1. - v[2] ) ;
}
Real derfct2_3_Q1_3D( const GeoVector& v )
{
    return -v[0] * ( 1. - v[1] );
}
Real derfct3_1_Q1_3D( const GeoVector& v )
{
    return v[1] * ( 1. - v[2] );
}
Real derfct3_2_Q1_3D( const GeoVector& v )
{
    return v[0] * ( 1. - v[2] );
}
Real derfct3_3_Q1_3D( const GeoVector& v )
{
    return -v[0] * v[1] ;
}
Real derfct4_1_Q1_3D( const GeoVector& v )
{
    return -v[1] * ( 1. - v[2] );
}
Real derfct4_2_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[2] );
}
Real derfct4_3_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] ) * v[1];
}
Real derfct5_1_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[1] ) * v[2];
}
Real derfct5_2_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] ) * v[2];
}
Real derfct5_3_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * ( 1. - v[1] );
}
Real derfct6_1_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[1] ) * v[2] ;
}
Real derfct6_2_Q1_3D( const GeoVector& v )
{
    return -v[0] * v[2];
}
Real derfct6_3_Q1_3D( const GeoVector& v )
{
    return v[0] * ( 1. - v[1] );
}
Real derfct7_1_Q1_3D( const GeoVector& v )
{
    return v[1] * v[2];
}
Real derfct7_2_Q1_3D( const GeoVector& v )
{
    return v[0] * v[2];
}
Real derfct7_3_Q1_3D( const GeoVector& v )
{
    return v[0] * v[1];
}
Real derfct8_1_Q1_3D( const GeoVector& v )
{
    return -v[1] * v[2];
}
Real derfct8_2_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[2];
}
Real derfct8_3_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[0] ) * v[1];
}

Real der2fct1_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct1_12_Q1_3D( const GeoVector& v )
{
    return 1. - v[2];
}
Real der2fct1_13_Q1_3D( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct1_21_Q1_3D( const GeoVector& v )
{
    return 1. - v[2];
}
Real der2fct1_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct1_23_Q1_3D( const GeoVector& v )
{
    return 1. - v[0];
}
Real der2fct1_31_Q1_3D( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct1_32_Q1_3D( const GeoVector& v )
{
    return 1. - v[0];
}
Real der2fct1_33_Q1_3D( const GeoVector& )
{
    return 0;
}

Real der2fct2_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_12_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[2] );
}
Real der2fct2_13_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[1] );
}
Real der2fct2_21_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[2] );
}
Real der2fct2_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct2_23_Q1_3D( const GeoVector& v )
{
    return v[0];
}
Real der2fct2_31_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[1] );
}
Real der2fct2_32_Q1_3D( const GeoVector& v )
{
    return v[0];
}
Real der2fct2_33_Q1_3D( const GeoVector& )
{
    return 0;
}

Real der2fct3_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_12_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[2] );
}
Real der2fct3_13_Q1_3D( const GeoVector& v )
{
    return -v[1];
}
Real der2fct3_21_Q1_3D( const GeoVector& v )
{
    return ( 1. - v[2] );
}
Real der2fct3_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct3_23_Q1_3D( const GeoVector& v )
{
    return -v[0];
}
Real der2fct3_31_Q1_3D( const GeoVector& v )
{
    return -v[1];
}
Real der2fct3_32_Q1_3D( const GeoVector& v )
{
    return -v[0];
}
Real der2fct3_33_Q1_3D( const GeoVector& )
{
    return 0;
}

Real der2fct4_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_12_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[2] );
}
Real der2fct4_13_Q1_3D( const GeoVector& v )
{
    return v[1];
}
Real der2fct4_21_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[2] );
}
Real der2fct4_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct4_23_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] );
}
Real der2fct4_31_Q1_3D( const GeoVector& v )
{
    return v[1];
}
Real der2fct4_32_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] );
}
Real der2fct4_33_Q1_3D( const GeoVector& )
{
    return 0;
}

Real der2fct5_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct5_12_Q1_3D( const GeoVector& v )
{
    return v[2];
}
Real der2fct5_13_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[1] );
}
Real der2fct5_21_Q1_3D( const GeoVector& v )
{
    return v[2];
}
Real der2fct5_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct5_23_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] );
}
Real der2fct5_31_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[1] );
}
Real der2fct5_32_Q1_3D( const GeoVector& v )
{
    return -( 1. - v[0] );
}
Real der2fct5_33_Q1_3D( const GeoVector& )
{
    return 0;
}

Real der2fct6_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_12_Q1_3D( const GeoVector& v )
{
    return -v[2];
}
Real der2fct6_13_Q1_3D( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct6_21_Q1_3D( const GeoVector& v )
{
    return -v[2];
}
Real der2fct6_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct6_23_Q1_3D( const GeoVector& v )
{
    return -v[0];
}
Real der2fct6_31_Q1_3D( const GeoVector& v )
{
    return 1. - v[1];
}
Real der2fct6_32_Q1_3D( const GeoVector& v )
{
    return -v[0];
}
Real der2fct6_33_Q1_3D( const GeoVector& )
{
    return 0;
}

Real der2fct7_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct7_12_Q1_3D( const GeoVector& v )
{
    return v[2];
}
Real der2fct7_13_Q1_3D( const GeoVector& v )
{
    return v[1];
}
Real der2fct7_21_Q1_3D( const GeoVector& v )
{
    return v[2];
}
Real der2fct7_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct7_23_Q1_3D( const GeoVector& v )
{
    return v[0];
}
Real der2fct7_31_Q1_3D( const GeoVector& v )
{
    return v[1];
}
Real der2fct7_32_Q1_3D( const GeoVector& v )
{
    return v[0];
}
Real der2fct7_33_Q1_3D( const GeoVector& )
{
    return 0;
}

Real der2fct8_11_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct8_12_Q1_3D( const GeoVector& v )
{
    return -v[2];
}
Real der2fct8_13_Q1_3D( const GeoVector& v )
{
    return -v[1];
}
Real der2fct8_21_Q1_3D( const GeoVector& v )
{
    return -v[2];
}
Real der2fct8_22_Q1_3D( const GeoVector& )
{
    return 0;
}
Real der2fct8_23_Q1_3D( const GeoVector& v )
{
    return 1. - v[0];
}
Real der2fct8_31_Q1_3D( const GeoVector& v )
{
    return -v[1];
}
Real der2fct8_32_Q1_3D( const GeoVector& v )
{
    return 1 -v[0];
}
Real der2fct8_33_Q1_3D( const GeoVector& )
{
    return 0;
}

//======================================================================
//
//                            RT0 Hexa  (3D)
//
//======================================================================
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

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8

*/

Real fct1_RT0_1_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct1_RT0_2_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct1_RT0_3_HEXA_3D( const GeoVector& v )
{
    return v[2] - 1.;
}

Real fct2_RT0_1_HEXA_3D( const GeoVector& v )
{
    return v[0] - 1.;
}
Real fct2_RT0_2_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct2_RT0_3_HEXA_3D( const GeoVector& )
{
    return 0.;
}

Real fct3_RT0_1_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct3_RT0_2_HEXA_3D( const GeoVector& v )
{
    return v[1] - 1.;
}
Real fct3_RT0_3_HEXA_3D( const GeoVector& )
{
    return 0.;
}

Real fct4_RT0_1_HEXA_3D( const GeoVector& v )
{
    return v[0];
}
Real fct4_RT0_2_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct4_RT0_3_HEXA_3D( const GeoVector& )
{
    return 0.;
}

Real fct5_RT0_1_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct5_RT0_2_HEXA_3D( const GeoVector& v )
{
    return v[1];
}
Real fct5_RT0_3_HEXA_3D( const GeoVector& )
{
    return 0.;
}

Real fct6_RT0_1_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct6_RT0_2_HEXA_3D( const GeoVector& )
{
    return 0.;
}
Real fct6_RT0_3_HEXA_3D( const GeoVector& v )
{
    return v[2];
}

Real fct1_DIV_RT0_HEXA_3D( const GeoVector& )
{
    return 1.;
}
Real fct2_DIV_RT0_HEXA_3D( const GeoVector& )
{
    return 1.;
}
Real fct3_DIV_RT0_HEXA_3D( const GeoVector& )
{
    return 1.;
}
Real fct4_DIV_RT0_HEXA_3D( const GeoVector& )
{
    return 1.;
}
Real fct5_DIV_RT0_HEXA_3D( const GeoVector& )
{
    return 1.;
}
Real fct6_DIV_RT0_HEXA_3D( const GeoVector& )
{
    return 1.;
}

//======================================================================
//
//                            RT0 Tetra  (3D)
//
//======================================================================

/*

                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2

SEE basisElSh.cc   for the ORIENTATION CONVENTIONS
   point 1: 0, 0, 0
   point 2: 1, 0, 0
   point 3: 0, 1, 0
   point 4: 0, 0, 1

   face 1: 2, 3, 4
   face 2: 1, 4, 3
   face 3: 1, 2, 4
   face 4: 1, 3, 2


*/

Real fct3_RT0_1_TETRA_3D( const GeoVector& v )
{
    return 2. * v[0];
}
Real fct3_RT0_2_TETRA_3D( const GeoVector& v )
{
    return 2. * v[1];
}
Real fct3_RT0_3_TETRA_3D( const GeoVector& v )
{
    return 2. * v[2];
}

Real fct4_RT0_1_TETRA_3D( const GeoVector& v )
{
    return 2. * v[0] - 2.;
}
Real fct4_RT0_2_TETRA_3D( const GeoVector& v )
{
    return 2. * v[1];
}
Real fct4_RT0_3_TETRA_3D( const GeoVector& v )
{
    return 2. * v[2];
}

Real fct2_RT0_1_TETRA_3D( const GeoVector& v )
{
    return 2. * v[0];
}
Real fct2_RT0_2_TETRA_3D( const GeoVector& v )
{
    return 2. * v[1] - 2.;
}
Real fct2_RT0_3_TETRA_3D( const GeoVector& v )
{
    return 2. * v[2];
}

Real fct1_RT0_1_TETRA_3D( const GeoVector& v )
{
    return 2. * v[0];
}
Real fct1_RT0_2_TETRA_3D( const GeoVector& v )
{
    return 2. * v[1];
}
Real fct1_RT0_3_TETRA_3D( const GeoVector& v )
{
    return 2. * v[2] - 2.;
}

Real fct1_DIV_RT0_TETRA_3D( const GeoVector& )
{
    return 6.;
}
Real fct2_DIV_RT0_TETRA_3D( const GeoVector& )
{
    return 6.;
}
Real fct3_DIV_RT0_TETRA_3D( const GeoVector& )
{
    return 6.;
}
Real fct4_DIV_RT0_TETRA_3D( const GeoVector& )
{
    return 6.;
}

// Transformation functions

std::vector<Real> lagrangianTransform(const std::vector<Real>& values)
{
    return values;
}

std::vector<Real> P1Bubble3DTransform(const std::vector<Real>& nodalValues)
{
    std::vector<Real> FEValues(nodalValues);
    FEValues[4]=256*nodalValues[4] - 64*(nodalValues[0]+nodalValues[1]+nodalValues[2]+nodalValues[3]);
    return FEValues;
}


//======================================================================
//
//                            P0  (0D)
//
//======================================================================
/*
                           1
*/

const RefFEScalar fePointP0( "Lagrange P0 on a point",
                             FE_P0_0D,
                             POINT,
                             1,                           // nb dof per vertex
                             0,                           // nb dof per edge
                             0,                           // nb dof per face
                             0,                           // nb dof per volume
                             1,                           // nb dof
                             1,                           // nb coor
                             fct_P0_0D,
                             derfct_P0_0D,
                             der2fct_P0_0D,
                             refcoor_P0_0D,
                             STANDARD_PATTERN,
                             ( RefFE* ) NULL,
                             &lagrangianTransform );

//======================================================================
//
//                            P1  (1D)
//
//======================================================================
/*
                           1-----2
*/

const RefFEScalar feSegP1( "Lagrange P1 on a segment", FE_P1_1D, LINE, 1, 0, 0, 0, 2, 1,
                           fct_P1_1D, derfct_P1_1D, der2fct_P1_1D, refcoor_P1_1D,
                           STANDARD_PATTERN, &fePointP0,&lagrangianTransform );

//======================================================================
//
//                            P2  (1D)
//
//======================================================================
/*
                           1--3--2
*/

const RefFEScalar feSegP2( "Lagrange P2 on a segment", FE_P2_1D, LINE, 1, 1, 0, 0, 3, 1,
                           fct_P2_1D, derfct_P2_1D, der2fct_P2_1D, refcoor_P2_1D,
                           STANDARD_PATTERN, &fePointP0,&lagrangianTransform );

//======================================================================
//
//                            P0  (2D)
//
//======================================================================
/*

                           |\
                           | \
                           | 1\
                            ---
*/

const RefFEScalar feTriaP0( "Lagrange P0 on a triangle", FE_P0_2D, TRIANGLE, 0, 0, 0, 1, 1, 2,
                            fct_P0_2D, derfct_P0_2D, der2fct_P0_2D, refcoor_P0_2D,
                            STANDARD_PATTERN, ( RefFE* ) NULL,&lagrangianTransform );

//======================================================================
//
//                            P1  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |  \
                           1---2
*/

const RefFEScalar feTriaP1( "Lagrange P1 on a triangle", FE_P1_2D, TRIANGLE, 1, 0, 0, 0, 3, 2,
                            fct_P1_2D, derfct_P1_2D, der2fct_P1_2D, refcoor_P1_2D,
                            STANDARD_PATTERN, &feSegP1,&lagrangianTransform );

//======================================================================
//
//                            P2  (2D)
//
//======================================================================
/*
                           3
                           |\
                           6 5
                           |  \
                           1-4-2
*/

const RefFEScalar feTriaP2( "Lagrange P2 on a triangle", FE_P2_2D, TRIANGLE, 1, 1, 0, 0, 6, 2,
                            fct_P2_2D, derfct_P2_2D, der2fct_P2_2D, refcoor_P2_2D,
                            STANDARD_PATTERN, &feSegP2,&lagrangianTransform );

//======================================================================
//
//                            Q0  (2D)
//
//======================================================================
/*
                            -------
                           |       |
                           |   1   |
                           |       |
                            -------
*/

const RefFEScalar feQuadQ0( "Lagrange Q0 on a quadrangle", FE_Q0_2D, QUAD, 0, 0, 1, 0, 1, 2,
                            fct_Q0_2D, derfct_Q0_2D, der2fct_Q0_2D, refcoor_Q0_2D,
                            STANDARD_PATTERN, ( RefFE* ) NULL,&lagrangianTransform );

//======================================================================
//
//                            Q1  (2D)
//
//======================================================================
/*
                           4-------3
                           |       |
                           |       |
                           |       |
                           1-------2
*/

const RefFEScalar feQuadQ1( "Lagrange Q1 on a quadrangle", FE_Q1_2D, QUAD, 1, 0, 0, 0, 4, 2,
                            fct_Q1_2D, derfct_Q1_2D, der2fct_Q1_2D, refcoor_Q1_2D,
                            STANDARD_PATTERN, &feSegP1,&lagrangianTransform );


//======================================================================
//
//                            Q2  (2D)
//
//======================================================================
/*
                           4---7---3
                           |       |
                           8   9   6
                           |       |
                           1---5---2
*/

const RefFEScalar feQuadQ2( "Lagrange Q2 on a quadrangle", FE_Q2_2D, QUAD, 1, 1, 1, 0, 9, 2,
                            fct_Q2_2D, derfct_Q2_2D, der2fct_Q2_2D, refcoor_Q2_2D,
                            STANDARD_PATTERN, &feSegP2,&lagrangianTransform );

//======================================================================
//
//                            P0  (3D)
//
//======================================================================
/*

               / .
              /  \.
             /  . \\
            / . 1  \\
           /.       \!
           ----------
*/
const RefFEScalar feTetraP0( "Lagrange P0 on a tetraedra", FE_P0_3D, TETRA, 0, 0, 0, 1, 1, 3,
                             fct_P0_3D, derfct_P0_3D, der2fct_P0_3D, refcoor_P0_3D,
                             STANDARD_PATTERN, &feTriaP0,&lagrangianTransform );

//======================================================================
//
//                            P1  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
const RefFEScalar feTetraP1( "Lagrange P1 on a tetraedra", FE_P1_3D, TETRA, 1, 0, 0, 0, 4, 3,
                             fct_P1_3D, derfct_P1_3D, der2fct_P1_3D, refcoor_P1_3D,
                             STANDARD_PATTERN, &feTriaP1,&lagrangianTransform );

//======================================================================
//
//                            P1bubble  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / . .5 \\
           /.       \!
         1 ----------2
*/
const RefFEScalar feTetraP1bubble( "Lagrange P1bubble on a tetraedra", FE_P1bubble_3D, TETRA, 1, 0, 0, 1, 5, 3,
                                   fct_P1bubble_3D, derfct_P1bubble_3D, der2fct_P1bubble_3D, refcoor_P1bubble_3D,
                                   STANDARD_PATTERN, &feTriaP1, &P1Bubble3DTransform );


//======================================================================
//
//                            P2  (3D)
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7    \6
           /.       \!
         1 -----5----2
*/
const RefFEScalar feTetraP2( "Lagrange P2 on a tetraedra", FE_P2_3D, TETRA, 1, 1, 0, 0, 10, 3,
                             fct_P2_3D, derfct_P2_3D, der2fct_P2_3D, refcoor_P2_3D,
                             STANDARD_PATTERN, &feTriaP2,&lagrangianTransform );
//======================================================================
//
//                            P2tilde  (3D)
// NAVIER-STOKES P2 Basis Oriented to the mass lumping
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7 .11 \6
           /.       \!
         1 -----5----2
*/
const RefFEScalar feTetraP2tilde( "Lagrange P2tilde on a tetraedra", FE_P2tilde_3D,
                                  TETRA, 1, 1, 0, 1, 11, 3, fct_P2tilde_3D,
                                  derfct_P2tilde_3D,
                                  der2fct_P2tilde_3D,
                                  refcoor_P2tilde_3D,
                                  STANDARD_PATTERN, &feTriaP2,&lagrangianTransform );

//======================================================================
//
//                            Q0  (3D)
//
//======================================================================
/*
        ________
       /.      /|
      / .     / |
     /_______/  |
     |  .  1 |  |
     |  .....|..|
     | .     | /
     |.      |/
     |_______|
*/
const RefFEScalar feHexaQ0( "Lagrange Q0 on a hexaedra", FE_Q0_3D, HEXA, 0, 0, 0, 1, 1, 3,
                            fct_Q0_3D, derfct_Q0_3D, der2fct_Q0_3D, refcoor_Q0_3D,
                            STANDARD_PATTERN, &feQuadQ0,&lagrangianTransform );

//======================================================================
//
//                            Q1  (3D)
//
//======================================================================
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
const RefFEScalar feHexaQ1( "Lagrange Q1 on a hexaedra", FE_Q1_3D, HEXA, 1, 0, 0, 0, 8, 3,
                            fct_Q1_3D, derfct_Q1_3D, der2fct_Q1_3D, refcoor_Q1_3D,
                            STANDARD_PATTERN, &feQuadQ1,&lagrangianTransform );

//======================================================================
//
//                            RT0 (3D)
//
//======================================================================
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

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8

*/
const RefFEHdiv feHexaRT0( "Lagrange RT0 on a hexaedra", FE_RT0_HEXA_3D, HEXA, 0, 0, 1, 0, 6, 3,
                           fct_RT0_HEXA_3D, fct_DIV_RT0_HEXA_3D, refcoor_RT0_HEXA_3D,
                           STANDARD_PATTERN, &feQuadQ0);

//======================================================================
//
//                            RT0 (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2


   face 1: 1, 3, 2
   face 2: 1, 2, 4
   face 3: 2, 3, 4
   face 4: 1, 4, 3
*/
const RefFEHdiv feTetraRT0( "Lagrange RT0 on a tetraedra", FE_RT0_TETRA_3D, TETRA, 0, 0, 1, 0, 4, 3,
                            fct_RT0_TETRA_3D, fct_DIV_RT0_TETRA_3D, refcoor_RT0_TETRA_3D,
                            STANDARD_PATTERN, &feTriaP0 );


//----------------------------------------------------------------------
//
//                       Mixed Hybrid FE
//
//----------------------------------------------------------------------

//======================================================================
//
//                           RT0 HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
//
//======================================================================
/*!

        8-------7
       /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8


*/

// N.B. : the hybrid classes and arrays depend on the quadrature rules,
//        geometrical mappings and other reference elements :
//        thus they must be defined AFTER the definitions of quadrule, geomap, refFE...

//! Total number of Boundary elements for the hybrid MFE for HEXA (= Number of faces, common for RT0,RT1...)
#define NB_BDFE_HYB_HEXA 6
static const StaticBdFE BdFE_RT0_HYB_HEXA_1( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                             refcoor_HYB_HEXA_FACE_1, 0 );
static const StaticBdFE BdFE_RT0_HYB_HEXA_2( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                             refcoor_HYB_HEXA_FACE_2, 1 );
static const StaticBdFE BdFE_RT0_HYB_HEXA_3( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                             refcoor_HYB_HEXA_FACE_3, 2 );
static const StaticBdFE BdFE_RT0_HYB_HEXA_4( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                             refcoor_HYB_HEXA_FACE_4, 3 );
static const StaticBdFE BdFE_RT0_HYB_HEXA_5( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                             refcoor_HYB_HEXA_FACE_5, 4 );
static const StaticBdFE BdFE_RT0_HYB_HEXA_6( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                             refcoor_HYB_HEXA_FACE_6, 5 );

static const StaticBdFE HybRT0HexaList[ NB_BDFE_HYB_HEXA ] =
{
    BdFE_RT0_HYB_HEXA_1, BdFE_RT0_HYB_HEXA_2,
    BdFE_RT0_HYB_HEXA_3, BdFE_RT0_HYB_HEXA_4,
    BdFE_RT0_HYB_HEXA_5, BdFE_RT0_HYB_HEXA_6
};

//const RefHybridFE feHexaRT0Hyb(NB_BDFE_HYB_HEXA,HybRT0HexaList,"Hybrid RT0 elements on a hexaedra",
//         FE_RT0_HYB_HEXA_3D, HEXA, 0,0,1,0,6,3,
//         refcoor_RT0HYB_HEXA,STANDARD_PATTERN);
static const StaticBdFE BdFE_RT0_HYB_HEXA_VdotN_1( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                                   refcoor_HYB_HEXA_FACE_1, 0, 1. );
static const StaticBdFE BdFE_RT0_HYB_HEXA_VdotN_2( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                                   refcoor_HYB_HEXA_FACE_2, 1, 1. );
static const StaticBdFE BdFE_RT0_HYB_HEXA_VdotN_3( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                                   refcoor_HYB_HEXA_FACE_3, 2, 1. );
static const StaticBdFE BdFE_RT0_HYB_HEXA_VdotN_4( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                                   refcoor_HYB_HEXA_FACE_4, 3, 1. );
static const StaticBdFE BdFE_RT0_HYB_HEXA_VdotN_5( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                                   refcoor_HYB_HEXA_FACE_5, 4, 1. );
static const StaticBdFE BdFE_RT0_HYB_HEXA_VdotN_6( feQuadQ0, geoBilinearQuad, quadRuleQuad4pt,
                                                   refcoor_HYB_HEXA_FACE_6, 5, 1. );

static const StaticBdFE HybRT0HexaVdotNList[ NB_BDFE_HYB_HEXA ] =
{
    BdFE_RT0_HYB_HEXA_VdotN_1, BdFE_RT0_HYB_HEXA_VdotN_2,
    BdFE_RT0_HYB_HEXA_VdotN_3, BdFE_RT0_HYB_HEXA_VdotN_4,
    BdFE_RT0_HYB_HEXA_VdotN_5, BdFE_RT0_HYB_HEXA_VdotN_6
};

const RefFEHybrid feHexaRT0Hyb( "Hybrid RT0 elements on a hexaedra", FE_RT0_HYB_HEXA_3D, HEXA,
                                0, 0, 1, 0, 6, 3, NB_BDFE_HYB_HEXA, HybRT0HexaList,
                                refcoor_RT0HYB_HEXA, STANDARD_PATTERN );

const RefFEHybrid feHexaRT0VdotNHyb( "Hybrid RT0 elements on a hexaedra", FE_RT0_HYB_HEXA_3D, HEXA,
                                     0, 0, 1, 0, 6, 3, NB_BDFE_HYB_HEXA, HybRT0HexaVdotNList,
                                     refcoor_RT0HYB_HEXA, STANDARD_PATTERN );


//======================================================================
//
//                           RT0 TETRA HYBRID (3D)
//                Element defined on FACES :  Q0 on each QUAD face.
//
//======================================================================
/*!

                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2

SEE basisElSh.cc   for the ORIENTATION CONVENTIONS
   point 1: 0, 0, 0
   point 2: 1, 0, 0
   point 3: 0, 1, 0
   point 4: 0, 0, 1

   face 1: 1, 3, 2
   face 2: 1, 2, 4
   face 3: 2, 3, 4
   face 4: 1, 4, 3


*/
// N.B. : the hybrid classes and arrays depend on the quadrature rules,
//        geometric mappings and other reference elements :
//        thus they must be defined AFTER the definitions of quadrule, geomap, refFE...

//! Total number of Boundary elements for the hybrid MFE for TETRA (= Number of faces. common for RT0,RT1...)
#define NB_BDFE_RT0_HYB_TETRA 4
static const StaticBdFE BdFE_RT0_HYB_TETRA_1( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                              refcoor_HYB_TETRA_FACE_1, 0 );
static const StaticBdFE BdFE_RT0_HYB_TETRA_2( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                              refcoor_HYB_TETRA_FACE_2, 1 );
static const StaticBdFE BdFE_RT0_HYB_TETRA_3( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                              refcoor_HYB_TETRA_FACE_3, 2 );
static const StaticBdFE BdFE_RT0_HYB_TETRA_4( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                              refcoor_HYB_TETRA_FACE_4, 3 );

static const StaticBdFE HybRT0TetraList[ NB_BDFE_RT0_HYB_TETRA ] =
{
    BdFE_RT0_HYB_TETRA_1, BdFE_RT0_HYB_TETRA_2,
    BdFE_RT0_HYB_TETRA_3, BdFE_RT0_HYB_TETRA_4
};

/*const RefHybridFE feTetraRT0Hyb (NB_BDFE_RT0_HYB_TETRA,HybRT0TetraList,"Hybrid RT0 elements on a tetrahedron",
    FE_RT0_HYB_TETRA_3D, TETRA, 0,0,1,0,4,3,
    refcoor_RT0HYB_TETRA,STANDARD_PATTERN);*/


static const StaticBdFE BdFE_RT0_HYB_TETRA_VdotN_1( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                                    refcoor_HYB_TETRA_FACE_1, 0, 2. );
static const StaticBdFE BdFE_RT0_HYB_TETRA_VdotN_2( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                                    refcoor_HYB_TETRA_FACE_2, 1, 2. );
static const StaticBdFE BdFE_RT0_HYB_TETRA_VdotN_3( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                                    refcoor_HYB_TETRA_FACE_3, 2, 2. / sqrt( 3. ) );
static const StaticBdFE BdFE_RT0_HYB_TETRA_VdotN_4( feTriaP0, geoLinearTria, quadRuleTria4pt,
                                                    refcoor_HYB_TETRA_FACE_4, 3, 2. );

static const StaticBdFE HybRT0TetraVdotNList[ NB_BDFE_RT0_HYB_TETRA ] =
{
    BdFE_RT0_HYB_TETRA_VdotN_1, BdFE_RT0_HYB_TETRA_VdotN_2,
    BdFE_RT0_HYB_TETRA_VdotN_3, BdFE_RT0_HYB_TETRA_VdotN_4
};

const RefFEHybrid feTetraRT0Hyb ( "Hybrid RT0 elements on a tetrahedron", FE_RT0_HYB_TETRA_3D, TETRA,
                                  0, 0, 1, 0, 4, 3, NB_BDFE_RT0_HYB_TETRA, HybRT0TetraList,
                                  refcoor_RT0HYB_TETRA, STANDARD_PATTERN );

const RefFEHybrid feTetraRT0VdotNHyb ( "Hybrid RT0 elements on a tetrahedron", FE_RT0_HYB_TETRA_3D, TETRA,
                                       0, 0, 1, 0, 4, 3, NB_BDFE_RT0_HYB_TETRA, HybRT0TetraVdotNList,
                                       refcoor_RT0HYB_TETRA, STANDARD_PATTERN );


}


