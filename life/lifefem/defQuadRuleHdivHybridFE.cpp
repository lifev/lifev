//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file defQuadRuleHybridFE.cpp
    @brief Hybrid FE rules.

    @author Simone Deparis <simone.deparis@epfl.ch>
    @date 09 Mar 2010

    The contents of this file were in  defQuadRuleFE.cpp
 */

#include <life/lifefem/refHdivFE.hpp>
#include <life/lifefem/refHybridFE.hpp>

namespace LifeV {

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
const RefHdivFE feHexaRT0( "RT0 on a hexaedra", FE_RT0_HEXA_3D, HEXA, 0, 0, 1, 0, 6, 3,
                           fct_RT0_3D, fct_DIV_RT0_3D, refcoor_RT0_3D,
                           allQuadRuleHexa, STANDARD_PATTERN );

//
// FARE I COMMENTI ESPLICATIVI
//
//
// Controllare argomenti
const RefHdivFE feTetraRT0( "RT0 on a tetra", FE_RT0_TETRA_3D, TETRA, 0, 0, 1, 0, 4, 3,
                            fct_RT0_3D_TETRA, fct_DIV_RT0_3D_TETRA, refcoor_RT0_3D,
                            allQuadRuleTetra, STANDARD_PATTERN );


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

const RefHybridFE feHexaRT0Hyb( NB_BDFE_HYB_HEXA, HybRT0HexaList, "Hybrid RT0 elements on a hexaedra",
                                FE_RT0_HYB_HEXA_3D, HEXA, 0, 0, 1, 0, 6, 3,
                                refcoor_RT0HYB_HEXA, STANDARD_PATTERN );

const RefHybridFE feHexaRT0VdotNHyb( NB_BDFE_HYB_HEXA, HybRT0HexaVdotNList, "Hybrid RT0 elements on a hexaedra",
                                     FE_RT0_HYB_HEXA_3D, HEXA, 0, 0, 1, 0, 6, 3,
                                     refcoor_RT0HYB_HEXA, STANDARD_PATTERN );




//======================================================================
//
//                           RT1 HYBRID (3D)
//                Element defined on FACES :  Q1 on each QUAD face.
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

SEE basisElSh.cc   for the ORIENTATION CONVENTIONS
   point 1: 0, 0, 0
   point 2: 1, 0, 0
   point 3: 1, 1, 0
   point 4: 1, 0, 0
   point 5: 0, 0, 1
   point 6: 1, 0, 1
   point 7: 1, 1, 1
   point 8: 1, 0, 1

   face 1: 1,4,3,2
   face 2: 1,5,8,4
   face 3: 1,2,6,5
   face 4: 2,3,7,6
   face 5: 3,4,8,7
   face 6: 5,6,7,8


*/
static const StaticBdFE BdFE_RT1_HYB_1( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
                                        refcoor_HYB_HEXA_FACE_1, 0 );
static const StaticBdFE BdFE_RT1_HYB_2( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
                                        refcoor_HYB_HEXA_FACE_2, 1 );
static const StaticBdFE BdFE_RT1_HYB_3( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
                                        refcoor_HYB_HEXA_FACE_3, 2 );
static const StaticBdFE BdFE_RT1_HYB_4( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
                                        refcoor_HYB_HEXA_FACE_4, 3 );
static const StaticBdFE BdFE_RT1_HYB_5( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
                                        refcoor_HYB_HEXA_FACE_5, 4 );
static const StaticBdFE BdFE_RT1_HYB_6( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
                                        refcoor_HYB_HEXA_FACE_6, 5 );

// NB_BDFE_HYB_HEXA previously defined (in RT0)
static const StaticBdFE HybRT1HexaList[ NB_BDFE_HYB_HEXA ] =
    {
        BdFE_RT1_HYB_1, BdFE_RT1_HYB_2,
        BdFE_RT1_HYB_3, BdFE_RT1_HYB_4,
        BdFE_RT1_HYB_5, BdFE_RT1_HYB_6
    };

/*const RefHybridFE feHexaRTVdotN1Hyb(NB_BDFE_HYB_HEXA,HybRT1HexaList,"Hybrid RT1 elements on a hexaedra",
          FE_RT1_HYB_HEXA_3D, HEXA, 0,0,4,0,24,3,
          refcoor_RT1HYB_HEXA,STANDARD_PATTERN);*/


static const StaticBdFE BdFE_RT1_HYB_VdotN_1( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
        refcoor_HYB_HEXA_FACE_1, 0, 1. );
static const StaticBdFE BdFE_RT1_HYB_VdotN_2( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
        refcoor_HYB_HEXA_FACE_2, 1, 1. );
static const StaticBdFE BdFE_RT1_HYB_VdotN_3( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
        refcoor_HYB_HEXA_FACE_3, 2, 1. );
static const StaticBdFE BdFE_RT1_HYB_VdotN_4( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
        refcoor_HYB_HEXA_FACE_4, 3, 1. );
static const StaticBdFE BdFE_RT1_HYB_VdotN_5( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
        refcoor_HYB_HEXA_FACE_5, 4, 1. );
static const StaticBdFE BdFE_RT1_HYB_VdotN_6( feQuadQ1, geoBilinearQuad, quadRuleQuad4pt,
        refcoor_HYB_HEXA_FACE_6, 5, 1. );

// NB_BDFE_HYB_HEXA previously defined (in RT0)
static const StaticBdFE HybRT1HexaVdotNList[ NB_BDFE_HYB_HEXA ] =
    {
        BdFE_RT1_HYB_VdotN_1, BdFE_RT1_HYB_VdotN_2,
        BdFE_RT1_HYB_VdotN_3, BdFE_RT1_HYB_VdotN_4,
        BdFE_RT1_HYB_VdotN_5, BdFE_RT1_HYB_VdotN_6
    };

const RefHybridFE feHexaRT1Hyb( NB_BDFE_HYB_HEXA, HybRT1HexaList, "Hybrid RT1 elements on a hexaedra",
                                FE_RT1_HYB_HEXA_3D, HEXA, 0, 0, 4, 0, 24, 3,
                                refcoor_RT1HYB_HEXA, STANDARD_PATTERN );
const RefHybridFE feHexaRT1VdotNHyb( NB_BDFE_HYB_HEXA, HybRT1HexaVdotNList, "Hybrid RT1 elements on a hexaedra",
                                     FE_RT1_HYB_HEXA_3D, HEXA, 0, 0, 4, 0, 24, 3,
                                     refcoor_RT1HYB_HEXA, STANDARD_PATTERN );

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
//        geometrical mappings and other reference elements :
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

/*const RefHybridFE feTetraRT0Hyb (NB_BDFE_RT0_HYB_TETRA,HybRT0TetraList,"Hybrid RT0 elements on a tetraedra",
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

const RefHybridFE feTetraRT0Hyb ( NB_BDFE_RT0_HYB_TETRA, HybRT0TetraList, "Hybrid RT0 elements on a tetraedra",
                                  FE_RT0_HYB_TETRA_3D, TETRA, 0, 0, 1, 0, 4, 3,
                                  refcoor_RT0HYB_TETRA, STANDARD_PATTERN );

const RefHybridFE feTetraRT0VdotNHyb ( NB_BDFE_RT0_HYB_TETRA, HybRT0TetraVdotNList, "Hybrid RT0 elements on a tetraedra",
                                       FE_RT0_HYB_TETRA_3D, TETRA, 0, 0, 1, 0, 4, 3,
                                       refcoor_RT0HYB_TETRA, STANDARD_PATTERN );


} // Namespace LifeV
