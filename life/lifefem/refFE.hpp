/*-*- mode: c++ -*-
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
#ifndef _REFFE_H
#define _REFFE_H

#include "lifeV.hpp"
#include "refEle.hpp"
#include "localDofPattern.hpp"

/*!
  \file refFE.h
  \brief Structure for a reference Lagrangian finite element
*/
namespace LifeV
{
/*!
  \class RefFE
  \brief The class for a reference Lagrangian finite element
  \author J.-F. Gerbeau
  \date 04/2002
 
  \par How to add a new finite element ?
 
  (1) in refEle.h : you declare the functions you need (fct1_Pipo_2D,
  derfct1_1_Pipo_2D, etc...), the static arrays containing these functions
  and the coordinates of the nodes on the reference element.
 
  \par
  (2) in defQuadRuleFE.cc : you define these functions (fct1_Pipo_2D, etc...)
 
  \par
  (3) in refFE.h, you declare your finite element:
  \code
  extern const RefFE fePipo;
  \endcode
 
  \par
  (4) in defQuadRuleFE.cc: you define your new element with a command like:
  \code
  const RefFE feTriaPipo("Pipo element on a triangle",TRIANGLE,1,0,0,0,3,2,
                         fct_Pipo_2D,derfct_Pipo_2D,der2fct_Pipo_2D,refcoor_Pipo_2D,allQuadRuleTria,STANDARD_PATTERN,&feSegP1);
  \endcode
  See documentation of RefFE::RefFE(...) for a precise description of all arguments
 
*/ 
/* Unique FE identifier*/

#define FE_P1_1D 1
#define FE_P2_1D 2

#define FE_P0_2D 10
#define FE_P1_2D 11
#define FE_P2_2D 12

#define FE_Q0_2D 13
#define FE_Q1_2D 14
#define FE_Q2_2D 15

#define FE_P0_3D 20
#define FE_P1_3D 21
#define FE_P1bubble_3D 22
#define FE_P2_3D 23
#define FE_P2tilde_3D 24

#define FE_Q0_3D 25
#define FE_Q1_3D 26
#define FE_Q2_3D 27

#define FE_RT0_HEXA_3D 31         //!< Vectorial space for Mixed FE
#define FE_RT0_TETRA_3D 32
#define FE_RT0_HYB_HEXA_3D 41     //!< for hybrid Mixed FE.
#define FE_RT1_HYB_HEXA_3D 42
#define FE_RT0_HYB_TETRA_3D 44


class RefFE:
            public RefEle,
            public LocalDofPattern
{
    const RefFE* _boundaryFE;
public:
    //! Type of finite element (FE_P1_2D, ..., see the #define at the beginning of refFE.h
    const int type;
    //! Constructor of a reference Lagrangian finite element.
    /*!
      Constructor of a reference finite element. The arguments are:
      \param _name  the name of the f.e.
      \param _type  the type of the f.e. (FE_P1_2D,... see the #define at the
      begining of refFE.h)
      \param _shape  the geometry belongs to enum ReferenceShapes {NONE, POINT,
      LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see basisElSh.h)
      \param _nbDofPerVertex  the number of degrees of freedom per vertex
      \param _nbDofPerEdge  the number of degrees of freedom per edge
      \param _nbDofPerFace  the number of degrees of freedom per face
      \param _nbDofPerVolume  the number of degrees of freedom per volume
      \param _nbDof  the total number of d.o.f ( = _nbDofPerVertex * nb vertex +
      _nbDofPerEdge * nb edges + etc...)
      \param _nbCoor  number of local coordinates
      \param phi  the static array containing the basis functions (defined in
      refEle.h)
      \param dPhi  the static array containing the derivatives of the basis
      functions (defined in refEle.h)
      \param d2Phi  the static array containing the second derivatives of the
      basis functions (defined in refEle.h)
      \param refCoor  the static array containing the coordinates of the nodes on
      the reference element (defined in refEle.h)
      \param sqr  a set of quadrature rule (defined in quadRule.cc)
      \param _patternType  in most of cases STANDARD_PATTERN, except for elements
      like P1isoP2 (to define a new pattern, add a new #define in refFE.h and
      code it in refFE.cc following the example of P1ISOP2_TRIA_PATTERN)
      \param bdRefFE  a pointer on the associated reference finite element on the boundary
    */
    RefFE( std::string _name, int _type, ReferenceShapes _shape,
           int _nbDofPerVertex, int _nbDofPerEdge, int _nbDofPerFace, int _nbDofPerVolume,
           int _nbDof, int _nbCoor, const Fct* phi, const Fct* dPhi, const Fct* d2Phi,
           const Real* _refCoor, const SetOfQuadRule& sqr, PatternType _patternType,
           const RefFE* bdRefFE );
    ~RefFE();
    friend std::ostream& operator << ( std::ostream& f, const RefFE& fe );
    //! return the natural reference finite element for the boundary
    inline const RefFE& boundaryFE() const
    {
        ASSERT_PRE( _boundaryFE , "No boundary FE defined" );
        return *_boundaryFE;
    }
};

//--------------------------------------------------

extern const RefFE feSegP1;
extern const RefFE feSegP2;

extern const RefFE feTriaP0;
extern const RefFE feTriaP1;
extern const RefFE feTriaP2;

extern const RefFE feQuadQ0;
extern const RefFE feQuadQ1;
extern const RefFE feQuadQ2;

extern const RefFE feTetraP0;
extern const RefFE feTetraP1;
extern const RefFE feTetraP1bubble;
extern const RefFE feTetraP2;
extern const RefFE feTetraP2tilde;

extern const RefFE feHexaQ0;
extern const RefFE feHexaQ1;
}
#endif
