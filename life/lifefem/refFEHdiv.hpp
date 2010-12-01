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
    @file
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10 May 2010

    A more detailed description of the file (if necessary)
 */

#ifndef REFFEHDIV_H
#define REFFEHDIV_H 1

#include<life/lifefem/refFE.hpp>

namespace LifeV
{

//! RefFEHdiv - Short description of the class

class RefFEHdiv
        : public RefFE
{
public:

    typedef RefFE::Fct Fct;

    //! Constructor of a reference Lagrangian finite element.
    /*!
      Constructor of a reference finite element. The arguments are:
      @param name  the name of the f.e.
      @param type  the type of the f.e. (FE_P1_2D,... see the #define at the
      begining of refFE.h)
      @param shape  the geometry belongs to enum ReferenceShapes {NONE, POINT,
      LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see basisElSh.h)
      @param nbDofPerVertex  the number of degrees of freedom per vertex
      @param nbDofPerEdge  the number of degrees of freedom per edge
      @param nbDofPerFace  the number of degrees of freedom per face
      @param nbDofPerVolume  the number of degrees of freedom per volume
      @param nbDof  the total number of d.o.f ( = nbDofPerVertex * nb vertex +
      nbDofPerEdge * nb edges + etc...)
      @param nbCoor  number of local coordinates
      @param phi  the static array containing the basis functions (defined in
      refEle.h)
      @param dPhi  the static array containing the derivatives of the basis
      functions (defined in refEle.h)
      @param d2Phi  the static array containing the second derivatives of the
      basis functions (defined in refEle.h)
      @param refCoor  the static array containing the coordinates of the nodes on
      the reference element (defined in refEle.h)
      @param patternType  in most of cases STANDARD_PATTERN, except for elements
      like P1isoP2 (to define a new pattern, add a new #define in refFE.h and
      code it in refFE.cc following the example of P1ISOP2_TRIA_PATTERN)
      @param bdRefFE  a pointer on the associated reference finite element on the boundary
    */
    RefFEHdiv( std::string          name,
               FE_TYPE              type,
               ReferenceShapes      shape,
               int                  nbDofPerVertex,
               int                  nbDofPerEdge,
               int                  nbDofPerFace,
               int                  nbDofPerVolume,
               int                  nbDof,
               int                  nbCoor,
               const Fct*           phi,
               const Fct*           divPhi,
               const Real*          refCoor,
               DofPatternType       patternType,
               const RefFE*         bdRefFE );

    ~RefFEHdiv()
    {
        DESTRUCTOR( "RefFEHdiv" );
    };

};

// Vector finite elements
extern const RefFEHdiv feHexaRT0;
extern const RefFEHdiv feTetraRT0;

} // Namespace LifeV

#endif /* REFFEHDIV_H */
