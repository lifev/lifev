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
    @brief Base structure for a reference finite element

    @author Jean-Frederic Gerbeau
    @date 00-04-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef _REFFE_H
#define _REFFE_H

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/ReferenceElement.hpp>
#include <lifev/core/fem/DOFLocalPattern.hpp>

namespace LifeV
{
/*!
  \class ReferenceFE
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
  extern const ReferenceFE fePipo;
  \endcode

  \par
  (4) in defQuadRuleFE.cc: you define your new element with a command like:
  \code
  const ReferenceFE feTriaPipo("Pipo element on a triangle",TRIANGLE,1,0,0,0,3,2,
  fct_Pipo_2D,derfct_Pipo_2D,der2fct_Pipo_2D,refcoor_Pipo_2D,allQuadRuleTria,STANDARD_PATTERN,&feSegP1);
  \endcode
  See documentation of ReferenceFE::ReferenceFE(...) for a precise description of all arguments

*/
/* Unique FE identifier*/
enum FE_TYPE
{
    FE_P0_0D = 1,
    FE_P0_1D,
    FE_P1_1D,
    FE_P2_1D,

    FE_P0_2D,
    FE_P1_2D,
    FE_P1bubble_2D,
    FE_P2_2D,

    FE_Q0_2D,
    FE_Q1_2D,
    FE_Q2_2D,

    FE_RT0_TRIA_2D,
    FE_RT0_HYB_TRIA_2D,

    FE_P0_3D,
    FE_P1_3D,
    FE_P1bubble_3D,
    FE_P2_3D,
    FE_P2tilde_3D,

    FE_Q0_3D,
    FE_Q1_3D,
    FE_Q2_3D,

    FE_RT0_HEXA_3D,         //!< Vectorial space for Mixed FE
    FE_RT0_TETRA_3D,
    FE_RT0_HYB_HEXA_3D,     //!< for hybrid Mixed FE.
    FE_RT1_HYB_HEXA_3D,
    FE_RT0_HYB_TETRA_3D
};


class ReferenceFE:
    public ReferenceElement,
    public DOFLocalPattern
{

public:

    //! @name Public Types
    //@{

    typedef ReferenceElement::function_Type function_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor of a reference Lagrangian finite element.
    /*!
      Constructor of a reference finite element. The arguments are:
      @param name  the name of the f.e.
      @param type  the type of the f.e. (FE_P1_2D,... see the #define at the
      begining of refFE.h)
      @param shape  the geometry belongs to enum ReferenceShapes {NONE, POINT,
      LINE, TRIANGLE, QUAD, HEXA, PRISM, TETRA}; (see ElementShapes.h)
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
      @param divPhi the static array containing the divergence of the basis function
      @param refCoor  the static array containing the coordinates of the nodes on
      the reference element (defined in refEle.h)
      @param patternType  in most of cases STANDARD_PATTERN, except for elements
      like P1isoP2 (to define a new pattern, add a new #define in refFE.h and
      code it in refFE.cc following the example of P1ISOP2_TRIA_PATTERN)
      @param bdRefFE  a pointer on the associated reference finite element on the boundary
    */
    ReferenceFE ( std::string          name,
                  FE_TYPE              type,
                  ReferenceShapes      shape,
                  Int                  nbDofPerVertex,
                  Int                  nbDofPerEdge,
                  Int                  nbDofPerFace,
                  Int                  nbDofPerVolume,
                  Int                  nbDof,
                  Int                  nbCoor,
                  Int                  FEDim,
                  const function_Type* phi,
                  const function_Type* dPhi,
                  const function_Type* d2Phi,
                  const function_Type* divPhi,
                  const Real*          refCoor,
                  DofPatternType       patternType,
                  const ReferenceFE*         bdRefFE );

    //! Destructor
    virtual ~ReferenceFE();

    //@}


    //! @name Get Methods
    //@{

    //! Check if the reference element has boundary elements
    bool hasBoundaryFE() const
    {
        return M_boundaryFE != NULL;
    }

    //! Getter for the boundary finite element
    /*
      The boundary of a finite element has to be understood in the sense of the trace. For example, the boundary finite element of the P0 finite element on a triangle is a P0 finite element on a segment, even if there is no degree of freedom located on the edges of the triangle for P0 finite element.
     */
    const ReferenceFE& boundaryFE() const
    {
        ASSERT ( M_boundaryFE , "No boundary FE defined" );
        return *M_boundaryFE;
    }

    //! Getter for the type of the finite element
    const FE_TYPE& type() const
    {
        return M_type;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    ReferenceFE();

    //! No copy constructor
    ReferenceFE (const ReferenceFE&);

    //@}


    //! Reference to the boundary finite element
    const ReferenceFE* M_boundaryFE;

    //! Type of finite element (FE_P1_2D, ..., see the #define at the beginning of refFE.h
    const FE_TYPE M_type;
};



}
#endif
