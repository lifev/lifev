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
#ifndef _REFFEDG_H
#define _REFFEDG_H

#include <life/lifecore/life.hpp>
#include <life/lifefem/refEleDG.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/localDofPattern.hpp>

/*!
  \file refFEDG.h
  \brief Structure for a discontinuous finite element
*/

namespace LifeV
{
/*!
  \class RefFEDG
  \brief The class for a reference discontinuous finite element
  \author D. A. Di Pietro
  \date 12/2003

  \par Remark
  This class differs from RefFE in that you must declare two local
  patterns before creating a new finite elements,referring to integrals
  on volumes and on faces respectively.

*/

// Unique DG FE identifiers

#define FE_DG_P1_1D 1
#define FE_DG_P2_1D 2

#define FE_DG_P0_2D 10
#define FE_DG_P1_2D 11
#define FE_DG_P2_2D 12

#define FE_DG_Q0_2D 13
#define FE_DG_Q1_2D 14
#define FE_DG_Q2_2D 15

#define FE_DG_P1_3D 21
#define FE_DG_P1bubble_3D 22
#define FE_DG_P2_3D 23
#define FE_DG_P2tilde_3D 24

#define FE_DG_Q0_3D 25
#define FE_DG_Q1_3D 26
#define FE_DG_Q2_3D 27

class RefFEDG:
public RefEleDG
{
  const RefFE* _boundaryFE;
 public:
  const LocalDofPattern& elPattern;
  const LocalDofPattern& facePattern;
 public:
  const int type;
  /*!
    RefFEDG is the standard constructor for the reference discontinuous element.
    Notice that geoMap property is used to map faces onto the reference element.
  */
  RefFEDG(std::string _name, int _type,
      ReferenceShapes _shape,
      int _nbDofPerVertex, int _nbDofPerEdge, int _nbDofPerFace, int _nbDofPerVolume, int _nbDof,
      int _nbCoor,
      const Fct* phi, const Fct* dPhi, const Fct* d2Phi,
      const Real* refCoor,
      const SetOfQuadRule& sqr, const LocalDofPattern& _elPattern, const RefFE* boundaryFE,
      ReferenceShapes _shapeFaces,
      int _nbFaces, int _nbGeoNodeFaces,
      const Real* refCoorFaces,
      const SetOfQuadRule& sqrFaces, const LocalDofPattern& _facePattern, const GeoMap _geoMap);
  ~RefFEDG();

  friend std::ostream& operator << (std::ostream& f, const RefFEDG& fe);

  inline const RefFE& boundaryFE() const
    {
      ASSERT_PRE(_boundaryFE, "No boundary FE defined");
      return *_boundaryFE;
    };
};
  //---------------------------------------------------------------

  extern const RefFEDG feDGSegP1;
  extern const RefFEDG feDGSegP2;

  extern const RefFEDG feDGTriaP0;
  extern const RefFEDG feDGTriaP1;
  extern const RefFEDG feDGTriaP2;

  extern const RefFEDG feDGQuadQ0;
  extern const RefFEDG feDGQuadQ1;
  extern const RefFEDG feDGQuadQ2;

  extern const RefFEDG feDGTetraP1;
  extern const RefFEDG feDGTetraP1bubble;
  extern const RefFEDG feDGTetraP2;
  extern const RefFEDG feDGTetraP2tilde;
}
#endif
