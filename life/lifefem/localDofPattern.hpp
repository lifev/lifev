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
#ifndef _LOCAL_DOF_PATTERN_HH
#define _LOCAL_DOF_PATTERN_HH

#include "lifeV.hpp"
#include <utility>

namespace LifeV
{
/*! Local pattern type
  This enum allows to distinguish the normal standard local pattern, which is a full pattern invoving all
  defgrees of freedom to special patterns. It is stared in LocalDofPattern for later use by Dof
*/
enum PatternType {STANDARD_PATTERN=1,P1ISOP2_TRIA_PATTERN=2};


class LocalDofPattern
{
protected:
  int* _patternFirst;//!< row index of the non zero terms of the element matrix
  int* _patternSecond;//!< column index of the non zero terms of the element matrix
  int _nbPattern; //!< Number of non-zero terms in the element matrix
  int _nbDiag;//!< Number of diagonal terms in the element matrix
  int _nbUpper;//!< Numer of upper terms in the element matrix
public:
  LocalDofPattern(int _nbLocalDof,int _nbDofPerVertex,int _nbDofPerEdge,int _nbDofPerFace,int _nbDofPerVolume,
		  PatternType _patternType);
  const int nbLocalDof;  //!< Total number of degrees of freedom (equal to refEle::nbDof)
  const int nbDofPerVertex;  //!< Number of degrees of freedom per vertex
  const int nbDofPerEdge;  //!< Number of degrees  of freedom per edge
  const int nbDofPerFace;  //!< Number of degrees  of freedom per face
  const int nbDofPerVolume;  //!< Number of degrees  of freedom per volume
  const PatternType patternType; //!< Type of pattern (STANDARD_PATTERN,P1ISOP2_TRIA_PATTERN,...)
  inline int nbPattern() const {return _nbPattern; }//!< Number of non-zero terms in the element matrix
  inline int nbDiag()  const {return _nbDiag; } //!< Number of diagonal terms in the element matrix
  inline int nbUpper() const {return _nbUpper; } //!< Numer of upper terms in the element matrix

  //!  patternFirst(i): row index in the element matrix of the i-th term of the pattern
  inline int patternFirst(int i) const
  {
    ASSERT_BD(i<_nbPattern)
      return _patternFirst[i];
  }
  //! patternSecond(i): column index in the element matrix of the i-th term of the pattern
  inline int patternSecond(int i) const
  {
    ASSERT_BD(i<_nbPattern)
      return _patternSecond[i];
  }
};
}
#endif

