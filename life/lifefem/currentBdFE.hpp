/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003 LifeV Team
  
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
#ifndef _CURRENTBDFE_H
#define _CURRENTBDFE_H

#include "lifeV.hpp"
#include "geoMap.hpp"
#include "refFE.hpp"
#include "geoMap.hpp"
#include "staticBdFE.hpp"
/*!
  \file currentBdFE.h
  \brief Structure for a current finite element on the boundary
*/

/*!
  \class CurrentBdFE
  \brief The class for a boundary finite element
  \author J.-F. Gerbeau
  \date 09/2002

  This class is used for a current boundary elements, i.e. the boundary of
  a CurrentFE. As for the CurrentFE (and contrarily to StaticBdFE) on must
  update this element with a geometrical element before using it. If only
  the measure is needed - to perform for example surface integrations -
  call updateMeas(...), which by the way compute also the tangent.  If, the
  normal at the integration point is needed, call updateMeasNormal(...).
  See the description of the base class StaticBdFE for further details.
*/

class CurrentBdFE:
  public StaticBdFE
{
public:
  CurrentBdFE(const RefFE& _refFE,const GeoMap& _geoMap);
  CurrentBdFE(const RefFE& _refFE,const GeoMap& _geoMap,const QuadRule& _qr);
  ~CurrentBdFE();
  /*!
    Compute only the coordinates of the nodes on the current boundary element
  */
  template<class GEOELE>
  void update(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasMeas       = false;
    _hasTangent    = false;
    _hasNormal     = false;
    _hasQuadPtCoor = false;
    _hasFirstDeriv = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points 
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
  };
  /*!
    Compute the arrays meas, weightMeas, tangent
    on the current boundary element
  */
  template<class GEOELE>
  void updateMeas(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasMeas       = true;
    _hasTangent    = true;
    _hasNormal     = false;
    _hasQuadPtCoor = false;
    _hasFirstDeriv = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points 
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the measure
    _comp_meas();
  };
  /*!
    Compute the arrays meas, weightMeas, tangent
    and quadrature points on the current boundary element
  */
  template<class GEOELE>
  void updateMeasQuadPt(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasMeas       = true;
    _hasTangent    = true;
    _hasNormal     = false;
    _hasQuadPtCoor = true;
    _hasFirstDeriv = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points 
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the measure
    _comp_meas();
    // compute the coordinates of the quad points
    _comp_quad_point_coor();
  };
  /*!
    Compute the arrays meas, weightMeas, tangent
    and normal on the current boundary element
  */
  template<class GEOELE>
  void updateMeasNormal(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasMeas = true;
    _hasTangent = true;
    _hasNormal = true;
    _hasQuadPtCoor = false;
    _hasFirstDeriv = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points 
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the measure and the normal
    _comp_meas_normal();    
  }
  /*!
    Compute the arrays meas, weightMeas, tangent,
    normal and quadrature points on the current boundary element
  */
  template<class GEOELE>
  void updateMeasNormalQuadPt(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasMeas = true;
    _hasTangent = true;
    _hasNormal = true;
    _hasQuadPtCoor = true;
    _hasFirstDeriv = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points 
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the measure and the normal
    _comp_meas_normal();
    // compute the coordinates of the quad points
    _comp_quad_point_coor();
  }
};
#endif
