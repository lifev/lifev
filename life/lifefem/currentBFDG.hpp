/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#ifndef _CURRENTBFDG_H
#define _CURRENTBFDG_H

#include "currentBdDG.hpp"
#include "bcHandler.hpp"

namespace LifeV
{

/*!
  \class CurrentBFDG
  \brief The class for a discontinuous finite element
  \author D. A. Di Pietro
  \date 12/2003
*/

class CurrentBFDG:public CurrentBdDG{
#ifdef TEST_PRE
 protected:
  bool _hasre;
#endif

 private:
  void _comp_re();
 public:
  BCType bcType;
  EntityFlag marker;

  KNMK<Real> re; //!< Contains the i-th coefficient of the piecewise polynomial approximation of the k-th component of r_e(jump(fi_j)) on the adjacent element.

 public:
  CurrentBFDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr);
  ~CurrentBFDG();

#ifdef TEST_PRE
  inline bool hasre() const {return _hasre;}
#endif

  template<class GEOELE>
    void updateBCType(const GEOELE& geoele, const BCHandler& BCh) {
    marker = EntityFlag(geoele.marker());
    bcType = BCh.boundaryType(marker);
  }
 //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent,
    normal, quadrature points, phiDer, mass and invMass
    on the current boundary element
  */
  template<class GEOELE, class GEOELEAD, class CURRFEAD>
    void updateMeasNormalQuadPtFirstDerivMassreAd(const GEOELE& geoele, const GEOELEAD& geoelead, CURRFEAD& currfead)
    {
#ifdef TEST_PRE
      _hasre = true;
#endif
      updateMeasNormalQuadPtFirstDerivMassAd(geoele, geoelead, currfead);
      _comp_re();
    }
};
}
#endif
