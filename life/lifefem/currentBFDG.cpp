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
#include <life/lifefem/currentBFDG.hpp>

namespace LifeV
{
CurrentBFDG::CurrentBFDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr):
  CurrentBdDG(_refFEDG, _refFE, _geoMap, _geoMapAd, _qr),

  re(_refFEDG.nbDof, _refFEDG.nbDof, _geoMapAd.nbCoor)
{
  CONSTRUCTOR("CurrentBFDG");
}

CurrentBFDG::~CurrentBFDG()
{
  DESTRUCTOR("CurrentBFDG");
}

void CurrentBFDG::_comp_re()
{
  KN<Real> rhs(nbNodeAd);

  for(int j = 0; j < nbNodeAd; j++){
    for(int icoor = 0; icoor < nbCoorAd; icoor++){
      rhs = 0.;

      // Compute right hand side vector
      for(int i = 0; i < nbNodeAd; i++){
    for(int ig = 0; ig < nbQuadPt; ig++){
      rhs(i) += - (phiAd(j, ig) * normal(icoor, ig)) * phiAd(i, ig) * weightMeas(ig);
    } // for ig
      } // for i

      // Solve the system
      rhs = invMassAd * rhs;

      // Store the solution
      for(int index = 0; index < nbNodeAd; index++){
    re(index, j, icoor) = rhs(index);
      }

    } // for icoor
  } // for j
}
}
