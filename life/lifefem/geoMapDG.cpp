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
#include <life/lifefem/geoMapDG.hpp>

namespace LifeV
{
  GeoMapDG::GeoMapDG(std::string _name,
         ReferenceShapes _shape,
         int _nbDof,
         int _nbCoor,
         const Fct* phi, const Fct* dPhi, const Fct* d2Phi, 
         const Real* refCoor, 
         const SetOfQuadRule& sqr,
         ReferenceShapes _shapeFaces, 
         int _nbFaces, int _nbGeoNodeFaces, 
         const Real* refCoorFaces, 
         const SetOfQuadRule& sqrFaces, const GeoMap& _geoMap, const GeoMap* bdMap):
  RefEleDG( _name, _shape, _nbDof, _nbCoor, phi, dPhi, d2Phi, refCoor, sqr, _shapeFaces, _nbFaces, _nbGeoNodeFaces, refCoorFaces, sqrFaces, _geoMap),
  _boundaryMap(bdMap)
{
  CONSTRUCTOR("GeoMap");
}

GeoMapDG::~GeoMapDG()
{
  DESTRUCTOR("GeoMap");
}

std::ostream& operator << (std::ostream& f,const GeoMapDG& geomap)
{
  f << "--------------------------------------------------------------------------------" << std::endl;
  f << "Geometrical mapping (discontinuous elements): " << geomap.name << std::endl;
  f << "--------------------------------------------------------------------------------" << std::endl;
  f << std::endl << "*** Shape : " << geomap.shape << std::endl;
  f << std::endl << "*** Local coordinate of the nodes :\n";
  for(int i=0;i<geomap.nbDof;i++){
    f << geomap.xi(i) << " " << geomap.eta(i) << " " << geomap.zeta(i) << std::endl;
  }
  for(int k=0;k<geomap._sqr->nbQuadRule;k++){
    const QuadRule& qr = geomap._sqr->quadRule(k);
    f << std::endl << "*** Quadrature rule : " << qr.name << std::endl;
    for(int ig=0;ig<qr.nbQuadPt;ig++){
      f << "    - Quadrature point : " << ig << std::endl;
      for(int i=0;i<geomap.nbDof;i++){
    f << "      Basif fct " << i << std::endl;
    f << "         Value = " << geomap.phi(i,ig,qr) << std::endl;
    f << "         Derivatives = " ;
    for(int icoor=0;icoor<geomap.nbCoor;icoor++) f << " " << geomap.dPhi(i,icoor,ig,qr);
    f << std::endl;
    f << "         Second derivatives = " ;
    for(int icoor=0;icoor<geomap.nbCoor;icoor++)
      for(int jcoor=0;jcoor<geomap.nbCoor;jcoor++)
      f << " " << geomap.d2Phi(i,icoor,jcoor,ig,qr);
    f << std::endl;
      }
    }
  }
  return f;
}
}
