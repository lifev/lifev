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
#include "refEleDG.hpp"
#include <set>

namespace LifeV
{

RefEleDG::RefEleDG(std::string _name, ReferenceShapes _shape,
		   int _nbDof, int _nbCoor, 
		   const Fct* phi, const Fct* dPhi, const Fct* d2Phi, 
		   const Real* refCoor, const SetOfQuadRule& sqr, 
		   ReferenceShapes _shapeFaces, 
		   int _nbFaces, int _nbGeoNodeFaces, 
		   const Real* refCoorFaces, const SetOfQuadRule& sqrFaces, const GeoMap& _geoMap):
  RefEle(_name, _shape, _nbDof, _nbCoor, phi, dPhi, d2Phi, refCoor, sqr),
  _sqrFaces(&sqrFaces),
  _refCoorFaces(refCoorFaces),
  shapeFaces(_shapeFaces),
  nbGeoNodeFaces(_nbGeoNodeFaces),
  nbFaces(_nbFaces),
  geoMap(_geoMap),
  _phiQuadFaces(_nbFaces, sqrFaces.totalNbQuadPoint() * _nbDof),
  _dPhiQuadFaces(_nbFaces, sqrFaces.totalNbQuadPoint() * _nbDof * _nbCoor),
  _d2PhiQuadFaces(_nbFaces, sqrFaces.totalNbQuadPoint() * _nbDof * _nbCoor * _nbCoor),
  _idxQuadFaces(sqrFaces.maxIdQuadRule() + 1),
  _idxDQuadFaces(sqrFaces.maxIdQuadRule() + 1),
  _idxD2QuadFaces(sqrFaces.maxIdQuadRule() + 1)
{
  CONSTRUCTOR("RefEleDG");

  Real xi = 0., eta = 0., zeta = 0.;
  int idx = 0, idxD = 0, idxD2 = 0;

  for(int iFace = 0; iFace < nbFaces; iFace++){
    idx = idxD = idxD2 = 0;
    for(int k = 0; k < sqrFaces.nbQuadRule; k++){
      _idxQuadFaces(sqrFaces.quadRule(k).id) = idx;
      _idxDQuadFaces(sqrFaces.quadRule(k).id) = idxD;
      _idxD2QuadFaces(sqrFaces.quadRule(k).id) = idxD2;

      for(int ig = 0; ig < sqrFaces.quadRule(k).nbQuadPt; ig++){
	const QuadPoint& pt = sqrFaces.quadRule(k).quadPoint(ig);
        FaceToElCoord(xi, eta, zeta, pt.x(), pt.y(), iFace);

	for(int i = 0; i < nbDof; i++){
	  _phiQuadFaces(iFace, idx) = this -> phi(i,xi,eta,zeta);
	  idx++;

	  for(int icoor = 0; icoor < nbCoor; icoor++){
	    _dPhiQuadFaces(iFace, idxD) = this -> dPhi(i, icoor, xi, eta, zeta);
	    idxD++;

	    for(int jcoor = 0; jcoor < nbCoor; jcoor++){
	      _d2PhiQuadFaces(iFace, idxD2) = this -> d2Phi(i, icoor, jcoor, xi, eta, zeta);
	      idxD2++;
	    } //for jcoor
	  } // for icoor
	} // for i
      } // for ig
    } // for k
  } // for iFace
 }

RefEleDG::~RefEleDG()
{
  DESTRUCTOR("RefEleDG");
}

void RefEleDG::check() const
{
  Real sumphi, sumdphi, sumd2phi;
  std::cout << "*** Check " << name << std::endl;

  for(int k = 0; k < _sqrFaces->nbQuadRule; k++){
    for(int iFace = 0; iFace < nbFaces; iFace++){
      sumphi = sumdphi = sumd2phi = 0.;
      const QuadRule& qr = _sqrFaces->quadRule(k);
      std::cout << std::endl << "    " << qr.name << std::endl;
      for(int ig = 0; ig < qr.nbQuadPt; ig++){
	for(int i = 0; i < nbDof; i++){
	  sumphi += phiFace(iFace, i, ig, qr) * qr.weight(ig);
	  for(int icoor = 0; icoor < nbCoor; icoor++){
	    sumdphi += dPhiFace(iFace, i, icoor, ig, qr) * qr.weight(ig);
	  } // for icoor
	  for(int icoor = 0; icoor < nbCoor; icoor++){
	    for(int jcoor = 0; jcoor < nbCoor; jcoor++){
	      sumd2phi += d2PhiFace(iFace, i, icoor, jcoor, ig, qr) * qr.weight(ig);
	    } // for jcoor
	  } // for icoor
	} // for i
      } // for ig

      std::cout << "    sum_i integral_Face" << iFace << " phi_i                                        = " << sumphi << std::endl;
      std::cout << "    sum_{i,icoor} integral_Face" << iFace << " dphi_i/dx_coor                       = " << sumdphi << std::endl;
      std::cout << "    sum_{i,icoor,jcoor} integral_Face" << iFace << " d2phi_i/dx_icoor dy_icoor      = " << sumd2phi << std::endl;
    } // for iFace
  } // for k
  std::cout << std::endl;
}

std::ostream& operator << (std::ostream& f, const RefEleDG& ele)
{
  // This part is common to continuous and discontinuous elements
  f << "-------------------------\n";
  f << "Reference Element: " << ele.name << std::endl;
  f << "-------------------------\n";
  f << std::endl << "*** Shape : " << ele.shape << std::endl;
  f << std::endl << "*** Local coordinates of the nodes :\n";
  for(int i=0;i<ele.nbDof;i++){
    f << ele.xi(i) << " " << ele.eta(i) << " " << ele.zeta(i) << std::endl;
  }
  for(int k = 0; k < ele._sqr -> nbQuadRule; k++){
    const QuadRule& qr = ele._sqr -> quadRule(k);
    f << std::endl << "*** Quadrature rule : " << qr.name << std::endl;
    for(int ig=0;ig<qr.nbQuadPt;ig++){
      f << "    - Quadrature point : " << ig << std::endl;
      //      f << "     number and values of basis functions = " << ele.phiQuadPt(ig,qr) << std::endl;
      for(int i=0;i<ele.nbDof;i++){
	f << "      Basif fct " << i << std::endl;
	f << "         Value = " << ele.phi(i,ig,qr) << std::endl;
	f << "         Derivatives = " ;
	for(int icoor=0;icoor<ele.nbCoor;icoor++) f << " " << ele.dPhi(i,icoor,ig,qr);
	f << std::endl;
	f << "         Second derivatives = " ;
	for(int icoor=0;icoor<ele.nbCoor;icoor++)
	  for(int jcoor=0;jcoor<ele.nbCoor;jcoor++)
	  f << " " << ele.d2Phi(i,icoor,jcoor,ig,qr);
	f << std::endl;
      } // for i
    } // for ig
  } // for k

  // This part is for discontinuous elements only
  f << std::endl << "*** Face shape: " << ele.shapeFaces << std::endl;
  f << std::endl << "*** Reference coordinates for the Faces  :" << std::endl;
  for(int iFace = 0; iFace < ele.nbFaces; iFace++){
    f << std::endl << "*** Face #" << iFace << std::endl;
    for(int i = 0; i < ele.nbGeoNodeFaces; i++){
      f << ele.xiFace(iFace,i) << " " << ele.etaFace(iFace,i) << " " << ele.zetaFace(iFace,i) << std::endl;
    }
  
    for(int k = 0; k < ele._sqrFaces -> nbQuadRule; k++){
      const QuadRule& qrFaces = ele._sqrFaces -> quadRule(k);
      f << std::endl << "*** Quadrature rule : " << qrFaces.name << std::endl;
      for(int ig = 0; ig < qrFaces.nbQuadPt; ig++){
	f << "   - Quadrature point : " << ig << std::endl;
	for(int i = 0; i < ele.nbDof; i++){
	  f << "      Basis fct " << i << std::endl;
	  f << "         Value = " << ele.phiFace(iFace, i, ig, qrFaces) << std::endl;
	  f << "         Derivatives = ";
	  for(int icoor = 0; icoor < ele.nbCoor; icoor++)
	    f << " " << ele.dPhiFace(iFace, i, icoor, ig, qrFaces);
	  f << std::endl;
	  f << "         Second derivatives = ";
	  for(int icoor = 0; icoor < ele.nbCoor; icoor++)
	    for(int jcoor = 0; jcoor < ele.nbCoor; jcoor++)
	      f << " " << ele.d2PhiFace(iFace, i, icoor, jcoor, ig, qrFaces);
	  f << std::endl;
	} // for i
      } // for ig
    } // for k
  } // for iFace
  return f;
}
}
