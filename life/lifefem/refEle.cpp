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
#include "refEle.hpp"
#include <set>

namespace LifeV
{

RefEle::RefEle(string _name,ReferenceShapes _shape,int _nbDof,int _nbCoor,const Fct* phi,const Fct* dPhi,
	       const Fct* d2Phi,const Real* refCoor,const SetOfQuadRule& sqr):
  _sqr(&sqr),
  _phi(phi),_dPhi(dPhi),_d2Phi(d2Phi),_refCoor(refCoor),name(_name),shape(_shape),nbDof(_nbDof),nbCoor(_nbCoor),
  _phiQuad(sqr.totalNbQuadPoint() * _nbDof ),
  _dPhiQuad(sqr.totalNbQuadPoint() * _nbDof * _nbCoor),
  _d2PhiQuad(sqr.totalNbQuadPoint() * _nbDof * _nbCoor * _nbCoor),
  _idxQuad(sqr.maxIdQuadRule()+1),
  _idxDQuad(sqr.maxIdQuadRule()+1),
  _idxD2Quad(sqr.maxIdQuadRule()+1)
{
  CONSTRUCTOR("RefEle");
  int idx=0,idxD=0,idxD2=0;
  for(int k=0;k<sqr.nbQuadRule;k++){
    _idxQuad( sqr.quadRule(k).id )  = idx;
    _idxDQuad( sqr.quadRule(k).id ) = idxD;
    _idxD2Quad( sqr.quadRule(k).id ) = idxD2;
    for(int ig=0;ig<sqr.quadRule(k).nbQuadPt;ig++){
      const QuadPoint& pt = sqr.quadRule(k).quadPoint(ig);
      for(int i=0;i<nbDof;i++){
	_phiQuad[idx] = _phi[i](pt.x(),pt.y(),pt.z());
	idx++;
	for(int icoor=0;icoor<nbCoor;icoor++){
	  _dPhiQuad[idxD] = _dPhi[i*nbCoor+icoor](pt.x(),pt.y(),pt.z());
	  idxD++;
	  for(int jcoor=0;jcoor<nbCoor;jcoor++){
	    _d2PhiQuad[idxD2] = _d2Phi[(i*nbCoor+icoor)*nbCoor+jcoor](pt.x(),pt.y(),pt.z());
	    idxD2++;
	  }
	}
      }
    }
  }
}

RefEle::~RefEle()
{
  DESTRUCTOR("RefEle");
}

void RefEle::check() const
{
  Real sumphi,sumdphi,sumd2phi;
  std::cout << "*** Check " << name << std::endl;
  for(int k=0;k<_sqr->nbQuadRule;k++){
    sumphi=sumdphi=sumd2phi=0.;
    const QuadRule& qr = _sqr->quadRule(k);
    std::cout << std::endl << "    " << qr.name << std::endl;
    for(int ig=0;ig<qr.nbQuadPt;ig++){
      for(int i=0;i<nbDof;i++){
	sumphi += phi(i,ig,qr) * qr.weight(ig);
	for(int icoor=0;icoor<nbCoor;icoor++) sumdphi += dPhi(i,icoor,ig,qr)*qr.weight(ig);
	for(int icoor=0;icoor<nbCoor;icoor++)
	  for(int jcoor=0;jcoor<nbCoor;jcoor++)
	    sumd2phi += d2Phi(i,icoor,jcoor,ig,qr)*qr.weight(ig);
      }
    }
    std::cout << "     sum_i integral phi_i                                   = " << sumphi << std::endl;
    std::cout << "     sum_{i,icoor} integral dphi_i/dx_icoor                 = " << sumdphi << std::endl;
    std::cout << "     sum_{i,icoor,jcoor} integral d2phi_i/dx_icoor dy_icoor = " << sumd2phi << std::endl;
  }
  std::cout << std::endl;
}

std::ostream& operator << (std::ostream& f,const RefEle& ele)
{
  f << "-------------------------\n";
  f << "Reference Element: " << ele.name << std::endl;
  f << "-------------------------\n";
  f << "*** Shape : " << ele.shape << std::endl;
  f << "*** Local coordinate of the nodes :\n";
  for(int i=0;i<ele.nbDof;i++){
    f << ele.xi(i) << " " << ele.eta(i) << " " << ele.zeta(i) << std::endl;
  }
  for(int k=0;k<ele._sqr->nbQuadRule;k++){
    const QuadRule& qr = ele._sqr->quadRule(k);
    f << "\n*** Quadrature rule : " << qr.name << std::endl;
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
      }
    }
  }
  return f;
}
}
