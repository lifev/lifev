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
#include "refHdivFE.hpp"
#include <set>

namespace LifeV
{

RefHdivFE::RefHdivFE(string _name, int _type, ReferenceShapes _shape,int _nbDofPerVertex,
		     int _nbDofPerEdge,int _nbDofPerFace,int _nbDofPerVolume,
		     int _nbDof,int _nbCoor,const FCT* fctphi,const FCT* fctdivPhi,
		     const Real* refCoor,const SetOfQuadRule& sqr,PatternType _patternType):
  LocalDofPattern(_nbDof,_nbDofPerVertex,_nbDofPerEdge,_nbDofPerFace,_nbDofPerVolume,_patternType),
  _sqr(&sqr),_phi(fctphi),_divPhi(fctdivPhi),
  _refCoor(refCoor),_phiQuad(sqr.totalNbQuadPoint() * _nbDof * _nbCoor ),
  _divPhiQuad(sqr.totalNbQuadPoint() * _nbDof),
  _idxQuad(sqr.maxIdQuadRule()+1),_idxDQuad(sqr.maxIdQuadRule()+1),
  name(_name),type(_type),shape(_shape),
  nbDof(_nbDof),nbCoor(_nbCoor)
{
  CONSTRUCTOR("RefHdivFE");

 int idx=0,idxD=0;
  for(int k=0;k<sqr.nbQuadRule;k++){
    _idxQuad( sqr.quadRule(k).id )  = idx;
    _idxDQuad( sqr.quadRule(k).id ) = idxD;
    for(int ig=0;ig<sqr.quadRule(k).nbQuadPt;ig++){
      const QuadPoint& pt = sqr.quadRule(k).quadPoint(ig);
      for(int i=0;i<nbDof;i++){
        for(int icomp=0;icomp<nbCoor;icomp++){
          _phiQuad[idx] = phi(i,icomp,pt.x(),pt.y(),pt.z());
          idx++;
        }
        _divPhiQuad[idxD] = _divPhi[i](pt.x(),pt.y(),pt.z());
        idxD++;
      }
    }
  }
}
RefHdivFE::~RefHdivFE()
{
  DESTRUCTOR("RefHdivFE");
}

void RefHdivFE::check() const
{
  /*
  Real sumphi,sumdivphi;
  cout << "*** Check " << name << endl;
  for(int k=0;k<_sqr->nbQuadRule;k++){
    sumphi=sumdivphi=0.;
    const QuadRule& qr = _sqr->quadRule(k);
    cout << endl << "    " << qr.name << endl;
    for(int ig=0;ig<qr.nbQuadPt;ig++){
      for(int i=0;i<nbDof;i++){
	sumphi += phi(i,ig,qr) * qr.weight(ig);
	for(int icoor=0;icoor<nbCoor;icoor++) sumdivphi += divPhi(i,icoor,ig,qr)*qr.weight(ig);
      }
    }
    cout << "     sum_i integral phi_i                                   = " << sumphi << endl;
    cout << "     sum_{i,icoor} integral divphi_i/dx_icoor               = " << sumdivphi << endl;
  }
  cout << endl;
  */
}

ostream& operator << (ostream& f,const RefHdivFE& fe)
{
  f << "-------------------------\n";
  f << "Reference Finite Element: " << fe.name << endl;
  f << "-------------------------\n";
  f << "*** Shape : " << fe.shape << endl;
  f << "*** Local coordinate of the nodes :\n";
  for(int i=0;i<fe.nbDof;i++){
    f << fe.xi(i) << " " << fe.eta(i) << " " << fe.zeta(i) << endl;
  }
  f << "*** Pattern :\n";
  for(int i=0;i<fe.nbPattern();i++) f << "(" << fe.patternFirst(i) << "," << fe.patternSecond(i) << ") \n";
  for(int k=0;k<fe._sqr->nbQuadRule;k++){
    const QuadRule& qr = fe._sqr->quadRule(k);
    f << "\n*** Quadrature rule : " << qr.name << endl;
    for(int ig=0;ig<qr.nbQuadPt;ig++){
      f << "    - Quadrature point : " << ig << endl;
      //      f << "     number and values of basis functions = " << fe.phiQuadPt(ig,qr) << endl;
      for(int i=0;i<fe.nbDof;i++){
	f << "      Basif fct " << i << endl;
        for(int icomp=0;icomp<fe.nbCoor;icomp++){
          f << "         Value = " << fe.phi(i,icomp,ig,qr) << endl;
        }
	f << "         Divergence = " ;
        f << " " << fe.divPhi(i,ig,qr);
        f << endl;
      }
    }
  }
  return f;
}
}
