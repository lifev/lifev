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
#include "quadRule.hpp"

QuadRule::QuadRule(const QuadPoint* pt,int _id,std::string _name,
		   ReferenceShapes _shape,int _nbQuadPt,int _degOfExact):
  _pt(pt),shape(_shape),id(_id),name(_name),
  nbQuadPt(_nbQuadPt),degOfExact(_degOfExact)
{
  CONSTRUCTOR("QuadRule");
}

QuadRule::~QuadRule(){
  DESTRUCTOR("QuadRule");
}

std::ostream& operator << (std::ostream& c,const QuadRule& qr)
{
  c << " name: " << qr.name << std::endl;
  c << " shape:" << (int)qr.shape << std::endl;
  c << " id: " << qr.id << std::endl;
  c << " nbQuadPt: " << qr.nbQuadPt << std::endl;
  c << " Points: \n";
  for(int i=0;i<qr.nbQuadPt;i++) c << qr._pt[i] << std::endl;
  return c;
}

SetOfQuadRule::SetOfQuadRule(const QuadRule* qr,int _nb)
  :_qr(qr),nbQuadRule(_nb)
{
  CONSTRUCTOR("SetOfQuadRule");
  _totalNbQuadPoint=0;
  _maxIdQuadRule=0;
  for(int i=0;i<nbQuadRule;i++){
    if(qr[i].shape != qr[0].shape ){
      std::cout << "Error in File : " << __FILE__ << " Line : " << __LINE__ << std::endl;
      ERROR_MSG("All quadrature rules of a set of quadrature rules \n should have the same geometric reference shape");
    }
    _totalNbQuadPoint += qr[i].nbQuadPt;
    if (qr[i].id > _maxIdQuadRule) _maxIdQuadRule = qr[i].id;
  }
  _posQuadRule = new int[ _maxIdQuadRule+1 ]; // don't forget the +1 !
  for(int i=0;i<nbQuadRule;i++){
    _posQuadRule[ qr[i].id ] = i;
  }
};

SetOfQuadRule::~SetOfQuadRule()
{
  DESTRUCTOR("SetOfQuadRule");
}

