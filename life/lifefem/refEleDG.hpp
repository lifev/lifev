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
#ifndef _REFELEDG_H_INCLUDE
#define _REFELEDG_H_INCLUDE


#include <life/lifefem/refEle.hpp>
#include <life/lifefem/geoMap.hpp>

/*!
  \file refEleDG.h
  \brief Structure for a discontinuous element
*/
namespace LifeV
{
/*!
  \class RefFEDG
  \brief The class for a reference discontinuous element
  \author D. A. Di Pietro
  \date 12/2003

*/

// Base class for Discontinuous Galerkin
class RefEleDG:public RefEle{
 protected:
  const SetOfQuadRule* _sqrFaces;
  const Real* _refCoorFaces;

 public:
  const ReferenceShapes shapeFaces;
  const int nbGeoNodeFaces;
  const int nbFaces; // WARNING: face numbering starts from 0
  const GeoMap& geoMap; // geoMap is used to map each face ((s, t) coordinates) in the reference element ((xi, eta, zeta) coordinates)

 private:
  KNM<Real> _phiQuadFaces;
  KNM<Real> _dPhiQuadFaces;
  KNM<Real> _d2PhiQuadFaces;
  KN<int> _idxQuadFaces;
  KN<int> _idxDQuadFaces;
  KN<int> _idxD2QuadFaces;

 public:

  RefEleDG(std::string _name, ReferenceShapes _shape,
       int _nbDof, int _nbCoor, 
       const Fct* phi, const Fct* dPhi, const Fct* d2Phi, 
       const Real* refCoor, const SetOfQuadRule& sqr, 
       ReferenceShapes _shapeFaces, 
       int _nbFaces, int _nbGeoNodeFaces, 
       const Real* refCoorFaces, const SetOfQuadRule& sqrFaces, const GeoMap& _geoMap);
  ~RefEleDG();

  inline Real xiFace(int iFace, int i) const{
    ASSERT_BD(iFace < nbFaces && i < nbGeoNodeFaces)
      return _refCoorFaces[iFace * nbGeoNodeFaces * 3 + 3 * i];
  }

  inline Real etaFace(int iFace, int i) const{
    ASSERT_BD(iFace < nbFaces && i < nbGeoNodeFaces)
      return _refCoorFaces[iFace * nbGeoNodeFaces * 3 + 3 * i + 1];
  }

  inline Real zetaFace(int iFace, int i) const{
    ASSERT_BD(iFace < nbFaces && i < nbGeoNodeFaces)
      return _refCoorFaces[iFace * nbGeoNodeFaces * 3 + 3 * i + 2];
  }

  inline Real refCoorFaces(int iFace, int i, int icoor){
    ASSERT_BD(iFace < nbFaces && i < nbGeoNodeFaces && icoor < nbCoor)
      return _refCoorFaces[iFace * nbGeoNodeFaces * 3 + 3 * i + icoor];
  }

  void FaceToElCoord(Real& xi, Real& eta, Real& zeta, Real s, Real t, int iFace){
    xi = eta = zeta = 0.;
    
    for(int i = 0; i < nbGeoNodeFaces; i++){
      xi += xiFace(iFace, i) * geoMap.phi(i, s, t, 0.);
      eta += etaFace(iFace, i) * geoMap.phi(i, s, t, 0.);
#if defined(THREEDIM)
      zeta += zetaFace(iFace, i) * geoMap.phi(i, s, t, 0.);
#endif
    
    }
  }
  // The qr argument is a quadrature rule on a face!!
  inline double phiFace(int iFace, int i, int ig, const QuadRule& qr) const{
    ASSERT_BD(iFace < nbFaces && i < nbDof && ig < qr.nbQuadPt)
      return _phiQuadFaces(iFace, _idxQuadFaces(qr.id) + ig * nbDof + i);
  }

  inline double dPhiFace(int iFace, int i, int icoor, int ig, const QuadRule& qr) const{
    ASSERT_BD(iFace < nbFaces && i < nbDof && icoor < nbCoor && ig < qr.nbQuadPt)
      return _dPhiQuadFaces(iFace, _idxDQuadFaces(qr.id) + (ig * nbDof + i) * nbCoor + icoor);
  }

  inline double d2PhiFace(int iFace, int i, int icoor, int jcoor, int ig, const QuadRule& qr) const{
    ASSERT_BD(iFace < nbFaces && i < nbDof && icoor < nbCoor && jcoor < nbCoor && ig < qr.nbQuadPt)
      return _d2PhiQuadFaces(iFace, _idxD2QuadFaces(qr.id) + ((ig * nbDof + i) * nbCoor + icoor) * nbCoor + jcoor);
  }

  void check() const;
  friend std::ostream& operator << (std::ostream& f, const RefEleDG& fe);

};

//======================================================================
//
//                            Discontinuous P1 (3D)
//
//======================================================================
/*                 
                4
               / .  
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/

// fToP for a linear tetrahedra (see basisElSh.cc):
// 1 3 2
// 1 2 4
// 2 3 4 
// 1 4 3
// with (see refEle.h):
// 1: 0. 0. 0.
// 2: 1. 0. 0.
// 3: 0. 1. 0.
// 4: 0. 0. 1.

//
static const Real refCoorFaces_P1_DG_3D[36] = {0., 0., 0.,
                           0., 1., 0.,
                           1., 0., 0.,

                           0., 0., 0.,
                           1., 0., 0.,
                           0., 0., 1.,

                           1., 0., 0.,
                           0., 1., 0.,
                           0., 0., 1.,

                           0., 0., 0.,
                           0., 0., 1.,
                           0., 1., 0.};
                        

/* static const GeoMap& geoMap_P1_DG_3D = geoLinearTria; */
}
#endif
                 
