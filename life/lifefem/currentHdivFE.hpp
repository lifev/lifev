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
#ifndef _CURRENTHDIVFE_H
#define _CURRENTHDIVFE_H

#include "lifeV.hpp"
#include "geoMap.hpp"
#include "refHdivFE.hpp"
/*!
  \file CurrentHdivFE.h
  \brief Structure for the current finite element
*/

/*!
  \class CurrentHdivFE
  \brief The class for a finite element 
  \author J.-F. Gerbeau & M. Belhadj & V. Martin
  \date 04/2002 - 08/2002
  
  
*/

class CurrentHdivFE{
private:
  void _comp_jacobian();
  UInt _currentId;
public:
  CurrentHdivFE(const RefHdivFE& _refHdivFE,const GeoMap& _geoMap,const QuadRule& _qr);
  ~CurrentHdivFE();
  const int nbGeoNode; //!< number of Geometric nodes (nodes for GeoMap)
  const int nbNode;    //!< number of dof for the element. Mean bugs may be based on a confusion between these two numbers...
  const int nbCoor;
  const int nbQuadPt;
  const int nbDiag;
  const int nbUpper;
  const int nbPattern;
  KNM<Real> point;     //!< The points that define the geometry of the current element.
  const RefHdivFE& refHdivFE;
  const GeoMap& geoMap;
  const QuadRule& qr;
  KNMK<Real> phi;     //!< Values of the basis vectorial functions in the current element on quad points.
  KNMK<Real> phiRef;  //!< Values of the basis vectorial functions in the reference element on quad points.
  KNM<Real> divPhi;   //!< Values of the basis divergence functions in the reference element on quad points.
  KNMK<Real> jacobian;   //!< Values of the jacobian on quad points.
  KNM<Real> phiGeo;   //!< Values of the basis geometrical functions on quad points.
  KNMK<Real> dPhiGeo; //!< Values of the derivatives of basis geo functions on quad points.
  KN<Real> weightDet;
  KN<Real> detJac;
  //--------------------------------------------------------------------------
  /*!
    return the id of the current element (updated with the update* functions)
   */
  inline UInt currentId() const {return _currentId;}
  /*!
    compute the coordinate (x,y,z)= F(xi,eta,zeta), (F: geo mappping)
    where (xi,eta,zeta) are the coor in the ref element
    and   (x,y,z) the coor in the current element
    (if the code is compiled in 2D mode then z=0 and zeta is disgarded)
  */
  void coorMap(Real& x,Real& y,Real& z,
	       const Real & xi,const Real & eta, const Real &
	       zeta) const;
  /*!
    return (x,y,z) = the global coordinates of the quadrature point ig
    in the current element
  */
  inline void coorQuadPt(Real& x,Real& y,Real& z,const int ig) const
  {
    coorMap(x,y,z,qr.quadPointCoor(ig,0),qr.quadPointCoor(ig,1),
	    qr.quadPointCoor(ig,2));
  }
  //!  patternFirst(i): row index in the element matrix of the i-th term of the pattern
  inline int patternFirst(int i) const{
      return refHdivFE.patternFirst(i);
  }
  //! patternSecond(i): column index in the element matrix of the i-th term of the pattern
  inline int patternSecond(int i) const{
      return refHdivFE.patternSecond(i);
  }
  template<class GEOELE>
  void updatePiola(const GEOELE& geoele)
  {
    // update the definition of the geo points 
    for(int i=0;i<nbGeoNode;i++){
      //      for(int icoor=0;icoor<nbCoor;icoor++)
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    _currentId = geoele.id();
    //! compute the jacobian...
    _comp_jacobian();
    /*!
      Computation of phi (in the current element) as a function of phiRef using the Piola transform.
       phi = 1/detJac * Jacobian * phiRef
    */
    Real mysum;
    for(Int ig = 0 ; ig < nbQuadPt ; ig ++ ){
      for(Int idof = 0; idof < nbNode ; idof ++ ){ 
	//! fixed a bug here : it was "nbGeoNode" before (instead of "nbNode"). V Martin 08/2002
	for(Int icoor = 0; icoor < nbCoor ; icoor ++){
	  mysum = 0.;
	  for(Int jcoor = 0 ; jcoor < nbCoor ; jcoor ++){
	    mysum += jacobian(icoor, jcoor, ig) * phiRef(idof, jcoor, ig);
	  }
	  phi(idof , icoor, ig) = mysum / detJac(ig); 
	}
      }
    }
  }

  /*!
    Return the measure of the current element
  */
  Real measure() const;

};
#endif
