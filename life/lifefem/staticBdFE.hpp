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
#ifndef _StaticBDFE_H
#define _StaticBDFE_H

#include "lifeV.hpp"
#include "geoMap.hpp"
#include "refFE.hpp"
#include "geoMap.hpp"
/*!
  \file staticBdFE.h
  \brief Structure for a Static boundary finite element
*/

/*!
  \class StaticBdFE
  \brief A class for static boundary finite element
  \author J.-F. Gerbeau & V. Martin
  \date 09/2002
  
  This class has two purposes:
  \par
  (1) it is a base class for standard boundary element (see CurrentBdFE.h)
  \par
  (2) it is used by refHybridFE as static boundary for a reference element
  
*/

class StaticBdFE{
protected:
  void _comp_meas();
  void _comp_meas_normal();
  void _comp_quad_point_coor();
  UInt _currentId;
#ifdef TEST_PRE
  bool _hasQR;
  bool _hasMeas;
  bool _hasFirstDeriv;
  bool _hasTangent;
  bool _hasNormal;
  bool _hasQuadPtCoor;
#endif
public:
  StaticBdFE(const RefFE& _refFE,const GeoMap& _geoMap);
  StaticBdFE(const RefFE& _refFE,const GeoMap& _geoMap,const QuadRule& _qr);
  StaticBdFE(const RefFE& _refFE,const GeoMap& _geoMap,const QuadRule& _qr,
	     const Real* refcoor, UInt currentid);
  ~StaticBdFE();

  const int nbGeoNode; //!< Number of geometrical nodes 
  const int nbNode;//!< Number of finite element node
  const int nbCoor;//!< Number of coordinates
  const int nbQuadPt;//!< Number of quadrature points
  KNM<Real> point;//!< The point that define the geometry
  const RefFE& refFE;//!< The reference finite element
  const GeoMap& geoMap;//!< The geometical mapping
  const QuadRule& qr;//!< The quadrature rule
  KNM<Real> phi;//!< Values of the basis functions on quadrature points
  KNMK<Real> dPhiRef;//! Values of the derivatives of the basis functions on quadrature points on the reference finite element
  KNMK<Real> dPhi; //!<Values of the derivatives of the basis functions on quadrature points on the current finite element \warning NOT YET IMPLEMENTED
  KNM<Real> phiGeo; //!<Values of the geometric basis functions on quadrature points 
  KNMK<Real> dPhiGeo;//!< Values of the derivatives of the geometric basis functions on quadrature points
  KN<Real> weightMeas; //!< Values of the weight times the measure on the quadrature points
  KN<Real> meas;//!< Values of the measures on the quadrature points
  KNM<Real> normal;//!< Values of the normal on the quadrature points
  KNMK<Real> tangent;//!< Values of the tangents on the quadrature points
  KNMK<Real> metric;//!< Metric tensor on the quadrature points
  KNM<Real> quadPt; //!< Coordinates of the quadrature points on the current element
  //--------------------------------------------------------------------------
#ifdef TEST_PRE
  /*!
    return true if a quadrature rule has been given 
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasQR() const {return _hasQR;}
  /*!
    return true if the measure has been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasMeas() const {return _hasMeas;}
  /*!
    return true if the first derivatives have been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasFirstDeriv() const {return _hasFirstDeriv;}
  /*!
    return true if the tangents have been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasTangent() const {return _hasTangent;}
  /*!
    return true if the normal has been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasNormal() const {return _hasNormal;}
  /*!
    return true if the coordinate of the quadrature points have been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasQuadPtCoor() const {return _hasQuadPtCoor;}
#endif
  /*!
    return the id of the current element (updated with the update* functions)
   */
  inline UInt currentId() const {return _currentId;}
  /*!
    return the coordinate (x,y,z)= F(xi,eta), (F: geo mappping)
    where (xi,eta) are the coor in the ref element
    and   (x,y,z) the coor in the current element
    (if the code is compiled in 2D mode then z=0 and eta is disgarded)
  */
  void coorMap(Real& x,Real& y,Real& z,
	       const Real & xi,const Real & eta) const;
  /*!
    return (x,y,z) = the global coordinates of the quadrature point ig
    in the current element. \warning this function is almost obsolete since if
    you call the function updateMeasQuadPt rather than updateMeas
    (for example), the coordinates of the quadrature points have already
    been computed and may be obtained via quadPt(ig,icoor). This is usually
    much less expensive since it avoids many calls to coorQuadPt
  */
  inline void coorQuadPt(Real& x,Real& y,Real& z,const int ig) const
  {
    ASSERT_PRE(_hasQR,"Needs a quadrature rule")
    coorMap(x,y,z,qr.quadPointCoor(ig,0),qr.quadPointCoor(ig,1) );
  }
  /*!
    Return the measure of the current element
    \warning either updateMeas(...) or updateMeasNormal(...)
    must have been called before
  */
  Real measure() const;//!< Return the measure of the current element
  /*!
    compute the integral over the current boundary element 
    \warning either updateMeas(...) or updateMeasNormal(...)
    must have been called before
  */
  template<typename functor>
  Real  integral(const functor & f) const
  {
    ASSERT_PRE(_hasMeas,"integral needs measure. Call an update function") 
    Real integ(0.0);
    Real x,y,z;
    for(int ig=0;ig<nbQuadPt;ig++){
      coorQuadPt(x,y,z,ig);
      integ += f(x,y,z)*weightMeas(ig);
    }
    return integ;
  }
  /*!
    compute the integral of f . n over the current boundary element 
    \warning  updateMeasNormal(...) must have been called before
  */
  template<typename functor>
    Real integral_n(const functor & f) const
    {
      ASSERT_PRE(_hasNormal,"integral_n needs measure and normal. Call the appropriate update function") 
      Real integ(0.0);
      Real x,y,z;
      Real ret[nbCoor+1];
      Real tmp;
      for(int ig=0;ig<nbQuadPt;++ig){
	coorQuadPt(x,y,z,ig);
	f(x,y,z,ret);
	tmp=0;
	for(int d=0;d<=nbCoor;++d){
	  tmp+=ret[d]*normal(d,ig);
	}
	integ += tmp* weightMeas(ig);
      }
      return integ;
    }
};
#endif
