#ifndef _CURRENTFEDG_H
#define _CURRENTFEDG_H

#include <life/lifecore/life.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/refFEDG.hpp>
#include <life/lifefem/geoMap.hpp>
/*!
  \file currentFEDG.h
  \brief Structure for the current discontinuous finite element
*/

namespace LifeV
{

/*!
  \class CurrentFEDG
  \brief The class for a discontinuous finite element
  \author D. A. Di Pietro
  \date 12/2003
*/

class CurrentFEDG{
private:
  void _comp_jacobian();
  void _comp_inv_jacobian();
  void _comp_quad_point_coor();
  void _comp_mass();
  void _comp_inv_mass();
  UInt _currentId;
#ifdef TEST_PRE
  bool _hasJac;
  bool _hasFirstDeriv;
  bool _hasSecondDeriv;
  bool _hasQuadPtCoor;
  bool _hasMass;
#endif
public:
  CurrentFEDG(const RefFEDG& _refFE,const GeoMap& _geoMap,const QuadRule& _qr);
  const int nbGeoNode;
  const int nbNode;
  const int nbCoor;
  const int nbQuadPt;
  const int nbDiag;
  const int nbUpper;
  const int nbPattern;
  KNM<Real> point;
  const RefFEDG& refFE;
  const GeoMap& geoMap;
  const QuadRule& qr;
  KNM<Real> phi;
  KNMK<Real> dPhiRef;
  KNMKL<Real> dPhiRef2;
  //
  KNMK<Real> phiDer;
  KNMK<Real> jacobian;
  KNMK<Real> tInvJac;
  KNM<Real> phiGeo;
  KNMK<Real> dPhiGeo;
  KN<Real> weightDet;
  KN<Real> detJac;
  KNM<Real> quadPt;//!< Coordinates of the quadrature points on the current element

  KNM<Real> mass;//!< Local mass matrix
  KNM<Real> invMass;//!< Inverse of the local mass matrix

  //second derivatives not yet done (four dimensions KNMKL ?)
  KNMKL<Real> phiDer2;
  //--------------------------------------------------------------------------
  /*!
    return the diameter of the element in the 1-norm
  */
  Real diameter() const;
  /*!
    return the id of the current element (updated with the update* functions)
   */
  inline UInt currentId() const {return _currentId;}
#ifdef TEST_PRE
  /*!
    return true if the determinant has been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasJac() const {return _hasJac;}
  /*!
    return true if the first derivatives have been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasFirstDeriv() const {return _hasFirstDeriv;}
  /*!
    return true if the second derivatives have been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasSecondDeriv() const {return _hasSecondDeriv;}
  /*!
    return true if the coordinate of the quadrature points have been updated
    (can ONLY be used if TEST_PRE is defined at compilation time)
   */
  inline bool hasQuadPtCoor() const {return _hasQuadPtCoor;}
  /*!
    return true if the local mass matrix has been updated
    (can ONLY be used if TEST_PRE is defined at compilation time
  */
  inline bool hasMass() const {return _hasMass;}
#endif
  /*!
    compute the coordinate (x,y,z)= F(xi,eta,zeta), (F: geo mapping)
    where (xi,eta,zeta) are the coor in the ref element
    and   (x,y,z) the coor in the current element
    (if the code is compiled in 2D mode then z=0 and zeta is disregarded)
  */
  void coorMap(Real& x,Real& y,Real& z,
           const Real & xi,const Real & eta, const Real &
           zeta) const;
  /*!  return (x,y,z) = the global coordinates of the quadrature point ig
    in the current element. \warning this function is almost obsolete since if
    you call the function updateFirstDerivQuadPt rather than updateFirstDeriv
    (for example), the coordinates of the quadrature points have already
    been computed and may be obtained via quadPt(ig,icoor). This is usually
    much less expensive since it avoids many calls to coorQuadPt
  */
  inline void coorQuadPt(Real& x,Real& y,Real& z,const int ig) const
  {
    coorMap(x,y,z,qr.quadPointCoor(ig,0),qr.quadPointCoor(ig,1),
        qr.quadPointCoor(ig,2));
  }
  //!  patternFirst(i): row index in the element matrix of the i-th term of the pattern
  inline int patternFirst(int i) const{
      return refFE.elPattern.patternFirst(i);
  }
  //! patternSecond(i): column index in the element matrix of the i-th term of the pattern
  inline int patternSecond(int i) const{
      return refFE.elPattern.patternSecond(i);
  }

  inline int patternFirstFaces(int i) const{
    return refFE.facePattern.patternFirst(i);
  }

  inline int patternSecondFaces(int i) const{
    return refFE.facePattern.patternSecond(i);
  }

  /*!
    minimal update: we just identify the id of the current element and
    of the geometrical points
  */
  template<class GEOELE>
  void update(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = false;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
  }
  /*!
    compute the arrays detJac, weightDet, jacobian on
    the current element
  */
  template<class GEOELE>
  void updateJac(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the jacobian...
    _comp_jacobian();
  }

  /*!
    compute the arrays detJac, weightDet, jacobian, mass and
    invMass on the current element
  */
template<class GEOELE>
  void updateJacMass(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the jacobian...
    _comp_jacobian();

    // the local mass matrix and its inverse
    _comp_inv_mass();
  }

  /*!
    compute the arrays detJac, weightDet, jacobian and quadPt
    on the current element
  */
  template<class GEOELE>
  void updateJacQuadPt(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = true;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the jacobian...
    _comp_jacobian();
    // and the coordinates of the quadrature points
    _comp_quad_point_coor();
  }

  /*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer on the current element
  */
  template<class GEOELE>
  void updateFirstDeriv(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      // for(int icoor=0;icoor<nbCoor;icoor++)
      // point(i,icoor) =  geoele.coor(i+1,icoor+1);
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();
    // product InvJac by dPhiRef:
    Real x;
    for(int ig=0;ig<nbQuadPt;ig++){
      for(int j=0;j<nbNode;j++){
    for(int icoor=0;icoor<nbCoor;icoor++){
      x = 0.;
      for(int jcoor=0;jcoor<nbCoor;jcoor++){
        x += tInvJac(icoor,jcoor,ig)*dPhiRef(j,jcoor,ig) ;
      }
      phiDer(j,icoor,ig)=x;
    }
      }
    }
  }

  /*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer and quadPt on the current element
  */
  template<class GEOELE>
  void updateFirstDerivQuadPt(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = true;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();
    // product InvJac by dPhiRef:
    Real x;
    for(int ig=0;ig<nbQuadPt;ig++){
      for(int j=0;j<nbNode;j++){
    for(int icoor=0;icoor<nbCoor;icoor++){
      x = 0.;
      for(int jcoor=0;jcoor<nbCoor;jcoor++){
        x += tInvJac(icoor,jcoor,ig)*dPhiRef(j,jcoor,ig) ;
      }
      phiDer(j,icoor,ig)=x;
    }
      }
    }
    // and the coordinates of the quadrature points
    _comp_quad_point_coor();
  }

  template<class GEOELE>
  void updateFirstDerivQuadPtMass(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = true;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();
    // product InvJac by dPhiRef:
    Real x;
    for(int ig=0;ig<nbQuadPt;ig++){
      for(int j=0;j<nbNode;j++){
    for(int icoor=0;icoor<nbCoor;icoor++){
      x = 0.;
      for(int jcoor=0;jcoor<nbCoor;jcoor++){
        x += tInvJac(icoor,jcoor,ig)*dPhiRef(j,jcoor,ig) ;
      }
      phiDer(j,icoor,ig)=x;
    }
      }
    }
    // and the coordinates of the quadrature points
    _comp_quad_point_coor();

    // the local mass matrix and its inverse
    _comp_inv_mass();
  }

  /*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer2 on the current element
  */
  template<class GEOELE>
  void updateSecondDeriv(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = false;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      // for(int icoor=0;icoor<nbCoor;icoor++)
      // point(i,icoor) =  geoele.coor(i+1,icoor+1);
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();

    Real x;
    for(int ig=0;ig<nbQuadPt;ig++){
     for(int j=0;j<nbNode;j++){
      for(int icoor=0;icoor<nbCoor;icoor++){
        for(int jcoor=0;jcoor<nbCoor;jcoor++){
           x=0.;
          for (int k1=0;k1<nbCoor;k1++){
           for (int k2=0;k2<nbCoor;k2++){
            x += tInvJac(icoor,k1,ig)*dPhiRef2(j,k1,k2,ig)*tInvJac(jcoor,k2,ig);
           }
          }
          phiDer2(j,icoor,jcoor,ig)=x;
         }
        }
       }
      }
  }

  /*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer2 on the current element
  */
  template<class GEOELE>
  void updateSecondDerivQuadPt(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = false;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = true;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      // for(int icoor=0;icoor<nbCoor;icoor++)
      // point(i,icoor) =  geoele.coor(i+1,icoor+1);
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();

    Real x;
    for(int ig=0;ig<nbQuadPt;ig++){
     for(int j=0;j<nbNode;j++){
      for(int icoor=0;icoor<nbCoor;icoor++){
        for(int jcoor=0;jcoor<nbCoor;jcoor++){
           x=0.;
          for (int k1=0;k1<nbCoor;k1++){
           for (int k2=0;k2<nbCoor;k2++){
            x += tInvJac(icoor,k1,ig)*dPhiRef2(j,k1,k2,ig)*tInvJac(jcoor,k2,ig);
           }
          }
          phiDer2(j,icoor,jcoor,ig)=x;
         }
        }
       }
      }
    // and the coordinates of the quadrature points
    _comp_quad_point_coor();
  }
  /*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer, phiDer2 on the current element
  */
  template<class GEOELE>
  void updateFirstSecondDeriv(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = false;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      // for(int icoor=0;icoor<nbCoor;icoor++)
      // point(i,icoor) =  geoele.coor(i+1,icoor+1);
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();

    Real x1,x2;
    for(int ig=0;ig<nbQuadPt;ig++){
     for(int j=0;j<nbNode;j++){
      for(int icoor=0;icoor<nbCoor;icoor++){
    x1=0.;
        for(int jcoor=0;jcoor<nbCoor;jcoor++){
          x2=0.;
          for (int k1=0;k1<nbCoor;k1++){
           for (int k2=0;k2<nbCoor;k2++){
            x2 += tInvJac(icoor,k1,ig)*dPhiRef2(j,k1,k2,ig)*tInvJac(jcoor,k2,ig);
           }
          }
          phiDer2(j,icoor,jcoor,ig)=x2;
         }
         phiDer(j,icoor,ig)=x1;
        }
       }
      }
  }
  /*!
    compute the arrays detJac, weightDet, jacobian,
    tInvJac, phiDer, phiDer2 on the current element
  */
  template<class GEOELE>
  void updateFirstSecondDerivQuadPt(const GEOELE& geoele)
  {
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = true;
    _hasMass = false;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      // for(int icoor=0;icoor<nbCoor;icoor++)
      // point(i,icoor) =  geoele.coor(i+1,icoor+1);
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();

    Real x1,x2;
    for(int ig=0;ig<nbQuadPt;ig++){
     for(int j=0;j<nbNode;j++){
      for(int icoor=0;icoor<nbCoor;icoor++){
    x1=0.;
        for(int jcoor=0;jcoor<nbCoor;jcoor++){
          x2=0.;
          for (int k1=0;k1<nbCoor;k1++){
           for (int k2=0;k2<nbCoor;k2++){
            x2 += tInvJac(icoor,k1,ig)*dPhiRef2(j,k1,k2,ig)*tInvJac(jcoor,k2,ig);
           }
          }
          phiDer2(j,icoor,jcoor,ig)=x2;
         }
         phiDer(j,icoor,ig)=x1;
        }
       }
      }
    // and the coordinates of the quadrature points
    _comp_quad_point_coor();
  }

  template<class GEOELE>
    void updateFirstSecondDerivQuadPtMass(const GEOELE& geoele)
{
#ifdef TEST_PRE
    _hasJac = true;
    _hasFirstDeriv = true;
    _hasSecondDeriv = true;
    _hasQuadPtCoor = true;
    _hasMass = true;
#endif
    _currentId = geoele.id();
    // update the definition of the geo points
    for(int i=0;i<nbGeoNode;i++){
      // for(int icoor=0;icoor<nbCoor;icoor++)
      // point(i,icoor) =  geoele.coor(i+1,icoor+1);
      point(i,0) = geoele.point(i+1).x();
      point(i,1) = geoele.point(i+1).y();
      point(i,2) = geoele.point(i+1).z();
    }
    // compute the inverse jacobian...
    _comp_inv_jacobian();

    Real x1,x2;
    for(int ig=0;ig<nbQuadPt;ig++){
     for(int j=0;j<nbNode;j++){
      for(int icoor=0;icoor<nbCoor;icoor++){
    x1=0.;
        for(int jcoor=0;jcoor<nbCoor;jcoor++){
          x2=0.;
          for (int k1=0;k1<nbCoor;k1++){
           for (int k2=0;k2<nbCoor;k2++){
            x2 += tInvJac(icoor,k1,ig)*dPhiRef2(j,k1,k2,ig)*tInvJac(jcoor,k2,ig);
           }
          }
          phiDer2(j,icoor,jcoor,ig)=x2;
         }
         phiDer(j,icoor,ig)=x1;
        }
       }
      }

    // the coordinates of the quadrature points...
    _comp_quad_point_coor();

    // the local mass matrix and its inverse
    _comp_inv_mass();
  }


  /*!
    Return the measure of the current element
  */
  Real measure() const;
};
}
#endif


//==============================================================================
// Useful stuff
//==============================================================================
template<typename P>
void Leverrier(LifeV::KNM<P>& A, LifeV::KNM<P>& invA);

template<typename MATRIX>
void MatMult(MATRIX& A, MATRIX& B, MATRIX& C);
