#ifndef _CURRENTBDDG_H
#define _CURRENTBDDG_H

#include <life/lifecore/life.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/geoMapDG.hpp>
#include <life/lifefem/refFEDG.hpp>
#include <life/lifefem/staticBdFE.hpp>
#include <life/lifefem/currentFEDG.hpp>

/*!
  \file currentBdDG.h
  \brief Structure for a boundary element (discontinuous elements)
*/

namespace LifeV
{
/*!
  \class CurrentBdDG
  \brief A class for boundary finite element (discontinuous elements)
  \author D. A. Di Pietro
  \date 12/2003

  This class provides basic data structures to handle face integration
  in discontinuous finite elements. Child classes are then provided to handle
  integration on internal and boundary faces.
*/

class CurrentBdDG:public StaticBdFE{
#ifdef TEST_PRE
 protected:
  bool _hasPointAd;
  bool _hasMassAd;
#endif

 private:
  /*!
    _comp_phi_ad_d_phi_face() computes the value of adjacent element's basis functions
    on face quadrature nodes
  */
  void _comp_phi_ad_d_phi_face();

  /*!
    _comp_inv_jacobian_ad() computes the inverse transposed jacobian of an adjacent element
    on face quad points once the pointAd array has been updated and stores it in tInvJacAd
    property
  */
  void _comp_inv_jacobian_ad();

 public:
  CurrentBdDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr);
  ~CurrentBdDG();

  int faceIDAd;
  const RefFEDG& refFEDG;

  const int nbNodeAd; //!< Number of local dof on each adjacent element
  const int nbGeoNodeAd; //!< Number of geo node for the adjacent elements
  const int nbCoorAd; //!< Number of coordinates for the adjacent elements

  KNM<Real> pointAd;
  const GeoMapDG& geoMapAd; //!< Geometrical mapping of the adjacent elements

  KNM<Real> phiAd; //!< Value of the basis functions of the adjacent elements on quadrature points

  KNMK<Real> dPhiRefAd; //!< Value of the derivatives of the basis functions of the first adjacent element with respect to local coordinates on face quadrature points

  KNMK<Real> jacobianAd; //!< Jacobian on face quadrature points
  KNMK<Real> tInvJacAd; //!< Transposed inverse Jacobian on face quadrature points
  KN<Real> detJacAd; //!< Determinant on face quad points
  KN<Real> weightDetAd; //!< Determinant times weight on face quad points

  KNMK<Real> phiDerAd; //!< Value of the derivatives of the basis functions of the first adjacent element on face quadrature points

  KNM<Real> phiGeoAd; //!< Value of the geo functions of an adjacent element on face quad points
  KNMK<Real> dPhiGeoAd; //!< Value of the derivatives of the geo functions of an adjacent element on face quad points

  KNM<Real> massAd; //!< Value of the local mass matrix of the adjacent element
  KNM<Real> invMassAd; //!< Value of the inverse of the local mass matrix of the adjacent element

#ifdef TEST_PRE
  inline bool hasPointAd() const{return _hasPointAd;}
  inline bool hasMassAd() const{return _hasMassAd;}
#endif

  /*!
    patternFirst(i): row index in the element matrix of the i-th term of the pattern
  */
  inline int patternFirst(int i) const{
    return refFEDG.facePattern.patternFirst(i);
  }

  /*!
    patternSecond(i): column index in the element matrix of the i-th term of the pattern
  */
  inline int patternSecond(int i) const{
    return refFEDG.facePattern.patternSecond(i);
  }

  //----------------------------------------------------------------------------
  /*!
    Compute only the coordinates of the nodes on the current boundary element
  */
  template<class GEOELE>
    void update(const GEOELE& geoele)
    {
#ifdef TEST_PRE
      _hasMeas       = false;
      _hasTangent    = false;
      _hasNormal     = false;
      _hasQuadPtCoor = false;
      _hasFirstDeriv = false;
      _hasPointAd    = false;
      _hasMassAd     = false;
#endif
      _currentId = geoele.id();
      // update the definition of the geo points
      for(int i=0;i<nbGeoNode;i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }
    }

  //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent
    on the current boundary element
  */
  template<class GEOELE>
    void updateMeas(const GEOELE& geoele)
    {
#ifdef TEST_PRE
      _hasMeas       = true;
      _hasTangent    = true;
      _hasNormal     = false;
      _hasQuadPtCoor = false;
      _hasFirstDeriv = false;
      _hasPointAd    = false;
      _hasMassAd     = false;
#endif
      _currentId = geoele.id();
      // update the definition of the geo points
      for(int i=0;i<nbGeoNode;i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }
      // compute the measure
      _comp_meas();
    };

  //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent
    and quadrature points on the current boundary element
  */
  template<class GEOELE>
    void updateMeasQuadPt(const GEOELE& geoele)
    {
#ifdef TEST_PRE
      _hasMeas       = true;
      _hasTangent    = true;
      _hasNormal     = false;
      _hasQuadPtCoor = true;
      _hasFirstDeriv = false;
      _hasPointAd    = false;
      _hasMassAd     = false;
#endif
      _currentId = geoele.id();
      // update the definition of the geo points
      for(int i=0;i<nbGeoNode;i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }
      // compute the measure
      _comp_meas();
      // compute the coordinates of the quad points
      _comp_quad_point_coor();
    };

  //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent
    and normal on the current boundary element
  */
  template<class GEOELE>
    void updateMeasNormal(const GEOELE& geoele)
    {
#ifdef TEST_PRE
      _hasMeas = true;
      _hasTangent = true;
      _hasNormal = true;
      _hasQuadPtCoor = false;
      _hasFirstDeriv = false;
      _hasPointAd    = false;
      _hasMassAd     = false;
#endif
      _currentId = geoele.id();
      // update the definition of the geo points
      for(int i=0;i<nbGeoNode;i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }
      // compute the measure and the normal
      _comp_meas_normal();
    }

  //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent,
    normal and quadrature points on the current boundary element
  */
  template<class GEOELE>
    void updateMeasNormalQuadPt(const GEOELE& geoele)
    {
#ifdef TEST_PRE
      _hasMeas = true;
      _hasTangent = true;
      _hasNormal = true;
      _hasQuadPtCoor = true;
      _hasFirstDeriv = false;
      _hasPointAd    = false;
      _hasMassAd     = false;
#endif
      _currentId = geoele.id();
      // update the definition of the geo points
      for(int i=0;i<nbGeoNode;i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }
      // compute the measure and the normal
      _comp_meas_normal();
      // compute the coordinates of the quad points
      _comp_quad_point_coor();
    }

  //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent,
    normal, quadrature points and phiDer on the current
    boundary element
  */
  template<class GEOELE, class GEOELEAD>
    void updateMeasNormalQuadPtFirstDerivAd(const GEOELE& geoele, const GEOELEAD& geoelead)
    {
#ifdef TEST_PRE
      _hasMeas = true;
      _hasTangent = true;
      _hasNormal = true;
      _hasQuadPtCoor = true;
      _hasFirstDeriv = true;
      _hasPointAd    = true ;
      _hasMassAd     = false;
#endif
      _currentId = geoele.id();
      faceIDAd = geoele.pos_first() - 1; // 0-base!!

      // update the definition of the geo points
      for(int i = 0; i < nbGeoNode; i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }

      // update the value of adjacent element's basis functions
      _comp_phi_ad_d_phi_face();

      // compute the measure and the normal
      _comp_meas_normal();

      // compute the coordinates of the quad points
      _comp_quad_point_coor();

      // update the definition of the geo points of the adjacent element
      for(int i = 0; i < nbGeoNodeAd; i++){
    pointAd(i, 0) = geoelead.point(i+1).x();
    pointAd(i, 1) = geoelead.point(i+1).y();
    pointAd(i, 2) = geoelead.point(i+1).z();
      }

      // compute the jacobian on quadrature points
      _comp_inv_jacobian_ad();

      // product tInvJacAd by dPhiRefAd
      Real x;
      for(int ig = 0; ig < nbQuadPt; ig++){
    for(int j = 0; j < nbNodeAd; j++){
      for(int icoor = 0; icoor < nbCoorAd; icoor++){
        x = 0.;
        for(int jcoor = 0; jcoor < nbCoorAd; jcoor++){
          x += tInvJacAd(icoor, jcoor, ig) * dPhiRefAd(j, jcoor, ig);
        } // for jcoor
        phiDerAd(j, icoor, ig) = x;
      } // for icoor
    } // for j
      } // for ig
    }

 //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent,
    normal, quadrature points, phiDer, mass and invMass
    on the current boundary element
  */
  template<class GEOELE, class GEOELEAD, class CURRFEAD>
    void updateMeasNormalQuadPtFirstDerivMassAd(const GEOELE& geoele, const GEOELEAD& geoelead, CURRFEAD& currfead)
    {
#ifdef TEST_PRE
      _hasMeas = true;
      _hasTangent = true;
      _hasNormal = true;
      _hasQuadPtCoor = true;
      _hasFirstDeriv = true;
      _hasPointAd    = true ;
      _hasMassAd     = false;
#endif
      _currentId = geoele.id();
      faceIDAd = geoele.pos_first() - 1; // 0-base!!

      // update the definition of the geo points
      for(int i = 0; i < nbGeoNode; i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }

      // update the value of adjacent element's basis functions
      _comp_phi_ad_d_phi_face();

      // compute the measure and the normal
      _comp_meas_normal();

      // compute the coordinates of the quad points
      _comp_quad_point_coor();

      // update the definition of the geo points of the adjacent element
      for(int i = 0; i < nbGeoNodeAd; i++){
    pointAd(i, 0) = geoelead.point(i+1).x();
    pointAd(i, 1) = geoelead.point(i+1).y();
    pointAd(i, 2) = geoelead.point(i+1).z();
      }

      // compute the jacobian on quadrature points
      _comp_inv_jacobian_ad();

      // product tInvJacAd by dPhiRefAd
      Real x;
      for(int ig = 0; ig < nbQuadPt; ig++){
    for(int j = 0; j < nbNodeAd; j++){
      for(int icoor = 0; icoor < nbCoorAd; icoor++){
        x = 0.;
        for(int jcoor = 0; jcoor < nbCoorAd; jcoor++){
          x += tInvJacAd(icoor, jcoor, ig) * dPhiRefAd(j, jcoor, ig);
        } // for jcoor
        phiDerAd(j, icoor, ig) = x;
      } // for icoor
    } // for j
      } // for ig

      // update massAd and invMassAd arrays
      currfead.updateJacMass(geoelead);

      massAd = currfead.mass;
      invMassAd = currfead.invMass;
    }
};
}
#endif
