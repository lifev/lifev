#ifndef _CURRENTIFDG_H
#define _CURRENTIFDG_H

#include <life/lifefem/currentBdDG.hpp>

namespace LifeV
{
class CurrentIFDG: public CurrentBdDG{
#ifdef TEST_PRE
 protected:
  bool _hasMassOp;
  bool _hasre;
#endif
 private:
  void _comp_corresp_node_op();
  void _comp_phi_dphi_op();
  void _comp_inv_jacobian_op();
/*   void _comp_re(); */
 public:
  CurrentIFDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr);
  ~CurrentIFDG();

  int faceIDOp;
  
  KNM<Real> pointOp; //!< Geometric nodes of the opposite element
  KNM<Real> refCoorFaceOp; //!< Local coordinates of the geometric nodes of the face (opposite element intrinsic coordinate system)
  KNM<Real> locQuadPtCoor; //!< Local coordinates of face quadrature points (opposite element intrinsic coordinate system)

  KNM<Real> phiOp; //!< Opposite element's basis functions on face quadrature nodes
  KNMK<Real> dPhiOp; //!< Opposite element's basis functions derivatives with respect to local coordinates on face quadrature nodes

  KNM<Real> phiGeoOp; //!< Opposite element's geometric basis functions on face quadrature nodes
  KNMK<Real> dPhiGeoOp; //!< Opposite element's geometric basis functions derivatives with respect to local coordinates on face quadrature nodes

  KNMK<Real> jacobianOp; //<! Jacobian of the opposite element on face quadrature nodes
  KNMK<Real> tInvJacOp; //<! Inverse transpost Jacobian of the opposite element on face quadrature nodes
  KN<Real> detJacOp; //<! Determinant of the Jacobian of the opposite element on face quadrature nodes
  KN<Real> weightDetOp; //<! Determinant of the Jacobian of the opposite element on face quadrature nodes times quadrature rule weight

  KNMK<Real> phiDerOp;

/*   KNM<Real> massOp; //!< Value of the local mass matrix of the opposite element */
/*   KNM<Real> invMassOp; //!< Value of the inverse of the local mass matrix of the opposite element */

/*   KNMK<Real> re; //!< Coefficients of the approximation of r_e */

/*   KNM<Real> massAdOp; */
/*   KNM<Real> invMassAdOp; */
/*   KNM<Real> phiAdOp; */

  KN<int> correspNodeOp;
 
#ifdef TEST_PRE
  inline bool hasMassOp() const {return _hasMassOp;}
  inline bool hasre() const {return _hasre;}
#endif

  void FaceToOpElCoor(Real &xi, Real &eta, Real &zeta, Real s, Real t) {
    xi = eta = zeta = 0.;

     for(int i = 0; i < nbGeoNode; i++){
      xi += refCoorFaceOp(i, 0) * geoMap.phi(i, s, t, 0.);
      eta += refCoorFaceOp(i, 1) * geoMap.phi(i, s, t, 0.);
#if defined(THREEDIM)
      zeta += refCoorFaceOp(i, 2) * geoMap.phi(i, s, t, 0.);
#endif
     }
  }

  template<class GEOELE, class GEOELEOP>
    void updateMeasNormalQuadPtFirstDerivOp(const GEOELE& geoele, const GEOELEOP& geoeleop)
    {
#ifdef TEST_PRE
      _hasMeas = true;
      _hasTangent = true;
      _hasNormal = true;
      _hasQuadPtCoor = true;
      _hasFirstDeriv = true;
      _hasMassOp = false;
      _hasre = false;
#endif
      _currentId = geoele.id();
      faceIDOp = geoele.pos_second() - 1;

      // update the definition of the geo points 
      for(int i = 0; i < nbGeoNode; i++){
    point(i,0) = geoele.point(i+1).x();
    point(i,1) = geoele.point(i+1).y();
    point(i,2) = geoele.point(i+1).z();
      }

      // update the definition of the geo points of the adjacent element
      for(int i = 0; i < nbGeoNodeAd; i++){
    pointOp(i, 0) = geoeleop.point(i+1).x();
    pointOp(i, 1) = geoeleop.point(i+1).y();
    pointOp(i, 2) = geoeleop.point(i+1).z();
      }

      // Compute the corresponding nodes on the opposite element
      _comp_corresp_node_op();

      // Opposite element basis functions on face quadrature points
      _comp_phi_dphi_op();

      // compute the jacobian on quadrature points
      _comp_inv_jacobian_op();

      // product tInvJacOp by dPhiRefAd
      Real x;
      for(int ig = 0; ig < nbQuadPt; ig++){
    for(int j = 0; j < nbNodeAd; j++){
      for(int icoor = 0; icoor < nbCoorAd; icoor++){
        x = 0.;
        for(int jcoor = 0; jcoor < nbCoorAd; jcoor++){
          x += tInvJacOp(icoor, jcoor, ig) * dPhiOp(j, jcoor, ig);
        } // for jcoor
        phiDerOp(j, icoor, ig) = x;
      } // for icoor
    } // for j
      } // for ig
    }

/*   template<class GEOELE, class GEOELEOP, class CURRFEOP> */
/*     void updateMeasNormalQuadPtFirstDerivMassOp(const GEOELE& geoele, const GEOELEOP& geoeleop, CURRFEOP& currfeop) */
/*     { */
/* #ifdef TEST_PRE */
/*       _hasMeas = true; */
/*       _hasTangent = true; */
/*       _hasNormal = true; */
/*       _hasQuadPtCoor = true; */
/*       _hasFirstDeriv = true; */
/*       _hasMassOp = true; */
/*       _hasre = false; */
/* #endif */
/*       _currentId = geoele.id(); */
/*       faceIDOp = geoele.pos_second() - 1; */

/*       // update the definition of the geo points  */
/*       for(int i = 0; i < nbGeoNode; i++){ */
/*     point(i,0) = geoele.point(i+1).x(); */
/*     point(i,1) = geoele.point(i+1).y(); */
/*     point(i,2) = geoele.point(i+1).z(); */
/*       } */

/*       // compute the measure and the normal */
/*       _comp_meas_normal(); */

/*       // compute the coordinates of the quad points */
/*       _comp_quad_point_coor(); */

/*       // update the definition of the geo points of the adjacent element */
/*       for(int i = 0; i < nbGeoNodeAd; i++){ */
/*     pointOp(i, 0) = geoeleop.point(i+1).x(); */
/*     pointOp(i, 1) = geoeleop.point(i+1).y(); */
/*     pointOp(i, 2) = geoeleop.point(i+1).z(); */
/*       } */

/*       // compute the jacobian on quadrature points */
/*       _comp_inv_jacobian_op(faceIDOp); */

/*       // product tInvJacOp by dPhiRefAd */
/*       Real x; */
/*       for(int ig = 0; ig < nbQuadPt; ig++){ */
/*     for(int j = 0; j < nbNodeAd; j++){ */
/*       for(int icoor = 0; icoor < nbCoorAd; icoor++){ */
/*         x = 0.; */
/*         for(int jcoor = 0; jcoor < nbCoorAd; jcoor++){ */
/*           x += tInvJacOp(icoor, jcoor, ig) * dPhiRefAd(faceIDOp, j, jcoor, ig); */
/*         } // for jcoor */
/*         phiDerOp(j, icoor, ig) = x; */
/*       } // for icoor */
/*     } // for j */
/*       } // for ig */

/*       // update massOp and invMassOp arrays */
/*       currfeop.updateJacMass(geoeleop); */

/*       massOp = currfeop.mass; */
/*       invMassOp = currfeop.invMass; */
/*     } */
/*   // updatere() assumes that massAd, invMassAd, massOp, invMassOp  */
/*   //have already been updated */
/*   void updatere() */
/*     { */
/* #ifdef TEST_PRE */
/*       ASSERT_PRE((hasMassAd()&&hasMassOp()), "Local mass matrices of adjacent elements must be updated") */
/*       _hasre = true; */
/* #endif */
/*       _comp_re(); */
/*     } */
};
}
#endif
