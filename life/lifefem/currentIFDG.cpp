#include <life/lifefem/currentIFDG.hpp>

namespace LifeV
{
CurrentIFDG::CurrentIFDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr):
  CurrentBdDG(_refFEDG, _refFE, _geoMap, _geoMapAd, _qr),
  pointOp(_geoMapAd.nbDof, _geoMapAd.nbCoor),
  refCoorFaceOp(_geoMap.nbDof, _geoMapAd.nbCoor),
  locQuadPtCoor(_qr.nbQuadPt, _geoMapAd.nbCoor),
  phiOp(_refFEDG.nbDof, _qr.nbQuadPt),
  dPhiOp(_refFEDG.nbDof, _geoMapAd.nbCoor, _qr.nbQuadPt),
  phiGeoOp(_geoMapAd.nbDof, _qr.nbQuadPt),
  dPhiGeoOp(_geoMapAd.nbDof, _geoMapAd.nbCoor, _qr.nbQuadPt),
  jacobianOp(_refFEDG.nbCoor, _refFEDG.nbCoor, _qr.nbQuadPt),
  tInvJacOp(_refFEDG.nbCoor, _refFEDG.nbCoor, _qr.nbQuadPt),
  detJacOp(_qr.nbQuadPt),
  weightDetOp(_qr.nbQuadPt),
  phiDerOp(_refFEDG.nbDof, _refFEDG.nbCoor, _qr.nbQuadPt),
//   massOp(_refFEDG.nbDof, _refFEDG.nbDof),
//   invMassOp(_refFEDG.nbDof, _refFEDG.nbDof),
//   re(2 * nbNodeAd, 2 * nbNodeAd, nbCoorAd),
//   massAdOp(2 * nbNodeAd, 2 * nbNodeAd),
//   invMassAdOp(2 * nbNodeAd, 2 * nbNodeAd),
//   phiAdOp(2 * nbNodeAd, _qr.nbQuadPt)
  correspNodeOp(_geoMap.nbDof)
{
  CONSTRUCTOR("CurrentIFDG");
}

CurrentIFDG::~CurrentIFDG()
{
  DESTRUCTOR("CurrentIFDG");
}

// Find the corresponding nodes on the opposite element

void CurrentIFDG::_comp_corresp_node_op() {
  bool same_point;
  Real xi, eta, zeta, s, t;

  for(int i = 0; i < nbGeoNode; i++){
    for(int j = 0; j < nbGeoNodeAd; j++){
 
      for(int icoor = 0; icoor < nbCoorAd; icoor++){
    same_point = false;

    if(point(i, icoor) == pointOp(j, icoor)) same_point = true;

    if(!same_point) break;
      }

      if(same_point){
    correspNodeOp(i) = j;
    break;
      }
    }
  }

  for(int i = 0; i < nbGeoNode; i++){
    for(int icoor = 0; icoor < nbCoorAd; icoor++){
      refCoorFaceOp(i, icoor) = refFEDG.refCoor(correspNodeOp(i), icoor);
    }
  }

  for(int i = 0; i < nbQuadPt; i++){
    s = qr.quadPointCoor(i, 0);
    t = qr.quadPointCoor(i, 1);

    FaceToOpElCoor(xi, eta, zeta, s, t);

    locQuadPtCoor(i, 0) = xi;
    locQuadPtCoor(i, 1) = eta;
    locQuadPtCoor(i, 2) = zeta;
  } 

}

// Compute the value of opposite element basis functions on face quadrature nodes
void CurrentIFDG::_comp_phi_dphi_op(){
  Real xi, eta, zeta;

  for(int i = 0; i < nbGeoNodeAd; i++){
    for(int ig = 0; ig < nbQuadPt; ig++){
      xi = locQuadPtCoor(ig, 0);
      eta = locQuadPtCoor(ig, 1);
      zeta = locQuadPtCoor(ig, 2);

      phiGeoOp(i, ig) = geoMapAd.phi(i, xi, eta, zeta);

      for(int icoor = 0; icoor < nbCoorAd; icoor++){
    dPhiGeoOp(i, icoor, ig) = geoMapAd.dPhi(i, icoor, xi, eta, zeta);
      }

    }
  }

  for(int i = 0; i < nbNodeAd; i++){
    for(int ig = 0; ig < nbQuadPt; ig++){
      xi = locQuadPtCoor(ig, 0);
      eta = locQuadPtCoor(ig, 1);
      zeta = locQuadPtCoor(ig, 2);

      phiOp(i, ig) = refFEDG.phi(i, xi, eta, zeta);

      for(int icoor = 0; icoor < nbCoorAd; icoor++){
    dPhiOp(i, icoor, ig) = refFEDG.dPhi(i, icoor, xi, eta, zeta);
      }

    }
  }
}

void CurrentIFDG::_comp_inv_jacobian_op()
{
  Real fctDer;

  // Derivatives of geo map basis functions of the opposite element
  for(int ig = 0; ig < nbQuadPt; ig++){
    for(int icoor = 0; icoor < nbCoorAd; icoor++){
      for(int jcoor = 0; jcoor < nbCoorAd; jcoor++){
    fctDer = 0.;

    for(int j = 0; j < nbGeoNodeAd; j++){
      fctDer += pointOp(j , icoor) * dPhiGeoOp(j, jcoor, ig);
    } // for j

    jacobianOp(icoor, jcoor, ig) = fctDer;
      } // for jcoor
    } // for icoor
  } // for ig

  // Determinant on face quad points and inverse transposed jacobian
#if defined(TWODIM)
  // *** 2D code ***
  for(int ig = 0; ig < nbQuadPt; ig++){
  a = jacobianOp(0, 0, ig);
  b = jacobianOp(0, 1, ig);
  c = jacobianOp(1, 0, ig);
  d = jacobianOp(1, 1, ig);
  det = a * d - b * c;
  weightDetOp(ig) = detJacOp(ig) * qr.weight(ig);
  tInvJacOp(0, 0, ig) = d / det;
  tInvJacOp(0, 1, ig) = - c / det;
  tInvJacOp(1, 0, ig) = - b / det;
  tInvJacOp(1, 1, ig) = a / det;
  }
#elif defined(THREEDIM)
  // *** 3D code ***
  Real a,b,c,d,e,f,g,h,i,ei,fh,bi,ch,bf,ce,det;
  for(int ig = 0; ig < nbQuadPt; ig++){
    a = jacobianOp(0,0,ig);
    b = jacobianOp(0,1,ig);
    c = jacobianOp(0,2,ig);
    d = jacobianOp(1,0,ig);
    e = jacobianOp(1,1,ig);
    f = jacobianOp(1,2,ig);
    g = jacobianOp(2,0,ig);
    h = jacobianOp(2,1,ig);
    i = jacobianOp(2,2,ig);
    ei=e*i;
    fh=f*h;
    bi=b*i;
    ch=c*h;
    bf=b*f;
    ce=c*e;
    det = a*(ei-fh) + d*(ch-bi) + g*(bf-ce);
    detJacOp(ig) = det;
    weightDetOp(ig) = detJacOp(ig) * qr.weight(ig);
    tInvJacOp(0,0,ig) = (  ei - fh )/det ;
    tInvJacOp(0,1,ig) = (-d*i + f*g)/det ;
    tInvJacOp(0,2,ig) = ( d*h - e*g)/det ;
      
    tInvJacOp(1,0,ig) = ( -bi + ch )/det ;
    tInvJacOp(1,1,ig) = ( a*i - c*g)/det ;
    tInvJacOp(1,2,ig) = (-a*h + b*g)/det ;
      
    tInvJacOp(2,0,ig) = (  bf - ce )/det ;
    tInvJacOp(2,1,ig) = (-a*f + c*d)/det ;
    tInvJacOp(2,2,ig) = ( a*e - b*d)/det ;
  }
#endif
}

// void CurrentIFDG::_comp_re()
// {
//   KN<Real> rhs(2 * nbNodeAd);

//   for(int i = 0; i < nbNodeAd; i++){
//     for(int j = 0; j < nbNodeAd; j++){
//       massAdOp(i, j) = massAd(i, j);
//       massAdOp(i + nbNodeAd, j + nbNodeAd) = massOp(i, j);

//       invMassAdOp(i, j) = invMassAd(i, j);
//       invMassAdOp(i + nbNodeAd, j + nbNodeAd) = invMassOp(i, j);
//     }
//   }

//   for(int i = 0; i < nbNodeAd; i++){
//     for(int ig = 0; ig < nbQuadPt; ig++){
//       phiAdOp(i, ig) = phiAd(faceIDAd, i, ig);
//       phiAdOp(i + nbNodeAd, ig) = phiAd(faceIDOp, i, ig);
//     }
//   }

//   Real n;

//     for(int j = 0; j < 2 * nbNodeAd; j++){
//       for(int icoor = 0; icoor < nbCoorAd; icoor++){
//     rhs = 0.;

//     // Compute right hand side vector
//       for(int i = 0; i < 2 * nbNodeAd; i++){
//         for(int ig = 0; ig < nbQuadPt; ig++){
//           if(j < nbNodeAd){
//         n = normal(icoor, ig);
//           }else{
//         n = - normal(icoor, ig);
//           }
//           rhs(i) += - phiAdOp(j, ig) * n * phiAdOp(i, ig) * weightMeas(ig);
//         } // for ig
//       } // for i

//     // Solve the system
//     rhs = invMassAdOp * rhs;

//     // Store the solution
//     for(int i = 0; i < 2 * nbNodeAd; i++){
//       re(i, j, icoor) = rhs(i);
//     }

//       } // for icoor
//     } // for j
// }
}
