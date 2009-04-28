#include <life/lifefem/currentBdDG.hpp>

namespace LifeV
{
CurrentBdDG::CurrentBdDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr):
  StaticBdFE(_refFE, _geoMap, _qr),
  refFEDG(_refFEDG),
  nbNodeAd(_refFEDG.nbDof),
  nbGeoNodeAd(_geoMapAd.nbDof),
  nbCoorAd(_geoMapAd.nbCoor),
  pointAd(_geoMapAd.nbDof, _refFEDG.nbCoor),
  geoMapAd(_geoMapAd),
  phiAd(_refFEDG.nbDof, _qr.nbQuadPt),
  dPhiRefAd(_refFEDG.nbDof, _refFEDG.nbCoor, _qr.nbQuadPt),
  jacobianAd(_refFEDG.nbCoor, _refFEDG.nbCoor, _qr.nbQuadPt),
  tInvJacAd(_refFEDG.nbCoor, _refFEDG.nbCoor, _qr.nbQuadPt),
  detJacAd(_qr.nbQuadPt),
  weightDetAd(_qr.nbQuadPt),
  phiDerAd(_refFEDG.nbDof, _geoMapAd.nbCoor, _qr.nbQuadPt),
  phiGeoAd(_geoMapAd.nbDof, _qr.nbQuadPt),
  dPhiGeoAd(_geoMapAd.nbDof, _refFEDG.nbCoor, _qr.nbQuadPt),
  massAd(_refFEDG.nbDof, _refFEDG.nbDof),
  invMassAd(_refFEDG.nbDof, _refFEDG.nbDof)
{
  CONSTRUCTOR("CurrentBdDG");

}

//------------------------------------------------------------------------------------

CurrentBdDG::~CurrentBdDG()
{
  DESTRUCTOR("CurrentBdDG");
}

//------------------------------------------------------------------------------------
void CurrentBdDG::_comp_phi_ad_d_phi_face() {
  for(int ig = 0; ig < qr.nbQuadPt; ig++){
    for(int i = 0; i < nbNodeAd; i++){
      phiAd(i, ig) = refFEDG.phiFace(faceIDAd, i, ig, qr);
      for(int icoor = 0; icoor < nbCoorAd; icoor++){
    dPhiRefAd(i, icoor, ig) = refFEDG.dPhiFace(faceIDAd, i, icoor, ig, qr);
      } // for icoor
    } // for i

    for(int k = 0; k < nbGeoNodeAd; k++){
      phiGeoAd(k, ig) = geoMapAd.phiFace(faceIDAd, k, ig, qr);
      for(int icoor = 0; icoor < nbCoorAd; icoor++){
    dPhiGeoAd(k, icoor, ig) = geoMapAd.dPhiFace(faceIDAd, k, icoor, ig, qr);
      } // for icoor
    } // for k

  } // for ig
}

//------------------------------------------------------------------------------------

void CurrentBdDG::_comp_inv_jacobian_ad()
{
  Real fctDer;

  // Derivatives of geo map of the adjacent element
  for(int ig = 0; ig < nbQuadPt; ig++){
    for(int icoor = 0; icoor < nbCoorAd; icoor++){
      for(int jcoor = 0; jcoor < nbCoorAd; jcoor++){
    fctDer = 0.;
    for(int j = 0; j < nbGeoNodeAd; j++){
      fctDer += pointAd(j, icoor) * dPhiGeoAd(j, jcoor, ig);
    } // for j
    jacobianAd(icoor, jcoor, ig) = fctDer;
      } // for jcoor
    } // for icoor
  } // for ig

  // Determinant on face quad points and inverse transposed jacobian
#if defined(TWODIM)
  // *** 2D code ***
  Real a,b,c,d,det;
  for(int ig = 0; ig < nbQuadPt; ig++){
  a = jacobianAd(0, 0, ig);
  b = jacobianAd(0, 1, ig);
  c = jacobianAd(1, 0, ig);
  d = jacobianAd(1, 1, ig);
  det = a * d - b * c;
  weightDetAd(ig) = detJacAd(ig) * qr.weight(ig);
  tInvJacAd(0, 0, ig) = d / det;
  tInvJacAd(0, 1, ig) = - c / det;
  tInvJacAd(1, 0, ig) = - b / det;
  tInvJacAd(1, 1, ig) = a / det;
  }
}
#elif defined(THREEDIM)
  // *** 3D code ***
  Real a,b,c,d,e,f,g,h,i,ei,fh,bi,ch,bf,ce,det;
  for(int ig = 0; ig < nbQuadPt; ig++){
    a = jacobianAd(0,0,ig);
    b = jacobianAd(0,1,ig);
    c = jacobianAd(0,2,ig);
    d = jacobianAd(1,0,ig);
    e = jacobianAd(1,1,ig);
    f = jacobianAd(1,2,ig);
    g = jacobianAd(2,0,ig);
    h = jacobianAd(2,1,ig);
    i = jacobianAd(2,2,ig);
    ei=e*i;
    fh=f*h;
    bi=b*i;
    ch=c*h;
    bf=b*f;
    ce=c*e;
    det = a*(ei-fh) + d*(ch-bi) + g*(bf-ce);
    detJacAd(ig) = det;
    weightDetAd(ig) = detJacAd(ig) * qr.weight(ig);
    tInvJacAd(0,0,ig) = (  ei - fh )/det;
    tInvJacAd(0,1,ig) = (-d*i + f*g)/det;
    tInvJacAd(0,2,ig) = ( d*h - e*g)/det;
      
    tInvJacAd(1,0,ig) = ( -bi + ch )/det;
    tInvJacAd(1,1,ig) = ( a*i - c*g)/det;
    tInvJacAd(1,2,ig) = (-a*h + b*g)/det;
      
    tInvJacAd(2,0,ig) = (  bf - ce )/det;
    tInvJacAd(2,1,ig) = (-a*f + c*d)/det;
    tInvJacAd(2,2,ig) = ( a*e - b*d)/det;
  }
}
#endif
}
