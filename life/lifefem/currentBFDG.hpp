#ifndef _CURRENTBFDG_H
#define _CURRENTBFDG_H

#include "currentBdDG.hpp"
#include "bcCond.hpp"

namespace LifeV
{
class CurrentBFDG:public CurrentBdDG{
#ifdef TEST_PRE
 protected:
  bool _hasre;
#endif

 private:
  void _comp_re();
 public:
  BCType bcType;
  EntityFlag marker;

  KNMK<Real> re; //!< Contains the i-th coefficient of the piecewise polynomial approximation of the k-th component of r_e(jump(fi_j)) on the adjacent element.

 public:
  CurrentBFDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr);
  ~CurrentBFDG();

#ifdef TEST_PRE
  inline bool hasre() const {return _hasre;}
#endif

  template<class GEOELE>
    void updateBCType(const GEOELE& geoele, const BC_Handler& BCh){
    marker = EntityFlag(geoele.marker());
    bcType = BCh.boundaryType(marker);
  }
 //----------------------------------------------------------------------------
  /*!
    Compute the arrays meas, weightMeas, tangent,
    normal, quadrature points, phiDer, mass and invMass
    on the current boundary element
  */
  template<class GEOELE, class GEOELEAD, class CURRFEAD>
    void updateMeasNormalQuadPtFirstDerivMassreAd(const GEOELE& geoele, const GEOELEAD& geoelead, CURRFEAD& currfead)
    {
#ifdef TEST_PRE
      _hasre = true;
#endif
      updateMeasNormalQuadPtFirstDerivMassAd(geoele, geoelead, currfead);
      _comp_re();
    }
};
}
#endif
