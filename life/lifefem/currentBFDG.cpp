#include "currentBFDG.hpp"

namespace LifeV
{
CurrentBFDG::CurrentBFDG(const RefFEDG& _refFEDG, const RefFE& _refFE, const GeoMap& _geoMap, const GeoMapDG& _geoMapAd, const QuadRule& _qr):
  CurrentBdDG(_refFEDG, _refFE, _geoMap, _geoMapAd, _qr),

  re(_refFEDG.nbDof, _refFEDG.nbDof, _geoMapAd.nbCoor)
{
  CONSTRUCTOR("CurrentBFDG");
}

CurrentBFDG::~CurrentBFDG()
{
  DESTRUCTOR("CurrentBFDG");
}

void CurrentBFDG::_comp_re()
{
  KN<Real> rhs(nbNodeAd);

  for(int j = 0; j < nbNodeAd; j++){
    for(int icoor = 0; icoor < nbCoorAd; icoor++){
      rhs = 0.;

      // Compute right hand side vector
      for(int i = 0; i < nbNodeAd; i++){
	for(int ig = 0; ig < nbQuadPt; ig++){
	  rhs(i) += - (phiAd(j, ig) * normal(icoor, ig)) * phiAd(i, ig) * weightMeas(ig);
	} // for ig
      } // for i

      // Solve the system
      rhs = invMassAd * rhs;

      // Store the solution
      for(int index = 0; index < nbNodeAd; index++){
	re(index, j, icoor) = rhs(index);
      }

    } // for icoor
  } // for j
}
}
