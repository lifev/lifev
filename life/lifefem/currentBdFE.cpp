#include "currentBdFE.hpp"

CurrentBdFE::CurrentBdFE(const RefFE& _refFE,const GeoMap& _geoMap,const QuadRule& _qr):
  StaticBdFE(_refFE,_geoMap,_qr)
{
  CONSTRUCTOR("CurrentBdFE");
}

CurrentBdFE::CurrentBdFE(const RefFE& _refFE,const GeoMap& _geoMap):
  StaticBdFE(_refFE,_geoMap)
{
  CONSTRUCTOR("CurrentBdFE (without quadrature rule)");
}

CurrentBdFE::~CurrentBdFE()
{
  DESTRUCTOR("CurrentBdFE")
}


