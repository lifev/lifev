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


