/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  
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
#include "dataNewton.hpp"

// Constructor
DataNewton::DataNewton(const GetPot& dfile, const std::string& section)
{
  _maxiter      =  dfile((section+"/maxiter").data(),100);
  _abstol       =  dfile((section+"/abstol").data(),0.0);
  _reltol       =  dfile((section+"/reltol").data(),0.0); 
  _etamax       =  dfile((section+"/etamax").data(),1.e-3);
  _linesearch   =  dfile((section+"/linesearch").data(),2);
}

// Destructor
DataNewton::~DataNewton() {}

 
// The max number of interations
UInt DataNewton::maxiter() const {
  return _maxiter;
}

// The absolute tolerance
Real DataNewton::abstol() const {
  return _abstol;
}

// The relative tolerance
Real DataNewton::reltol() const {
  return _reltol;
}

// The relative tolerance
Real DataNewton::etamax() const {
  return _etamax;
}

// The linesearch option
UInt DataNewton::linesearch() const {
  return _linesearch;
}

// Output
void DataNewton::showMe( std::ostream& c) const
{
  // 
  c << "maxiter        = " << _maxiter << std::endl; 
  c << "abstol         = " << _abstol << std::endl; 
  c << "reltol         = " << _reltol << std::endl; 
  c << "etamax         = " << _reltol << std::endl; 
  c << "linesearch     = " << _linesearch << std::endl; 
}

