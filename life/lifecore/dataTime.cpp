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
#include "dataTime.hpp"

// Constructor
DataTime::DataTime(const GetPot& dfile, const std::string& section)
{
  _dt      =  dfile((section+"/timestep").data(),1.); 
  _order_bdf = dfile((section+"/order_bdf").data(),1);
}

// Destructor
DataTime::~DataTime() {}

 
// The time step
Real DataTime::timestep() const {
  return _dt;
}

// Order of the bdf formula
unsigned int DataTime::order_bdf() const {
	return _order_bdf;
}	

// Output
void DataTime::showMe( std::ostream& c ) const
{
  // time step
  c << "timestep  = " << _dt << std::endl; 
  c << "order bdf = " << _order_bdf << std::endl; 
}

