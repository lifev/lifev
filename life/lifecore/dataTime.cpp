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

