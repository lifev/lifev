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

