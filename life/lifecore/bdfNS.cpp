#include "bdfNS.hpp"

BdfNS::BdfNS(const UInt n):_bdf_u(n),_bdf_p(max((UInt)1,n-1))
{}

Bdf& BdfNS::bdf_u() 
{
  return _bdf_u;
}

Bdf& BdfNS::bdf_p() 
{
  return _bdf_p;
}


BdfNS::~BdfNS() {};

