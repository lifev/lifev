#ifndef __ud_functions_H
#define __ud_functions_H 1

#include <cstdlib>

#include "lifeV.hpp"


Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
 switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break; 
 default:
   ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
   break;
  }   
}

Real fZero(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break;  
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }  
}

Real g1(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break; 
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break; 
  }  
}

// Initial displacement and velocity 
Real d0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
  switch(i) {
  case 1:
        return -z*(z-10)*x/50;
    return 0;
    break;
  case 2:
       return -z*(z-10)*y/50;
    return 0;
    break;
  case 3:
    return 0.0;
    break;  
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }  
}

Real w0(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i)
{
 
  switch(i) {
  case 1:
    return 0.0;
    break;
  case 2:
    return 0.0;
    break;
  case 3:
    return 0.0;
    break;   
  default:
    ERROR_MSG("This entrie is not allowed: ud_functions.hpp");
    break;
  }  
}


#endif /* __ud_functions_H */
