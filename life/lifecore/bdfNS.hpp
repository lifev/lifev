/*!
  \file bdfNS.h
  \author A. Veneziani
  \date 04/2003 
  \version 1.0

  \brief File containing a class for an easy handling of different order time discretizations/extrapolations
         BDF based SPECIFIC FOR THE NAVIER-STOKES PROBLEM 
         The idea is to couple a BDF of order q with a pressure incremental approach of order q-1 (see Van Kan, Prohl, Guermond, ecc.)
         If q=1, we still have an incremental pressure treatment (see Guermond)
         REM: At the moment, the couple BDF of order q + extrapolation of order q seems unstable (why ?)
*/
#ifndef _BDF_NS_H
#define _BDF_NS_H
#include <string>
#include <iostream>
#include <algorithm>
#include "GetPot.hpp"
#include "lifeV.hpp"
#include "vecUnknown.hpp"
#include "bdf.hpp"


class BdfNS
{
 public:
  // ! Constructor
  BdfNS(const UInt n);

  ~BdfNS();

  Bdf& bdf_u();
  Bdf& bdf_p();


 private:
  Bdf _bdf_u;
  Bdf _bdf_p;
};



#endif
