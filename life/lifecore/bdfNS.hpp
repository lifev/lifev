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
