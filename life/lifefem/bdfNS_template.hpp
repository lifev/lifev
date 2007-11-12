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
#ifndef _BDF_NS_TEMPLATE_H
#define _BDF_NS_TEMPLATE_H
#include <string>
#include <iostream>
#include <algorithm>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifearray/vecUnknown.hpp>
#include <life/lifefem/bdf_template.hpp>

namespace LifeV
{
template<typename VectorType = EpetraVector<double> >
class BdfTNS
{
public:
    // ! Constructor
    BdfTNS( const UInt n );

    ~BdfTNS() {}

    inline BdfT<VectorType>& bdf_u() { return _bdf_u; }
    inline BdfT<VectorType>& bdf_p() { return _bdf_p; }



private:
    BdfT<VectorType> _bdf_u;
    BdfT<VectorType> _bdf_p;
};


template<typename VectorType>
BdfTNS<VectorType>::BdfTNS( const UInt n )
    :
    _bdf_u( n ),
    _bdf_p( std::max( UInt( 1 ), n - 1 ) )
{}


}
#endif
