/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
#ifndef _SELECTMARKER_HH_
#define _SELECTMARKER_HH_

#include <life/lifecore/debug.hpp>
#include <life/lifemesh/regionMesh3D.hpp>

namespace LifeV
{
//! \file selectMarker.hpp
//! This file contains the standard selector for internal entities


//! Functor class that tells whether an entity flag corresponds to an internal face
class InternalEntitySelector
{
public:
    //! The default watermark used when standard contructor is adopted
    static const EntityFlag defMarkFlag;
    InternalEntitySelector();
    InternalEntitySelector(const EntityFlag & w);
    /*! operator returning true if flag corresponds to internal entity
      If the EntityFlag is greater that the watermark then it is assumed that
      the associated geometry entity is internal
    */
    bool operator()(EntityFlag const &) const;
private:
    //! The current watermark
    EntityFlag waterMarkFlag;
};
}


#endif
