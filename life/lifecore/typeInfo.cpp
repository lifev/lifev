/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-10-13

  Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
/**
   \file typeInfo.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-10-13
 */
#include <cassert>

#include <life/lifecore/life.hpp>
#include <life/lifecore/typeInfo.hpp>

namespace LifeV
{

// ===============================
// Constructors and Destructor
// ===============================

TypeInfo::TypeInfo()
{
    class Nil {};
    _M_info = &typeid(Nil);
    assert( _M_info != 0 );

}

TypeInfo::TypeInfo(const std::type_info& ti)
        :
        _M_info(&ti)
{
    assert( _M_info != 0 );
}

TypeInfo::TypeInfo( TypeInfo const& ti )
        :
        _M_info( ti._M_info )
{
    assert( _M_info != 0 );
}

TypeInfo::~TypeInfo()
{
}

// ===============================
// Public methods
// ===============================

bool
TypeInfo::before(const TypeInfo& rhs) const
{
    assert( _M_info != 0 );
    return _M_info->before(*rhs._M_info);
}

// ======================
// Get methods
// ======================

const std::type_info&
TypeInfo::typeInfo() const
{
    assert( _M_info != 0 );
    return *_M_info;
}

const char*
TypeInfo::name() const
{
    assert( _M_info != 0 );
    return _M_info->name();
}

}
