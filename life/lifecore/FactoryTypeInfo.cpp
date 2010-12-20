//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

This file is part of LifeV.

LifeV is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LifeV is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER
/*!
  @file
  @brief Type information class

  @date 13-10-2004
  @author  Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#include <cassert>

#include <life/lifecore/life.hpp>
#include <life/lifecore/FactoryTypeInfo.hpp>

namespace LifeV
{

// ===============================
// Constructors and Destructor
// ===============================

FactoryTypeInfo::FactoryTypeInfo()
{
    class Nil {};
    M_info = &typeid(Nil);
    assert( M_info != 0 );

}

FactoryTypeInfo::FactoryTypeInfo(const std::type_info& ti):
    M_info(&ti)
{
    assert( M_info != 0 );
}

FactoryTypeInfo::FactoryTypeInfo( FactoryTypeInfo const& ti ):
    M_info( ti.M_info )
{
    assert( M_info != 0 );
}

FactoryTypeInfo::~FactoryTypeInfo()
{
}

// ===============================
// Public methods
// ===============================

bool
FactoryTypeInfo::before(const FactoryTypeInfo& rhs) const
{
    assert( M_info != 0 );
    return M_info->before(*rhs.M_info);
}

// ======================
// Get methods
// ======================

const std::type_info&
FactoryTypeInfo::typeInfo() const
{
    assert( M_info != 0 );
    return *M_info;
}

const char*
FactoryTypeInfo::name() const
{
    assert( M_info != 0 );
    return M_info->name();
}

}
