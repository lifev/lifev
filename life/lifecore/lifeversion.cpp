/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-02-19

  Copyright (C) 2005 EPFL

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
   \file lifeversion.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-02-19
 */
#include <life/lifecore/life.hpp>

namespace LifeV
{
unsigned int version()
{
    return LIFEV_VERSION;
}

unsigned int versionMajor()
{
    return LIFEV_VERSION_MAJOR;
}

unsigned int versionMinor()
{
    return LIFEV_VERSION_MINOR;
}

unsigned int versionMicro()
{
    return LIFEV_VERSION_MICRO;
}

char const* versionString()
{
    return LIFEV_VERSION_STRING;
}
}

