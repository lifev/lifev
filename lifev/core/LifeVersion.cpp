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
    @brief  The file was created from KDE/kdelibs/kdecore/kdeversion.hpp and
    accomodated to LifeV needs.

    @date 2005-02-19
    @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
    @mantainer Simone Deparis <simone.deparis@epfl.ch>

 */

#include <lifev/core/LifeV.hpp>
#include <lifev/core/LifeVersion.hpp>

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

