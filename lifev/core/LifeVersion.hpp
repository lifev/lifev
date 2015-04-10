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

#ifndef _LIFEV_VERSION_H_
#define _LIFEV_VERSION_H_

//#include <lifev/core/LifeV.hpp>

#define LIFEV_MAKE_VERSION( a,b,c ) (((a) << 16) | ((b) << 8) | (c))

#if !defined(LIFEV_VERSION_MAJOR)
#define LIFEV_VERSION_MAJOR 3
#endif

#if !defined(LIFEV_VERSION_MINOR)
#define LIFEV_VERSION_MINOR 8
#endif

#if !defined(LIFEV_VERSION_MICRO)
#define LIFEV_VERSION_MICRO 8
#endif

#if !defined(LIFEV_VERSION)
#define LIFEV_VERSION \
  LIFEV_MAKE_VERSION(LIFEV_VERSION_MAJOR,LIFEV_VERSION_MINOR,LIFEV_VERSION_MICRO)
#endif

#define LIFEV_IS_VERSION(a,b,c) ( LIFEV_VERSION >= LIFEV_MAKE_VERSION(a,b,c) )

/**
 * Namespace for general LIFEV functions.
 */
namespace LifeV
{
/**
 * Returns the encoded number of LIFEV's version, see the LIFEV_VERSION macro.
 * In contrary to that macro this function returns the number of the actually
 * installed LIFEV version, not the number of the LIFEV version that was
 * installed when the program was compiled.
 * @return the version number, encoded in a single uint
 * @since 0.7
 */
unsigned int version();

/**
 * Returns the major number of LIFEV's version, e.g.
 * 0 for LIFEV 0.7
 * @return the major version number
 * @since 0.7
 */
unsigned int versionMajor();

/**
 * Returns the minor number of LIFEV's version, e.g.
 * 7 for LIFEV 0.7.0
 * @return the minor version number
 * @since 0.7
 */
unsigned int versionMinor();

/**
 * Returns the micro number of LIFEV's version, e.g.
 * 0 for LIFEV 0.7.0
 * @return the extra information
 * @since 0.7
 */
unsigned int versionMicro();

/**
 * Returns the LIFEV version as string, e.g. "0.7.0".
 * @return the LIFEV version. You can keep the string forever
 * @since 0.7
 */
char const* versionString();
}

#endif // _LIFEV_VERSION_H_
