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
    @brief LifeV main header file


    @author Luca Formaggia
    @author Jean-Fred Gerbeau.
    @author Christophe Prud'homme

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>

    @date 01-01-2004, 01-12-2010

    Defines typedefs and macros common to ALL lifeV.h software
    it must be includes in all translation units.
*/


/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 */

#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include <lifeconfig.h>

#ifndef __cplusplus
#error You must use C++ for LifeV
#endif

#ifndef _LIFEV_HH_
#define _LIFEV_HH_

#include <stdint.h>

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numeric>
#include <iosfwd>
#include <string>

#include <life/lifecore/lifeassert.hpp>

namespace LifeV
{

/*!
  @page types_page LifeV Types
  @section types Types
  @subsection real Real Numbers

  LifeV defines a number of types that are used in the library.

  -# \c Real 64 bits real number type

  @section ints Integers

  LifeV defines a number of integer type that have controlled bit
  size. These types are constructed automatically by LifeV in order to have
  platform independant integer types.

  Here is the list of signed integers:

  -# \c int1_type  a 1 bit signed integer
  -# \c int8_type a 8 bit signed integer
  -# \c int16_type a 16 bit signed integer
  -# \c int32_type a 32 bit signed integer
  -# \c int64_type a 64 bit signed integer

  Here is the list of unsigned integers:

  -# \c uint1_type a 1 bit unsigned integer
  -# \c uint8_type a 8 bit unsigned integer
  -# \c uint16_type a 16 bit unsigned integer
  -# \c uint32_type a 32 bit unsigned integer
  -# \c uint64_type a 64 bit unsigned integer

  LifeV defines a number of useful aliases for integers
  -# \c Int an alias to int32_type
  -# \c UInt an alias to uint32_type used for adressing
  -# \c ID an alias to id_type used to identify local numbering or components
  -# \c size_type an alias to size_t used as indices for arrays, vectors or matrices

*/

// Create type that are machine independent

//! Generic real data
typedef double Real;

typedef int8_t  int8_type;
typedef int16_t int16_type;
typedef int32_t int32_type;
typedef int64_t int64_type;

typedef uint8_t  uint8_type;
typedef uint16_t uint16_type;
typedef uint32_t uint32_type;
typedef uint64_t uint64_type;

//! Generic integer data
typedef int32_type  Int;

//! generic unsigned integer (used mainly for addressing)
typedef uint32_type UInt;

//! IDs (which starts ALWAYS from 1)
typedef uint32_type ID;

// For now only 3 dimensional problems.
extern const UInt nDimensions;
#define NDIM 3

} // end namespace LifeV

#endif
