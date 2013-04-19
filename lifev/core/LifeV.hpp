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

#include <lifev/core/Core_config.h>
#include <LifeV_config.h>

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
#include <limits>
#include <set>
#include <list>
#include <map>
#include <vector>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/function.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// deprecated attribute for LifeV functions
// the macro is needed to avoid problems with compilers other than gcc
// other compiler specific implementation of the deprecated attribute can be
// added with #elif defined macros
#ifdef __GNUC__
#define LIFEV_DEPRECATED( func ) func __attribute__ ((deprecated))
#else
#define LIFEV_DEPRECATED( func ) func
#endif

// macro to avoid warning in opt mode for variables that are only needed for
// dbg mode
// note: these problems arise when there is no clear separation between a
// function that does some work and return some state, so it should be
// avoided as much as possible.
#define LIFEV_UNUSED(x) ((void)x)

#include <lifev/core/util/LifeAssert.hpp>

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
  -# \c flag_Type an alias to uint32_type used for boolean flags

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

//! IDs
typedef uint32_type ID;

//! bit-flag with up to 32 different flags
typedef uint32_type flag_Type;

//! flag related free functions and functors
namespace Flag
{
/** @name FlagFunctions
 *  They implement basic operations on boolean flags
 */
//@{
//! It returns true if all bit-flags common set in refFlag are also set in inputFlag
inline bool testAllSet ( flag_Type const& inputFlag, flag_Type const& refFlag )
{
    return ( inputFlag  & refFlag ) == refFlag;
}

//! returns true if at least one flag set in refFlag is set in inputFlag
inline bool testOneSet ( flag_Type const& inputFlag, flag_Type const& refFlag )
{
    return inputFlag  & refFlag;
}

//! returns false if at least one flag set in refFlag is set in inputFlag
inline bool testOneNotSet ( flag_Type const& inputFlag, flag_Type const& refFlag )
{
    return ! (inputFlag & refFlag);
}

//! turns on the refFlag active bits in inputFlag
inline flag_Type turnOn  ( flag_Type const& inputFlag, flag_Type const& refFlag )
{
    return inputFlag  | refFlag;
}

//! turns off the refFlag active bits in inputFlag
inline flag_Type turnOff ( flag_Type const& inputFlag, flag_Type const& refFlag )
{
    return inputFlag  & ~refFlag;
}

//! switches the refFlag active bits in inputFlag
inline flag_Type change ( flag_Type const& inputFlag, flag_Type const& refFlag )
{
    return inputFlag  ^ refFlag;
}

//! replaces the given flag with the reference one. This method is introduced with the same
//! signature of the other methods in order to be used as a policy
inline flag_Type replaceFlag  ( flag_Type const& /*inputFlag*/, flag_Type const& refFlag )
{
    return refFlag;
}

//! showMe method to print out flag status
//! the flag is converted to its binary form ( right -> left corresponds to first -> last flag )
void showMe ( flag_Type const& flag, std::ostream& out = std::cout );
//@}

//end namespace Flag
}
// For now only 3 dimensional problems.
extern const UInt nDimensions;

// used to denote an ID not set. We use max of Int (instead of UInt) since
// MPI does not have the UInt type.
const ID NotAnId = std::numeric_limits<Int>::max();
#define NDIM 3

//! clearVector
/*!
 * This is a general purpose utility that clears up a std::vector<T>
 * making sure that it does not uses up memory after the call
 * Useful when memory is an issue, since clear() does not free memory
 * */
template<typename T>
void clearVector (T& stdVector)
{
    stdVector.clear();
    T().swap (stdVector);
}
//! resizeVector
/*!
 * This is a general purpose utility that resizes up a std::vector<T>
 * making sure that it does not uses up more memory after the call
 * Useful when memory is an issue, since resize() does not free memory
 */
template<typename T>
void resizeVector (T& stdVector, UInt const& newsize)
{
    stdVector.resize (newsize);
    if (stdVector.capacity() > stdVector.size() )
    {
        T (stdVector).swap (stdVector);
    }
}

} // end namespace LifeV

#endif
