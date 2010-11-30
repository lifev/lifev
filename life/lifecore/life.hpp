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
/*!
  \file life.hpp

  LifeV main header file

  \author Luca Formaggia
  \author Jean-Fred Gerbeau.
  \author Christophe Prud'homme

  #Purposes Defines typedefs and macros common to ALL lifeV.h software
  it must be includes in all translation units.
*/
# ifndef __cplusplus
# error You must use C++ for LifeV
# endif

# ifndef _LIFEV_HH_
# define _LIFEV_HH_



#include <boost/mpl/multiplies.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/lower_bound.hpp>
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/sizeof.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/base.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/begin_end.hpp>

// only available in boost 1.32
#include <boost/mpl/eval_if.hpp>

# include <cstdlib>

# include <iostream>
# include <cmath>
# include <numeric>
# include <iosfwd>
# include <string>

#include <life/lifecore/lifemacros.hpp>
#include <life/lifecore/lifeversion.hpp>
#include <life/lifecore/lifeassert.hpp>
#include <life/lifecore/debug.hpp>


namespace LifeV
{

/*!  \page types_page LifeV Types
  \section types Types
  \subsection real Real Numbers

  LifeV defines a number of types that are used in the library.

  -# \c Real 64 bits real number type

  \section ints Integers

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
  -# \c SInt an alias to int16_type
  -# \c Int an alias to int32_type
  -# \c UInt an alias to uint32_type used for adressing
  -# \c USInt an alias to uint16_type
  -# \c id_type an alias to uint32_type used to identify local numbering or components
  -# \c ID an alias to id_type used to identify local numbering or components
  -# \c Index_t an alias to id_type used to identify local numbering or components
  -# \c dim_type an alias to uint8_type used to identify dimensions
  -# \c size_type an alias to size_t used as indices for arrays, vectors or matrices

*/

typedef double Real;

//
// Create type that are machine independent
// \warning should test here the boost version
//


/*! \namespace detail
  \internal
*/
namespace detail
{
template<int bit_size>
class no_int
{
private:
    no_int();
};

template< int bit_size >
struct integer
{
    typedef boost::mpl::list<signed char,signed short, signed int, signed long int, signed long long> builtins_;
    typedef typename boost::mpl::base< typename boost::mpl::lower_bound<
            boost::mpl::transform_view< builtins_, boost::mpl::multiplies< boost::mpl::sizeof_<boost::mpl::placeholders::_1>, boost::mpl::int_<8> >
            >
        , boost::mpl::integral_c<size_t, bit_size>
        >::type >::type iter_;

    typedef typename boost::mpl::end<builtins_>::type last_;
    typedef typename boost::mpl::eval_if<
          boost::is_same<iter_,last_>
        , boost::mpl::identity< no_int<bit_size> >
        , boost::mpl::deref<iter_>
        >::type type;
};
}

typedef detail::integer<1>::type  int1_type;

typedef detail::integer<8>::type  int8_type;
typedef detail::integer<16>::type int16_type;
typedef detail::integer<32>::type int32_type;
typedef detail::integer<64>::type int64_type;

/*! \namespace detail
  \internal
*/
namespace detail
{
template< int bit_size >
struct unsigned_integer
{
    typedef boost::mpl::list<unsigned char,unsigned short,unsigned int,unsigned long int, unsigned long long> builtins_;
    typedef typename boost::mpl::base< typename boost::mpl::lower_bound<
        boost::mpl::transform_view< builtins_
      , boost::mpl::multiplies< boost::mpl::sizeof_<boost::mpl::placeholders::_1>, boost::mpl::int_<8> >
    >
        , boost::mpl::integral_c<size_t, bit_size>
    >::type >::type iter_;

    typedef typename boost::mpl::end<builtins_>::type last_;
    typedef typename boost::mpl::eval_if<
        boost::is_same<iter_,last_>
        , boost::mpl::identity< no_int<bit_size> >
        , boost::mpl::deref<iter_>
    >::type type;
};
}
typedef detail::unsigned_integer<1>::type  uint1_type;
typedef detail::unsigned_integer<8>::type  uint8_type;
typedef detail::unsigned_integer<16>::type uint16_type;
typedef detail::unsigned_integer<32>::type uint32_type;
typedef detail::unsigned_integer<64>::type uint64_type;

//! Generic integer data
typedef int32_type  Int;
typedef int16_type  SInt;
//! generic unsigned integer (used mainly for addressing)
typedef uint32_type UInt;
typedef uint16_type USInt;

//! IDs (which starts ALWAYS from 1)
typedef uint32_type id_type;
typedef id_type ID;
typedef id_type Index_t;

//! dimension type
typedef uint8_type dim_type;

//! Indices (starting from 0)
typedef size_t size_type;

// For now only 3 dimensional problems.

extern const UInt nDimensions;
#define NDIM 3


} // end namespace LifeV

#endif

