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
  @brief Assert macros for LifeV

  @date 19-02-2005
  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef LIFEASSERT_HPP
#define LIFEASSERT_HPP 1

/*!  \page macros_page LifeV Macros
  \section assert_macros Assertion Macros

  \subsection assertions Assertions

  LifeV defines a few macros to test pre and post conditions. To
  enable them all you define the preprocessing variable \c
  LIFEV_CHECK_ALL . This is done easily by configuring LifeV with the
  flag \c --enable-debug which will include the option
  \c -DLIFEV_CHECK_ALL

  -# #LIFEV_ASSERT

  \subsection hints LifeV C++ Compiler Hints

  -# #INLINE
  -# #LIFEV_RESTRICT

  \subsection attribute_macro LifeV Attribute Macros

  -# #LIFEV_EXPORT and #LIFEV_NO_EXPORT
  -# #LIFEV_PACKED
  -# #LIFEV_DEPRECATED
  -# #LIFEV_ISLIKELY and #LIFEV_ISUNLIKELY
*/

#include <life/lifecore/SmartAssert.hpp>

#define ABORT() std::abort()

# define ERROR_MSG(A)  \
   do { std::cerr << std::endl << std::endl << A << std::endl << std::endl ; ABORT() ; } while (0)

# define ASSERT0(X,A) if ( !(X) ) \
ERROR_MSG(A << std::endl << "Error in file" << __FILE__ << " line " << __LINE__) ;

# define ASSERT_PRE0(X,A) if ( !(X) ) \
ERROR_MSG(A << std::endl << "Precondition Error " << "in file " << __FILE__ \
     << " line " << __LINE__) ;

# define ASSERT_POS0(X,A) if ( !(X) ) \
ERROR_MSG(A << std::endl <<"Postcondition Error " << "in file " << __FILE__ \
     << " line " << __LINE__) ;

# define ASSERT_INV0(X,A)  if ( !(X) ) \
ERROR_MSG(A <<std::endl <<  "Invariant Error " << "in file " << __FILE__  \
   << " line " << __LINE__) ;

# define ASSERT_BD0(X)  if ( !(X) ) \
ERROR_MSG("Array bound error " << "in file " << __FILE__  \
   << " line " << __LINE__) ;

#ifdef  LIFEV_CHECK_ALL
#define CHECK_KN
#define TEST_PRE
#define TEST_POS
#define TEST_INV
#define TEST_BOUNDS
#define NOINLINE
#undef  NDEBUG
#endif /* LIFEV_CHECK_ALL */

#ifdef NDEBUG
#define ASSERT(X,A)
#else
#define ASSERT(X,A) ASSERT0(X,A)
#endif

#ifdef TEST_PRE
#define ASSERT_PRE(X,A) ASSERT_PRE0(X,A)
#else
#define ASSERT_PRE(X,A)
#endif

#ifdef TEST_POS
#define ASSERT_POS(X,A) ASSERT_POS0(X,A)
#else
#define ASSERT_POS(X,A)
#endif

#ifdef TEST_INV
#define ASSERT_INV(X,A) ASSERT_INV0(X,A)
#else
#define ASSERT_INV(X,A)
#endif

#ifdef TEST_BOUNDS
#define ASSERT_BD(X) ASSERT_BD0(X)
#else
#define ASSERT_BD(X)
#endif

#endif /* LIFEASSERT_HPP */

