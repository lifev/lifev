/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
/*
  \file lifeV.hpp

  LifeV main header file

  \author Luca Formaggia
  \author Jean-Fred Gerbeau.

  #Purposes Defines typedefs and macros common to ALL lifeV.h software
  it must be includes in all translation units.
*/
# ifndef __cplusplus
# error You must use C++ for LifeV
# endif

# ifndef _LIFEV_HH_
# define _LIFEV_HH_

# include <cstdlib>

# include <iostream>
# include <cmath>
# include <numeric>
# include <iosfwd>
# include <string>

#define ABORT() std::abort()


// debugging macros
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

//
// switch all debugging on if LIFEV_CHECK_ALL is set
//
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

#ifdef NOINLINE
#define INLINE
#else
#define INLINE inline
#endif

// If there is restricted use it!

#ifdef HAS_RESTRICT
#define RESTRICT restrict
#else
#define RESTRICT
#endif
namespace LifeV
{
// Typedefs

//! Real data
typedef double Real;
//! Generic integer data
typedef int Int;
//! generic unsigned integer (used mainly for addressing)
typedef uint  UInt;
typedef unsigned short int USInt;

//! IDs (which starts ALWAYS from 1)
typedef UInt ID;
//! specialised version of ID for Degree of Freedom ID
typedef ID ID_Dof;
//! Indices (starting from 0)
typedef UInt Index_t;

// typedef for indices

#ifdef INT_BCNAME
typedef int  BCName;
//! nullBCName identify a NULL Bundary condition
const BCName nullBCName=0;
#else
typedef std::string BCName;
const BCName nullBCName; // The empty string!
#endif

// Marker is used for storing attributes normally related to boundary condition
// application. At the moment is just an integer

//typedef int Marker;
//const Marker NULLMARKER=0; // 0 is reserved for the NULLMARKER
// global variables
/* -------------------------------------------------------
 I want to know wheter I am working in 2 or 3 Dimensions
 I will set both a constant global variable (nDimensions)
 and a cpp macro (NDIM) to the correct values, for later use.
---------------------------------------------------------*/
#if defined(TWODIM) && defined(THREEDIM)
#error Dont be silly. Either Two or Three dimensions.
#endif
#if defined(TWODIM)
#define NDIM 2
const UInt nDimensions=2;
#elif defined(THREEDIM)
#define NDIM 3
const UInt nDimensions=3;
#else
#error You MUST compile with either -DTWODIM of -DTHREEDIM set, sorry.
#endif

#undef VERSION_2D
# define MSG(A)  \
   do { std::cout << A << std::endl ; } while (0)

// Trace the constructors and destructors
#ifdef TRACE_CONSTRUCTOR
#define CONSTRUCTOR(A) MSG("Constructor of "<<A<<" (file "<<__FILE__<< ", line " << __LINE__ << ")");
#else
#define CONSTRUCTOR(A)
#endif

#ifdef TRACE_DESTRUCTOR
#define DESTRUCTOR(A) MSG("Destructor of "<<A<<" (file "<<__FILE__<< ", line " << __LINE__ << ")");
#else
#define DESTRUCTOR(A)
#endif
} // end namespace LifeV

#endif

