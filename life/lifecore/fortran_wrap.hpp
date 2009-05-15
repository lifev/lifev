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
#ifndef FORTRAN_FROM_CPLUSPLUS
#define FORTRAN_FROM_CPLUSPLUS

#include <cassert>
#include <cstddef>
#include <cstring>

namespace LifeV
{
//! \file fortran_wrap.h
/*
  Useful utilities for interfacing FORTRAN77 subroutines with C++
  code. Mainly based on the utilitieas developed by Carsten A. Arnholm
  (http://www.lctn.u-nancy.fr/design/hcomp.html), slightly simplified by
  Luca Formaggia (March 2002).
*/

typedef int I_F77; // INTEGER 4 bytes
typedef float R4_F77; // REAL 4 bytes
typedef double R8_F77; // DOUBLE PRECISION 8 bytes
typedef int L_F77; // LOGICAL 4 bytes

#define FALSE_F77 0
#define TRUE_F77 1

#define SUBROUTINE_F77 extern "C" void
#define INTEGER_FUNCTION_F77 extern "C" int
#define REAL_FUNCTION_F77 extern "C" float
#define LOGICAL_FUNCTION_F77 extern "C" int
#define DOUBLE_PRECISION_FUNCTION_F77 extern "C" double
#define F77NAME(X) X ## _

//!  A matrix class for interfacing with fortran
/*!
 
  A minimal class used when passing multi-dimensional array arguments from
  C++ to FORTRAN 77 (received as FORTRAN arrays), and subsequently returned
  back to C++ as properly aranged C++ arrays.
 
  \paragraph Problem
 
  FORTRAN organises data in a "column-first" order,
  while C++ organises data in a "row-first" order.
 
  \paragraph Solution
 
  (1) The FMATRIX class can take a C++ array as a constructor
  parameter. A FORTRAN compatible copy of the array is
  then made. The destructor will then copy the result back
  to the original c++ array.
 
 
  (2) The FMATRIX class provides "subscript operators" allowing
  the programmer to read and write from the array, using
  FORTRAN-like syntax and indexing semantics.
 
  /author Carsten A. Arnholm, 04-MAR-1996 (Modified by L. Formaggia MAR 2002)
 
*/



template <class T>
class FMATRIX
{
public:
    FMATRIX( size_t dim1, size_t dim2 = 1 );
    FMATRIX( T* cpparr, size_t dim1, size_t dim2 = 1 );
    operator T*();
    T& operator() ( size_t index1, size_t index2 = 0 ); // numbering from 0
    ~FMATRIX();
public:
    size_t dim[ 7 ]; // size of each dimension
    T* cpprep; // original c++ array
    T* f77rep; // array used by FORTRAN
    const size_t ndim; // number of array dimensions
};


//! F77 compatible character class
/*!
 *
 A minimal class used when passing string arguments from C++ to FORTRAN 77
 (received as FORTRAN 77 CHARACTER strings), and subsequently returned back
 to C++ as properly zero terminated strings.
 
 \paragraph Method used for zero-termination:
 
 When the CHARACTER destructor is activated the zero-termination of the
 c-string is automatically managed. Zero termination is also done each time
 a string array is subscripted using
 
 \verbatim
     CHARACTER::operator()(size_t index)
     \endverbatim
 
     \paragraph FORTRAN Assumptions:
 
     (1) F77 truncates strings when CHARACTER variable is short
 
     (2) F77 pads variable with blanks when assigned string is short
 
     (3) F77 represents a string as a pointer followed by a length
 
     (4) A string array is stored in contiguous memory
 
     \author: Carsten A. Arnholm, 20-AUG-1995
 
*/


class CHARACTER
{
public:
    CHARACTER( char* cstring );
    CHARACTER( char* cstring, const size_t lstr );
    ~CHARACTER();
    CHARACTER operator() ( size_t index );
    void pad( size_t first,size_t howmany = 1 );
    void operator=( char* str );
    operator char*();
public:
    char* rep; // Actual string
    size_t len; // String length
};

/*=====================================================================
  IMPLEMENTATIONS
  ====================================================================*/

// ============== FMATRIX

template <class T>
FMATRIX<T>::FMATRIX( size_t dim1, size_t dim2 )
        : cpprep( NULL ), f77rep( new T[ dim1*dim2 ] ), ndim( 2 )
{
    dim[ 0 ] = dim1;
    dim[ 1 ] = dim2;
    dim[ 2 ] = 0;
    dim[ 3 ] = 0;
    dim[ 4 ] = 0;
    dim[ 5 ] = 0;
    dim[ 6 ] = 0;
}

template <class T>
FMATRIX<T>::FMATRIX( T* cpparr, size_t dim1, size_t dim2 )
        : cpprep( cpparr ), f77rep( new T[ dim1*dim2 ] ), ndim( 2 )
{
    dim[ 0 ] = dim1;
    dim[ 1 ] = dim2;
    dim[ 2 ] = 0;
    dim[ 3 ] = 0;
    dim[ 4 ] = 0;
    dim[ 5 ] = 0;
    dim[ 6 ] = 0;

    // make a FORTRAN-compatible copy of the array
    size_t index_cpp = 0;
    size_t index_f77;
    for ( size_t i = 0;i < dim[ 0 ];i++ )
    {
        for ( size_t j = 0;j < dim[ 1 ];j++ )
        {
            index_f77 = j * dim[ 0 ] + i;
            f77rep[ index_f77 ] = cpprep[ index_cpp++ ];
        }
    }
}

template <class T>
FMATRIX<T>::operator T*()
{
    // Pass the FORTRAN representation when calling a function
    return f77rep;
}

template <class T>
T& FMATRIX<T>::operator() ( size_t index1, size_t index2 )
{
    assert( ndim == 2 ); // only 2d arrays supported (so far)
    // indexing according to F77 conventions
    size_t index_f77 = index2 * dim[ 0 ] + index1;
    // return a reference to the array element
    return *( f77rep + index_f77 );
}

template <class T>
FMATRIX<T>::~FMATRIX()
{
    if ( cpprep )
    {
        assert( ndim == 2 ); // only 2d arrays supported (so far)
        // copy back from FORTRAN to C++ array
        size_t index_cpp;
        size_t index_f77 = 0;
        for ( size_t j = 0;j < dim[ 1 ];j++ )
        {
            for ( size_t i = 0;i < dim[ 0 ];i++ )
            {
                index_cpp = i * dim[ 1 ] + j;
                cpprep[ index_cpp ] = f77rep[ index_f77++ ];
            }
        }
    }
    // delete the FORTRAN copy of the arry
    delete[] f77rep;
}

// ============== CHARACTER

inline CHARACTER::CHARACTER( char* cstring )
        : rep( cstring ), len( strlen( cstring ) )
{}

inline CHARACTER::CHARACTER( char* cstring, const size_t lstr )
        : rep( cstring ), len( lstr )
{
    // find position from where to start padding
    size_t slen = strlen( rep ); // upper limit
    size_t actual = ( slen < len ) ? slen : len; // actual <= len.
    for ( size_t i = actual;i < len;i++ )
        rep[ i ] = ' '; // Do the padding.
}

inline CHARACTER::~CHARACTER()
{
    if ( rep[ len ] == '\0' )
        return ; // catches string constants
    for ( int i = len - 1;i >= 0;i-- )
    {
        if ( rep[ i ] == '\0' )
            break; // already zero terminated
        if ( rep[ i ] != ' ' )
        { // non-blank discovered, so
            rep[ i + 1 ] = '\0'; // zero-terminate and jump out
            break;
        }
    }
}

inline CHARACTER CHARACTER::operator() ( size_t index )
{
    // Construct a temporary CHARACTER object for the array element
    // identified by "index" in order to zero-terminate that element
    size_t pos = index * len; // start pos of array element
    CHARACTER element( rep + pos, len ); // construct new CHARACTER.
    return element; // destructor called here.
}

inline void CHARACTER::pad( size_t first, size_t howmany )
{
    size_t pos = 0, i = 0, stop = first + howmany - 1;
    for ( size_t index = first; index <= stop; index++ )
    {
        pos = index * len;
        size_t slen = strlen( rep + pos ); // upper limit
        size_t actual = ( slen < len ) ? slen : len;
        for ( i = pos + actual;i < pos + len;i++ )
            rep[ i ] = ' '; // Do the padding.
    }
}

inline void CHARACTER::operator=( char* str )
{
    strncpy( rep, str, len ); // this will copy a zero if str < rep
    rep[ len - 1 ] = '\0'; // zero terminate in case strncpy did not
    size_t slen = strlen( rep ); // upper limit
    size_t actual = ( slen < len ) ? slen : len; // actual <= len.
    for ( size_t i = actual;i < len;i++ )
        rep[ i ] = ' '; // Do the padding.
}

inline CHARACTER::operator char*()
{
    return rep;
}
}
#endif
