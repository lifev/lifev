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
  @brief Fortran wrapper for C++

  @date 1-03-2002
  @author Luca Formaggia

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
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
  code. Mainly based on the utilities developed by Carsten A. Arnholm
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
template <class ScalarType>
class FMATRIX
{
public:
    //! @name Public typedefs
    //@{
    typedef ScalarType scalar_Type;
    //! @name Constructors and destructor
    //@{
    FMATRIX() {}
    FMATRIX( size_t dim1, size_t dim2 = 1 );
    FMATRIX( ScalarType* cppArray, size_t dim1, size_t dim2 = 1 );
    virtual ~FMATRIX();
    //@}

    //! @name Operators
    //@{
    operator ScalarType*();
    ScalarType& operator() ( size_t index1, size_t index2 = 0 ); // numbering from 0
    //@}

public:
    size_t M_arrayDimensions[ 7 ]; // size of each dimension
    ScalarType* M_cppRepresentation; // original c++ array
    ScalarType* M_fortranRepresentation; // array used by FORTRAN
    const size_t M_numDimensions; // number of array dimensions
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
    //! @name Constructors and destructor
    //@{
    CHARACTER() {}
    CHARACTER( char* cstring );
    CHARACTER( char* cstring, const size_t stringLength );
    virtual ~CHARACTER();
    //@}

    //! @name Operators
    //@{
    CHARACTER operator() ( size_t index );
    void operator=( char* str );
    operator char*();
    //@}

    //! @name Methods
    //@{
    void pad( size_t first, size_t paddingSize = 1 );
    //@}

public:
    char* M_representation; // Actual string
    size_t M_length; // String length
};

// ======================
// IMPLEMENTATIONS
// ======================

// =======================================
// FMATRIX Constructors and destructor
// =======================================

template <class ScalarType>
FMATRIX<ScalarType>::FMATRIX( size_t dim1, size_t dim2 ):
    M_cppRepresentation( NULL ),
    M_fortranRepresentation( new ScalarType[ dim1*dim2 ] ),
    M_numDimensions( 2 )
{
    M_arrayDimensions[ 0 ] = dim1;
    M_arrayDimensions[ 1 ] = dim2;
    M_arrayDimensions[ 2 ] = 0;
    M_arrayDimensions[ 3 ] = 0;
    M_arrayDimensions[ 4 ] = 0;
    M_arrayDimensions[ 5 ] = 0;
    M_arrayDimensions[ 6 ] = 0;
}

template <class ScalarType>
FMATRIX<ScalarType>::FMATRIX( ScalarType* cppArray, size_t dim1, size_t dim2 ):
    M_cppRepresentation( cppArray ),
    M_fortranRepresentation( new ScalarType[ dim1*dim2 ] ),
    M_numDimensions( 2 )
{
    M_arrayDimensions[ 0 ] = dim1;
    M_arrayDimensions[ 1 ] = dim2;
    M_arrayDimensions[ 2 ] = 0;
    M_arrayDimensions[ 3 ] = 0;
    M_arrayDimensions[ 4 ] = 0;
    M_arrayDimensions[ 5 ] = 0;
    M_arrayDimensions[ 6 ] = 0;

    // make a FORTRAN-compatible copy of the array
    size_t index_cpp = 0;
    size_t index_f77;
    for ( size_t i = 0; i < M_arrayDimensions[ 0 ]; i++ )
    {
        for ( size_t j = 0; j < M_arrayDimensions[ 1 ]; j++ )
        {
            index_f77 = j * M_arrayDimensions[ 0 ] + i;
            M_fortranRepresentation[ index_f77 ] = M_cppRepresentation[ index_cpp++ ];
        }
    }
}

template <class ScalarType>
FMATRIX<ScalarType>::~FMATRIX()
{
    if ( M_cppRepresentation )
    {
        assert( M_numDimensions == 2 ); // only 2d arrays supported (so far)
        // copy back from FORTRAN to C++ array
        size_t index_cpp;
        size_t index_f77 = 0;
        for ( size_t j = 0; j < M_arrayDimensions[ 1 ]; j++ )
        {
            for ( size_t i = 0; i < M_arrayDimensions[ 0 ]; i++ )
            {
                index_cpp = i * M_arrayDimensions[ 1 ] + j;
                M_cppRepresentation[ index_cpp ] = M_fortranRepresentation[ index_f77++ ];
            }
        }
    }
    // delete the FORTRAN copy of the arry
    delete[] M_fortranRepresentation;
}

// ============================
// FMATRIX Operators
// ============================

template <class ScalarType>
FMATRIX<ScalarType>::operator ScalarType*()
{
    // Pass the FORTRAN representation when calling a function
    return M_fortranRepresentation;
}

template <class ScalarType>
ScalarType& FMATRIX<ScalarType>::operator() ( size_t index1, size_t index2 )
{
    assert( M_numDimensions == 2 ); // only 2d arrays supported (so far)
    // indexing according to F77 conventions
    size_t index_f77 = index2 * M_arrayDimensions[ 0 ] + index1;
    // return a reference to the array element
    return *( M_fortranRepresentation + index_f77 );
}

// =======================================
// CHARACTER Constructors and destructor
// =======================================

inline CHARACTER::CHARACTER( char* cstring ):
    M_representation( cstring ),
    M_length( strlen( cstring ) )
{}

inline CHARACTER::CHARACTER( char* cstring, const size_t stringLength ):
    M_representation( cstring ),
    M_length( stringLength )
{
    // find position from where to start padding
    size_t slen = strlen( M_representation ); // upper limit
    size_t actual = ( slen < M_length ) ? slen : M_length; // actual <= M_length.
    for ( size_t i = actual; i < M_length; i++ )
        M_representation[ i ] = ' '; // Do the padding.
}

inline CHARACTER::~CHARACTER()
{
    if ( M_representation[ M_length ] == '\0' )
        return ; // catches string constants
    for ( int i = M_length - 1; i >= 0; i-- )
    {
        if ( M_representation[ i ] == '\0' )
            break; // already zero terminated
        if ( M_representation[ i ] != ' ' )
        { // non-blank discovered, so
            M_representation[ i + 1 ] = '\0'; // zero-terminate and jump out
            break;
        }
    }
}

// ================================
// CHARACTER Operators
// ================================

inline CHARACTER CHARACTER::operator() ( size_t index )
{
    // Construct a temporary CHARACTER object for the array element
    // identified by "index" in order to zero-terminate that element
    size_t pos = index * M_length; // start pos of array element
    CHARACTER element( M_representation + pos, M_length ); // construct new CHARACTER.
    return element; // destructor called here.
}

inline void CHARACTER::operator=( char* str )
{
    strncpy( M_representation, str, M_length ); // this will copy a zero if str < rep
    M_representation[ M_length - 1 ] = '\0'; // zero terminate in case strncpy did not
    size_t slen = strlen( M_representation ); // upper limit
    size_t actual = ( slen < M_length ) ? slen : M_length; // actual <= M_length.
    for ( size_t i = actual; i < M_length; i++ )
        M_representation[ i ] = ' '; // Do the padding.
}

inline CHARACTER::operator char*()
{
    return M_representation;
}

// ================================
// CHARACTER Methods
// ================================

inline void CHARACTER::pad( size_t first, size_t paddingSize )
{
    size_t pos = 0, i = 0, stop = first + paddingSize - 1;
    for ( size_t index = first; index <= stop; index++ )
    {
        pos = index * M_length;
        size_t slen = strlen( M_representation + pos ); // upper limit
        size_t actual = ( slen < M_length ) ? slen : M_length;
        for ( i = pos + actual; i < pos + M_length; i++ )
            M_representation[ i ] = ' '; // Do the padding.
    }
}

} // Namespace LifeV

#endif // FORTRAN_FROM_CPLUSPLUS
