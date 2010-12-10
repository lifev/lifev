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
  @brief Collection of basic vector and matrix operations

  @date 10-10.2010
  @author

  @contributor
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef TAB_H
#define TAB_H

#include <cmath>

#undef max
#undef min

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifearray/RNM.hpp>
#include <life/lifecore/life.hpp>

namespace LifeV
{
/*!
  @brief Basic vector and matrix operations

  A collection of typedefs, functions (and functor object) and class
  that provide basic vector and matrix operations.
*/

//! @name Constants
//@{
const Real Pi = 3.14159265358979323846264338328;
const Real TGV = 1e20;
//@}

//! @name Using clauses
//@{
/* reduce operations */
using boost::numeric::ublas::index_norm_inf; // TODO: REMOVE UNUSED
using boost::numeric::ublas::norm_1; // TODO: REMOVE UNUSED
using boost::numeric::ublas::norm_2;
using boost::numeric::ublas::norm_inf;
using boost::numeric::ublas::sum;

/* binary operations */
using boost::numeric::ublas::inner_prod; // TODO: REMOVE UNUSED
using boost::numeric::ublas::outer_prod; // TODO: REMOVE UNUSED
//@}

//! @name Typedefs // TODO: AFTER REMOVING R2 and R3 check unused and remove
//@{
typedef KN<Real> RN;
typedef KNM<Real> RNM;
typedef KNM_<Real> RNM_;
typedef KN_<Real> RN_;
typedef boost::numeric::ublas::matrix<Real> Tab2d;
typedef Tab2d Matrix;
typedef boost::numeric::ublas::matrix_range<Tab2d> Tab2dView;
typedef boost::numeric::ublas::range TabRange;
typedef boost::numeric::ublas::zero_matrix<Real> ZeroMatrix;
typedef boost::numeric::ublas::scalar_matrix<Real> ScalarMatrix; // TODO: REMOVE UNUSED
typedef boost::numeric::ublas::identity_matrix<Real> IdentityMatrix; // TODO: REMOVE UNUSED
typedef boost::numeric::ublas::vector<Real> Tab1d;
typedef boost::numeric::ublas::vector_range<Tab1d> Tab1dView;
typedef boost::numeric::ublas::vector<Int> Tab1dInt; // TODO: REMOVE UNUSED
typedef boost::numeric::ublas::vector<Int> Tab1dViewInt; // TODO: REMOVE UNUSED
typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::unit_vector<Real> UnitVector; // TODO: REMOVE UNUSED
typedef boost::numeric::ublas::scalar_vector<Real> ScalarVector;
typedef boost::numeric::ublas::zero_vector<Real> ZeroVector;
// Type for storing the geometric vectors (instead of x,y,z)
typedef boost::numeric::ublas::vector<Real> GeoVector;
//@}

//! @name 2D vector class // TODO: REMOVE UNUSED
class R2
{
public:
//! @name Constructors and destructor
//@{
    //! Default constructor
    R2 () : x( 0 ), y( 0 ) {}
    //! Constructor
    /*!
      @param a, b - Real
    */
    R2 ( Real a, Real b ) : x( a ), y( b ) {}

    //! Constructor
    /*!
      @param A, B - R2 objects
    */
    R2 ( R2 A, R2 B ) : x( B.x - A.x ), y( B.y - A.y ) {}
//@}

//! @name Operators
//@{
    //! Addition operator
    R2 operator+( R2 P ) const
    {
        return R2( x + P.x, y + P.y );
    }

    //! Addition-update operator
    R2 operator+=( R2 P )
    {
        x += P.x;
        y += P.y;
        return *this;
    }

    //! Substraction operator
    R2 operator-( R2 P ) const
    {
        return R2( x -P.x, y - P.y );
    }

    //! Substration-update operator
    R2 operator-=( R2 P )
    {
        x -= P.x;
        y -= P.y;
        return *this;
    }

    //! Opposite operator
    R2 operator-() const
    {
        return R2( -x, -y );
    }

    //! Should be removed
    R2 operator+() const
    {
        return * this;
    }

    //! Comma operator - Scalar product
    Real operator,( R2 P ) const
    {
        return x * P.x + y * P.y;
    }

    //! Hat operator - Mixt product
    Real operator^( R2 P ) const
    {
        return x * P.y - y * P.x;
    }

    //! Multiplication operator
    R2 operator*( Real c ) const
    {
        return R2( x * c, y * c );
    }

    //! Division operator
    R2 operator/( Real c ) const
    {
        return R2( x / c, y / c );
    }

    //! Index operator
    Real & operator[] ( Int i )
    {
        return ( &x ) [ i ];
    }
//@}

//! @name Methods
//@{
    //! Computes the perpendicular of a vector
    R2 perp()
    {
        return R2( -y, x );
    }
//@}

//! Friend functions
//@{
    //! Multiplication operator
    friend R2 operator*( Real c, R2 P )
    {
        return P * c;
    }

    //! Division operator
    friend R2 operator/( Real c, R2 P )
    {
        return P / c;
    }

    //! Stream insertion operator
    friend std::ostream& operator << ( std::ostream& f, const R2 & P )
    {
        f << P.x << ' ' << P.y ;
        return f;
    }

    //! Stream extraction operator
    friend std::istream& operator >>( std::istream& f, R2 & P )
    {
        f >> P.x >> P.y ;
        return f;
    }
//@}

//! Data members
//@{
    Real x, y;
//@}
};

//! @name 3D vector class // TODO: REMOVE UNUSED
class R3: public R2
{
public:
//! @name Constructors and destructor
//@{
    //! Default constructor
    R3 () : z( 0 ) {}
    //! Constructor
    /*!
      @param a, b, c - Real
    */
    R3 ( Real a, Real b, Real c ) : R2( a, b ), z( c ) {}
    //! Constructor
    /*!
      @param A, B - R3 object
    */
    R3 ( R3 A, R3 B ) : R2( B.x - A.x, B.y - A.y ), z( B.z - A.z ) {}
//@}

//! @name Operators
//@{
    //! Addition operator
    R3 operator+( R3 P ) const
    {
        return R3( x + P.x, y + P.y, z + P.z );
    }

    //! Addition-update operator
    R3 operator+=( R3 P )
    {
        x += P.x;
        y += P.y;
        z += P.z;
        return *this;
    }

    //! Substraction operator
    R3 operator-( R3 P ) const
    {
        return R3( x -P.x, y - P.y, z - P.z );
    }

    //! Substration-update operator
    R3 operator-=( R3 P )
    {
        x -= P.x;
        y -= P.y;
        z -= P.z;
        return *this;
    }

    //! Opposite operator
    R3 operator-() const
    {
        return R3( -x, -y, -z );
    }

    //! Should be removed
    R3 operator+() const
    {
        return * this;
    }

    //! Comma operator - Scalar product
    Real operator,( R3 P ) const
    {
        return x * P.x + y * P.y + z * P.z;
    }

    //! Hat operator - Mixt product
    R3 operator^( R3 P ) const
    {
        return R3( y * P.z - z * P.y , P.x * z - x * P.z, x * P.y - y * P.x );
    }

    //! Multiplication operator
    R3 operator*( Real c ) const
    {
        return R3( x * c, y * c, z * c );
    }

    //! Division operator
    R3 operator/( Real c ) const
    {
        return R3( x / c, y / c, z / c );
    }

    //! Index operator
    Real & operator[] ( Int i )
    {
        return ( &x ) [ i ];
    }
//@}

//! @name Friend functions
//@{
    //! Multiplication operator
    friend R3 operator*( Real c, R3 P )
    {
        return P * c;
    }

    //! Division operator
    friend R3 operator/( Real c, R3 P )
    {
        return P / c;
    }

    //! Stream insertion operator
    friend std::ostream& operator <<( std::ostream& f, const R3 & P )
    {
        f << P.x << ' ' << P.y << ' ' << P.z ;
        return f;
    }

    //! Stream extraction operator
    friend std::istream& operator >>( std::istream& f, R3 & P )
    {
        f >> P.x >> P.y >> P.z ;
        return f;
    }
//@}

//! @name Data members
//@{
    Real z;
//@}
};

//! @name Functor that returns the infinity norm of a templated value
struct norm_inf_adaptor // TODO: REMOVE UNUSED
{
    template<typename T>
    Real operator()( T const& v ) const
    {
        return  v.NormInf();
    }
};

//! @name Functions
//@{
//! Templated function that returns the minimum of two values
/*!
  @param a, b - const reference to T
  @return - T
*/
template <class T>
inline T Min ( const T &a, const T &b )
{
    return a < b ? a : b;
}

//! Templated function that returns the maximum of two values
/*!
  @param a, b - const reference to T
  @return - T
*/
template <class T>
inline T Max ( const T &a, const T & b )
{
    return a > b ? a : b;
}

//! Templated function that return the absolute value
/*!
  @param b - const reference to T
  @return - T
*/
template <class T>
inline T Abs ( const T &a )
{
    return a < 0 ? -a : a;
}

//! Templated function that exchanges two values
/*!
  @param a, b - reference to T
*/
template <class T>
inline void Exchange ( T& a, T& b )
{
    T c = a;
    a = b;
    b = c;
}

//! Templated function that returns the maximum of three values
/*!
  @param a, b, c - const reference to T
  @return - T
*/
template <class T>
inline T Max ( const T &a, const T & b, const T & c )
{
    return Max( Max( a, b ), c );
}

//! Templated function that returns the minimum of three values
/*!
  @param a, b, c - const reference to T
  @return - T
*/
template <class T>
inline T Min ( const T &a, const T & b, const T & c )
{
    return Min( Min( a, b ), c );
}

//! Function that computes the square of the area designated by three 2D vectors
/*!
  @param A, B, C  - const R2 objects
  @return Real
*/
inline Real Area2( const R2 A, const R2 B, const R2 C ) // TODO: REMOVE UNUSED
{
    return ( B - A ) ^ ( C - A );
}

//! Function that computes the square of the Euclidian norm of a 2D vector
/*!
  @param A - const reference to R2 object
  @return Real
*/
inline Real Norme2_2( const R2 & A ) // TODO: REMOVE UNUSED
{
    return ( A, A );
}

//! Function that computed the Euclidian norm of a 2D vector
/*!
  @param A - const reference to R2 object
  @return Real
*/
inline Real Norme2( const R2 & A ) // TODO: REMOVE UNUSED
{
    return sqrt( ( A, A ) );
}

//! Function that computes the square of the Euclidian norm of a 3D vector
/*!
  @param A - const reference to R3 object
  @return - Real
*/
inline Real Norme2_2( const R3 & A ) // TODO: REMOVE UNUSED
{
    return ( A, A );
}

//! Function that computed the Euclidian norm of a 3D vector
/*!
  @param A - const reference to R3 object
  @return - Real
*/
inline Real Norme2( const R3 & A ) // TODO: REMOVE UNUSED
{
    return sqrt( ( A, A ) );
}

//! Function that computes the infinity norm of a 2D vector
/*!
  @param A - const reference to R2 object
  @return - Real
*/
inline Real Norme_infty( const R2 & A ) // TODO: REMOVE UNUSED
{
    return Max( Abs( A.x ), Abs( A.y ) );
}

//! Function that computes the infinity norm of a 3D vector
/*!
  @param A - const reference to R2 object
  @return - Real
*/
inline Real Norme_infty( const R3 & A ) // TODO: REMOVE UNUSED
{
    return Max( Abs( A.x ), Abs( A.y ), Abs( A.z ) );
}

//! Function that computes the scalar product of two arbitrarily sized vectors
/*!
  @param ex_v1, ex_v2 - const references to Vector objects
  @return - Real
*/
inline Real dot( Vector const &ex_v1, Vector const &ex_v2 )
{
    return boost::numeric::ublas::inner_prod( ex_v1, ex_v2 );
}

//! Overload of the stream insertion operator for printing arbitrarily sized vectors
/*!
  @param s - reference to ostream
  @param x - const reference to Vector object
*/
inline std::ostream & operator << (std::ostream & s, const Vector& x)
{
    for (UInt iz = 0; iz < x.size() ; iz ++ )
        s << x[iz] << " " ;

    return s;
}
//@}
}

#endif // TAB_H
