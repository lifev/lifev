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
#ifndef _TAB_H_INCLUDED
#define _TAB_H_INCLUDED


#include <cmath>

#undef max
#undef min

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <lifeV.hpp>

#include <RNM.hpp>


namespace LifeV
{
const double Pi = 3.14159265358979323846264338328;
const double TGV = 1e20;

template <class T>
inline T Min ( const T &a, const T &b )
{
    return a < b ? a : b;
}
template <class T>
inline T Max ( const T &a, const T & b )
{
    return a > b ? a : b;
}
template <class T>
inline T Abs ( const T &a )
{
    return a < 0 ? -a : a;
}

template <class T>
inline void Exchange ( T& a, T& b )
{
    T c = a;
    a = b;
    b = c;
}
template <class T>
inline T Max ( const T &a, const T & b, const T & c )
{
    return Max( Max( a, b ), c );
}
template <class T>
inline T Min ( const T &a, const T & b, const T & c )
{
    return Min( Min( a, b ), c );
}

// The reals R

typedef double R;

// The class R2

class R2
{
public:
    R x, y;
    R2 () : x( 0 ), y( 0 )
    {}
    ;
    R2 ( R a, R b ) : x( a ), y( b )
    {}
    R2 ( R2 A, R2 B ) : x( B.x - A.x ), y( B.y - A.y )
    {}
    R2 operator+( R2 P ) const
    {
        return R2( x + P.x, y + P.y );
    }
    R2 operator+=( R2 P )
    {
        x += P.x;
        y += P.y;
        return *this;
    }
    R2 operator-( R2 P ) const
    {
        return R2( x -P.x, y - P.y );
    }
    R2 operator-=( R2 P )
    {
        x -= P.x;
        y -= P.y;
        return *this;
    }
    R2 operator-() const
    {
        return R2( -x, -y );
    }
    R2 operator+() const
    {
        return * this;
    }
    R operator,( R2 P ) const
    {
        return x * P.x + y * P.y;
    } // produit scalaire
    R operator^( R2 P ) const
    {
        return x * P.y - y * P.x;
    } // produit mixte
    R2 operator*( R c ) const
    {
        return R2( x * c, y * c );
    }
    R2 operator/( R c ) const
    {
        return R2( x / c, y / c );
    }
    R & operator[] ( int i )
    {
        return ( &x ) [ i ];
    }
    R2 perp()
    {
        return R2( -y, x );
    } // the perpendicular
    friend R2 operator*( R c, R2 P )
    {
        return P * c;
    }
    friend R2 operator/( R c, R2 P )
    {
        return P / c;
    }

    friend std::ostream& operator << ( std::ostream& f, const R2 & P )
    {
        f << P.x << ' ' << P.y ;
        return f;
    }
    friend std::istream& operator >>( std::istream& f, R2 & P )
    {
        f >> P.x >> P.y ;
        return f;
    }
};

// The class R3
class R3: public R2
{
public:
    R z;
    R3 () : z( 0 )
    {}
    ;
    R3 ( R a, R b, R c ) : R2( a, b ), z( c )
    {}
    R3 ( R3 A, R3 B ) : R2( B.x - A.x, B.y - A.y ), z( B.z - A.z )
    {}
    R3 operator+( R3 P ) const
    {
        return R3( x + P.x, y + P.y, z + P.z );
    }
    R3 operator+=( R3 P )
    {
        x += P.x;
        y += P.y;
        z += P.z;
        return *this;
    }
    R3 operator-( R3 P ) const
    {
        return R3( x -P.x, y - P.y, z - P.z );
    }
    R3 operator-=( R3 P )
    {
        x -= P.x;
        y -= P.y;
        z -= P.z;
        return *this;
    }
    R3 operator-() const
    {
        return R3( -x, -y, -z );
    }
    R3 operator+() const
    {
        return * this;
    }
    R operator,( R3 P ) const
    {
        return x * P.x + y * P.y + z * P.z;
    } // produit scalaire
    R3 operator^( R3 P ) const
    {
        return R3( y * P.z - z * P.y , P.x * z - x * P.z, x * P.y - y * P.x );
    } // produit mixte
    R3 operator*( R c ) const
    {
        return R3( x * c, y * c, z * c );
    }
    R3 operator/( R c ) const
    {
        return R3( x / c, y / c, z / c );
    }
    R & operator[] ( int i )
    {
        return ( &x ) [ i ];
    }
    friend R3 operator*( R c, R3 P )
    {
        return P * c;
    }
    friend R3 operator/( R c, R3 P )
    {
        return P / c;
    }
    friend std::ostream& operator <<( std::ostream& f, const R3 & P )
    {
        f << P.x << ' ' << P.y << ' ' << P.z ;
        return f;
    }
    friend std::istream& operator >>( std::istream& f, R3 & P )
    {
        f >> P.x >> P.y >> P.z ;
        return f;
    }
};

inline R Area2( const R2 A, const R2 B, const R2 C )
{
    return ( B -A ) ^ ( C - A );
}
inline R Norme2_2( const R2 & A )
{
    return ( A, A );
}
inline R Norme2( const R2 & A )
{
    return sqrt( ( A, A ) );
}
inline R Norme2_2( const R3 & A )
{
    return ( A, A );
}
inline R Norme2( const R3 & A )
{
    return sqrt( ( A, A ) );
}
inline R Norme_infty( const R2 & A )
{
    return Max( Abs( A.x ), Abs( A.y ) );
}
inline R Norme_infty( const R3 & A )
{
    return Max( Abs( A.x ), Abs( A.y ), Abs( A.z ) );
}
// inline R Theta(R2 P){ return atan2(P.y,P.x);}
}

//# include "RNM.hpp"
//# include "RNM.hpp"

namespace LifeV
{
typedef KNM<R> RNM;
typedef KNM_<R> RNM_;
typedef KN<R> RN;
typedef KN_<R> RN_;

#if 0
typedef boost::numeric::ublas::matrix<R> RNM;
typedef boost::numeric::ublas::matrix<R> RNM_;
typedef boost::numeric::ublas::vector<R> RN;
typedef boost::numeric::ublas::vector<R> RN_;
#endif

typedef boost::numeric::ublas::matrix<double> Tab2d;
typedef Tab2d Matrix;
typedef boost::numeric::ublas::matrix_range<Tab2d> Tab2dView;
typedef boost::numeric::ublas::range TabRange;

typedef boost::numeric::ublas::zero_matrix<double> ZeroMatrix;
typedef boost::numeric::ublas::scalar_matrix<double> ScalarMatrix;
typedef boost::numeric::ublas::identity_matrix<double> IdentityMatrix;


typedef boost::numeric::ublas::vector<double> Tab1d;
typedef boost::numeric::ublas::vector_range<Tab1d> Tab1dView;
typedef boost::numeric::ublas::vector<int> Tab1dInt;
typedef boost::numeric::ublas::vector<int> Tab1dViewInt;

typedef boost::numeric::ublas::vector<double> Vector;
typedef boost::numeric::ublas::unit_vector<double> UnitVector;
typedef boost::numeric::ublas::scalar_vector<double> ScalarVector;
typedef boost::numeric::ublas::zero_vector<double> ZeroVector;

inline double
dot( Vector const &ex_v1, Vector const &ex_v2 )
{
    return boost::numeric::ublas::inner_prod( ex_v1, ex_v2 );
}

/* reduce operations */
using boost::numeric::ublas::sum;
using boost::numeric::ublas::norm_1;
using boost::numeric::ublas::norm_2;
using boost::numeric::ublas::norm_inf;
using boost::numeric::ublas::index_norm_inf;

/* binary operations */
using boost::numeric::ublas::outer_prod;
using boost::numeric::ublas::inner_prod;

struct norm_inf_adaptor
{
    template<typename E>
    Real operator()( E const& __v ) const
        {
            return norm_inf( __v );
        }
};
}
#endif
