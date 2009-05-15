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
namespace LifeV
{
template <class R>
const KN_<R>& KN_<R>::operator oper ( const Mul_KNM_KN_<R> & u )
{
    K_assert ( SameShape( u.A.shapei ) && !constant() );
    R * l( v );
    KN_<const_R> li( u.A( 0, '.' ) ); //  first line
    for ( int i = 0;i < n;i++, l += step, ++li )
        * l oper ( li, u.b );
    return *this;
}


template <class R>
const KN_<R>& KN_<R>::operator oper ( const Add_KN_<R> & u )
{
    K_assert( u.a.N() == N() );
    int stepa( u.a.step ), stepb( u.b.step );
    R * l( v );
    const_R *aa( u.a ), *bb( u.b );
    for ( int i = 0;i < n;i++, l += step, aa += stepa, bb += stepb )
        * l oper * aa + *bb;
    return *this;
}

template <class R>
const KN_<R>& KN_<R>::operator oper ( const Sub_KN_<R> & u )
{
    K_assert( u.a.N() == N() );
    int stepa( u.a.step ), stepb( u.b.step );
    R * l( v );
    const_R *aa( u.a ), *bb( u.b );
    for ( int i = 0;i < n;i++, l += step, aa += stepa, bb += stepb )
        * l oper * aa - *bb;
    return *this;
}

template <class R>
const KN_<R>& KN_<R>::operator oper ( const Mulc_KN_<R> & u )
{
    K_assert( u.a.N() == N() );
    int stepa( u.a.step );
    R * l( v );
    const_R *aa( u.a ), bb( u.b ) ;
    for ( int i = 0;i < n;i++, l += step, aa += stepa )
        * l oper * aa * bb;
    return *this;
}

template <class R>
const KN_<R>& KN_<R>::operator oper ( const Add_Mulc_KN_<R> & u )
{
    K_assert( u.a.N() == N() );
    const int stepa( u.a.step ), stepb( u.b.step );
    const R ca( u.ca ), cb( u.cb );
    R * l( v );
    const R *aa( u.a ), *bb( u.b );
    for ( int i = 0;i < n;i++, l += step, aa += stepa, bb += stepb )
        * l oper * aa * ca + *bb * cb;
    return *this;
}
}
