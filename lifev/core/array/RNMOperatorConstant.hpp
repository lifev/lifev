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
    @brief Operator implementation of RNM.hpp

    @contributor Miguel Fernandez
                 Matteo Astorino <matteo.astorino@epfl.ch>
    @mantainer Matteo Astorino <matteo.astorino@epfl.ch>

 */
// Modif Miguel 02/12/2002
// "operator=" was changed to "operator oper" when
// tab are consecutive in memory
namespace LifeV
{
template <class R>
const KN_<R>& KN_<R>::operator oper ( const_R a )
{
    R * l( v );
    for ( int i = 0; i < n; i++, l += step )
        * l oper a;
    return *this;
}

template <class R>
inline const KNM_<R> & KNM_<R>::operator oper ( const_R a )
{
    if ( IsVector1() )
        KN_<R>::operator oper ( a );
    else
    {
        KN_<R> lj( operator() ( '.',0 ) ); //  (.,.,O)
        for ( int j = 0; j < M(); ++j, ++lj )
            lj oper a;
    }
    return *this;
}

template <class R>
inline const KNMK_<R> & KNMK_<R>::operator oper ( const_R a )
{
    if ( IsVector1() )
        KN_<R>::operator oper ( a );
    else
    {
        KNM_<R> lj( operator() ( '.','.', 0 ) ); //  (.,.,O)
        int j = K();
        while ( j-- )
        {
            lj oper a;
            ++lj;
        }
    }
    return *this;
}

template <class R>
const KN_<R>& KN_<R>::operator oper ( const KN_<const_R> & u )
{
    K_assert( u.n == n );
    R * l( v );
    const R *r( u );
    for ( int i = 0; i < n; i++, l += step, r += u.step )
        * l oper * r;
    return *this;
}

template <class R>
inline const KNM_<R> & KNM_<R>::operator oper ( const KNM_<const_R> & u )
{
    if ( IsVector1() && u.IsVector1() )
    {
        KN_<R>::operator oper ( u );
    }
    else
    {
        KN_<R> lj( operator() ( '.',0 ) ); //  (.,O)
        KN_<const_R> uj( u( '.', 0 ) );
        int j = M();
        while ( j-- )
        {
            lj oper uj;
            ++lj;
            ++uj;
        }
    }
    return *this;
}


template <class R>
inline const KNMK_<R> & KNMK_<R>::operator oper ( const KNMK_<const_R> & u )
{
    if ( IsVector1() && u.IsVector1() )
        KN_<R>::operator oper ( u );
    else
    {
        K_assert( K() == u.K() );
        KNM_<R> lj( operator() ( '.','.', 0 ) ); //  (.,O)
        KNM_<const_R> uj( u( '.', '.', 0 ) );
        int j = K();
        while ( j-- )
        {
            lj oper uj;
            ++lj;
            ++uj;
        }
    }
    return *this;
}

#undef oper
}
