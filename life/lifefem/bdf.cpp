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
#include <stdexcept>
#include <sstream>

#include "bdf.hpp"

namespace LifeV
{
Bdf::Bdf( const UInt n )
    :
    _M_order( n ),
    _M_size( 0 ),
    _M_alpha( n + 1 ),
    _M_beta( n )
{
    if ( n <= 0 || n > BDF_MAX_ORDER )
    {
        std::ostringstream __ex;
        __ex << "Error: wrong BDF order\n"
        << " you want to use BDF order " << n << "\n"
        << " we support BDF order from 1 to " << BDF_MAX_ORDER << "\n";
        throw std::invalid_argument( __ex.str() );
    }
    switch ( n )
    {
    case 1:
        _M_alpha[ 0 ] = 1.; // Backward Euler
        _M_alpha[ 1 ] = 1.;
        _M_beta[ 0 ] = 1.; // u^{n+1} \approx u^n
        break;
    case 2:
        _M_alpha[ 0 ] = 3. / 2.;
        _M_alpha[ 1 ] = 2.;
        _M_alpha[ 2 ] = -1. / 2.;
        _M_beta[ 0 ] = 2.;
        _M_beta[ 1 ] = -1.;
        break;
    case 3:
        _M_alpha[ 0 ] = 11. / 6.;
        _M_alpha[ 1 ] = 3.;
        _M_alpha[ 2 ] = -3. / 2.;
        _M_alpha[ 3 ] = 1. / 3.;
        _M_beta[ 0 ] = 3.;
        _M_beta[ 1 ] = -3.;
        _M_beta[ 2 ] = 1.;
        break;
    }
    _M_unknowns.resize( n );
}



Bdf::~Bdf()
{}


void Bdf::initialize_unk( Vector u0 )
{
    std::vector< Vector >::iterator iter = _M_unknowns.begin();
    std::vector< Vector >::iterator iter_end = _M_unknowns.end();

    _M_size = u0.size();

    for ( iter = _M_unknowns.begin() ; iter != iter_end; iter++ )
    {
        *iter = u0;
    }

    return ;
}


void Bdf::initialize_unk( std::vector<Vector> uv0 )
{
    std::vector< Vector >::iterator iter = _M_unknowns.begin();
    std::vector< Vector >::iterator iter_end = _M_unknowns.end();

    std::vector< Vector >::iterator iter0 = uv0.begin();

    _M_size = iter0->size();

    UInt n0 = uv0.size();

    // Check if uv0 has the right dimensions
    ASSERT( n0 >= _M_order, "Initial data are not enough for the selected BDF" )

    // if n0>n, only the first n inital data will be considered
    if ( n0 > _M_order )
    {
        std::cout << "The initial data set is larger than needed by the BDF."
        << std::endl;
        std::cout << "Only the first " << _M_order << " data will be considered. "
        << std::endl;
    }
    for ( iter = _M_unknowns.begin() ; iter != iter_end; iter++ )
    {
        *iter = *iter0;
        iter0++;
    }

    return ;
}

//
//
//





const
std::vector<Vector>& Bdf::unk() const
{
    return _M_unknowns;
}


double
Bdf::coeff_der( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 & i < _M_order + 1,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_alpha[ i ];
}

double
Bdf::coeff_ext( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 & i < _M_order,
            "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _M_beta[ i ];
}

void
Bdf::showMe() const
{
    std::cout << "*** BDF Time discretization of order " << _M_order << " ***"
              << std::endl;
    std::cout << "    Coefficients: " << std::endl;
    for ( UInt i = 0;i < _M_order + 1;++i )
        std::cout << "       alpha(" << i << ") = " << _M_alpha[ i ]
                  << std::endl;
    for ( UInt i = 0;i < _M_order;++i )
        std::cout << "       beta (" << i << ") = " << _M_beta[ i ]
                  << std::endl;

    std::cout << "    " << _M_unknowns.size() << " unknown vector"
              << (_M_unknowns.size() == 1 ? "" : "s") << " of length "
              << _M_size << std::endl;

    /*   std::vector< Vector >::iterator itg=_u_prec.begin(); */
    /*   std::vector< Vector >::iterator itge=_u_prec.end(); */
    /*   for ( ; itg!=itge;++itg){ */
    /*     std::cout << "u_prec: " << itg-_u_prec.begin() << std::endl; */
    /*     Vector::iterator itl=itg->begin(); */
    /*     Vector::iterator itle=itg->end(); */
    /*     std::cout << "[ "; */
    /*     for ( ;itl!=itle; ++itl) */
    /*       std::cout << *itl << ", " ;   */

    /*     std::cout << "]" << std::endl; */
    /*   }  */


    return ;
}

void
Bdf::shift_right( Vector const& unk_curr )
{
#if 1
    std::vector< Vector >::iterator it = _M_unknowns.end() - 1;
    std::vector< Vector >::iterator itm1 = _M_unknowns.end() - 1;
    std::vector< Vector >::iterator itb = _M_unknowns.begin();

    for ( ; it != itb; --it )
    {
        itm1--;
        *it = *itm1;
    }
    *itb = unk_curr;

#else
    std::vector< Vector >::reverse_iterator it = _M_unknowns.rbegin();
    std::vector< Vector >::reverse_iterator itm1 = _M_unknowns.rbegin();
    std::vector< Vector >::reverse_iterator itb = _M_unknowns.rend();

    for ( ; it != itb; ++it )
    {
        ++itm1;
        *it = *itm1;
    }
    _M_unknowns.front() = unk_curr;
#endif

}


Vector
Bdf::time_der( Real dt ) const
{
    Vector ut( _M_size );

    for ( UInt j = 0;j < _M_size;++j )
        ut[ j ] = 0.;

    for ( UInt j = 0;j < _M_size;++j )
    {
        for ( UInt i = 0;i < _M_order;++i )
            ut[ j ] = ut[ j ] + _M_alpha[ i + 1 ] / dt * _M_unknowns[ i ][ j ];
    }

    return ut;
}

Vector
Bdf::time_der() const
{
    Vector ut( _M_size );

    for ( UInt j = 0;j < _M_size;++j )
        ut[ j ] = 0.;

    for ( UInt j = 0;j < _M_size;++j )
    {
        for ( UInt i = 0;i < _M_order;++i )
            ut[ j ] = ut[ j ] + _M_alpha[ i + 1 ] * _M_unknowns[ i ][ j ];
    }

    return ut;
}

Vector
Bdf::extrap() const
{
    Vector ue( _M_size );
    for ( UInt j = 0;j < _M_size;++j )
        ue[ j ] = 0.;

    for ( UInt j = 0;j < _M_size;++j )
    {
        for ( UInt i = 0;i < _M_order;++i )
            ue[ j ] = ue[ j ] + _M_beta[ i ] * _M_unknowns[ i ][ j ];
    }

    return ue;
}
}
