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
Bdf::Bdf( const UInt n ) :
        _n( n ), _s( 0 ), _alpha( n + 1 ), _beta( n )
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
        _alpha[ 0 ] = 1.; // Backward Euler
        _alpha[ 1 ] = 1.;
        _beta[ 0 ] = 1.; // u^{n+1} \approx u^n
        break;
    case 2:
        _alpha[ 0 ] = 3. / 2.;
        _alpha[ 1 ] = 2.;
        _alpha[ 2 ] = -1. / 2.;
        _beta[ 0 ] = 2.;
        _beta[ 1 ] = -1.;
        break;
    case 3:
        _alpha[ 0 ] = 11. / 6.;
        _alpha[ 1 ] = 3.;
        _alpha[ 2 ] = -3. / 2.;
        _alpha[ 3 ] = 1. / 3.;
        _beta[ 0 ] = 3.;
        _beta[ 1 ] = -3.;
        _beta[ 2 ] = 1.;
        break;
    }
    _unk.resize( n );
}



Bdf::~Bdf()
{}


void Bdf::initialize_unk( Vector u0 )
{
    std::vector< Vector >::iterator iter = _unk.begin();
    std::vector< Vector >::iterator iter_end = _unk.end();

    _s = u0.size();

    for ( iter = _unk.begin() ; iter != iter_end; iter++ )
    {
        *iter = u0;
    }

    return ;
}


void Bdf::initialize_unk( std::vector<Vector> uv0 )
{
    std::vector< Vector >::iterator iter = _unk.begin();
    std::vector< Vector >::iterator iter_end = _unk.end();

    std::vector< Vector >::iterator iter0 = uv0.begin();

    _s = iter0->size();

    UInt n0 = uv0.size();

    // Check if uv0 has the right dimensions
    ASSERT( n0 >= _n, "Initial data are not enough for the selected BDF" )

    // if n0>n, only the first n inital data will be considered
    if ( n0 > _n )
    {
        std::cout << "The initial data set is larger than needed by the BDF."
        << std::endl;
        std::cout << "Only the first " << _n << " data will be considered. "
        << std::endl;
    }
    for ( iter = _unk.begin() ; iter != iter_end; iter++ )
    {
        *iter = *iter0;
        iter0++;
    }

    return ;
}

//
//
//





const std::vector<Vector>& Bdf::unk() const
{
    return _unk;
}


double Bdf::coeff_der( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 & i < _n + 1, "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _alpha[ i ];
}

double Bdf::coeff_ext( UInt i ) const
{
    // Pay attention: i is c-based indexed
    ASSERT( i >= 0 & i < _n, "Error in specification of the time derivative coefficient for the BDF formula (out of range error)" );
    return _beta[ i ];
}

void Bdf::showMe() const
{
    std::cout << "*** BDF Time discretization of order " << _n << " ***"
              << std::endl;
    std::cout << "    Coefficients: " << std::endl;
    for ( UInt i = 0;i < _n + 1;++i )
        std::cout << "       alpha(" << i << ") = " << _alpha[ i ]
                  << std::endl;
    for ( UInt i = 0;i < _n;++i )
        std::cout << "       beta (" << i << ") = " << _beta[ i ]
                  << std::endl;

    std::cout << "    " << _unk.size() << " unknown vector" 
              << (_unk.size() == 1 ? "" : "s") << " of length "
              << _s << std::endl;

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

void Bdf::shift_right( Vector unk_curr )
{
    std::vector< Vector >::iterator it = _unk.end() - 1;
    std::vector< Vector >::iterator itm1 = _unk.end() - 1;
    std::vector< Vector >::iterator itb = _unk.begin();

    for ( ; it != itb; --it )
    {
        itm1--;
        *it = *itm1;
    }
    *itb = unk_curr;

    return ;
}


Vector Bdf::time_der( Real dt ) const
{
    Vector ut( _s );

    for ( UInt j = 0;j < _s;++j )
        ut[ j ] = 0.;

    for ( UInt j = 0;j < _s;++j )
    {
        for ( UInt i = 0;i < _n;++i )
            ut[ j ] = ut[ j ] + _alpha[ i + 1 ] / dt * _unk[ i ][ j ];
    }

    return ut;
}

Vector Bdf::time_der() const
{
    Vector ut( _s );

    for ( UInt j = 0;j < _s;++j )
        ut[ j ] = 0.;

    for ( UInt j = 0;j < _s;++j )
    {
        for ( UInt i = 0;i < _n;++i )
            ut[ j ] = ut[ j ] + _alpha[ i + 1 ] * _unk[ i ][ j ];
    }

    return ut;
}

Vector Bdf::extrap() const
{
    Vector ue( _s );
    for ( UInt j = 0;j < _s;++j )
        ue[ j ] = 0.;

    for ( UInt j = 0;j < _s;++j )
    {
        for ( UInt i = 0;i < _n;++i )
            ue[ j ] = ue[ j ] + _beta[ i ] * _unk[ i ][ j ];
    }

    return ue;
}
}
