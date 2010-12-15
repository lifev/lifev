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
  @brief Switches class

  @date 13-12-2010
  @author

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#include <life/lifecore/life.hpp>
#include <life/lifecore/switch.hpp>

namespace LifeV
{

// =======================
// Public methods
// =======================

bool Switch::set ( std::string const & a )
{
    iterator i = find( a );
    if ( i == end() )
    {
        return false;
    }

    else
    {
        i->second = true;
        return true;
    }
}

bool Switch::set ( const char * a )
{
    std::string temp( a );
    return set
           ( temp );
}

bool Switch::unset( std::string const & a )
{
    iterator i = find( a );
    if ( i == end() )
    {
        return false;
    }

    else
    {
        i->second = false;
        return true;
    }
}

bool Switch::unset( const char * a )
{
    std::string temp( a );
    return unset( temp );

}

bool Switch::toggle( std::string const & a )
{
    iterator i = find( a );
    if ( i == end() )
    {
        return false;
    }

    else
    {
        i->second = ! ( i->second );
        return true;
    }
}

bool Switch::toggle( const char * a )
{
    std::string temp( a );
    return toggle( temp );

}

void Switch::create( std::string const & a, bool status )
{
    iterator i = find( a );
    if ( i == end() )
    {
        insert( std::make_pair( a, status ) );
    }

    else
    {
        i->second = status;
    }
}

void Switch::create( const char * a, bool status )
{
    std::string temp( a );
    create( temp, status );
}


std::pair<bool, bool> Switch::status( std::string const & a ) const
{
    const_iterator i = find( a );
    if ( i == end() )
    {
        return std::make_pair( false, false );
    }

    else
    {
        return std::make_pair( true, i->second );
    }
}

std::pair<bool, bool> Switch::status( const char * a ) const
{
    std::string temp( a );
    return status( temp );
}


bool
Switch::test( std::string const & a ) const
{
    const_iterator i = find( a );
    if ( i == end() )
    {
        return false;
    }
    else
    {
        return i->second;
    }
}

bool Switch::test( const char * a ) const
{
    std::string temp( a );
    return test( temp );
}

std::ostream & Switch::showMe( bool verbose, std::ostream & out ) const
{
    out << std::endl << " Status of switches" << std::endl;
    for ( const_iterator i = begin(); i != end(); ++i )
    {
        out << "Switch named: " << i->first << " Value= " << i->second << std::endl;
    }
    return out;
}

}
