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
#include "life.hpp" // only for ASSERTs

#include "switch.hpp"

namespace LifeV
{
bool
Switch::set
    ( std::string const & a )
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

bool
Switch::set
    ( const char * a )
{
    std::string temp( a );
    return set
               ( temp );
}




bool
Switch::unset( std::string const & a )
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

bool
Switch::unset( const char * a )
{
    std::string temp( a );
    return unset( temp );

}


bool
Switch::toggle( std::string const & a )
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

bool
Switch::toggle( const char * a )
{
    std::string temp( a );
    return toggle( temp );

}


void
Switch::create( std::string const & a, bool status )
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

void
Switch::create( const char * a, bool status )
{
    std::string temp( a );
    create( temp, status );
}


std::pair<bool, bool>
Switch::status( std::string const & a ) const
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

std::pair<bool, bool>
Switch::status( const char * a ) const
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

bool
Switch::test( const char * a ) const
{
    std::string temp( a );
    return test( temp );
}

std::ostream & Switch::showMe( bool /*verbose*/, std::ostream & out ) const
{
    out << std::endl << " Status of switches" << std::endl;
    for ( const_iterator i = begin(); i != end(); ++i )
    {
        out << "Switch named: " << i->first << " Value= " << i->second << std::endl;
    }
    return out;

}
}
