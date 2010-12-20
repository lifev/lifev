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
  @brief String utilities

  @date 13-12-2010
  @author

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#include <life/lifecore/StringUtility.hpp>

namespace LifeV
{
std::istream & eatLine( std::istream & s )    // eat a whole line from std::istream
{
    while ( s.get() != '\n' && s . good() )
        {}
    return s ;
}

std::istream & eatComments( std::istream & s )    //eat lines starting with '!%#;$'
{
    char c = 'a';
    s.get( c ) ;
    while ( c == '!' ||
            c == '%' ||
            c == '#' ||
            c == ';' ||
            c == '$' )
    {
        s >> eatLine ;
        s.get( c ) ;
    }
    return s.putback( c ) ;
}

std::istream & nextGoodLine( std::istream & s, std::string & line )
{
    s >> eatComments;
    getline( s, line );
    return s;
}

std::string& setStringLength( std::string& s, unsigned int len, char c )
{
    /*
      always return a std::string with len characters
        - if the s has more than len characters : keep only the first len
        - if the s has less than len characters : complete with c until len
    */
    std::string stmp( len, c );
    if ( s.length() > len )
    {
        s.erase( len, s.length() );
    }
    stmp.replace( 0, s.length(), s );
    s = stmp;
    return s;
}

int atoi( const std::string & s )
{
    return ::atoi( s.c_str() );
}

std::string operator+( const std::string & str, const int i )
{
    int digits = abs( i ) / 10;
    char* str_i = new char[ digits ];
    sprintf( str_i, "%i", i );
    std::string str2 = str + str_i;
    delete[] str_i;
    return str2;
}

std::string operator+( const std::string & str, const long i )
{
    int digits = abs( i ) / 10;
    char* str_i = new char[ digits ];
    sprintf( str_i, "%ld", i );
    std::string str2 = str + str_i;
    delete[] str_i;
    return str2;
}

std::string operator+( const std::string & str, const unsigned int i )
{
    int digits = i / 10;
    char* str_i = new char[ digits ];
    sprintf( str_i, "%u", i );
    std::string str2 = str + str_i;
    delete[] str_i;
    return str2;
}
}
