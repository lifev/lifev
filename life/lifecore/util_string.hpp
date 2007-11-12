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
#ifndef UTIL_STRING_H
#define UTIL_STRING_H
# include <iosfwd>
# include <iostream>

# include <string>
# include <cstring>
# include <cstdlib>
# include <cstdio>
# include <list>
# include <sstream>

namespace LifeV
{
/*! \file util_string.h
\brief Special structures for handling mesh faces and sides
\version 0.0 Experimental   5/2/00. Luca Formaggia
Some utilities for handling ascii files
*/

/*! It gets a the next line from std::istream
*/
std::istream & eatline( std::istream & s );
//!skip lines starting with '!%#;$'
std::istream & eat_comments( std::istream & s );
//!  gets next uncommented line
std::istream & next_good_line( std::istream & s, std::string & line );
/*!
    always return a std::string with len characters
      - if the s has more than len characters : keep only the first len
      - if the s has less than len characters : complete with c until len
*/
std::string& setStringLength( std::string& s, unsigned int len, char c );


//! extends atoi to STL std::strings (from Stroustrup)
int atoi( const std::string & s );

std::string operator+( const std::string & str, const int i );
std::string operator+( const std::string & str, const long int i );
std::string operator+( const std::string & str, const unsigned int i );

template <typename T>
void parseList( const std::string& slist, std::list<T>& list )
{
    std::string stringList = slist;
//     std::set<T> setList;
    if ( slist == "" )
    {
        return;
    }

    int commaPos = 0;

    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );

        std::stringstream stream;
        stream <<  stringList.substr( 0, commaPos ).c_str();

        T var;
        stream >> var;
        list.push_back( var );

        stringList = stringList.substr( commaPos + 1 );
    }

}

}
#endif
