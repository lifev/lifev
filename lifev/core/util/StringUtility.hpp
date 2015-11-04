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

#ifndef STRING_UTILITY_H
#define STRING_UTILITY_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iosfwd>
#include <sstream>


#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include <lifev/core/LifeV.hpp>

namespace LifeV
{
/*! \file util_string.h
\brief Special structures for handling mesh faces and sides
\version 0.0 Experimental   5/2/00. Luca Formaggia
Some utilities for handling ascii files
*/

/*! It gets a the next line from std::istream
*/
std::istream& eatLine ( std::istream& s );
//!skip lines starting with '!%#;$'
std::istream& eatComments ( std::istream& s );
//!  gets next uncommented line
std::istream& nextGoodLine ( std::istream& s, std::string& line );
/*!
    always return a std::string with len characters
      - if the s has more than len characters : keep only the first len
      - if the s has less than len characters : complete with c until len
*/
std::string& setStringLength ( std::string& s, unsigned int len, char c );


//! extends atoi to STL std::strings (from Stroustrup)
int atoi ( const std::string& s );

std::string operator+ ( const std::string& str, const int i );
std::string operator+ ( const std::string& str, const long int i );
std::string operator+ ( const std::string& str, const unsigned int i );

template <typename EntryType>
void parseList ( const std::string& slist, std::list<EntryType>& list )
{
    std::string stringList = slist;
    if ( slist == "" )
    {
        return;
    }

    int commaPos = 0;

    while ( commaPos != (int) std::string::npos )
    {
        commaPos = stringList.find ( "," );

        std::stringstream stream;
        stream <<  stringList.substr ( 0, commaPos ).c_str();

        EntryType var;
        stream >> var;
        list.push_back ( var );

        stringList = stringList.substr ( commaPos + 1 );
    }
}

// @author Cristiano Malossi
// Convert a std::string to a number ( Int, bool, Real, ... )
inline Real string2number ( const std::string& s )
{
    std::stringstream out;
    out << s;

    Real n;
    out >> n;

    return n;

    // Temporary disabled
    //return std::to_string( s );
}

// @author Cristiano Malossi
// Convert a number ( Int, bool, Real, ... ) to a std::string
template <typename NumberType>
inline std::string number2string ( const NumberType& n )
{
    return std::to_string ( n );
}

// @author Cristiano Malossi
// Convert an Enum to a std::string using a map as a library for conversion
template < typename EnumeratorType >
inline std::string enum2String ( const EnumeratorType& Enum,
                                 const std::map < std::string,
                                 EnumeratorType > & Map )
{
    for ( typename std::map<std::string, EnumeratorType>::const_iterator j = Map.begin(); j != Map.end() ; ++j )
        if ( j->second == Enum )
        {
            return j->first;
        }

    return "NO_TYPE_FOUND";
}

// @author Cristiano Malossi
// Convert a string made by NumberTypes separated by commas, to a vector of numbers
template< typename NumberType >
void string2numbersVector ( const std::string& string,
                            std::vector< NumberType >& numberVector )
{
    //Split the string
    std::vector< std::string > stringVector;
    boost::split ( stringVector, string, boost::is_any_of ( "," ) );

    //Convert to the right type
    for ( UInt i ( 0 ); i < static_cast<UInt> ( stringVector.size() ); ++i )
    {
        numberVector.push_back ( static_cast<NumberType> ( std::atoi ( stringVector[i].c_str() ) ) );
    }
}

} // Namespace LifeV

#endif //  STRING_UTILITY_H
