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
  @brief Small classes to manage list data strings

  @date 1-11-2002
  @author J. F. Gerbeau

  @contributor
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#include <stdexcept>
#include <sstream>

#include <lifev/core/util/StringData.hpp>
#include <lifev/core/util/LifeDebug.hpp>

namespace LifeV
{

// ===============
// Constructors
// ===============

StringData::StringData ( std::string str, Int val, std::string help ) :
    M_string ( str ), M_value ( val ), M_help ( help )
{}

StringDataList::StringDataList ( std::string title ) :
    M_title ( title )
{}

// ===============
// Public methods
// ===============

void StringDataList::add ( std::string str, Int val, std::string help )
{
    M_list.push_back ( StringData ( str, val, help ) );
}

void StringDataList::showMe ( std::ostream& c, bool val ) const
{
    c << M_title << " : " << std::endl;
    for ( std::vector<StringData>::const_iterator ds = M_list.begin();
            ds != M_list.end(); ++ds )
    {
        c << "   " << ds->string() << " : " << ds->help();
        if ( val )
        {
            c << " (" << ds->value() << ")";
        }
        c << std::endl;
    }
}

Int StringDataList::value ( const std::string& str ) const
{
    std::vector<StringData>::const_iterator ds = M_list.begin();
    while ( ds != M_list.end() )
    {
        if ( ds->string() == str )
        {
            return ds->value();
        }
        ++ds;
    };
    std::ostringstream exception;
    exception << "Error in " << LIFEV_FUNCINFO  << ": "
              << str << " is not in the list of possible choices for '"
              << M_title;
    throw std::invalid_argument ( exception.str() );
}

}
