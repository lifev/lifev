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

#include "dataString.hpp"

namespace LifeV
{
DataString::DataString(std::string str,int val,std::string help):
  _str(str),_val(val),_help(help)
{}

DataStringList::DataStringList(std::string title):
  _title(title)
{}

void DataStringList::add(std::string str,int val,std::string help)
{
  _list.push_back(DataString(str,val,help));
}

void DataStringList::showMe(std::ostream& c,bool val) const
{
    c << _title << " : " << std::endl;
    for( std::vector<DataString>::const_iterator ds = _list.begin();
         ds != _list.end();ds++)
    {
        c << "   "  << ds->str() <<  " : " << ds->help();
        if(val) c << " (" << ds->val() << ")";
        c << std::endl;
    }
}

int DataStringList::value(const std::string& str) const
{
    std::vector<DataString>::const_iterator ds = _list.begin();
    while(ds != _list.end())
    {
        if(ds->str() == str) return ds->val();
        ds++;
    };
    std::ostringstream __ex;
    __ex << "Error in " << __PRETTY_FUNCTION__ << ": " << str << " is not in the list of possible choices for '" << _title;
    throw std::invalid_argument( __ex.str() );
}
}
