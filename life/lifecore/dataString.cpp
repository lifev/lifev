/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003 LifeV Team
  
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
#include "dataString.hpp"

DataString::DataString(string str,int val,string help):
  _str(str),_val(val),_help(help)
{}

DataStringList::DataStringList(string title):
  _title(title)
{}

void DataStringList::add(string str,int val,string help)
{
  _list.push_back(DataString(str,val,help));
}

void DataStringList::showMe(ostream& c,bool val) const
{
  c << _title << " : " << endl;
  for(vector<DataString>::const_iterator ds = _list.begin();ds != _list.end();ds++){
    c << "   "  << ds->str() <<  " : " << ds->help();
    if(val) c << " (" << ds->val() << ")";
    c << endl;
  }  
}

int DataStringList::value(const string& str) const
{
  vector<DataString>::const_iterator ds = _list.begin();
  while(ds != _list.end()){
    if(ds->str() == str) return ds->val();
    ds++;
  };
  cout << "Error : " << str <<
    " is not in the list of possible choices for '" << _title << "'\n";
  showMe();
  exit(1);
  return 0;
}

