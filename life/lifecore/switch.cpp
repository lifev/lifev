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
#include "switch.hpp"

bool
Switch::set(string const & a)
{
  iterator i=find(a);
  if (i == end())
    {
      return false;
    }
  
  else{
    i->second=true;
    return true;
  }
}

bool
Switch::set(const char *  a)
{
  string temp(a);
  return set(temp);
}

	      


bool
Switch::unset(string const & a)
{
  iterator i=find(a);
  if (i == end())
    {
      return false;
    }
  
  else{
    i->second=false;
    return true;
  }
}

bool
Switch::unset(const char *  a)
{
  string temp(a);
  return unset(temp);
  
}


bool
Switch::toggle(string const & a)
{
  iterator i=find(a);
  if (i == end())
    {
      return false;
    }
  
  else{
    i->second=! (i->second);
    return true;
  }
}

bool
Switch::toggle(const char *  a)
{
  string temp(a);
  return toggle(temp);
  
}


void
Switch::create(string const & a, bool status)
{
  iterator i=find(a);
  if (i == end())
    {
      insert(make_pair(a,status));
    }
  
  else{
    i->second=status;
  }
}

void
Switch::create(const char *  a, bool status)
{
  string temp(a);
  create(temp,status);
}


pair<bool,bool>
Switch::status(string const & a) const
{
  const_iterator i=find(a);
  if (i == end())
    {
      return make_pair(false,false);
    }
  
  else{
    return make_pair(true,i->second);
  }
}

pair<bool,bool>
Switch::status(const char *  a) const
{
  string temp(a);
  return status(temp);
}

  
bool
Switch::test(string const & a) const
{
  const_iterator i=find(a);
  if (i == end())
    {
      return false;
    }
  else{
    return i->second;
  }
}

  bool
Switch::test(const char *  a) const
{
  string temp(a);
  return test(temp);
}

ostream & Switch::showMe(bool verbose, ostream & out) const
{
  out<<endl<< " Status of switches"<<endl;
  for (const_iterator i=begin(); i != end(); ++i) {
    out<< "Switch named: " <<i->first<<" Value= "<<i->second<<endl;
  }
  return out;
  
}
