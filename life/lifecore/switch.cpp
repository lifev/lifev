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
