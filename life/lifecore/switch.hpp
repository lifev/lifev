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
#ifndef _NEW_SWITCHES_H_
#define _NEW_SWITCHES_H_
#include<string>
#include<iostream>
#include<map>

#include "util_string.hpp"

namespace LifeV
{
/*! A better Switches class */
//! I use standard constructor/destructors
class Switch: public std::map<std::string,bool>
{

public:
    typedef  std::map<std::string,bool>::iterator iterator;

    //! It sets the switch  It does NOT create a new switch (use create);
    //! Returns true if the switch s existed
    bool set(std::string const & a);
    bool set(const char *  a);
    //! It unsets the switch. It does NOT create a new switch (use create)
    //! Returns true if the switch s existed
    bool unset(std::string const & a);
    bool unset(const char *  a);
    //! It toggles the switch
    //! Returns true if the switch s existed
    bool toggle(std::string const & a);
    bool toggle(const char *  a);
    //! Creates a switch (default is OFF)
    //! It sets the switch if it is there
    void create(std::string const & a, bool status=false);
    void create(const char *  a, bool status=false);

    std::ostream & showMe(bool verbose=false, std::ostream & out = std::cout) const;
    /*! Returns the status of the switch:
      true-> set  false-> unset or non existent.
      It does not distinguish between not existent and unset switches: use status if you want
      the full information */
    bool test(std::string const & a) const;
    bool test(const char *  a) const;
    /*! Returns a std::pair of bools : the first refers to the
      existence in the list of switches, the second on the status (true or false)
      If the first bool is false also the second is false */
    std::pair<bool,bool> status(std::string const & a) const;
    std::pair<bool,bool> status(const char *  a) const;
private:
};
}
#endif
