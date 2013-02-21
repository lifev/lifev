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
  @brief Switches class

  @date 13-12-2010
  @author

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef SWITCH_H
#define SWITCH_H

#include<iostream>
#include<map>
#include<string>

#include <lifev/core/util/StringUtility.hpp>

namespace LifeV
{

/*! A better Switches class */
//! I use standard constructor/destructors
class Switch: public std::map<std::string, bool>
{

public:
    //! @name Public typedefs
    //@{
    typedef std::map<std::string, bool>::iterator iterator_Type;
    typedef std::map<std::string, bool>::const_iterator iteratorConst_Type;
    //@}

    //! @name Public methods
    //@{
    //! It sets the switch  It does NOT create a new switch (use create);
    //! Returns true if the switch s existed
    bool set ( std::string const& a );

    bool set ( const char* a );

    //! It unsets the switch. It does NOT create a new switch (use create)
    //! Returns true if the switch s existed
    bool unset ( std::string const& a );

    bool unset ( const char* a );

    //! It toggles the switch
    //! Returns true if the switch s existed
    bool toggle ( std::string const& a );

    bool toggle ( const char* a );

    //! Creates a switch (default is OFF)
    //! It sets the switch if it is there
    void create ( std::string const& a, bool status = false );
    void create ( const char* a, bool status = false );

    std::ostream& showMe ( bool verbose = false, std::ostream& out = std::cout ) const;

    /*! Returns the status of the switch:
      true-> set  false-> unset or non existent.
      It does not distinguish between not existent and unset switches: use status if you want
      the full information */
    bool test ( std::string const& a ) const;

    bool test ( const char* a ) const;

    /*! Returns a std::pair of bools : the first refers to the
      existence in the list of switches, the second on the status (true or false)
      If the first bool is false also the second is false */
    std::pair<bool, bool> status ( std::string const& a ) const;

    std::pair<bool, bool> status ( const char* a ) const;
    //@}
};

} // Namespace LifeV

#endif // SWITCH_H
