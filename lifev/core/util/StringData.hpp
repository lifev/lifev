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

#ifndef STRING_DATA_H
#define STRING_DATA_H

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

/*!
  @class StringData

  A data string contains a string (the name used
  in a data file), a integer value (the value used inside
  the code), a small help text
*/
class StringData
{
public:
    //! @name Constructor
    //@{
    StringData() {}
    StringData ( std::string str, Int val, std::string help );

    virtual ~StringData() {}
    //@}

    //! @name Get methods
    //@{
    inline const std::string& help() const
    {
        return M_help;
    }

    inline std::string& help()
    {
        return M_help;
    }

    inline const std::string& string () const
    {
        return M_string;
    }

    inline std::string& string ()
    {
        return M_string;
    }

    inline const Int& value () const
    {
        return M_value;
    }

    inline Int& value ()
    {
        return M_value;
    }
    //@}

private:
    //! @name Private members
    //@{
    std::string M_string;
    Int M_value;
    std::string M_help;
    //@}
};

/*!
  @class StringDataList

  To build a list of data string.

  Example:

  You want to propose to the user three solvers: cg, gmres, cgs
  inside the code these solvers are identified by the values 1, 2, 3

  You first create the list:

  StringDataList solver_list("My solvers");
  solver_list.add("cg",1,"preconditioned conjugate gradient method");
  solver_list.add("gmres",2,"preconditioned gmres method");
  solver_list.add("cgs",2,"preconditioned cg squared method");

  Next, you read (for example) a GetPot data file, where the user
  gives you a string string_string. To get the value corresponding
  the string provided by the user:

  solver = aztec_solver_list.value(user_string);

  if the "user_string" belongs to the list, then solver has the corresponding
  value,
  else, the list of possible choices (with help) is given and
  the code stop with an error.
*/
class StringDataList
{
public:
    //! @name Constructor
    //@{
    StringDataList() {}
    StringDataList ( std::string title );

    virtual ~StringDataList() {}
    //@}

    //! @name Public methods
    //@{
    void add ( std::string str, Int val, std::string help );

    //! State printing function (val=true:the values are shown)
    void showMe ( std::ostream& c = std::cout, bool val = false ) const;

    Int value ( const std::string& str ) const;
    //@}

private:
    //! @name Private data members
    //@{
    std::string M_title;
    std::vector<StringData> M_list;
    //@}
};
}
#endif // STRING_DATA_H
