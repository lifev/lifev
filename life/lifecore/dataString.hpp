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
#ifndef _DATA_STRING_H_
#define _DATA_STRING_H_
#include <iostream>
#include <vector>
#include <string>

namespace LifeV
{
/*!
  \file dataSting.h
  \author J.-F. Gerbeau
  \date 11/2002
  \brief Small classes to manage list data strings
*/

/*!
  \class DataString

  A data string contains a string (the name used
  in a data file), a integer value (the value used inside
  the code), a small help text
*/
class DataString
{
  std::string _str;
  int    _val;
  std::string _help;
public:
  DataString(std::string str,int val,std::string help);
  inline const std::string& str () const {return _str;};
  inline const int&    val () const {return _val;};
  inline const std::string& help() const {return _help;};
  inline  std::string& str ()  {return _str;};
  inline  int&    val ()  {return _val;};
  inline  std::string& help()  {return _help;};
};

/*!
  \class DataStringList

  To build a list of data string.

  Example:

  You want to propose to the user three solvers: cg, gmres, cgs
  inside the code these solvers are identified by the values 1, 2, 3

  You first create the list:

  DataStringList solver_list("My solvers");
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
class DataStringList
{
    std::string _title;
    std::vector<DataString> _list;
public:
    DataStringList(std::string title);
    void add(std::string str,int val,std::string help);
    void showMe(std::ostream& c= std::cout,bool val=false) const;//!<val=true:the values are shown
    int  value(const std::string& str) const;
};
}
#endif
