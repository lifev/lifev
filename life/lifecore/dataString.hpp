#ifndef _DATA_STRING_H_
#define _DATA_STRING_H_
#include <iostream>
#include <vector>
#include <string>

using namespace std;
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
  string _str;
  int    _val;
  string _help;
public:
  DataString(string str,int val,string help);
  inline const string& str () const {return _str;};
  inline const int&    val () const {return _val;};
  inline const string& help() const {return _help;};
  inline  string& str ()  {return _str;};
  inline  int&    val ()  {return _val;};
  inline  string& help()  {return _help;};
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
  string _title;
  vector<DataString> _list;
public:
  DataStringList::DataStringList(string title);
  void add(string str,int val,string help);
  void showMe(ostream& c=cout,bool val=false) const;//!<val=true:the values are shown
  int  value(const string& str) const;
};

#endif
