#ifndef _NEW_SWITCHES_H_
#define _NEW_SWITCHES_H_
#include "lifeV.hpp" // only for ASSERTs
#include<string>
#include <iostream>
#include<map>
#include "util_string.hpp"
/*! A better Switches class */
//! I use standard constructor/destructors
class Switch: public std::map<string,bool>
{
  
public:
  typedef  map<string,bool>::iterator iterator;
  
  //! It sets the switch  It does NOT create a new switch (use create);
  //! Returns true if the switch s existed
  bool set(string const & a);
  bool set(const char *  a);
  //! It unsets the switch. It does NOT create a new switch (use create)
  //! Returns true if the switch s existed
  bool unset(string const & a);
  bool unset(const char *  a);
  //! It toggles the switch 
  //! Returns true if the switch s existed
  bool toggle(string const & a);
  bool toggle(const char *  a);
  //! Creates a switch (default is OFF)
  //! It sets the switch if it is there
  void create(string const & a, bool status=false);
  void create(const char *  a, bool status=false);

  ostream & showMe(bool verbose=false,ostream & out=cout) const;
  /*! Returns the status of the switch:
    true-> set  false-> unset or non existent.
    It does not distinguish between not existent and unset switches: use status if you want
    the full information */
  bool test(string const & a) const;
  bool test(const char *  a) const;
  /*! Returns a pair of bools : the first refers to the
    existence in the list of switches, the second on the status (true or false)
    If the first bool is false also the second is false */
  pair<bool,bool> status(string const & a) const;
  pair<bool,bool> status(const char *  a) const;
private:
};
#endif
