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
#ifndef UTIL_STRING_H
#define UTIL_STRING_H
# include <iosfwd> 
# include <iostream> 

#if defined(__Linux)
# include <cstdlib> 
# include <cstring>
#elif defined(__OSF1)
# include <stdlib.h>
# include <string.h>
#endif

# include <string>
# include <stdlib.h>
# include <stdio.h>

using namespace std;

/*! \file util_string.h
\brief Special structures for handling mesh faces and sides
\version 0.0 Experimental   5/2/00. Luca Formaggia
Some utilities for handling ascii files
*/

/*! It gets a the next line from istream
*/
istream & eatline(istream & s);
//!skip lines starting with '!%#;$'
istream & eat_comments(istream & s);
//!  gets next uncommented line
istream & next_good_line(istream & s, string & line);
/*!
    always return a string with len characters
      - if the s has more than len characters : keep only the first len
      - if the s has less than len characters : complete with c until len
*/
string& setStringLength(string& s,unsigned int len,char c);


//! extends atoi to STL strings (from Stroustrup)
int atoi(const string & s);

string operator+(const string & str, const int i);
string operator+(const string & str, const long int i);
string operator+(const string & str, const unsigned int i);

#endif
