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
#include "util_string.hpp"

istream & eatline(istream & s)  // eat a whole line from istream
{
  while ( s . get() != '\n' && s . good() ) {}
  return s ;
}

istream & eat_comments(istream & s)  //eat lines starting with '!%#;$'
{
  char c ;
  s . get(c) ;
  while ( c == '!' ||
	  c == '%' ||
	  c == '#' ||
	  c == ';' ||
	  c == '$' ) {
    s >> eatline ;
    s . get(c) ; 
  }
  return s . putback(c) ;
}

string& setStringLength(string& s,unsigned int len,char c)
{
  /*
    always return a string with len characters
      - if the s has more than len characters : keep only the first len
      - if the s has less than len characters : complete with c until len
  */
  string stmp(len,c);
  if(s.length()>len){
    s.erase(len,s.length());
  }
  stmp.replace(0,s.length(),s);
  s = stmp;
  return s;
}

int atoi(const string & s)
{
  return atoi(s.c_str());
}

istream & next_good_line(istream & s, string & line) 
{
  s>>eat_comments;
  getline(s,line);
  return s;
}

string operator+(const string & str, const int i)
{
  int digits = abs(i)/10;
  char* str_i=new char[digits];
  sprintf(str_i,"%i",i);
  string str2=str+str_i;
  delete[] str_i;
  return str2;
}

string operator+(const string & str, const long i)
{
  int digits = abs(i)/10;
  char* str_i=new char[digits];
  sprintf(str_i,"%ld",i);
  string str2=str+str_i;
  delete[] str_i;
  return str2;
}

string operator+(const string & str, const unsigned int i)
{
  int digits = i/10;
  char* str_i=new char[digits];
  sprintf(str_i,"%u",i);
  string str2=str+str_i;
  delete[] str_i;
  return str2;
}
