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
#include "util_string.hpp"

namespace LifeV
{
std::istream & eatline(std::istream & s)  // eat a whole line from std::istream
{
  while ( s . get() != '\n' && s . good() ) {}
  return s ;
}

std::istream & eat_comments(std::istream & s)  //eat lines starting with '!%#;$'
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

std::string& setStringLength(std::string& s,unsigned int len,char c)
{
  /*
    always return a std::string with len characters
      - if the s has more than len characters : keep only the first len
      - if the s has less than len characters : complete with c until len
  */
  std::string stmp(len,c);
  if(s.length()>len){
    s.erase(len,s.length());
  }
  stmp.replace(0,s.length(),s);
  s = stmp;
  return s;
}

int atoi(const std::string & s)
{
  return ::atoi(s.c_str());
}

std::istream & next_good_line(std::istream & s, std::string & line)
{
  s>>eat_comments;
  getline(s,line);
  return s;
}

std::string operator+(const std::string & str, const int i)
{
  int digits = abs(i)/10;
  char* str_i=new char[digits];
  sprintf(str_i,"%i",i);
  std::string str2=str+str_i;
  delete[] str_i;
  return str2;
}

std::string operator+(const std::string & str, const long i)
{
  int digits = abs(i)/10;
  char* str_i=new char[digits];
  sprintf(str_i,"%ld",i);
  std::string str2=str+str_i;
  delete[] str_i;
  return str2;
}

std::string operator+(const std::string & str, const unsigned int i)
{
  int digits = i/10;
  char* str_i=new char[digits];
  sprintf(str_i,"%u",i);
  std::string str2=str+str_i;
  delete[] str_i;
  return str2;
}
}
