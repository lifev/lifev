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
