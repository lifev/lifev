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
//  -*- c++ -*- 
//  GetPot Version 1.0                                        Sept/13/2002
//
//  WEBSITE: http://getpot.sourceforge.net
//  
//  This library is  free software; you can redistribute  it and/or modify
//  it  under  the terms  of  the GNU  Lesser  General  Public License  as
//  published by the  Free Software Foundation; either version  2.1 of the
//  License, or (at your option) any later version.
//  
//  This library  is distributed in the  hope that it will  be useful, but
//  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
//  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
//  Lesser General Public License for more details.
//  
//  You  should have  received a  copy of  the GNU  Lesser  General Public
//  License along  with this library; if  not, write to  the Free Software
//  Foundation, Inc.,  59 Temple Place,  Suite 330, Boston,  MA 02111-1307
//  USA
//  
//  (C) 2001-2002 Frank R. Schaefer  
//==========================================================================
#include <cctype> 
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <cstring>

#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>

#include "GetPot.hpp"
using getpot::StringVector;


#define victorate(TYPE, VARIABLE, ITERATOR) \
	for(std::vector<TYPE>::const_iterator (ITERATOR) = (VARIABLE).begin(); \
	    (ITERATOR) != (VARIABLE).end(); ++(ITERATOR))

// Quick fix.
#define inline

///////////////////////////////////////////////////////////////////////////////
// (*) constructors, destructor, assignment operator
//.............................................................................
//
inline 
GetPot::GetPot()
  : cursor(0), search_loop_f(true),  search_failed_f(true), nominus_cursor(-1)
{ }

inline 
GetPot::GetPot(int argc_, char *argv_[])
  : cursor(0), search_loop_f(true),  search_failed_f(true), nominus_cursor(-1)
{ 
    // make an internal copy of the argument list:
    StringVector __argv;
    __argv.push_back(std::string(argv_[0]));
    for(int i=1; i<argc_; i++) {
	std::string tmp(argv_[i]); // recall the problem with temporaries,
	__argv.push_back(tmp);  // reference counting in arguement lists ...
    }
    __parse_argument_vector(__argv); 
}


inline 
GetPot::GetPot(const char* FileName)
  : cursor(0), search_loop_f(true),  search_failed_f(true), nominus_cursor(-1)
{ 
    StringVector __argv;
    __argv.push_back(FileName); 
    StringVector args = read_in_file(FileName);
    __argv.insert(__argv.begin()+1, args.begin(), args.end());
    __parse_argument_vector(__argv); 
}
// begin luca modif 29 Dec 2002
inline 
GetPot::GetPot(const StringVector & FileNameList)
  : cursor(0), search_loop_f(true),  search_failed_f(true), nominus_cursor(-1)
{ 
    StringVector __argv;
    __argv.push_back(FileNameList[0]);
    for (StringVector::const_iterator s=FileNameList.begin(); s != FileNameList.end(); ++s)
      {
	StringVector args = read_in_file(s->c_str());
	__argv.insert(__argv.end(), args.begin(), args.end());
      }
    __parse_argument_vector(__argv); 
}
// end luca modif
inline 
GetPot::GetPot(const GetPot& Other)
{ GetPot::operator=(Other); }

inline 
GetPot::~GetPot()
{ }

inline GetPot&
GetPot::operator=(const GetPot& Other)
{
    if (&Other != this) {
	argv            = Other.argv;
	variables       = Other.variables;
	cursor          = Other.cursor;
	nominus_cursor  = Other.nominus_cursor;
	search_failed_f = Other.search_failed_f;
	idx_nominus     = Other.idx_nominus;
	search_loop_f   = Other.search_loop_f;
    }
    return *this;
}


inline void
GetPot::__parse_argument_vector(const StringVector& ARGV)
{
    // build internal databases:
    //   1) array with no-minus arguments (usually used as filenames)
    //   2) variable assignments:
    //             'variable name' '=' number | string
    section = "";
    StringVector  section_stack;
    unsigned i=0;
    StringVector::const_iterator it = ARGV.begin();
    argv.push_back(*it);
    for(it++; it != ARGV.end(); it++, i++) {
	std::string arg = *it;
	if( arg.length() == 0 ) continue;

	// -- [section] labels
	if( arg.length() > 1 && arg[0] == '[' && arg[arg.length()-1] == ']' ) {
	    const std::string Name = DBE_expand_string(arg.substr(1, arg.length()-2));
	    section = process_section_label(Name, section_stack);
	    // new section --> append to list of sections
	    if( find(section_list.begin(), section_list.end(), section) == section_list.end() )
		if( section.length() != 0 ) section_list.push_back(section);
	    argv.push_back(arg);
	}
	else {
	    arg = section + DBE_expand_string(arg);
	    argv.push_back(arg);
	}
	    
	// -- separate array for nominus arguments
	if( arg[0] != '-' ) idx_nominus.push_back(unsigned(i));

	// -- variables: does arg contain a '=' operator ?
	for(const char* p = arg.c_str(); *p ; p++) {
	    if( *p == '=' ) {
		// set terminating 'zero' to treat first part as single string
		// => arg (from start to 'p') = Name of variable
		//    p+1     (until terminating zero) = value of variable
		char* o = const_cast<char*>(p++);
		*o = '\0';
		const GetPot::variable* Var = find_variable(arg.c_str());
		if( Var == 0 ) variables.push_back(variable(arg.c_str(), p)); 
		else           (const_cast<GetPot::variable*>(Var))->take(p);
		*o = '=';
		break;
	    }
	}
    }
}
inline StringVector
GetPot::read_in_file(const char* FileName)
{
    std::ifstream  i(FileName);
    if( ! i ) return StringVector();
    // argv[0] == the filename of the file that was read in
    return read_in_stream(i);
}

inline StringVector
GetPot::read_in_stream(std::istream& istr)
{
    StringVector  brute_tokens;
    while(istr) {
	__skip_whitespace(istr);
	const std::string Token = __get_next_token(istr);
	// g++ 2.96 does not have template char_trait<char>
	//
#ifdef __GNUG__ 
	int tmpfix = Token[0];
	if( Token.length() == 0 || tmpfix == EOF) break;
#else
	if( Token.length() == 0 ||
	    std::char_traits<char>::eq_int_type(Token[0],  std::char_traits<char>::eof()) ) break;
#endif
	brute_tokens.push_back(Token);
    }

    // -- reduce expressions of token1'='token2 to a single 
    //    string 'token1=token2'
    // -- copy everything into 'argv'
    // -- arguments preceded by something like '[' name ']' (section)
    //    produce a second copy of each argument with a prefix '[name]argument'
    unsigned i1 = 0;
    unsigned i2 = 1;
    unsigned i3 = 2;
    
    StringVector  arglist;
    while( i1 < brute_tokens.size() ) {
	const std::string& SRef = brute_tokens[i1];
	// 1) concatinate 'abcdef' '=' 'efgasdef' to 'abcdef=efgasdef'
	// note: java.lang.String: substring(a,b) = from a to b-1
	//        C++ string:      substr(a,b)    = from a to a + b
	if( i2 < brute_tokens.size() && brute_tokens[i2] == "=" ) {
	    if( i3 >= brute_tokens.size() )
		arglist.push_back(brute_tokens[i1] + brute_tokens[i2]);
	    else
		arglist.push_back(brute_tokens[i1] + brute_tokens[i2] + brute_tokens[i3]);
	    i1 = i3+1; i2 = i3+2; i3 = i3+3;
	    continue;
	}
	else {
	    arglist.push_back(SRef);
	    i1=i2; i2=i3; i3++;
	}
    }
    return arglist;
}

inline void
GetPot::__skip_whitespace(std::istream& istr)
  // find next non-whitespace while deleting comments
{ 
  do {
    char tmp;
    istr >> tmp;    
    // found a non whitespace 
    if( tmp == '#' ) {
      // comment => skip until end of line
      while( tmp != '\n' ) {
	if( ! istr ) { istr.putback(tmp); return; }
	istr.get(tmp);
      }
      continue; 
    }
    else {
      istr.putback(tmp);
      return;
    }
  } while( istr );
}

inline const std::string
GetPot::__get_next_token(std::istream& istr)
    // get next concatinates string token. consider quotes that embrace
    // whitespaces
{
    std::string token;
    int    tmp = 0; 
    int    last_letter = 0;
    while(1+1 == 2) {
	last_letter = tmp; tmp = istr.get();
	if( tmp == EOF 
	    || ((tmp == ' ' || tmp == '\t' || tmp == '\n') && last_letter != '\\') ) {
	    return token;
	}
	else if( tmp == '\'' && last_letter != '\\' ) {
           // QUOTES: un-backslashed quotes => it's a string
	    token += __get_string(istr);
	    continue;
	}
	else if( tmp == '{' && last_letter == '$') {
	    token += '{' + __get_until_closing_bracket(istr);
	    continue;
	}
	else if( tmp == '$' && last_letter == '\\') {
	    token += tmp; tmp = 0;  //  so that last_letter will become = 0, not '$';
	    continue;
	}
	else if( tmp == '\\' && last_letter != '\\') 
	    continue;              // don't append un-backslashed backslashes
	token += tmp;
    }
}

inline const std::string
GetPot::__get_string(std::istream& istr)
    // parse input until next matching '
{
    std::string str;
    int    tmp = 0;
    int    last_letter = 0;
    while(1 + 1 == 2) {
	last_letter = tmp; tmp = istr.get();
	if( tmp == EOF)  return str;
	// un-backslashed quotes => it's the end of the string
	else if( tmp == '\'' && last_letter != '\\')  return str;
	else if( tmp == '\\' && last_letter != '\\')  continue; // don't append 

	str += tmp;
    }
}

inline const std::string
GetPot::__get_until_closing_bracket(std::istream& istr)
    // parse input until next matching }
{
    std::string str = ""; 
    int    tmp = 0;
    int    last_letter = 0;
    int    brackets = 1;
    while(1 + 1 == 2) {
	last_letter = tmp; tmp = istr.get();
	if( tmp == EOF) return str;
	else if( tmp == '{' && last_letter == '$') brackets += 1;
	else if( tmp == '}') {
	    brackets -= 1;
            // un-backslashed brackets => it's the end of the string
	    if( brackets == 0) return str + '}';
	    else if( tmp == '\\' && last_letter != '\\') 
                continue;  // do not append an unbackslashed backslash
	}
	str += tmp;
    }
}

inline std::string
GetPot::process_section_label(const std::string& Section, 
			      StringVector& section_stack)
{
    std::string sname = Section;
    //  1) subsection of actual section ('./' prefix)
    if( sname.length() >= 2 && sname.substr(0, 2) == "./" ) {
	sname = sname.substr(2);
    }
    //  2) subsection of parent section ('../' prefix)
    else if( sname.length() >= 3 && sname.substr(0, 3) == "../" ) {
	do {
	    if( section_stack.end() != section_stack.begin() ) 
		section_stack.pop_back();
	    sname = sname.substr(3);
	} while( sname.substr(0, 3) == "../" );
    }
    // 3) subsection of the root-section
    else {
	section_stack.erase(section_stack.begin(), section_stack.end());
	// [] => back to root section
    }

    if( sname != "" ) {
	// parse section name for 'slashes'
	unsigned i=0;
	while( i < sname.length() ) {
	    if( sname[i] == '/' ) {
		section_stack.push_back(sname.substr(0,i));
		if( i+1 < sname.length()-1 )
		    sname = sname.substr(i+1);
		i = 0;
	    }
	    else 
		i++;
	}	   
	section_stack.push_back(sname);
    } 
    std::string section = "";
    if( section_stack.size() != 0 )
	victorate(std::string, section_stack, it) {
	    section += *it + "/";
	}
    return section;
}
// convert string to DOUBLE, if not possible return Default
inline double   
GetPot::__convert_to_type(const std::string& String, double Default) const
{
    double tmp;
    if( sscanf(String.c_str(),"%lf", &tmp) != 1 ) return Default;
    return tmp;
}

// convert string to INT, if not possible return Default
inline int   
GetPot::__convert_to_type(const std::string& String, int Default) const
{
    int tmp;
    if( sscanf(String.c_str(),"%i", &tmp) != 1 ) return Default;
    return tmp;    
}

//////////////////////////////////////////////////////////////////////////////
// (*) cursor oriented functions
//.............................................................................
const std::string  
GetPot::__get_remaining_string(const std::string& String, const std::string& Start) const
    // Checks if 'String' begins with 'Start' and returns the remaining String.
    // Returns None if String does not begin with Start.
{
    if( Start == "" ) return String;
    // note: java.lang.String: substring(a,b) = from a to b-1
    //        C++ string:      substr(a,b)    = from a to a + b
    if( String.find(Start) == 0 ) return String.substr(Start.length());
    else                          return "";
}

//     -- search for a certain argument and set cursor to position
inline bool
GetPot::search(const char* Option)
{ 
    unsigned      OldCursor = cursor;
    const std::string  SearchTerm = prefix + Option;

    if( OldCursor >= argv.size() ) OldCursor = argv.size() - 1;
    search_failed_f = true;
    
    // (*) first loop from cursor position until end
    unsigned  c = cursor;
    for(; c < argv.size(); c++) {
	if( argv[c] == SearchTerm ) 
	{ cursor = c; search_failed_f = false; return true; }
    }
    if( ! search_loop_f ) return false;
	
    // (*) second loop from 0 to old cursor position
    for(c = 1; c < OldCursor; c++) {
	if( argv[c] == SearchTerm ) 
	{ cursor = c; search_failed_f = false; return true; }
    }
    // in case nothing is found the cursor stays where it was
    return false;
}


inline bool
GetPot::search(unsigned No, const char* P, ...)
{
    if( No == 0 ) return false;

    // search for the first argument
    if( search(P) == true ) return true;
    
    // start interpreting variable argument list
    va_list ap;
    va_start(ap, P);
    for(unsigned i=1; i<No; i++) {
	char* Opt = va_arg(ap, char *);
	if( search(Opt) == true ) 
	{ va_end(ap); return true; }
    }
    va_end(ap);
    return false;
}

inline void  
GetPot::reset_cursor()
{ search_failed_f = false; cursor = 0; }

inline void         
GetPot::init_multiple_occurrence()
{ disable_loop(); reset_cursor(); }
///////////////////////////////////////////////////////////////////////////////
// (*) direct access to command line arguments
//.............................................................................
//
inline const char*   
GetPot::operator[](unsigned idx) const
{ return idx<argv.size() ? argv[idx].c_str() : 0; }

inline int   
GetPot::get(unsigned Idx, int Default) const
{ 
    if( Idx >= argv.size() ) return Default;
    return __convert_to_type(argv[Idx], Default);
}

inline double   
GetPot::get(unsigned Idx, const double& Default) const
{ 
    if( Idx >= argv.size() ) return Default;
    return __convert_to_type(argv[Idx], Default);
}

inline const char*   
GetPot::get(unsigned Idx, const char* Default) const
{ 
    if( Idx >= argv.size() ) return Default;
    else                     return argv[Idx].c_str();
}

inline unsigned 
GetPot::size() const
{ return argv.size(); }
//     -- next() function group
inline int  
GetPot::next(int Default)    
{ 
    if( search_failed_f ) return Default;
    cursor++; 
    if( cursor >= argv.size() ) 
    { cursor = argv.size(); return Default; }

    const std::string Remain = __get_remaining_string(argv[cursor], prefix);
	
    return Remain != "" ? __convert_to_type(Remain, Default) : Default;
}

inline double 
GetPot::next(const double& Default) 
{ 
    if( search_failed_f ) return Default;
    cursor++;

    if( cursor >= argv.size() ) 
    { cursor = argv.size(); return Default; }

    std::string Remain = __get_remaining_string(argv[cursor], prefix);
	
    return Remain != "" ? __convert_to_type(Remain, Default) : Default;
}

inline const char*   
GetPot::next(const char* Default) 
{ 
    if( search_failed_f ) return Default;
    cursor++; 

    if( cursor >= argv.size() )
    { cursor = argv.size(); return Default; }

    const std::string Remain = __get_remaining_string(argv[cursor], prefix);
	
    return Remain != "" ? Remain.c_str() : Default;
}

//     -- follow() function group 
//        distinct option to be searched for
inline int 
GetPot::follow(int Default, const char* Option) 
{ 
    if( search(Option) == false ) return Default;
    return next(Default);
}

inline double  
GetPot::follow(const double& Default, const char* Option)
{ 
    if( search(Option) == false ) return Default;
    return next(Default);
}

inline const char* 
GetPot::follow(const char* Default, const char* Option)
{ 
    if( search(Option) == false ) return Default;
    return next(Default);
}

//     -- second follow() function group 
//        multiple option to be searched for
inline int 
GetPot::follow(int Default, unsigned No, const char* P, ...) 
{ 
    if( No == 0 ) return Default;
    if( search(P) == true ) return next(Default);
    
    va_list ap;
    va_start(ap, P);
    for(unsigned i=1; i<No; i++) {
	char* Opt = va_arg(ap, char *);
	if( search(Opt) == true ) { 
	    va_end(ap); 
	    return next(Default); 
	}
    }
    va_end(ap); 
    return Default;
}

inline double  
GetPot::follow(const double& Default, unsigned No, const char* P, ...)
{
    if( No == 0 ) return Default;
    if( search(P) == true ) return next(Default);
    
    va_list ap;
    va_start(ap, P);
    for(unsigned i=1; i<No; i++) {
	char* Opt = va_arg(ap, char *);
	if( search(Opt) == true ) { 
	    va_end(ap); 
	    return next(Default); 
	}
    }
    va_end(ap); 
    return Default;
}

inline const char* 
GetPot::follow(const char* Default, unsigned No, const char* P, ...)
{
    if( No == 0 ) return Default;
    if( search(P) == true ) return next(Default);
    
    va_list ap;
    va_start(ap, P);
    for(unsigned i=1; i<No; i++) {
	char* Opt = va_arg(ap, char *);
	if( search(Opt) == true ) { 
	    va_end(ap); 
	    return next(Default); 
	}
    }
    va_end(ap); 
    return Default;
}
///////////////////////////////////////////////////////////////////////////////
// (*) directly connected options
//.............................................................................
//
inline int
GetPot::direct_follow(int Default, const char* Option) 
{
    const char* FollowStr = match_starting_string(Option);
    if( FollowStr == 0 )                    return Default;
    if( ++cursor >= argv.size() ) cursor = argv.size();
    return __convert_to_type(FollowStr, Default);
}

inline double 
GetPot::direct_follow(const double& Default, const char* Option) 
{
    const char* FollowStr = match_starting_string(Option);
    if( FollowStr == 0 )                   return Default;
    if( ++cursor >= argv.size() ) cursor = argv.size();
    return __convert_to_type(FollowStr, Default);
}

inline const char* 
GetPot::direct_follow(const char* Default, const char* Option)
{
    if( search_failed_f ) return Default;
    const char* FollowStr = match_starting_string(Option);
    if( FollowStr == 0 )          return Default;
    if( ++cursor >= argv.size() ) cursor = argv.size();
    return FollowStr;
}

inline const char*
GetPot::match_starting_string(const char* StartString) 
    // pointer  to the place where the string after 
    //          the match inside the found argument starts.
    // 0        no argument matches the starting string.
{
    const unsigned N = strlen(StartString);
    unsigned OldCursor = cursor;
    if( OldCursor >= argv.size() ) OldCursor = argv.size() - 1;
    search_failed_f = true;

    // (*) first loop from cursor position until end
    unsigned c = cursor;
    for(; c < argv.size(); c++) {
	if( strncmp(StartString, argv[c].c_str(), N) == 0)
	{ cursor = c; search_failed_f = false; return &(argv[c].c_str()[N]); }
    }

    if( ! search_loop_f ) return false;

    // (*) second loop from 0 to old cursor position
    for(c = 1; c < OldCursor; c++) {
	if( strncmp(StartString, argv[c].c_str(), N) == 0)
	{ cursor = c; search_failed_f = false; return &(argv[c].c_str()[N]); }
    }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
// (*) search for flags
//.............................................................................
//
inline bool
GetPot::options_contain(const char* FlagList) const
{
    // go through all arguments that start with a '-' (but not '--')
    std::string str;
    for(StringVector::const_iterator it = argv.begin();
	it != argv.end();
	it++) {
	str = __get_remaining_string(*it, prefix);

	if( str.length() >= 2 && str[0] == '-' && str[1] != '-' )
	    if( check_flags(str, FlagList) ) return true;
    }
    return false;
}

inline bool
GetPot::argument_contains(unsigned Idx, const char* FlagList) const
{
    if( Idx >= argv.size() ) return false;

    if( prefix == "" )
        // search argument for any flag in flag list
	return check_flags(argv[Idx], FlagList);

    // if a prefix is set, then the argument index is the index
    //   inside the 'namespace'
    // => only check list of arguments that start with prefix
    unsigned no_matches = 0;
    for(unsigned i=0; i<argv.size(); i++) {
	const std::string Remain = __get_remaining_string(argv[i], prefix);
	if( Remain != "") {
	    no_matches += 1;
	    if( no_matches == Idx)
		return check_flags(Remain, FlagList);
	}
    }
    // no argument in this namespace
    return false;
}

inline bool
GetPot::check_flags(const std::string& Str, const char* FlagList) const
{
    for(const char* p=FlagList; *p != '\0' ; p++)
	if( Str.find(*p) != std::string::npos ) return true; // found something
    return false;
}
///////////////////////////////////////////////////////////////////////////////
// (*) nominus arguments
inline
void 
GetPot::reset_nominus_cursor()
{ nominus_cursor = -1; }

inline
StringVector  
GetPot::nominus_vector() const
    // return vector of nominus arguments
{
    StringVector nv;
    for(std::vector<unsigned>::const_iterator it = idx_nominus.begin();
	it != idx_nominus.end();
	it++) {
	nv.push_back(argv[*it]);
    }
    return nv;
}

inline
unsigned
GetPot::nominus_size() const
{ return idx_nominus.size(); }

inline
const char* 
GetPot::next_nominus()
{
  if( nominus_cursor < int(idx_nominus.size()) - 1 )
    return argv[idx_nominus[++nominus_cursor]].c_str();
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// (*) variables
//.............................................................................
//
inline int 
GetPot::operator()(const char* VarName, int Default) const
{
    const variable*  sv = find_variable(VarName);
    if( sv == 0 ) return Default;
    return __convert_to_type(sv->original, Default);
}

inline double
GetPot::operator()(const char* VarName, const double& Default) const
{
    const variable*  sv = find_variable(VarName);
    if( sv == 0 ) return Default;
    return __convert_to_type(sv->original, Default);
}

inline const char*
GetPot::operator()(const char* VarName, const char* Default) const
{
    const variable*  sv = find_variable(VarName);
    if( sv == 0 ) return Default;
    return sv->original.c_str();
}

inline int 
GetPot::operator()(const char* VarName, int Default, unsigned Idx) const
{
    const variable* sv = find_variable(VarName);
    if( sv == 0 ) return Default;
    const std::string*  element = sv->get_element(Idx);
    if( element == 0 ) return Default;
    return __convert_to_type(*element, Default);
}

inline double
GetPot::operator()(const char* VarName, const double& Default, unsigned Idx) const
{
    const variable* sv = find_variable(VarName);
    if( sv == 0 ) return Default;
    const std::string*  element = sv->get_element(Idx);
    if( element == 0 ) return Default;
    return __convert_to_type(*element, Default);
}

inline const char*
GetPot::operator()(const char* VarName, const char* Default, unsigned Idx) const
{
    const variable*  sv = find_variable(VarName);
    if( sv == 0 ) return Default;
    const std::string* element = sv->get_element(Idx);     
    if( element == 0 )  return Default;
    return element->c_str();
}

inline unsigned
GetPot::vector_variable_size(const char* VarName) const
{
    const variable*  sv = find_variable(VarName);
    if( sv == 0 ) return 0;
    return sv->value.size();
}

inline StringVector
GetPot::get_variable_names() const
{
    StringVector  result;
    for(std::vector<GetPot::variable>::const_iterator it = variables.begin();
	it != variables.end(); it++) {
	const std::string Tmp = __get_remaining_string((*it).name, prefix);
	if( Tmp != "" ) result.push_back(Tmp);
    }
    return result;
}

// (*) cursor oriented functions -------------------------------------------
void
GetPot::set_prefix(const char* Prefix)
{ prefix = std::string(Prefix); }

bool
GetPot::search_failed() const
{ return search_failed_f; }

//     -- enable/disable search for an option in loop
void
GetPot::disable_loop()
{ search_loop_f = false; }

void
GetPot::enable_loop()
{ search_loop_f = true; }

inline StringVector
GetPot::get_section_names() const
{ return section_list; }

inline const GetPot::variable*
GetPot::find_variable(const char* Name) const
{
    std::string name = prefix + Name;
    for(std::vector<variable>::const_iterator it = variables.begin();
	it != variables.end();
	it++)
	if( (*it).name == name ) return &(*it);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// (*) ouput (basically for debugging reasons
//.............................................................................
//
inline int 
GetPot::print() const
{
    std::cout << "argc = " << argv.size() << std::endl;
    for(StringVector::const_iterator it = argv.begin(); it != argv.end(); it++)
	std::cout << *it << std::endl;
    std::cout << std::endl;
    return 1;
}

  // (*) dollar bracket expressions (DBEs) ------------------------------------
  //
  //     1) Entry Function: DBE_expand_string()
  //        Takes a string such as
  //
  //          "${+ ${x} ${y}}   Subject-${& ${section} ${subsection}}:   ${title}"
  //
  //        calls DBE_expand() for each of the expressions
  //
  //           ${+ ${x} ${y}}
  //           ${& ${section} ${subsection}}
  //           ${Title}
  //
  //        and returns the string
  //
  //          "4711 Subject-1.01:   Mit den Clowns kamen die Schwaene"
  //
  //        assuming that
  //            x          = "4699"
  //            y          = "12"
  //            section    = "1."
  //            subsection = "01"
  //            title      = "Mit den Clowns kamen die Schwaene"
  //
  //      2) DBE_expand():
  //
  //           checks for the command, i.e. the 'sign' that follows '${'
  //           divides the argument list into sub-expressions using
  //           DBE_get_expr_list()
  //
  //           ${+ ${x} ${y}}                 -> "${x}"  "${y}"
  //           ${& ${section} ${subsection}}  -> "${section}" "${subsection}"
  //           ${Title}                       -> Nothing, variable expansion
  //
  //      3) DBE_expression_list():
  //
  //           builds a vector of unbracketed whitespace separated strings, i.e.
  //
  //           "  ${Number}.a ${: Das Marmorbild} AB-${& Author= ${Eichendorf}-1870}"
  //
  //           is split into a vector
  //
  //              [0] ${Number}.a
  //              [1] ${: Das Marmorbild}
  //              [2] AB-${& Author= ${Eichendorf}}-1870
  //
  //           Each sub-expression is expanded using expand(). 
  //---------------------------------------------------------------------------    
std::string 
GetPot::DBE_expand_string(const std::string str) 
{
    // Parses for closing operators '${ }' and expands them letting
    // white spaces and other letters as they are.
    std::string   new_string = "";
    unsigned open_brackets = 0;
    unsigned first = 0;
    for(unsigned i = 0;  i<str.size(); i++) {
	if( i < str.size() - 2 && str.substr(i, 2) == "${" ) {
	    if( open_brackets == 0 ) first = i+2;
	    open_brackets++;
	}
	else if( str[i] == '}' && open_brackets > 0) {
	    open_brackets -= 1;
	    if( open_brackets == 0 ) {
		const std::string Replacement = DBE_expand(str.substr(first, i - first));
		new_string += Replacement;
	    }
	}
	else if( open_brackets == 0 )
	    new_string += str[i];
    }
    return new_string;
}

StringVector
GetPot::DBE_get_expr_list(const std::string str_, const unsigned ExpectedNumber)
    // ensures that the resulting vector has the expected number
    // of arguments, but they may contain an error message
{
    std::string str = str_;
    // Separates expressions by non-bracketed whitespaces, expands them
    // and puts them into a list.

    // (1) eat initial whitespaces
    std::string::size_type i = str.find_first_not_of(" \f\n\r\t\v");
//     unsigned i=0;
//     for(; i < str.size(); ++i)
// 	if( ! isspace(str[i]) ) break;

    StringVector   expr_list;
    unsigned         open_brackets = 0;
    std::vector<unsigned> start_idx;
    unsigned         start_new_string = i;
    unsigned         l = str.size();

    // (2) search for ${ } expressions ...
    while( i < l ) {
 	const char letter = str[i];
	// whitespace -> end of expression
	if( isspace(letter)  &&  open_brackets == 0) {
	    expr_list.push_back(str.substr(start_new_string, i - start_new_string));
	    bool no_breakout_f = true;
	    for(i++; i < l ; i++) {
		if( ! isspace(str[i]) ) 
		{ no_breakout_f = false; start_new_string = i; break; }
	    }
	    if( no_breakout_f ) {
		// end of expression list
		if( expr_list.size() < ExpectedNumber ) {
		    const std::string   pre_tmp("<< ${ }: missing arguments>>");
		    StringVector tmp(ExpectedNumber - expr_list.size(), pre_tmp); 
		    expr_list.insert(expr_list.end(), tmp.begin(), tmp.end());
		}
		return expr_list;
	    }
	}
      
	// dollar-bracket expression
	if( str.length() >= i+2 && str.substr(i, 2) == "${" ) {
	    open_brackets++;
	    start_idx.push_back(i+2);
	}
	else if( letter == '}' && open_brackets > 0) {
	    unsigned start = start_idx[start_idx.size()-1];
	    start_idx.pop_back();
	    const std::string Replacement = DBE_expand(str.substr(start, i-start));
	    if( start < 3 )
		str = Replacement + str.substr(i+1);
	    else
		str = str.substr(0, start-2) + Replacement + str.substr(i+1);
	    l = str.size();              
	    i = start + Replacement.size() - 3;
	    open_brackets--;
	}
	i++;
    }

    // end of expression list
    expr_list.push_back(str.substr(start_new_string, i-start_new_string));

    if( expr_list.size() < ExpectedNumber ) {
	const std::string   pre_tmp("<< ${ }: missing arguments>>");
	StringVector tmp(ExpectedNumber - expr_list.size(), pre_tmp); 
	expr_list.insert(expr_list.end(), tmp.begin(), tmp.end());
    }

    return expr_list;
}

const GetPot::variable*
GetPot::DBE_get_variable(std::string VarName)
{
    static GetPot::variable ev;
    std::string secure_Prefix = prefix;
    
    prefix = section;
    // (1) first search in currently active section
    const GetPot::variable* var = find_variable(VarName.c_str());
    if( var != 0 ) { prefix = secure_Prefix; return var; }
    
    // (2) search in root name space
    prefix = "";
    var = find_variable(VarName.c_str());
    if( var != 0 ) { prefix = secure_Prefix; return var; }

    prefix = secure_Prefix;
    
    // error occured => variable name == ""
    char* tmp = new char[VarName.length() + 25];
    sprintf(tmp, "<<${ } variable '%s' undefined>>", VarName.c_str());
    ev.name = "";
    ev.original = std::string(tmp);
    delete [] tmp;
    return &ev;
}
  
std::string 
GetPot::DBE_expand(const std::string expr)
{
    // ${: } pure text
    if( expr[0] == ':' )
	return expr.substr(1);
      
    // ${& expr expr ... } text concatination
    else if( expr[0] == '&' ) {
	const StringVector A = DBE_get_expr_list(expr.substr(1), 1);

	StringVector::const_iterator it = A.begin();
	std::string result = *it++;
	for(; it != A.end(); it++) result += *it;
	
	return result;
    }
      
    // ${<-> expr expr expr} text replacement
    else if( expr.length() >= 3 && expr.substr(0, 3) == "<->" ) {
	StringVector A = DBE_get_expr_list(expr.substr(3), 3);
	unsigned tmp = 0;
	const unsigned L = A[1].length();
	while( (tmp = A[0].find(A[1])) != std::string::npos ) {
	    A[0].replace(tmp, L, A[2]);
	}
	return A[0];
    }
    // ${+ ...}, ${- ...}, ${* ...}, ${/ ...} expressions
    else if( expr[0] == '+' ) {
	StringVector A = DBE_get_expr_list(expr.substr(1), 2);
	StringVector::const_iterator it = A.begin();
	double result = __convert_to_type(*it++, 0.0);
	for(; it != A.end(); it++)
	    result += __convert_to_type(*it, 0.0);

	return double2string(result);
    }
    else if( expr[0] == '-' ) {
	StringVector A = DBE_get_expr_list(expr.substr(1), 2);
	StringVector::const_iterator it = A.begin();
	double result = __convert_to_type(*it++, 0.0);
	for(; it != A.end(); it++)
	    result -= __convert_to_type(*it, 0.0);

	return double2string(result);
    }
    else if( expr[0] == '*' ) {
	StringVector A = DBE_get_expr_list(expr.substr(1), 2);
	StringVector::const_iterator it = A.begin();
	double result = __convert_to_type(*it++, 0.0);
	for(; it != A.end(); it++)
	    result *= __convert_to_type(*it, 0.0);

	return double2string(result);
    }
    else if( expr[0] == '/' ) {
	StringVector A = DBE_get_expr_list(expr.substr(1), 2);
	StringVector::const_iterator it = A.begin();
	double result = __convert_to_type(*it++, 0.0);
	if( result == 0 ) return "0.0";
	for(it++; it != A.end(); it++) {
	    const double Q = __convert_to_type(*it, 0.0);
	    if( Q == 0.0 ) return "0.0";	    
	    result /= Q;
	}
	return double2string(result);
    }

    // ${^ ... } power expressions
    else if( expr[0] == '^' ) {
	StringVector A = DBE_get_expr_list(expr.substr(1), 2);
	StringVector::const_iterator it = A.begin();
	double result = __convert_to_type(*it++, 0.0);
	for(; it != A.end(); it++)
	    result = pow(result, __convert_to_type(*it, 0.0));
	return double2string(result);
    }
    
    // ${==  } ${<=  } ${>= } comparisons (return the number of the first 'match'
    else if( expr.length() >= 2 && 
	     ( expr.substr(0,2) == "==" || expr.substr(0,2) == ">=" || 
	       expr.substr(0,2) == "<=" || expr[0] == '>'           || expr[0] == '<')) {
	// differentiate between two and one sign operators
	unsigned op = 0;
	enum { EQ, GEQ, LEQ, GT, LT };
	if      ( expr.substr(0, 2) == "==" ) op = EQ;
	else if ( expr.substr(0, 2) == ">=" ) op = GEQ;
	else if ( expr.substr(0, 2) == "<=" ) op = LEQ;
	else if ( expr[0] == '>' )            op = GT;
	else    /*                     "<" */ op = LT;
	
	StringVector a;
	if ( op == GT || op == LT ) a = DBE_get_expr_list(expr.substr(1), 2);
	else                        a = DBE_get_expr_list(expr.substr(2), 2);

	std::string   x_orig = a[0];          
	double   x = __convert_to_type(x_orig, 1e37);
	unsigned i = 1;

	StringVector::const_iterator y_orig = a.begin();
	for(y_orig++; y_orig != a.end(); y_orig++) { 
	    double y = __convert_to_type(*y_orig, 1e37);

	    // set the strings as reference if one wasn't a number
	    if ( x == 1e37 || y == 1e37 ) {
		// it's a string comparison
		if( (op == EQ  && x_orig == *y_orig) || (op == GEQ && x_orig >= *y_orig) ||
		    (op == LEQ && x_orig <= *y_orig) || (op == GT  && x_orig >  *y_orig) ||
		    (op == LT  && x_orig <  *y_orig) )
		    return int2string(i);
	    }
	    else {
		// it's a number comparison
		if( (op == EQ  && x == y) || (op == GEQ && x >= y) ||
		    (op == LEQ && x <= y) || (op == GT  && x >  y) || 
		    (op == LT  && x <  y) )
		    return int2string(i);
	    }
	    i++;
	}
	
	// nothing fulfills the condition => return 0
	return "0";
    }
    // ${?? expr expr} select 
    else if( expr.length() >= 2 && expr.substr(0, 2) == "??" ) {
      StringVector a = DBE_get_expr_list(expr.substr(2), 2);
      double x = __convert_to_type(a[0], 1e37);
      // last element is always the default argument
      if( x == 1e37 || x < 0 || x >= a.size() - 1 ) return a[a.size()-1];

      // round x to closest integer
      return a[int(x+0.5)];
    }
    // ${? expr expr expr} if then else conditions
    else if( expr[0] == '?' ) {
	StringVector a = DBE_get_expr_list(expr.substr(1), 2);
	if( __convert_to_type(a[0], 0.0) == 1.0 ) return a[1];
	else if( a.size() > 2 )	return a[2];
    }
    // ${! expr} maxro expansion 
    else if( expr[0] == '!' ) {
	const GetPot::variable* Var = DBE_get_variable(expr.substr(1));
	// error
	if( Var->name == "" ) return std::string(Var->original);

	const StringVector A = DBE_get_expr_list(Var->original, 2);
	return A[0];
    }
    // ${@: } - string subscription
    else if( expr.length() >= 2 && expr.substr(0,2) == "@:" ) {
	const StringVector A = DBE_get_expr_list(expr.substr(2), 2);
	double x = __convert_to_type(A[1], 1e37);
	
	// last element is always the default argument
	if( x == 1e37 || x < 0 || x >= A[0].size() - 1)
	    return "<<1st index out of range>>";
	
	if( A.size() > 2 ) {
	    double y = __convert_to_type(A[2], 1e37);
	    if ( y != 1e37 && y > 0 && y <= A[0].size() - 1 && y > x )
		return A[0].substr(int(x+0.5), int(y+1.5) - int(x+0.5));
	    else if( y == -1 )
		return A[0].substr(int(x+0.5));
	    return "<<2nd index out of range>>";
	}
	else {
	    char* tmp = new char[2];
	    tmp[0] = A[0][int(x+0.5)]; tmp[1] = '\0';
	    std::string result(tmp);
	    delete [] tmp;
	    return result;
	}
    }
    // ${@ } - vector subscription
    else if( expr[0] == '@' ) {
	StringVector          A   = DBE_get_expr_list(expr.substr(1), 2);
	const GetPot::variable* Var = DBE_get_variable(A[0]);
	// error
	if( Var->name == "" ) {    
	    // make a copy of the string if an error occured
	    // (since the error variable is a static variable inside get_variable())
	    return std::string(Var->original);
	}

	double x = __convert_to_type(A[1], 1e37);
      
	// last element is always the default argument
	if (x == 1e37 || x < 0 || x >= Var->value.size() )
	    return "<<1st index out of range>>";

	if ( A.size() > 2) {
	    double y = __convert_to_type(A[2], 1e37);
	    int    begin = int(x+0.5);
	    int    end = 0;
	    if ( y != 1e37 && y > 0 && y <= Var->value.size() && y > x)
		end = int(y+1.5);
	    else if( y == -1 )
		end = Var->value.size();
	    else
		return "<<2nd index out of range>>";	    	    

	    std::string result = *(Var->get_element(begin));
	    for(int i = begin+1; i < end; i++) 
		result += std::string(" ") + *(Var->get_element(i));
	    return result;
	}
	else
	    return *(Var->get_element(int(x+0.5)));
    }
    
    const StringVector    A = DBE_get_expr_list(expr, 1);
    const GetPot::variable* B = DBE_get_variable(A[0]);
    
    // make a copy of the string if an error occured
    // (since the error variable is a static variable inside get_variable())
    if( B->name == "" ) return std::string(B->original);
    else                return B->original;
}


///////////////////////////////////////////////////////////////////////////////
// (*) unidentified flying objects
//.............................................................................
//
inline bool
GetPot::__search_string_vector(const StringVector& VecStr, const std::string& Str) const
{
    victorate(std::string, VecStr, itk) {
	if( *itk == Str ) return true;
    }
    return false;
}

inline StringVector 
GetPot::unidentified_arguments(unsigned Number, 
			       const char* KnownArgument1, ...) const
{
    StringVector known_arguments;

    // (1) create a vector of known arguments
    if( Number == 0 ) return StringVector();

    va_list ap;
    va_start(ap, KnownArgument1);
    known_arguments.push_back(std::string(KnownArgument1));
    for(unsigned i=1; i<Number; i++)
	known_arguments.push_back(std::string(va_arg(ap, char *)));
    va_end(ap);
    
    return unidentified_arguments(known_arguments);
}

inline StringVector 
GetPot::unidentified_arguments(const StringVector& Knowns) const
{
    StringVector ufos;
    StringVector::const_iterator it = argv.begin();
    it++; // forget about argv[0] (application or filename)
    for(; it != argv.end(); it++) {
	// -- argument belongs to prefixed section ?
	const std::string arg = __get_remaining_string(*it, prefix);
	if( arg == "" ) continue;

	// -- check if in list
	if( __search_string_vector(Knowns, arg) == false)
	    ufos.push_back(*it);
    }
    return ufos;
}


inline StringVector 
GetPot::unidentified_options(unsigned Number, 
			     const char* KnownOption1, ...) const
{
    StringVector known_options;

    // (1) create a vector of known arguments
    if( Number == 0 ) return StringVector();

    va_list ap;
    va_start(ap, KnownOption1);
    known_options.push_back(std::string(KnownOption1));
    for(unsigned i=1; i<Number; i++)
	known_options.push_back(std::string(va_arg(ap, char *)));
    va_end(ap);
    
    return unidentified_options(known_options);
}

inline StringVector 
GetPot::unidentified_options(const StringVector& Knowns) const
{
    StringVector ufos;
    StringVector::const_iterator it = argv.begin();
    it++; // forget about argv[0] (application or filename)
    for(; it != argv.end(); it++) {
	// -- argument belongs to prefixed section ?
	const std::string arg = __get_remaining_string(*it, prefix);
	if( arg == "" ) continue;

	// is argument really an option (starting with '-') ?
	if( arg.length() < 1 || arg[0] != '-' ) continue;

	if( __search_string_vector(Knowns, arg) == false)
	    ufos.push_back(*it);
    }
    
    return ufos;
}

inline std::string
GetPot::unidentified_flags(const char* KnownFlagList,
			   int ArgumentNumber) const
{
    std::string         ufos;
    StringVector known_arguments;
    std::string         KFL(KnownFlagList);

    // (2) iteration over '-' arguments (options)
    if( ArgumentNumber == -1 ) {
	StringVector::const_iterator it = argv.begin();
	it++; // forget about argv[0] (application or filename)
	for(; it != argv.end(); it++) {
	    // -- argument belongs to prefixed section ?
	    const std::string arg = __get_remaining_string(*it, prefix);
	    if( arg == "" ) continue;
	    
	    // -- does arguments start with '-' (but not '--')
	    if     ( arg.length() < 2 ) continue;
	    else if( arg[0] != '-' )    continue;
	    else if( arg[1] == '-' )    continue;
	    
	    // -- check out if flags inside option are contained in KnownFlagList
	    const char* p=arg.c_str();
	    p++; // skip starting minus
	    for(; *p != '\0' ; p++)
		if( KFL.find(*p) == std::string::npos ) ufos += *p;
	}
    }
    // (1) check specific argument
    else {
	// -- only check arguments that start with prefix
	int no_matches = 0;
	for(unsigned i=1; i<argv.size(); i++) {
	    const std::string Remain = __get_remaining_string(argv[i], prefix);
	    if( Remain != "") {
		no_matches++;
		if( no_matches == ArgumentNumber) {
		    // -- the right argument number inside the section is found
		    // => check it for flags
		    const char* p = Remain.c_str();
		    p++; // skip starting minus
		    for(; *p != '\0' ; p++)
			if( KFL.find(*p) == std::string::npos ) ufos += *p;
		    return ufos;
		}
	    }
	}
    }
    return ufos;
}

inline StringVector 
GetPot::unidentified_variables(unsigned Number, 
			       const char* KnownVariable1, ...) const
{
    StringVector known_variables;

    // create vector of known arguments
    if( Number == 0 ) return StringVector();

    va_list ap;
    va_start(ap, KnownVariable1);
    known_variables.push_back(std::string(KnownVariable1));
    for(unsigned i=1; i<Number; i++)
	known_variables.push_back(std::string(va_arg(ap, char *)));
    va_end(ap);
    
    return unidentified_variables(known_variables);
}

inline StringVector 
GetPot::unidentified_variables(const StringVector& Knowns) const
{
    StringVector ufos;

    victorate(GetPot::variable, variables, it) {
	// -- check if variable has specific prefix
	const std::string var_name = __get_remaining_string((*it).name, prefix);
	if( var_name == "" ) continue;

	// -- check if variable is known
	if( __search_string_vector(Knowns, var_name) == false)
	    ufos.push_back((*it).name);
    }
    return ufos;    
}

inline StringVector 
GetPot::unidentified_sections(unsigned Number, 
			      const char* KnownSection1, ...) const
{
    StringVector known_sections;

    // (1) create a vector of known arguments
    if( Number == 0 ) return StringVector();

    va_list ap;
    va_start(ap, KnownSection1);
    known_sections.push_back(std::string(KnownSection1));
    for(unsigned i=1; i<Number; i++) {
	std::string tmp = std::string(va_arg(ap, char *));
	if( tmp.length() == 0 ) continue;
	if( tmp[tmp.length()-1] != '/' ) tmp += '/';
	known_sections.push_back(tmp);
    }
    va_end(ap);

    return unidentified_sections(known_sections);
}

inline StringVector 
GetPot::unidentified_sections(const StringVector& Knowns) const
{
    StringVector ufos;

    victorate(std::string, section_list, it) {
	// -- check if section conform to prefix
	const std::string sec_name = __get_remaining_string(*it, prefix);
	if( sec_name == "" ) continue;

	// -- check if section is known
	if( __search_string_vector(Knowns, sec_name) == false )
	    ufos.push_back(*it);
    }
    
    return ufos;    
}

inline StringVector 
GetPot::unidentified_nominuses(unsigned Number, const char* Known, ...) const
{
    StringVector known_nominuses;

    // create vector of known arguments
    if( Number == 0 ) return StringVector();

    va_list ap;
    va_start(ap, Known);
    known_nominuses.push_back(std::string(Known));
    for(unsigned i=1; i<Number; i++) {
	std::string tmp = std::string(va_arg(ap, char *));
	if( tmp.length() == 0 ) continue;
	known_nominuses.push_back(tmp);
    }
    va_end(ap);
    
    return unidentified_nominuses(known_nominuses);
}

inline StringVector 
GetPot::unidentified_nominuses(const StringVector& Knowns) const
{
    StringVector ufos;

    // (2) iterate over all arguments 
    StringVector::const_iterator it = argv.begin();
    it++; // forget about argv[0] (application or filename)
    for(; it != argv.end(); it++) {
	// -- check if nominus part of prefix
	const std::string arg = __get_remaining_string(*it, prefix);
	if( arg == "" ) continue;

	if( arg.length() < 1 )                                  continue;
	// option ? --> not a nomius 
	if( arg[0] == '-' )                                     continue;
	// section ? --> not a real nominus
	if( arg[0] == '[' && arg[arg.length()-1] == ']' ) continue;
	// variable definition ? --> not a real nominus
	bool continue_f = false;
	for(unsigned i=0; i<arg.length() ; i++)
	    if( arg[i] == '=' ) { continue_f = true; break; }
	if( continue_f )                                           continue;

	// real nominuses are compared with the given list
	if( __search_string_vector(Knowns, arg) == false )
	    ufos.push_back(*it);
    }    
    return ufos;
}


///////////////////////////////////////////////////////////////////////////////
// (*) variable class
//.............................................................................
//
inline
GetPot::variable::~variable()
{}

inline
GetPot::variable::variable()
{}

inline
GetPot::variable::variable(const variable& Other)
  : name(Other.name), value(Other.value), original(Other.original)
{ 
}

inline
GetPot::variable::variable(const char* Name, const char* Value)
    : name(Name) 
{ 
    // make a copy of the 'Value'
    take(Value);
}

inline const std::string*
GetPot::variable::get_element(unsigned Idx) const
{ if( Idx >= value.size() ) return 0; else return &(value[Idx]); }

inline void
GetPot::variable::Swap(variable& Other)
{
  name.swap(Other.name);
  value.swap(Other.value);
  original.swap(Other.original);
}

inline void
GetPot::variable::take(const char* Value)
{ 
  // Orignial used strtok_r and was not exception safe due
  // to use of naked pointer.  Man page of strtok says:
  // "Never use these functions...", so I changed the
  // implementation. (2002-09-20, hedman@pdc.kth.se)
  //
  original = std::string(Value);
  value.clear();

  std::istringstream is(original);
  std::string slask;
  while ( is >> slask )
    value.push_back(slask);
}

// Simplify conversion.  Could be template instead.
std::string
GetPot::double2string(const double& Value) const { 
  std::ostringstream result;
  result << Value;
  return result.str();
}

// Simplify conversion.  Could be template instead.
std::string
GetPot::int2string(const int& Value) const { 
  std::ostringstream result;
  result << Value;
  return result.str();
}

#undef victorate
#undef reduce
#undef inline
