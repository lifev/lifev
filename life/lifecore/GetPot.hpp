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
//! \file util_string.h
//  -*- c++ -*-
//  GetPot Version 1.0-LifeV             Sept/13/2002
//
//
//
/*!
  WEBSITE: http://getpot.sourceforge.net
 
  This library is  free software; you can redistribute  it and/or modify
  it  under  the terms  of  the GNU  Lesser  General  Public License  as
  published by the  Free Software Foundation; either version  2.1 of the
  License, or (at your option) any later version.
  
  This library  is distributed in the  hope that it will  be useful, but
  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
  Lesser General Public License for more details.
  
  You  should have  received a  copy of  the GNU  Lesser  General Public
  License along  with this library; if  not, write to  the Free Software
  Foundation, Inc.,  59 Temple Place,  Suite 330, Boston,  MA 02111-1307
  USA
  
  (C) 2001-2002 Frank R. Schaefer  
*/ 
//==========================================================================
#ifndef __GETPOT_H__
#define __GETPOT_H__

#include <string>
#include <vector>

namespace getpot
{
typedef std::vector<std::string> StringVector;
}

//! Handles command line options and files
/*
\version 1.0-LifeV
This class handles the parsing of command line and options file according to the
rules set up in the getpot manual.
 
This is a version modified for the LifeV library, based on the official release 1.0
 
<ul>
<li>It has been split into two files: one for the declarations (GetPot.h) and
one for the definitions (GetPot.cc) in order to comply with the LifeV
general organization.</li>
<li>A hacking has been made in order to accept as input a list of filenames:
the data base will be constructed as the union of the data contained in all
the files in the list</li>
</ul>
 
\note: in order to upgrade from a newer official version the same
modifications have to be carried out.
 
*/
class GetPot
{
    //--------
public:
    // (*) constructors, destructor, assignment operator -----------------------
    GetPot();
    GetPot( const GetPot& );
    GetPot( int argc_, char* argv_[] );
    GetPot( const char* FileName );
    // begin luca 29/12/2002
    GetPot( const getpot::StringVector & FileNameList );
    // end luca 29/12/2002
    ~GetPot();
    GetPot& operator=( const GetPot& );

    // (*) direct access to command line arguments -----------------------------
    const char* operator[] ( unsigned Idx ) const;
    int get
        ( unsigned Idx, int Default ) const;
    double get
        ( unsigned Idx, const double& Default ) const;
    const char* get
        ( unsigned Idx, const char* Default ) const;
    unsigned size() const;

    // (*) flags ---------------------------------------------------------------
    bool options_contain( const char* FlagList ) const;
    bool argument_contains( unsigned Idx, const char* FlagList ) const;

    // (*) variables -----------------------------------------------------------
    //     -- scalar values
    int operator() ( const char* VarName, int Default ) const;
    double operator() ( const char* VarName, const double& Default ) const;
    const char* operator() ( const char* VarName, const char* Default ) const;
    //     -- vectors
    int operator() ( const char* VarName, int Default, unsigned Idx ) const;
    double operator() ( const char* VarName, const double& Default, unsigned Idx ) const;
    const char* operator() ( const char* VarName, const char* Default, unsigned Idx ) const;
    unsigned vector_variable_size( const char* VarName ) const;
    getpot::StringVector get_variable_names() const;
    getpot::StringVector get_section_names() const;


    // (*) cursor oriented functions -------------------------------------------
    void set_prefix( const char* Prefix );
    bool search_failed() const;

    //     -- enable/disable search for an option in loop
    void disable_loop();
    void enable_loop();

    //     -- reset cursor to position '1'
    void reset_cursor();
    void init_multiple_occurrence();

    //     -- search for a certain option and set cursor to position
    bool search( const char* option );
    bool search( unsigned No, const char* P, ... );
    //     -- get argument at cursor++
    int next( int Default );
    double next( const double& Default );
    const char* next( const char* Default );
    //     -- search for option and get argument at cursor++
    int follow( int Default, const char* Option );
    double follow( const double& Default, const char* Option );
    const char* follow( const char* Default, const char* Option );
    //     -- search for one of the given options and get argument that follows it
    int follow( int Default, unsigned No, const char* Option, ... );
    double follow( const double& Default, unsigned No, const char* Option, ... );
    const char* follow( const char* Default, unsigned No, const char* Option, ... );
    //     -- directly followed arguments
    int direct_follow( int Default, const char* Option );
    double direct_follow( const double& Default, const char* Option );
    const char* direct_follow( const char* Default, const char* Option );

    // (*) nominus arguments ---------------------------------------------------
    void reset_nominus_cursor();
    getpot::StringVector nominus_vector() const;
    unsigned nominus_size() const;
    const char* next_nominus();

    // (*) unidentified flying objects -----------------------------------------
    getpot::StringVector unidentified_arguments( unsigned Number, const char* Known, ... ) const;
    getpot::StringVector unidentified_arguments( const getpot::StringVector& Knowns ) const;

    getpot::StringVector unidentified_options( unsigned Number, const char* Known, ... ) const;
    getpot::StringVector unidentified_options( const getpot::StringVector& Knowns ) const;

    // Two modes:
    //  ArgumentNumber >= 0 check specific argument
    //  ArgumentNumber == -1 check all options starting with one '-' for flags
    std::string unidentified_flags( const char* Known,
                                    int ArgumentNumber = -1 ) const;

    getpot::StringVector unidentified_variables( unsigned Number, const char* Known, ... ) const;
    getpot::StringVector unidentified_variables( const getpot::StringVector& Knowns ) const;

    getpot::StringVector unidentified_sections( unsigned Number, const char* Known, ... ) const;
    getpot::StringVector unidentified_sections( const getpot::StringVector& Knowns ) const;

    getpot::StringVector unidentified_nominuses( unsigned Number, const char* Known, ... ) const;
    getpot::StringVector unidentified_nominuses( const getpot::StringVector& Knowns ) const;

    // (*) output --------------------------------------------------------------
    int print() const;

private:
    // (*) Type Declaration ----------------------------------------------------
    struct variable
    {
        //-----------
        // Variable to be specified on the command line or in input files.
        // (i.e. of the form var='12 312 341')

        // -- constructors, destructors, assignment operator
        ~variable();
        variable();
        variable( const variable& );
        variable( const char* Name, const char* Value );
        variable& operator=( const variable& Other )
        {
            variable temp( Other );
            Swap( temp );
            return *this;
        }

        void Swap( variable& Other );
        void take( const char* Value );

        // -- get a specific element in the string vector
        //    (return 0 if not present)
        const std::string* get_element( unsigned Idx ) const;


        // -- data memebers
        std::string name;      // identifier of variable
        getpot::StringVector value;     // value of variable stored in vector
        std::string original;  // value of variable as given on command line
    };

    // (*) variables -----------------------------------------------------------
    std::string prefix;          // prefix automatically added in queries
    std::string section;         // (for dollar bracket parsing)
    getpot::StringVector section_list;    // list of all parsed sections
    //     -- argurment vector
    getpot::StringVector argv;            // vector of command line arguments stored as strings
    unsigned cursor;          // cursor for argv
    bool search_loop_f;   // shall search start at beginning after
    //                               // reaching end of arg array ?
    bool search_failed_f; // flag indicating a failed search() operation
    //                               // (e.g. next() functions react with 'missed')

    //     --  nominus vector
    int nominus_cursor; // cursor for nominus_pointers
    std::vector<unsigned> idx_nominus;     // indecies of 'no minus' arguments

    //    -- intern variables
    //       (arguments of the form "variable=value")
    std::vector<variable> variables;

    // (*) helper functions ----------------------------------------------------
    //     -- produce three basic data vectors:
    //          - argument vector
    //          - nominus vector
    //          - variable dictionary
    void __parse_argument_vector( const getpot::StringVector& ARGV );

    //     -- helpers for argument list processing
    //        * search for a variable in 'variables' array
    const variable* find_variable( const char* ) const;
    //        * support finding directly followed arguments
    const char* match_starting_string( const char* StartString );
    //        * support search for flags in a specific argument
    bool check_flags( const std::string& Str, const char* FlagList ) const;
    //        * type conversion if possible
    int __convert_to_type( const std::string& String, int Default ) const;
    double __convert_to_type( const std::string& String, double Default ) const;
    //        * prefix extraction
    const std::string __get_remaining_string( const std::string& String, const std::string& Start ) const;
    //        * search for a specific string
    bool __search_string_vector( const getpot::StringVector& Vec,
                                 const std::string& Str ) const;

    //     -- helpers to parse input files
    //        create an argument vector based on data found in an input file, i.e.:
    //           1) delete '#'-comments
    //           2) contract assignment expressions, such as
    //                   my-variable   =    '007 J. B.'
    //             into
    //                   my-variable='007 J. B.'
    //           3) interprete sections like '[../my-section]' etc.
    void __skip_whitespace( std::istream& istr );
    const std::string __get_next_token( std::istream& istr );
    const std::string __get_string( std::istream& istr );
    const std::string __get_until_closing_bracket( std::istream& istr );

    getpot::StringVector read_in_stream( std::istream& istr );
    getpot::StringVector read_in_file( const char* FileName );
    std::string process_section_label( const std::string& Section,
                                       getpot::StringVector& section_stack );

    //      -- dollar bracket expressions
    std::string DBE_expand_string( const std::string str );
    std::string DBE_expand( const std::string str );
    const GetPot::variable* DBE_get_variable( const std::string str );
    getpot::StringVector DBE_get_expr_list( const std::string str, const unsigned ExpectedNumber );

    std::string double2string( const double& Value ) const;
    std::string int2string( const int& Value ) const;
};

#endif // __GETPOT_H__


