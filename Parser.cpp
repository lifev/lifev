//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 * @file
 * @brief Parser
 *
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 * @author Gilles Fourestey  <gilles.fourestey@epfl.ch>
 * @date 07-04-2009
 */

#include <lifemc/lifecore/Parser.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
Parser::Parser( const bool& applyRules ) :
    M_strings       (),
    M_variables     (),
    M_results       (),
    M_nResults      ( 0 ),
    M_calculator    ( M_variables, M_results, M_nResults ),
    M_applyRules    ( applyRules )
{

#ifdef DEBUG
    Debug( 5030 ) << "Parser::Parser"<< "\n";
#endif

    // Set default variables
    setDefaultVariables();
}

Parser::Parser( const std::string& string, const bool& applyRules ) :
    M_strings       (),
    M_variables     (),
    M_results       (),
    M_nResults      ( 0 ),
    M_calculator    ( M_variables, M_results,M_nResults ),
    M_applyRules    ( applyRules )
{

#ifdef DEBUG
    Debug( 5030 ) << "Parser::Parser"<< "\n";
#endif

    // Set default variables
    setDefaultVariables();

    // Set the string (it is not necessary to empty M_results)
    setString( string );
}

Parser::Parser( const Parser& parser ) :
    M_strings       ( parser.M_strings ),
    M_variables     ( parser.M_variables ),
    M_results       ( parser.M_results ),
    M_nResults      ( parser.M_nResults ),
    M_calculator    ( parser.M_calculator ),
    M_applyRules    ( parser.M_applyRules )
{
}

// ===================================================
// Operators
// ===================================================
Parser&
Parser::operator=( const Parser& parser )
{
    if ( this != &parser )
    {
        M_strings    = parser.M_strings;
        M_variables  = parser.M_variables;
        M_results    = parser.M_results;
        M_nResults   = parser.M_nResults;
        M_calculator = parser.M_calculator;
        M_applyRules = parser.M_applyRules;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
const Real&
Parser::evaluate( const UInt& ID )
{
    for ( UInt i = 0; i < M_strings.size(); ++i )
        spirit::parse( M_strings[i].begin(), M_strings[i].end(), M_calculator, spirit::space_p );

#ifdef DEBUG
    Debug( 5030 ) << "Parser::evaluate          results[ "<< (ID - 1) << "]: " << M_results[ID - 1] << "\n";
#endif

    M_nResults = 0; //Reset for next evaluation

    return M_results[ID - 1];
}

UInt
Parser::countSubstring( const std::string& substring )
{
    UInt count( 0 );
    std::string::size_type position( 0 );

    for ( ;; )
    {
        position = M_strings.back().find( substring, position );

        if ( position == std::string::npos )
            break;

        ++count;
        position += substring.length(); // start next search after this substring
    }

    return count;
}

// ===================================================
// Set Methods
// ===================================================
void
Parser::setString( const std::string& string, const std::string& stringSeparator )
{

#ifdef DEBUG
    Debug( 5030 ) << "Parser::setString:          strings: " << string << "\n";
#endif

    M_strings.clear();
    boost::split( M_strings, string, boost::is_any_of( stringSeparator ) );

#ifdef DEBUG
    Debug( 5030 ) << "                               applyRules: " << M_applyRules << "\n";
#endif

    if ( M_applyRules )
        for ( UInt i = 0; i < M_strings.size(); ++i )
            ruleTheString( M_strings[i] );

    setupResults();
}

void
Parser::setVariable( const std::string& name, const Real& value )
{

#ifdef DEBUG
    Debug( 5030 ) << "Parser::setVariable    variables[" << name << "]: " << value << "\n";
#endif

    M_variables[name] = value;
}

// ===================================================
// Get Methods
// ===================================================
const Real&
Parser::getVariable( const std::string& name )
{

#ifdef DEBUG
    Debug( 5030 ) << "Parser::getVariable    variables[" << name << "]: " << M_variables[name] << "\n";
#endif

    return M_variables[name];
}

// ===================================================
// Private functions
// ===================================================
inline void
Parser::setDefaultVariables()
{
    //setVariable( "pi", Pi ); //Using the tab of LifeV
    setVariable( "pi", 3.141592653589793 );
    //setVariable( "pi", 3.1415926535897932384626433832795 ); //Better only with long Real!
    setVariable( "e", 2.718281828459046 );
}

inline void
Parser::setupResults()
{
    M_nResults = 0;

    //Reserve the space for results
    M_results.clear();
    M_results.reserve( countSubstring( "," ) + 1 );

#ifdef DEBUG
    Debug( 5030 ) << "Parser::setupResults      dimension: " << countSubstring( "," ) + 1 << "\n";
#endif

    //std::vector<std::string> tempVectorString;
    //boost::split( tempVectorString, M_strings.back(), boost::is_any_of(stringSeparator) );
    //M_results.reserve( tempVectorString.size() );
}

inline void
Parser::ruleTheString( std::string& string )
{

#ifdef DEBUG
    Debug( 5030 ) << "Parser::ruleTheString" << "\n";
    Debug( 5030 ) << "                        (before) - string: " << string << "\n";
#endif

    // Remove spaces from the string
    boost::replace_all( string, " ", "" );

    // Solve the problem of signed expressions
    boost::replace_all( string, "(-", "(-1*" );
    boost::replace_all( string, ",-", ",-1*" );
    boost::replace_all( string, "=-", "=-1*" );
    if ( string[0] == '-' )
        string.insert( 1, "1*" );

    // Convert the string to lower case
    boost::to_lower( string );

    // Apply the rules for the parser
    boost::replace_all( string, "sqrt", "Q" );
    boost::replace_all( string, "exp", "E" );
    boost::replace_all( string, "log", "L" );
    boost::replace_all( string, "log10", "M" );
    boost::replace_all( string, "sin", "S" );
    boost::replace_all( string, "cos", "C" );
    boost::replace_all( string, "tan", "T" );

    //Apply exceptions for the parser
    boost::replace_all( string, "visCity", "viscosity" );

#ifdef DEBUG
    Debug( 5030 ) << "                         (after) - string: " << string << "\n";
#endif

}

} // Namespace LifeV
