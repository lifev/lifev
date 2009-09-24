/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>,
			 Gilles Fourestey  <gilles.fourestey@epfl.ch>
       Date: 2009-04-07

  Copyright (C) 2009 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file SpiritParser.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \author Gilles Fourestey  <gilles.fourestey@epfl.ch>
   \date 2009-04-07
 */

#include <lifemc/lifecore/SpiritParser.hpp>

namespace LifeV {

// ===================================================
//! Constructor & Destructor
// ===================================================
SpiritParser::SpiritParser( const bool& applyRules ) :
	M_strings 					( ),
	M_variables					( ),
	M_results					( ),
	M_nResults					( 0 ),
	M_calculator				( M_variables, M_results, M_nResults ),
	M_applyRules				( applyRules )
{

#ifdef DEBUG
    Debug( 5030 ) << "SpiritParser::SpiritParser"<< "\n";
#endif

	// Set default variables
	setDefaultVariables();
}



SpiritParser::SpiritParser( const std::string& string, const bool& applyRules ) :
	M_strings 					( ),
	M_variables					( ),
	M_results					( ),
	M_nResults					( 0 ),
	M_calculator				( M_variables, M_results, M_nResults ),
	M_applyRules				( applyRules )
{

#ifdef DEBUG
    Debug( 5030 ) << "SpiritParser::SpiritParser"<< "\n";
#endif

	// Set default variables
	setDefaultVariables();

	// Set the string (it is not necessary to empty M_results)
	setString( string );
}

SpiritParser::SpiritParser( const SpiritParser& parser ) :
	M_strings 					( parser.M_strings ),
	M_variables					( parser.M_variables ),
	M_results					( parser.M_results ),
	M_nResults					( parser.M_nResults ),
	M_calculator				( parser.M_calculator ),
	M_applyRules				( parser.M_applyRules )
{
}



// ===================================================
//! Methods
// ===================================================
SpiritParser&
SpiritParser::operator=( const SpiritParser& parser )
{
    if ( this != &parser )
    {
    	M_strings 		= parser.M_strings;
    	M_variables		= parser.M_variables;
    	M_results 		= parser.M_results;
    	M_nResults 		= parser.M_nResults;
    	M_calculator 	= parser.M_calculator;
    	M_applyRules 	= parser.M_applyRules;
    }

	return *this;
}

void
SpiritParser::setString( const std::string& string, const std::string& stringSeparator )
{
#ifdef DEBUG
    Debug( 5030 ) << "SpiritParser::setString:          strings: " << string 		<< "\n";
#endif

    M_strings.clear();
    boost::split( M_strings, string, boost::is_any_of(stringSeparator) );
#ifdef DEBUG
    Debug( 5030 ) << "                               applyRules: " << M_applyRules 	<< "\n";
#endif

	if ( M_applyRules )
		for ( UInt i = 0; i < M_strings.size(); ++i )
			ruleTheString( M_strings[i] );

	setupResults();
}

void
SpiritParser::setVariable( const std::string& name, const Real& value )
{
#ifdef DEBUG
	Debug( 5030 ) << "SpiritParser::setVariable    variables[" << name << "]: " << value << "\n";
#endif

	M_variables[name] = value;
}

const Real&
SpiritParser::getVariable( const std::string& name )
{
#ifdef DEBUG
	Debug( 5030 ) << "SpiritParser::getVariable    variables[" << name << "]: " << M_variables[name] << "\n";
#endif

	return M_variables[name];
}

const Real&
SpiritParser::evaluate( const UInt& ID )
{
	for ( UInt i = 0; i < M_strings.size(); ++i )
		boost::spirit::parse(M_strings[i].begin(), M_strings[i].end(), M_calculator, boost::spirit::space_p);

#ifdef DEBUG
    Debug( 5030 ) << "SpiritParser::evaluate          results[ "<< (ID - 1) << "]: " << M_results[ID - 1] << "\n";
#endif

    M_nResults = 0; //Reset for next evaluation

	return M_results[ID - 1];
}

UInt
SpiritParser::countSubstring( const std::string& substring )
{
	UInt count( 0 );
	std::string::size_type position( 0 );

	for ( ; ; )
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
//! Private functions
// ===================================================
inline void
SpiritParser::setDefaultVariables( void )
{
	setVariable( "pi", Pi );
	//setVariable( "pi", 3.141592653589793 );
	//setVariable( "pi", 3.1415926535897932384626433832795 ); //Better only with long double!
	setVariable( "e" , 2.718281828459046 );
}



inline void
SpiritParser::setupResults( void )
{
	M_nResults = 0;

	//Reserve the space for results
	M_results.clear();
	M_results.reserve( countSubstring( "," ) + 1 );

#ifdef DEBUG
	Debug( 5030 ) << "SpiritParser::setupResults      dimension: " << countSubstring( "," ) + 1 << "\n";
#endif

	//std::vector<std::string> tempVectorString;
	//boost::split( tempVectorString, M_strings.back(), boost::is_any_of(stringSeparator) );
	//M_results.reserve( tempVectorString.size() );
}

inline void
SpiritParser::ruleTheString( std::string& string )
{
#ifdef DEBUG
	Debug( 5030 ) << "SpiritParser::ruleTheString" << "\n";
	Debug( 5030 ) << "                        (before) - string: " << string 	<< "\n";
#endif

	// Remove spaces from the string
	boost::replace_all( string, " ",  "" );

	// Solve the problem of signed expressions
	boost::replace_all( string, "(-",  "(-1*" );
	boost::replace_all( string, ",-",  ",-1*" );
	boost::replace_all( string, "=-",  "=-1*" );
	if ( string[0] == '-' )
		string.insert(1,"1*");

	// Convert the string to lower case
	boost::to_lower( string );

	// Apply the rules for the parser
	boost::replace_all( string, "sqrt",  "Q" );
	boost::replace_all( string, "exp",   "E" );
	boost::replace_all( string, "log",   "L" );
	boost::replace_all( string, "log10", "M" );
	boost::replace_all( string, "sin",   "S" );
	boost::replace_all( string, "cos",   "C" );
	boost::replace_all( string, "tan",   "T" );

	//Apply exceptions for the parser
	boost::replace_all( string, "visCity",   "viscosity" );

#ifdef DEBUG
	Debug( 5030 ) << "                         (after) - string: " << string 	<< "\n";
#endif
}

} // Namespace LifeV
