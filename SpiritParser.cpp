/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
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
   \date 2009-04-07
 */

#include <lifemc/lifecore/SpiritParser.hpp>





// ===================================================
//! Constructor & Destructor
// ===================================================
SpiritParser::SpiritParser( const bool& applyRules ) :
	M_strings 					( ),
	M_result					( 0. ),
	M_calculator				( M_variables, M_result ),
	M_applyRules				( applyRules )
{
	// Set default variables
	setDefaultVariables();
}



SpiritParser::SpiritParser( const std::string& string, const bool& applyRules ) :
	M_strings 					( ),
	M_result					( 0. ),
	M_calculator				( M_variables, M_result ),
	M_applyRules				( applyRules )
{
	// Set default variables
	setDefaultVariables();

	// Set the string
	setString( string );
}



SpiritParser::SpiritParser( const SpiritParser& parser ) :
	M_strings 					( parser.M_strings ),
	M_variables					( parser.M_variables ),
	M_result					( parser.M_result ),
	M_calculator				( parser.M_calculator ),
	M_applyRules				( parser.M_applyRules )
{
}



SpiritParser&
SpiritParser::operator=( const SpiritParser& parser )
{
    if ( this != &parser )
    {
    	M_strings 		= parser.M_strings;
    	M_variables		= parser.M_variables;
    	M_result 		= parser.M_result;
    	M_calculator 	= parser.M_calculator;
    	M_applyRules 	= parser.M_applyRules;
    }

	return *this;
}





// ===================================================
//! Methods
// ===================================================
void
SpiritParser::setString( const std::string& string, const std::string& stringSeparator )
{
    Debug( 5030 ) << "SpiritParser::setString:        M_strings: " << string 		<< "\n";
    Debug( 5030 ) << "                             M_applyRules: " << M_applyRules 	<< "\n";

    boost::split( M_strings, string, boost::is_any_of(stringSeparator) );

	if ( M_applyRules )
		for (UInt i = 0; i < M_strings.size(); ++i)
			ruleTheString( M_strings[i] );
}



void
SpiritParser::setVariable( const std::string& name, const Real& value )
{
	M_variables[name] = value;

    Debug( 5030 ) << "SpiritParser::setVariable: M_variables[" << name << "]: " << value << "\n";
}



Real&
SpiritParser::evaluate( void )
{
	for (UInt i = 0; i < M_strings.size(); ++i)
		boost::spirit::parse(M_strings[i].begin(), M_strings[i].end(), M_calculator, boost::spirit::space_p);

    Debug( 5030 ) << "SpiritParser::evaluate:          M_result: " << M_result << "\n";

	return M_result;
}





// ===================================================
//! Private functions
// ===================================================
inline void
SpiritParser::setDefaultVariables( void )
{
	// Set default variables
	M_variables["pi"] = 3.141592653589792;
	M_variables["e"]  = 2.718281828459046;
}



inline void
SpiritParser::ruleTheString( std::string& string )
{
	Debug( 5030 ) << "SpiritParser::ruleTheString: " << "\n";
	Debug( 5030 ) << "                     (before) - M_strings: " << string 	<< "\n";

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

	Debug( 5030 ) << "                       (after) - M_string: " << string 	<< "\n";
}
