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

#include <lifemc/lifefem/SpiritParser.hpp>





// ===================================================
//! Constructor & Destructor
// ===================================================
SpiritParser::SpiritParser( const bool& applyRules ) :
	M_string 					( "0" ),
	M_result					( 0. ),
	M_calculator				( M_variables, M_result ),
	M_applyRules				( applyRules )
{
	// Set default variables
	M_variables["pi"] = 3.141592653589792;
	M_variables["e"]  = 2.718281828459046;
}



SpiritParser::SpiritParser( const std::string& string, const bool& applyRules ) :
	M_result					( 0. ),
	M_calculator				( M_variables, M_result ),
	M_applyRules				( applyRules )
{
	// Set the string
	setString( string );

	// Set default variables
	M_variables["pi"] = 3.141592653589792;
	M_variables["e"]  = 2.718281828459046;
}





// ===================================================
//! Methods
// ===================================================
void
SpiritParser::setString( const std::string& string )
{
	M_string = string;
	if ( M_applyRules )
		ruleTheString( M_string );
}



void
SpiritParser::setVariable( const std::string& name, const Real& value )
{
	M_variables[name] = value;
}



Real&
SpiritParser::evaluate( void )
{
    boost::spirit::parse(M_string.begin(), M_string.end(), M_calculator, boost::spirit::space_p);

	return M_result;
}





// ===================================================
//! Private functions
// ===================================================
inline void
SpiritParser::ruleTheString( std::string& string )
{
	// Convert the string to lower case
	stringToLowerCase( string );

	// Apply the rules for the parser
	replaceContent( string, "sqrt",  "Q" );
	replaceContent( string, "exp",   "E" );
	replaceContent( string, "log",   "L" );
	replaceContent( string, "log10", "M" );
	replaceContent( string, "sin",   "S" );
	replaceContent( string, "cos",   "C" );
	replaceContent( string, "tan",   "T" );
}

inline void
SpiritParser::stringToLowerCase( std::string& string )
{
	for (std::string::iterator IT = string.begin() ; IT != string.end() ; ++IT)
		*IT = std::tolower((unsigned char)*IT);
}



inline void
SpiritParser::replaceContent( std::string& string, const std::string& content, const std::string& replace )
{
	std::size_t dimension = content.size();
	std::size_t position;

	for (UInt i(0) ; ; ++i)
	{
		position = string.find( content );
		if ( position == std::string::npos )
			break;
		else
			string.replace( position, dimension, replace );
	}
}
