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
	//M_variables					( new variables_type ),
	//M_results					( new results_type ),
	M_variables					( ),
	M_results					( ),
	M_nResults					( 0 ),
	//M_calculator				( *M_variables, *M_results, M_nResults ),
	M_calculator				( M_variables, M_results, M_nResults ),
	M_applyRules				( applyRules )
{
	// Set default variables
	setDefaultVariables();
}



SpiritParser::SpiritParser( const std::string& string, const bool& applyRules ) :
	M_strings 					( ),
	//M_variables					( new variables_type ),
	//M_results					( new results_type ),
	M_variables					( ),
	M_results					( ),
	M_nResults					( M_nResults ),
	//M_calculator				( *M_variables, *M_results, M_nResults ),
	M_calculator				( M_variables, M_results, M_nResults ),
	M_applyRules				( applyRules )
{
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





// ===================================================
//! Methods
// ===================================================
void
SpiritParser::setString( const std::string& string, const std::string& stringSeparator )
{
    Debug( 5030 ) << "SpiritParser::setString:          strings: " << string 		<< "\n";

    boost::split( M_strings, string, boost::is_any_of(stringSeparator) );

    Debug( 5030 ) << "                             M_applyRules: " << M_applyRules 	<< "\n";
	if ( M_applyRules )
		for ( UInt i = 0; i < M_strings.size(); ++i )
			ruleTheString( M_strings[i] );

	setupResults();
}



void
SpiritParser::setVariable( const std::string& name, const Real& value )
{
	Debug( 5030 ) << "SpiritParser::setVariable: M_variables[" << name << "]: " << value << "\n";

	//M_variables->operator[](name) = value;
	M_variables[name] = value;
}



Real&
SpiritParser::evaluate( const UInt& ID )
{
	for (UInt i = 0; i < M_strings.size(); ++i)
		boost::spirit::parse(M_strings[i].begin(), M_strings[i].end(), M_calculator, boost::spirit::space_p);

    //Debug( 5030 ) << "SpiritParser::evaluate:       M_results[ "<< (ID - 1) << "]: " << M_results->operator[](ID - 1) << "\n";
    Debug( 5030 ) << "SpiritParser::evaluate:       M_results[ "<< (ID - 1) << "]: " << M_results[ID - 1] << "\n";

    M_nResults = 0; //Reset for next evaluation
	//return M_results->operator[](ID - 1);
	return M_results[ID - 1];
}





// ===================================================
//! Private functions
// ===================================================
inline void
SpiritParser::setDefaultVariables( void )
{
	setVariable( "pi", 3.141592653589792 );
	setVariable( "e" , 2.718281828459046 );
}



inline void
SpiritParser::setupResults( const std::string& stringSeparator )
{
	M_nResults = 0;

	//Reserve the space for results
	std::vector<std::string> tempVectorString;
	boost::split( tempVectorString, M_strings.back(), boost::is_any_of(stringSeparator) );
	//(*M_results).reserve( tempVectorString.size() );
	M_results.reserve( tempVectorString.size() );
}



inline void
SpiritParser::ruleTheString( std::string& string )
{
	Debug( 5030 ) << "SpiritParser::ruleTheString: " << "\n";
	Debug( 5030 ) << "                      (before) - M_string: " << string 	<< "\n";

	// Remove spaces from the string
	boost::replace_all( string, " ",  "" );

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

} // Namespace LifeV
