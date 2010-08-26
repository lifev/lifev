//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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

#include "Parser.hpp"

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
Parser::Parser() :
    M_strings       (),
    M_results       (),
    M_calculator    (),
    M_evaluate      ( true )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5030 ) << "Parser::Parser"<< "\n";
#endif

    M_calculator.SetDefaultVariables();
}

Parser::Parser( const std::string& String ) :
    M_strings       (),
    M_results       (),
    M_calculator    (),
    M_evaluate      ( true )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5030 ) << "Parser::Parser( string, applyRules )"<< "\n";
#endif

    M_calculator.SetDefaultVariables();
    SetString( String );
}

Parser::Parser( const Parser& parser ) :
    M_strings       ( parser.M_strings ),
    M_results       ( parser.M_results ),
    M_calculator    ( parser.M_calculator ),
    M_evaluate      ( parser.M_evaluate )
{
}

// ===================================================
// Operators
// ===================================================
Parser&
Parser::operator=( const Parser& Parser )
{
    if ( this != &Parser )
    {
        M_strings    = Parser.M_strings;
        M_results    = Parser.M_results;
        //M_calculator = Parser.M_calculator; //NOT WORKING!!!
        M_evaluate   = Parser.M_evaluate;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
const Real&
Parser::Evaluate( const UInt& ID )
{
    if ( M_evaluate )
    {
        M_results.clear();
        String_Iterator start, end;

        for ( UInt i = 0; i < M_strings.size(); ++i )
        {
            start = M_strings[i].begin();
            end   = M_strings[i].end();
#ifdef HAVE_BOOST_SPIRIT_QI
            qi::phrase_parse( start, end, M_calculator, ascii::space, M_results);
#else
            std::cout << "!!! ERROR: Boost version < 1.41 !!!" << std::endl;
            // This generate an error ---------
            int *a = new int(0); int *b; b = a;
            delete a; delete b;
            // --------------------------------
#endif
        }

        M_evaluate = false;
    }

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5030 ) << "Parser::evaluate          results[ "<< (ID - 1) << "]: " << M_results[ID - 1] << "\n";
#endif

    return M_results[ID - 1];
}

UInt
Parser::CountSubstring( const std::string& Substring )
{
    UInt count( 0 );
    std::string::size_type position( 0 );

    for ( ;; )
    {
        position = M_strings.back().find( Substring, position );

        if ( position == std::string::npos )
            break;

        ++count;
        position += Substring.length(); // start next search after this substring
    }

    return count;
}

void
Parser::ClearVariables()
{
    M_calculator.ClearVariables();
    M_evaluate = true;
}

// ===================================================
// Set Methods
// ===================================================
void
Parser::SetString( const std::string& String, const std::string& StringSeparator )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5030 ) << "Parser::setString:          strings: " << String << "\n";
#endif

    M_strings.clear();
    boost::split( M_strings, String, boost::is_any_of( StringSeparator ) );

    //Remove white space to speed up the parser
    for ( UInt i = 0; i < M_strings.size(); ++i )
        boost::replace_all( M_strings[i], " ", "" );

    //Reserve the space for results
    M_results.clear();
    M_results.reserve( CountSubstring( "," ) + 1 );

    M_evaluate = true;
}

void
Parser::SetVariable( const std::string& Name, const Real& Value )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5030 ) << "Parser_Utility::SetVariable    variables[" << Name << "]: " << Value << "\n";
#endif

    M_calculator.SetVariable( Name, Value);

    M_evaluate = true;
}

// ===================================================
// Get Methods
// ===================================================
const Real&
Parser::GetVariable( const std::string& Name )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5030 ) << "Parser_Utility::GetVariable    variables[" << Name << "]: " << M_calculator.GetVariable( Name ) << "\n";
#endif

    return M_calculator.GetVariable( Name );
}

} // Namespace LifeV
