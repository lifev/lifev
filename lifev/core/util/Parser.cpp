//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the Parser interface
 *
 *  @date 07-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/core/util/Parser.hpp>

namespace LifeV
{

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
    debugStream ( 5030 ) << "Parser::Parser" << "\n";
#endif

    M_calculator.setDefaultVariables();
}

Parser::Parser ( const std::string& string ) :
    M_strings       (),
    M_results       (),
    M_calculator    (),
    M_evaluate      ( true )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5030 ) << "Parser::Parser( string )" << "\n";
#endif

    M_calculator.setDefaultVariables();
    setString ( string );
}

Parser::Parser ( const Parser& parser ) :
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
Parser::operator= ( const Parser& parser )
{
    if ( this != &parser )
    {
        std::cerr << "!!! ERROR: Operator= not working !!!" << std::endl;
        std::exit ( EXIT_FAILURE );

        M_strings    = parser.M_strings;
        M_results    = parser.M_results;
        //M_calculator = parser.M_calculator; //NOT WORKING!!!
        M_evaluate   = parser.M_evaluate;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
const Real&
Parser::evaluate ( const ID& id )
{
    if ( M_evaluate )
    {
        M_results.clear();
        stringIterator_Type start, end;

        for ( UInt i (0); i < M_strings.size(); ++i )
        {
            start = M_strings[i].begin();
            end   = M_strings[i].end();
#ifdef HAVE_BOOST_SPIRIT_QI
#ifdef ENABLE_SPIRIT_PARSER
            qi::phrase_parse ( start, end, M_calculator, ascii::space, M_results );
#else
            std::cerr << "!!! ERROR: The Boost Spirit parser has been disabled !!!" << std::endl;
            std::exit ( EXIT_FAILURE );
#endif /* ENABLE_SPIRIT_PARSER */
#else
            std::cerr << "!!! ERROR: Boost version < 1.41 !!!" << std::endl;
            std::exit ( EXIT_FAILURE );
#endif
        }
        M_evaluate = false;
    }

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5030 ) << "Parser::evaluate          results[ " << id << "]: " << M_results[id] << "\n";
#endif

    return M_results[id];
}

UInt
Parser::countSubstring ( const std::string& substring ) const
{
    UInt counter ( 0 );
    std::string::size_type position ( 0 );

    for ( ;; )
    {
        position = M_strings.back().find ( substring, position );

        if ( position == std::string::npos )
        {
            break;
        }

        ++counter;
        position += substring.length(); // start next search after this substring
    }

    return counter;
}

void
Parser::clearVariables()
{
    M_calculator.clearVariables();
    M_evaluate = true;
}

// ===================================================
// Set Methods
// ===================================================
void
Parser::setString ( const std::string& string, const std::string& stringSeparator )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5030 ) << "Parser::setString         strings: " << string << "\n";
#endif

    M_strings.clear();
    boost::split ( M_strings, string, boost::is_any_of ( stringSeparator ) );

    //Remove white space to speed up the parser
    for ( UInt i = 0; i < M_strings.size(); ++i )
    {
        boost::replace_all ( M_strings[i], " ", "" );
    }

    //Reserve the space for results
    M_results.clear();
    M_results.reserve ( countSubstring ( "," ) + 1 );

    M_evaluate = true;
}

void
Parser::setVariable ( const std::string& name, const Real& value )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5030 ) << "Parser::setVariable       variables[" << name << "]: " << value << "\n";
#endif

    M_calculator.setVariable ( name, value);

    M_evaluate = true;
}

// ===================================================
// Get Methods
// ===================================================
const Real&
Parser::variable ( const std::string& name )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5030 ) << "Parser::variable          variables[" << name << "]: " << M_calculator.variable ( name ) << "\n";
#endif

    return M_calculator.variable ( name );
}

} // Namespace LifeV
