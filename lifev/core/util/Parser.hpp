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

#ifndef Parser_H
#define Parser_H 1

#include <lifev/core/util/LifeDebug.hpp>
#include <lifev/core/util/ParserSpiritGrammar.hpp>

namespace LifeV
{

//! Parser - A string parser for algebraic expressions
/*!
 *  @author Cristiano Malossi
 *
 *  \c Parser is a general interface class for \c LifeV algebraic parsers.
 *  At the present time it works only with \c boost::spirit::qi.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *
 *  The syntax is very simple:
 *
 *  <CODE>
 *  Parser parser;<BR>
 *  parser.setString( "-sqrt(4)+1*2" );<BR>
 *  Real result = parser.evaluate();<BR>
 *  </CODE>
 *
 *  You can change the string at any time:
 *
 *  <CODE>
 *  parser.setString( "c=2; [0., c, c*c, c*c*c]" );<BR>
 *  Real result0 = parser.evaluate(0); // 0<BR>
 *  Real result1 = parser.evaluate(1); // c<BR>
 *  Real result2 = parser.evaluate(2); // c*c<BR>
 *  Real result3 = parser.evaluate(3); // c*c*c<BR>
 *  </CODE>
 *
 *  See \c ParserSpiritGrammar class for more details on the expression syntax.
 */
class Parser
{
public:

    //! @name Public Types
    //@{

    /*! @typedef stringsVector_Type */
    //! Type definition for the vector containing the string segments
    typedef std::vector< std::string >                       stringsVector_Type;

    /*! @typedef stringIterator_Type */
    //! Type definition for the iterator over the strings
    typedef std::string::const_iterator                      stringIterator_Type;

    /*! @typedef calculator_Type */
    //!Type definition for the parser interpreter
    typedef ParserSpiritGrammar< stringIterator_Type >       calculator_Type;

    /*! @typedef results_Type */
    //! Type definition for the results
    typedef calculator_Type::results_Type                    results_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor (it needs a manual call to setString)
    explicit Parser();

    //! Constructor
    /*!
     * @param string expression to parse
     */
    explicit Parser ( const std::string& string );

    //! Copy constructor
    /*!
     * @param parser Parser
     */
    explicit Parser ( const Parser& parser );

    //! Destructor
    virtual ~Parser() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * <b> Important note: this method is not working right now. </b>
     *
     * @param parser Parser
     * @return reference to a copy of the class
     */
    Parser& operator= ( const Parser& parser );

    //@}


    //! @name Methods
    //@{

    //! Evaluate the expression
    /*!
     * @param id expression index (starting from 0)
     * @return computed value
     */
    const Real& evaluate ( const ID& id = 0 );

    //! Count how many substrings are present in the string (utility for BCInterfaceFunctionParser)
    /*!
     * @param substring string to find
     * @return number of substring
     */
    UInt countSubstring ( const std::string& substring ) const;

    //! Clear all the variables.
    void clearVariables();

    //@}


    //! @name Set Methods
    //@{

    //! Set string function
    /*!
     * @param string Expression to evaluate
     * @param stringSeparator Separator identifier (default -> ";")
     */
    void setString ( const std::string& string, const std::string& stringSeparator = ";" );

    //! Set/replace a variable
    /*!
     * @param name name of the parameter
     * @param value value of the parameter
     */
    void setVariable ( const std::string& name, const Real& value );

    //@}


    //! @name Get Methods
    //@{

    //! Get variable
    /*!
     * @param name name of the parameter
     * @return value of the variable
     */
    const Real& variable ( const std::string& name );

    //@}

private:

    stringsVector_Type  M_strings;

    results_Type        M_results;

    calculator_Type     M_calculator;

    bool                M_evaluate;
};

} // Namespace LifeV

#endif /* Parser_H */
