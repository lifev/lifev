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
 *  @contributor Gilles Fourestey  <gilles.fourestey@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef Parser_H
#define Parser_H 1

#include <life/lifecore/debug.hpp>
#include <lifemc/lifecore/ParserDefinitions.hpp>
#include <lifemc/lifecore/ParserSpiritGrammar.hpp>

namespace LifeV
{

//! Parser - A string parser for algebraic expressions
/*!
 *  @author(s) Cristiano Malossi, Gilles Fourestey
 *
 *  See \c ParserSpiritGrammar class for more details.
 */
class Parser
{
public:

    //! @name Public Types
    //@{

    typedef std::vector< std::string >                       stringsVector_Type;
    typedef std::string::const_iterator                      stringIterator_Type;
    typedef ParserSpiritGrammar< stringIterator_Type >       calculator_Type;
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
    explicit Parser( const std::string& string );

    //! Copy constructor
    /*!
     * @param parser Parser
     */
    explicit Parser( const Parser& parser );

    //! Destructor
    virtual ~Parser() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param parser Parser
     * @return reference to a copy of the class
     */
    Parser& operator=( const Parser& parser );

    //@}


    //! @name Methods
    //@{

    /*! Evaluate the expression
     * @param id expression index (starting from 1)
     * @return computed value
     */
    const Real& evaluate( const ID& id = 1 );

    /*! Count how many times a substring is present in the string (utility for BCInterfaceFunction)
     *
     * @param substring string to find
     * @return number of substring
     */
    UInt countSubstring( const std::string& substring );

    //! Clear all the variables.
    void clearVariables();

    //@}


    //! @name Set Methods
    //@{

    /*! Set string function
     *
     * @param string Expression to evaluate
     * @param stringSeparator Separator identifier (default -> ";")
     */
    void setString( const std::string& string, const std::string& stringSeparator = ";" );

    /*! Set/replace a variable
     *
     * @param name name of the parameter
     * @param value value of the parameter
     */
    void setVariable( const std::string& name, const Real& value );

    //@}


    //! @name Get Methods
    //@{

    /*! Get variable
     *
     * @param name name of the parameter
     * @return value of the variable
     */
    const Real& variable( const std::string& name );

    //@}

private:

    stringsVector_Type  M_strings;

    results_Type        M_results;

    calculator_Type     M_calculator;

    bool                M_evaluate;
};

} // Namespace LifeV

#endif /* Parser_H */
