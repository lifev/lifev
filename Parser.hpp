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

#ifndef Parser_H
#define Parser_H 1

#include <lifemc/lifecore/Parser_Definitions.hpp>
#include <lifemc/lifecore/Parser_SpiritGrammar.hpp>

namespace LifeV {

//! Parser - A string parser for algebraic expressions
/*!
 *  @author(s) Cristiano Malossi, Gilles Fourestey
 */
class Parser
{
public:

    typedef Parser_SpiritGrammar<String_Iterator>        Calculator_Type;

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor (it needs a manual call to setString)
    Parser();

    //! Constructor
    /*!
     * @param String expression to parse
     */
    Parser( const std::string& String );

    //! Copy constructor
    /*!
     * @param Parser Parser
     */
    Parser( const Parser& Parser );

    //! Destructor
    ~Parser() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param Parser Parser
     * @return reference to a copy of the class
     */
    Parser& operator=( const Parser& Parser );

    //@}


    //! @name Methods
    //@{

    /*! Evaluate the expression
     * @param ID expression index (starting from 1)
     * @return computed value
     */
    const Real& Evaluate( const UInt& ID = 1 );

    /*! Count how many times a substring is present in the string (utility for BCInterfaceFunction)
     *
     * @param Substring string to find
     * @return number of substring
     */
    UInt CountSubstring( const std::string& Substring );

    //! Clear all the variables.
    void ClearVariables();

    //@}

    //! @name Set Methods
    //@{

    /*! Set string function
     *
     * @param String Expression to evaluate
     * @param StringSeparator Separator identifier (default -> ";")
     */
    void SetString( const std::string& String, const std::string& StringSeparator = ";" );

    /*! Set/replace a variable
     *
     * @param Name name of the parameter
     * @param Value value of the parameter
     */
    void SetVariable( const std::string& Name, const Real& Value );

    //@}


    //! @name Get Methods
    //@{

    /*! Get variable
     *
     * @param Name name of the parameter
     * @return value of the variable
     */
    const Real& GetVariable( const std::string& Name );

    //@}
private:

    Strings_Type        M_strings;
    Results_Type        M_results;

    Calculator_Type     M_calculator;

    bool                M_evaluate;

};

} // Namespace LifeV

#endif /* Parser_H */
