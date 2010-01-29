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
#ifdef HAVE_BOOST_SPIRIT_CLASSIC
    #include <lifemc/lifecore/Parser_SpiritClassic.hpp>
#elif  HAVE_BOOST_SPIRIT_QI
    //#include <lifemc/lifecore/Parser_SpiritQI.hpp>
#endif

namespace LifeV {

//! Parser - A string parser for algebraic expressions
/*!
 *  @author(s) Cristiano Malossi, Gilles Fourestey
 */
class Parser
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Empty string constructor (it needs a manual call to setString)
    /*!
     * @param applyRules use rules for notation
     */
    //! Constructor
    Parser( const bool& applyRules = true );

    //! Constructor
    /*!
     * @param string expression to parse
     * @param applyRules use rules for notation
     */
    Parser( const std::string& string, const bool& applyRules = true );

    //! Copy constructor
    /*!
     * @param Parser Parser
     */
    Parser( const Parser& parser );

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
    Parser& operator=( const Parser& parser );

    //@}


    //! @name Methods
    //@{

    /*! Evaluate the expression
     * @param ID expression index (starting from 1)
     * @return computed value
     */
    const Real& evaluate( const UInt& ID = 1 );

    /*! Count how many times a substring is present in the string (utility for BCInterfaceFunction)
     *
     * @param substring string to find
     * @return number of substring
     */
    UInt countSubstring( const std::string& substring );

    void showMeVariables() const
    {
        M_calculator.showMeVariables();
    }

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
    const Real& getVariable( const std::string& name );

    //@}
private:

    Strings_Type     M_strings;
    Variables_Type   M_variables;
    Results_Type     M_results;
    UInt             M_nResults;
    SpiritCalculator M_calculator;
    bool             M_applyRules;

    //! @name Private functions
    //@{

    //! Setup results
    inline void setupResults();

    //! Set default variables
    inline void setDefaultVariables();

    //! Apply rules to the string
    inline void ruleTheString( std::string& string );

    //@}

};

} // Namespace LifeV

#endif /* Parser_H */
