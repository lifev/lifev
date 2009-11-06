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
 \file SpiritParser.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \author Gilles Fourestey  <gilles.fourestey@epfl.ch>
 \date 2009-04-07
 */

#ifndef __SpiritParser_H
#define __SpiritParser_H 1

// SpiritParser need optimization level at least = -O1
//#pragma GCC optimization_level 2
//#pragma GCC optimize ("-O2")     // Should work with GCC v. 4.4+

#include <boost/algorithm/string.hpp>
#include <boost/spirit.hpp>
#include <boost/spirit/phoenix/binders.hpp>

#include <map>
#include <iomanip>
#include <string>

namespace LifeV {

//! SpiritCalculator - A string parser calculator based on \c boost::spirit
/*!
 *  @author(s) Cristiano Malossi, Gilles Fourestey
 *
 *  \c SpiritCalculator is a \c boost::spirit based class to perform
 *  evaluation of \c std::string expressions.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *	Let's look at an example: suppose that you have this function:
 *
 *	[u,v,w] = f(x,y,z,t)
 *
 *	where
 *
 *	u(x)   = a*b*x
 *	v(x,y) = a/b*sqrt(x^2 + y^2)
 *	w(t)   = b*t;
 *
 *	with "a" and "b" constants such that a=5.12345, b=9.999999.
 *	You can use this syntax to implement it inside the spiritParser:
 *
 *	string = "a=5.12345 ; b=9.999999 ; (a*b*x, a/b*sqrt(x^2 + y^2), b*t)"
 *
 *	where semicolons (";") separate constants and commas (",") separate output functions.
 *
 *	NOTE:
 *  Currently SpiritParser works with the following operators:
 *  \verbatim
 *  +, -, *, /, ^, sqrt(), sin(), cos(), tan(), exp(), log(), log10(), >, <.
 *  \endverbatim
 *
 *  SpiritParser need optimization level at least -O1
 *
 *  \TODO Fix a problem when calling destructors - NOW FIXED USING REFERENCE INSTEAD OF SHARED_PTR;
 *  \TODO Find a better way to manage results (M_results, M_nResults, setResult(), ...)
 *  \TODO Avoid the use of ruleTheString (if possible)!
 *
 */
struct SpiritCalculator: boost::spirit::grammar< SpiritCalculator >
{
public:

    typedef boost::spirit::grammar< SpiritCalculator > super;

    typedef std::map< std::string, double >            variables_type;
    typedef std::vector< double >                      results_type;
    typedef std::vector< std::string >                 string_type;

    SpiritCalculator( variables_type& variables, results_type& results, unsigned int& nResults ) :
        super       (),
        M_variables ( variables ),
        M_results   ( results ),
        M_nResults  ( nResults )
    {
    }

    SpiritCalculator( const SpiritCalculator& calculator ) :
        super       ( calculator ),
        M_variables ( calculator.M_variables ),
        M_results   ( calculator.M_results ),
        M_nResults  ( calculator.M_nResults )
    {
    }

    SpiritCalculator& operator=( const SpiritCalculator& calculator )
    {
        if ( this != &calculator )
        {
            super::operator=( calculator );
            M_variables = calculator.M_variables;
            M_results   = calculator.M_results;
            M_nResults  = calculator.M_nResults;
        }

        return *this;
    }

    ~SpiritCalculator() {}

    // Closures
    struct assignment_closure: boost::spirit::closure< assignment_closure, std::string, double >
    {
        member1 name;
        member2 value;
    };

    struct value_closure: boost::spirit::closure< value_closure, double >
    {
        member1 value;
    };

    struct string_closure: boost::spirit::closure< string_closure, std::string >
    {
        member1 name;
    };

    template< typename ScannerT >
    struct definition
    {
    public:

        definition( const SpiritCalculator& self )
        {
            using namespace boost::spirit;
            using namespace phoenix;

            statement = ( assignment
                        | command
                        | operation
                        ) >> (end_p | ';');

            assignment = identifier[assignment.name = arg1] >> '=' >> expression[assignment.value = arg1]
                       [ bind(&SpiritCalculator::define)(self, assignment.name, assignment.value) ];

            command = as_lower_d["showme"][bind(&SpiritCalculator::showMeVariables)(self)];

            operation = ( expression[bind(&SpiritCalculator::setResult)(self, arg1)]
                        | expression_vector
                        );

            expression_vector = '(' >> ( expression[bind(&SpiritCalculator::setResult)(self, arg1)]
                                         >> *(',' >> expression[bind(&SpiritCalculator::setResult)(self, arg1)])
                                       ) >> ')';

            expression = term[expression.value = arg1] >> *( lowOperations[expression.value = arg1] );

            term = subterm[term.value = arg1] >> *( midOperations[term.value = arg1] );

            subterm = factor[subterm.value = arg1] >> *( topOperations[subterm.value = arg1] );

            lowOperations = ('+' >> term[lowOperations.value += arg1])
                          | ('-' >> term[lowOperations.value -= arg1]);

            midOperations = ('*' >> subterm[midOperations.value *= arg1])
                          | ('/' >> subterm[midOperations.value /= arg1]);

            topOperations =  ('^' >> factor[topOperations.value = bind(&SpiritCalculator::pow)(self, topOperations.value, arg1)])
                           | ('>' >> factor[topOperations.value = bind(&SpiritCalculator::more)(self, topOperations.value, arg1)])
                           | ('<' >> factor[topOperations.value = bind(&SpiritCalculator::less)(self, topOperations.value, arg1)]);

            factor = ( literal[factor.value = arg1]
                     | function[factor.value = arg1]
                     | identifier[factor.value = bind(&SpiritCalculator::lookup)(self, arg1)]
                     | group[factor.value = arg1]
                     );

            literal = longest_d[ int_p[literal.value = arg1] | real_p[literal.value = arg1] ];

            function = ( ("S" >> group[function.value = bind(&SpiritCalculator::sin)(self, arg1)])
                       | ('C' >> group[function.value = bind(&SpiritCalculator::cos)(self, arg1)])
                       | ('T' >> group[function.value = bind(&SpiritCalculator::tan)(self, arg1)])
                       | ('Q' >> group[function.value = bind(&SpiritCalculator::sqrt)(self, arg1)])
                       | ('E' >> group[function.value = bind(&SpiritCalculator::exp)(self, arg1)])
                       | ('L' >> group[function.value = bind(&SpiritCalculator::log)(self, arg1)])
                       | ('M' >> group[function.value = bind(&SpiritCalculator::log10)(self, arg1)])
                       );
            identifier = lexeme_d[ ( alpha_p | '_') >> *( alnum_p | '_') ]
                                             [identifier.name = construct_<std::string>(arg1, arg2)];

            group = '(' >> expression[group.value = arg1] >> ')';
    	}

        boost::spirit::rule< ScannerT > const&
        start() const
        {
            return statement;
        }

        boost::spirit::rule< ScannerT > statement, command, operation, expression_vector;
        boost::spirit::rule< ScannerT, assignment_closure::context_t > assignment;
        boost::spirit::rule< ScannerT, string_closure::context_t > identifier;
        boost::spirit::rule< ScannerT, value_closure::context_t > expression, topOperations,
                                                                  midOperations, lowOperations,
                                                                  term, subterm, factor, function,
                                                                  literal, group;
    };

    // Member functions that are called in semantic actions.
    inline void define( const std::string& name, const double& value ) const
    {
        M_variables[name] = value;
    }

    inline void showMeVariables() const
    {
        variables_type::const_iterator it;

        std::cout << "SpiritParser showMe: " << std::setprecision( 30 ) << std::endl;

        for ( it = M_variables.begin(); it != M_variables.end(); ++it )
            std::cout << "                     " << it->first << " = " << it->second << std::endl;
        std::cout << std::endl;
    }

    inline double lookup( const std::string& name ) const
    {
        if ( M_variables.find( name ) == M_variables.end() )
        {
            std::cerr << "!!! Warning: SpiritParser has undefined name " << name << " !!!" << std::endl;
            return 0.0;
        }
        else
            return M_variables.find( name )->second;
    }

    inline double pow( const double& a, const double& b ) const
    {
        return std::pow( a, b );
    }

    inline double more( const double& a, const double& b ) const
    {
        return a > b;
    }

    inline double less( const double& a, const double& b ) const
    {
        return a < b;
    }

    inline double sqrt( const double& a ) const
    {
        return std::sqrt( a );
    }

    inline double exp( const double& a ) const
    {
        return std::exp( a );
    }

    inline double log( const double& a ) const
    {
        return std::log( a );
    }

    inline double log10( const double& a ) const
    {
        return std::log10( a );
    }

    inline double sin( const double& a ) const
    {
        return std::sin( a );
    }

    inline double cos( const double& a ) const
    {
        return std::cos( a );
    }

    inline double tan( const double& a ) const
    {
        return std::tan( a );
    }

    inline void setResult( const double& result ) const
    {
        M_results[M_nResults] = result;
        M_nResults++;
    }

private:

    variables_type& M_variables;
    results_type&   M_results;
    unsigned int&   M_nResults;
};
//! SpiritParser - A LifeV interface for SpiritCalculator
/*!
 *  @author(s) Cristiano Malossi, Gilles Fourestey
 */
class SpiritParser
{
public:

    typedef SpiritCalculator::variables_type    variables_type;
    typedef SpiritCalculator::results_type      results_type;
    typedef SpiritCalculator::string_type       string_type;

    //! @name Constructors & Destructor
    //@{

    //! Empty string constructor (it needs a manual call to setString)
    /*!
     * \param applyRules		- use rules for notation
     */
    //! Constructor
    SpiritParser( const bool& applyRules = true );

    //! Constructor
    /*!
     * \param string			- expression to parse
     * \param applyRules		- use rules for notation
     */
    SpiritParser( const std::string& string, const bool& applyRules = true );

    //! Copy constructor
    /*!
     * \param Parser			- SpiritParser
     */
    SpiritParser( const SpiritParser& parser );

    //! Destructor
    ~SpiritParser() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * \param Parser - SpiritParser
     */
    SpiritParser& operator=( const SpiritParser& parser );

    /*! Set string function
     *
     * \param string - Expression to evaluate
     * \param stringSeparator - Separator identifier (default -> ";")
     */
    void setString( const std::string& string, const std::string& stringSeparator = ";" );

    /*! Set/replace a variable
     *
     * \param name - name of the parameter
     * \param value - value of the parameter
     */
    void setVariable( const std::string& name, const double& value );

    /*! Get variable
     *
     * \param name - name of the parameter
     */
    const double& getVariable( const std::string& name );

    /*! Evaluate the expression
     */
    const double& evaluate( const unsigned int& ID = 1 );

    /*! Count how many times a substring is present in the string (utility for BCInterfaceFunction)
     *
     * \param substring - string to find
     */
    unsigned int countSubstring( const std::string& substring );

    void showMeVariables( void ) const
    {
        M_calculator.showMeVariables();
    }

    //@}

private:

    string_type      M_strings;
    variables_type   M_variables;
    results_type     M_results;
    unsigned int     M_nResults;
    SpiritCalculator M_calculator;
    bool             M_applyRules;

    //! @name Private functions
    //@{

    //! Setup results
    inline void setupResults( void );

    //! Set default variables
    inline void setDefaultVariables( void );

    //! Apply rules to the string
    inline void ruleTheString( std::string& string );

    //@}

};

} // Namespace LifeV

#endif /* __SpiritParser_H */
