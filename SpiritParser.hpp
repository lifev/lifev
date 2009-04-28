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
   \file SpiritParser.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-07
 */

#ifndef __SpiritParser_H
#define __SpiritParser_H 1





// ===================================================
//! Include
// ===================================================
#include <life/lifecore/life.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/spirit.hpp>
#include <boost/spirit/phoenix/binders.hpp>





// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;





/*!
 * \class SpiritCalculator
 * \brief A string parser calculator based on \c boost::spirit
 *
 *  @author Cristiano Malossi
 *  @see
 *
 *  \c SpiritCalculator is a \c boost::spirit based class to perform
 *  evaluation of \c std::string expressions.
 *
 *  Currently it works with the following operators:
 *  \verbatim
 *  +, -, *, /, ^, sqrt(), sin(), cos(), tan(), exp(), log(), log10(), >, <.
 *  \endverbatim
 *
 */
struct SpiritCalculator : boost::spirit::grammar<SpiritCalculator>
//     :
//     public LifeV::Application
{
public:

	SpiritCalculator( std::map<std::string, Real>& variables, Real& result ) :
		M_variables	( variables ),
		M_result 	( result )
	{}

	SpiritCalculator( const SpiritCalculator& calculator ) :
		M_variables ( calculator.M_variables ),
		M_result 	( calculator.M_result )
	{}

	SpiritCalculator& operator=( const SpiritCalculator& calculator )
	{
	    if ( this != &calculator )
	    {
	    	M_variables		= calculator.M_variables;
	    	M_result 		= calculator.M_result;
	    }

		return *this;
	}

	// Closures
    struct assignment_closure : boost::spirit::closure<assignment_closure, std::string, Real>
    {
    	member1 name;
    	member2 value;
    };

    struct value_closure : boost::spirit::closure<value_closure, Real>
    {
        member1 value;
    };

    struct string_closure : boost::spirit::closure<string_closure, std::string>
    {
        member1 name;
    };

    template <typename ScannerT>
    struct definition
    {
    public:

        definition(SpiritCalculator const& self)
        {
    		using namespace boost::spirit;
    		using namespace phoenix;

    		identifier
    			= lexeme_d
    			[
    				( alpha_p | '_')
    			>> *( alnum_p | '_')
    			][identifier.name = construct_<std::string>(arg1, arg2)]
    			;

    		group
    			=  '('
    			>> expression[group.value = arg1]
    			>> ')'
    			;

    		assignment
				= identifier[assignment.name = arg1]
				>> '='
				>> expression[assignment.value = arg1][bind(&SpiritCalculator::define)(self, assignment.name, assignment.value)]
				;

    		statement
    			= (   assignment
    				| expression[bind(&SpiritCalculator::setResult)(self, arg1)]
    			  )
    			  >> (end_p | ',')
    			;

    		literal
    			= longest_d
    			[
    			  int_p[literal.value = arg1]
    			| real_p[literal.value = arg1]
    			]
    			;

    		function
				= ('Q' >> group[function.value = bind(&SpiritCalculator::sqrt)(self, arg1)])
				| ('E' >> group[function.value = bind(&SpiritCalculator::exp)(self, arg1)])
				| ('L' >> group[function.value = bind(&SpiritCalculator::log)(self, arg1)])
				| ('M' >> group[function.value = bind(&SpiritCalculator::log10)(self, arg1)])
				| ("S" >> group[function.value = bind(&SpiritCalculator::sin)(self, arg1)])
				| ('C' >> group[function.value = bind(&SpiritCalculator::cos)(self, arg1)])
				| ('T' >> group[function.value = bind(&SpiritCalculator::tan)(self, arg1)])
				;

    		factor
    			= literal[factor.value = arg1]
    			| group[factor.value = arg1]
    			| function[factor.value = arg1]
    			| identifier[factor.value = bind(&SpiritCalculator::lookup)(self, arg1)]
    			;

    		term
    			= factor[term.value = arg1]
    			>> *( ('*' >> factor[term.value *= arg1])
    				| ('/' >> factor[term.value /= arg1])
    				| ('^' >> factor[term.value = bind(&SpiritCalculator::pow)(self, term.value, arg1)])
					| ('>' >> factor[term.value = bind(&SpiritCalculator::more)(self, term.value, arg1)])
					| ('<' >> factor[term.value = bind(&SpiritCalculator::less)(self, term.value, arg1)])
    				)
    			;

    		expression
    			= term[expression.value = arg1]
    			>> *( ('+' >> term[expression.value += arg1])
    				| ('-' >> term[expression.value -= arg1])
    				)
    			;
    	}

        boost::spirit::rule<ScannerT> const&
        start() const { return statement; }

        boost::spirit::rule<ScannerT> statement;
        boost::spirit::rule<ScannerT, assignment_closure::context_t> assignment;
        boost::spirit::rule<ScannerT, string_closure::context_t> identifier;
        boost::spirit::rule<ScannerT,  value_closure::context_t> expression, term, factor, function, literal, group;
    };

    // Member functions that are called in semantic actions.
    void define(const std::string& name, const Real value) const
    {
    	M_variables[name] = value;
    }

    Real lookup(const std::string& name) const
    {
    	if (M_variables.find(name) == M_variables.end())
    	{
    	    std::cerr << "undefined name: " << name << '\n';
    	    return 0.0;
    	}
    	else
    	   return (*M_variables.find(name)).second;
    }

    Real pow( const Real a, const Real b ) const
    {
		return std::pow(a,b);
    }

    Real more( const Real a, const Real b ) const
    {
		return a > b;
    }

    Real less( const Real a, const Real b ) const
    {
		return a < b;
    }

    Real sqrt( const Real a ) const
    {
    	return std::sqrt(a);
    }

    Real exp( const Real a ) const
    {
    	return std::exp(a);
    }

    Real log( const Real a ) const
    {
    	return std::log(a);
    }

    Real log10( const Real a ) const
    {
    	return std::log10(a);
    }

    Real sin( const Real a ) const
    {
    	return std::sin(a);
    }

    Real cos( const Real a ) const
    {
    	return std::cos(a);
    }

    Real tan( const Real a ) const
    {
    	return std::tan(a);
    }

    void setResult( const Real R) const
    {
    	M_result = R;
    }

private:

    std::map<std::string, Real>&						M_variables;
	Real&												M_result;

};





/*!
 * \class SpiritParser
 * \brief A LifeV interface for SpiritCalculator
 *
 *  @author Cristiano Malossi
 *  @see
 */
class SpiritParser
//     :
//     public LifeV::Application
{
public:

	// ===================================================
	//! Public functions
	// ===================================================

    /** @name Constructors & Destructor
     */
    //@{

	//! Empty string constructor (need a manual call to setString)
	/*!
	 * \param applyRules - use rules for notation
     */
    //! Constructor
	SpiritParser( const bool& applyRules=true );

	//! Constructor
	/*!
	 * \param string - expression to parse
	 * \param applyRules - use rules for notation
     */
	SpiritParser( const std::string& string, const bool& applyRules=true );

	//! Copy constructor
	/*!
	 * \param Parser - SpiritParser
     */
	SpiritParser( const SpiritParser& parser );

	//! Operator =
	/*!
	 * \param Parser - SpiritParser
     */
	SpiritParser& operator=( const SpiritParser& parser );

    //! Destructor
    ~SpiritParser() {}

    //@}



    /** @name Method
     */

    //@{
    /*! Set string function
     *
     * \param string          - Expression to evaluate
     * \param stringSeparator - Separator identifier (default -> ",")
     */
    void setString( const std::string& string,  const std::string& stringSeparator="," );

    /*! Set/replace a variable
     *
     * \param name  - name of the parameter
     * \param value - value of the parameter
     */
    void setVariable( const std::string& name, const Real& value );

    /*! Evaluate the expression
     */
    Real& evaluate( void );

    //@}

private:

	// ===================================================
	//! Private variables
	// ===================================================

	std::vector<std::string>							M_strings;
    std::map<std::string, Real>							M_variables;
	Real												M_result;
	SpiritCalculator									M_calculator;
	bool												M_applyRules;



	// ===================================================
	//! Private functions
	// ===================================================

    /** @name Private functions
     */
    //@{

	//! Set default variables
	inline void setDefaultVariables( void );

	//! Apply rules to the string
	inline void ruleTheString( std::string& string );

    //@}

};

#endif /* __SpiritParser_H */
