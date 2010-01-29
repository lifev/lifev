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
 * @brief Parser_SpiritClassic
 *
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 * @author Gilles Fourestey  <gilles.fourestey@epfl.ch>
 * @date 07-04-2009
 */

#ifndef Parser_SpiritClassic_H
#define Parser_SpiritClassic_H 1

#include <lifemc/lifecore/Parser_Definitions.hpp>

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
 *	You can use this syntax to implement it inside the Parser_SpiritClassic:
 *
 *	string = "a=5.12345 ; b=9.999999 ; (a*b*x, a/b*sqrt(x^2 + y^2), b*t)"
 *
 *	where semicolons (";") separate constants and commas (",") separate output functions.
 *
 *	NOTE:
 *  Currently Parser_SpiritClassic works with the following operators:
 *  \verbatim
 *  +, -, *, /, ^, sqrt(), sin(), cos(), tan(), exp(), log(), log10(), >, <.
 *  \endverbatim
 *
 *  Parser_SpiritClassic need optimization level at least -O1
 *
 *  \TODO Fix a problem when calling destructors - NOW FIXED USING REFERENCE INSTEAD OF SHARED_PTR;
 *  \TODO Find a better way to manage results (M_results, M_nResults, setResult(), ...)
 *  \TODO Avoid the use of ruleTheString (if possible)!
 *
 */

class SpiritCalculator : public spirit::grammar< SpiritCalculator >
{
public:

    //! @name Constructors & Destructor
    //@{

    SpiritCalculator( Variables_Type& variables, Results_Type& results, UInt& nResults ) :
        spirit::grammar< SpiritCalculator > (),
        M_variables                         ( variables ),
        M_results                           ( results ),
        M_nResults                          ( nResults )
    {
    }

    SpiritCalculator( const SpiritCalculator& calculator ) :
        spirit::grammar< SpiritCalculator > ( calculator ),
        M_variables                         ( calculator.M_variables ),
        M_results                           ( calculator.M_results ),
        M_nResults                          ( calculator.M_nResults )
    {
    }

    SpiritCalculator& operator=( const SpiritCalculator& calculator )
    {
        if ( this != &calculator )
        {
            spirit::grammar< SpiritCalculator >::operator=( calculator );
            M_variables = calculator.M_variables;
            M_results   = calculator.M_results;
            M_nResults  = calculator.M_nResults;
        }

        return *this;
    }

    ~SpiritCalculator() {}

    //@}

    // Closures
    struct assignment_closure: spirit::closure< assignment_closure, std::string, Real >
    {
        member1 name;
        member2 value;
    };

    struct value_closure: spirit::closure< value_closure, Real >
    {
        member1 value;
    };

    struct string_closure: spirit::closure< string_closure, std::string >
    {
        member1 name;
    };

    template< typename ScannerT >
    struct definition
    {
    public:

        definition( const SpiritCalculator& self )
        {
            using namespace boost::spirit::classic;
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

        //boost::spirit::rule< ScannerT > const&
        spirit::rule< ScannerT > const&
        start() const
        {
            return statement;
        }

        spirit::rule< ScannerT > statement, command, operation, expression_vector;
        spirit::rule< ScannerT, assignment_closure::context_t > assignment;
        spirit::rule< ScannerT, string_closure::context_t > identifier;
        spirit::rule< ScannerT, value_closure::context_t > expression, topOperations,
                                                           midOperations, lowOperations,
                                                           term, subterm, factor, function,
                                                           literal, group;
    };

    //! @name Methods
    //@{

    // Member functions that are called in semantic actions.
    inline void define( const std::string& name, const Real& value ) const
    {
        M_variables[name] = value;
    }

    inline void showMeVariables() const
    {
        Variables_Type::const_iterator it;

        std::cout << "Parser_SpiritClassic showMe: " << std::setprecision( 30 ) << std::endl;

        for ( it = M_variables.begin(); it != M_variables.end(); ++it )
            std::cout << "                     " << it->first << " = " << it->second << std::endl;
        std::cout << std::endl;
    }

    inline Real lookup( const std::string& name ) const
    {
        if ( M_variables.find( name ) == M_variables.end() )
        {
            std::cerr << "!!! Warning: Parser_SpiritClassic has undefined name " << name << " !!!" << std::endl;
            return 0.0;
        }
        else
            return M_variables.find( name )->second;
    }

    inline Real pow( const Real& a, const Real& b ) const
    {
        return std::pow( a, b );
    }

    inline Real more( const Real& a, const Real& b ) const
    {
        return a > b;
    }

    inline Real less( const Real& a, const Real& b ) const
    {
        return a < b;
    }

    inline Real sqrt( const Real& a ) const
    {
        return std::sqrt( a );
    }

    inline Real exp( const Real& a ) const
    {
        return std::exp( a );
    }

    inline Real log( const Real& a ) const
    {
        return std::log( a );
    }

    inline Real log10( const Real& a ) const
    {
        return std::log10( a );
    }

    inline Real sin( const Real& a ) const
    {
        return std::sin( a );
    }

    inline Real cos( const Real& a ) const
    {
        return std::cos( a );
    }

    inline Real tan( const Real& a ) const
    {
        return std::tan( a );
    }

    inline void setResult( const Real& result ) const
    {
        M_results[M_nResults] = result;
        M_nResults++;
    }

    //@}

private:

    Variables_Type& M_variables;
    Results_Type&   M_results;
    UInt&           M_nResults;
};

} // Namespace LifeV

#endif /* Parser_SpiritClassic_H */
