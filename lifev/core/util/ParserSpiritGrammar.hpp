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
 *  @brief File containing the Boost Spirit parser grammar
 *
 *  @date 05-02-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef Parser_SpiritGrammar_H
#define Parser_SpiritGrammar_H 1

#include <lifev/core/util/ParserDefinitions.hpp>

namespace LifeV
{

#if ( !defined(HAVE_BOOST_SPIRIT_QI) || !defined(ENABLE_SPIRIT_PARSER) )

/// @cond
//! ParserSpiritGrammar - An empty implementation for boost version < 1.41
template < typename IteratorType = std::string::const_iterator, typename ResultsType = std::vector < Real > >
class ParserSpiritGrammar
{
public:

    typedef IteratorType                                        iterator_Type;
    typedef boost::iterator_range< iterator_Type >              iteratorRange_Type;
    typedef ResultsType                                         results_Type;

    ParserSpiritGrammar() : M_real (0.) {}
    ParserSpiritGrammar ( const ParserSpiritGrammar& ) : M_real (0.) {}
    ~ParserSpiritGrammar() {}

    ParserSpiritGrammar& operator= ( const ParserSpiritGrammar& )
    {
        return *this;
    }

    void clearVariables() {}

    void setDefaultVariables() {}
    void setVariable ( const std::string&, const Real& ) {}

    Real& variable ( const std::string& )
    {
        return M_real;
    }

private:

    Real M_real;
};
/// @endcond

#else

//! ParserSpiritGrammar - A string parser grammar based on \c boost::spirit::qi
/*!
 *  @author(s) Cristiano Malossi
 *
 *  \c ParserSpiritGrammar is a \c boost::spirit::qi based class to perform
 *  evaluation of \c std::string expressions.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *  Let's consider the following example: suppose that we have this function:
 *
 *  [u,v,w] = f(x,y,z,t)
 *
 *  where
 *
 *  u(x)   = a*b*x
 *  v(x,y) = a/b*sqrt(x^2 + y^2)
 *  w(t)   = b*t;
 *
 *  with "a" and "b" constants such that a=5.12345, b=9.999999.
 *
 *  To evaluate function f(x,y,z,t), we use this syntax:
 *
 *  <CODE>
 *  string = "a=5.12345 ; b=9.999999 ; (a*b*x, a/b*sqrt(x^2 + y^2), b*t)"
 *  </CODE>
 *
 *  where semicolons (";") separate constants and commas (",") separate output functions.
 *
 *  NOTE:
 *  Currently ParserSpiritGrammar works with the following operators:
 *  \verbatim
 *  +, -, *, /, ^, sqrt(), sin(), cos(), tan(), exp(), log(), log10(), >, <.
 *  \endverbatim
 *
 */
template < typename IteratorType = std::string::const_iterator, typename ResultsType = std::vector < Real > >
class ParserSpiritGrammar : public qi::grammar< IteratorType, ResultsType(), ascii::space_type >
{
public:

    //! @name Public Types
    //@{

    /*! @typedef iterator_Type */
    //! Type definition for the iterator
    typedef IteratorType                                                 iterator_Type;

    /*! @typedef iteratorRange_Type */
    //!Type definition for the iterator range
    typedef boost::iterator_range< iterator_Type >                       iteratorRange_Type;

    /*! @typedef results_Type */
    //! Type definition for the results
    typedef ResultsType                                                  results_Type;

    /*! @typedef qiRuleVoid_Type */
    //! Type definition for the void rule of the parser
    typedef qi::rule< iterator_Type, void(), ascii::space_type >         qiRuleVoid_Type;

    /*! @typedef qiRuleReal_Type */
    //! Type definition for the Real rule of the parser
    typedef qi::rule< iterator_Type, Real(), ascii::space_type >         qiRuleReal_Type;

    /*! @typedef qiRuleResults_Type */
    //! Type definition for the result rule of the parser
    typedef qi::rule< iterator_Type, results_Type(), ascii::space_type > qiRuleResults_Type;

    /*! @typedef qiSymbolReal_Type */
    //! Type definition for the Real symbol of the parser
    typedef qi::symbols<char, Real >                                     qiSymbolReal_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ParserSpiritGrammar();

    //! Copy constructor
    /*!
     * @param ParserSpiritGrammar ParserSpiritGrammar
     */
    explicit ParserSpiritGrammar ( const ParserSpiritGrammar& spiritGrammar );

    //! Destructor
    virtual ~ParserSpiritGrammar() {}

    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param SpiritGrammar ParserSpiritGrammar
     * @return reference to a copy of the class
     */
    ParserSpiritGrammar& operator= ( const ParserSpiritGrammar& spiritGrammar );

    //@}


    //! @name Methods
    //@{

    //! Assign a variable using a \c boost::iterator_range
    /*!
     * @param stringIteratorRange name of the parameter
     * @param value value of the parameter
     */
    void assignVariable ( const iteratorRange_Type& stringIteratorRange, const Real& value )
    {
        setVariable ( std::string ( stringIteratorRange.begin(), stringIteratorRange.end() ), value );
    }

    //! Clear all the variables.
    void clearVariables()
    {
        M_variable.clear();
    }

    /* TODO Implement the showMe and use the M_command to enable it.
     * Note that the M_command is a rule to execute commands given in the strings, e.g.:
     * string = "a=1; b=2; showMe; a*t+b" means that first we set the variables "a" and "b",
     * then we call the command showMe (to display the variables) and finally we evaluate
     * the function, i.e., "f=a*t+b", where t is the time.
    */
    //    //! Show all the variables
    //    void ShowMe() { std::cout << "To be implemented" << std::endl; }

    //@}


    //! @name Set Methods
    //@{

    //! Set default variables
    void setDefaultVariables();

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
    Real& variable ( const std::string& name )
    {
        return M_variable.at ( name );
    }

    //@}

private:

    //! @name Private Methods
    //@

    //! Phoenix wrapper for \c std::sin
    /*!
     * @param value input value
     * @return sin( value )
     */
    Real sin ( const Real& value ) const
    {
        return std::sin ( value );
    }

    //! Phoenix wrapper for \c std::cos
    /*!
     * @param value input value
     * @return cos( value )
     */
    Real cos ( const Real& value ) const
    {
        return std::cos ( value );
    }

    //! Phoenix wrapper for \c std::tan
    /*!
     * @param value input value
     * @return tan( value )
     */
    Real tan ( const Real& value ) const
    {
        return std::tan ( value );
    }

    //! Phoenix wrapper for \c std::pow
    /*!
     * @param base input base
     * @param exponent input exponent
     * @return pow( base, exponent )
     */
    Real pow ( const Real& Base, const Real& Exponent ) const
    {
        return std::pow ( Base, Exponent );
    }

    //! Phoenix wrapper for \c std::sqrt
    /*!
     * @param value input value
     * @return sqrt( value )
     */
    Real sqrt ( const Real& value ) const
    {
        return std::sqrt ( value );
    }

    //! Phoenix wrapper for \c std::exp
    /*!
     * @param value input value
     * @return exp( value )
     */
    Real exp ( const Real& value ) const
    {
        return std::exp ( value );
    }

    //! Phoenix wrapper for \c std::log
    /*!
     * @param value input value
     * @return log( value )
     */
    Real log ( const Real& value ) const
    {
        return std::log ( value );
    }

    //! Phoenix wrapper for \c std::log10
    /*!
     * @param value input value
     * @return log10( value )
     */
    Real log10 ( const Real& value ) const
    {
        return std::log10 ( value );
    }

    //@}

    qiRuleResults_Type      M_start;

    qiRuleVoid_Type         M_assignment;
    //    qiRuleVoid_Type         M_command;

    qiRuleReal_Type         M_expression;
    qiRuleReal_Type         M_compare;
    qiRuleReal_Type         M_plusMinus;
    qiRuleReal_Type         M_multiplyDivide;
    qiRuleReal_Type         M_elevate;
    qiRuleReal_Type         M_element;
    qiRuleReal_Type         M_number;
    qiRuleReal_Type         M_function;
    qiRuleReal_Type         M_group;

    qiSymbolReal_Type       M_variable;
};



// ===================================================
// Constructors & Destructor
// ===================================================
template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >::ParserSpiritGrammar() :
    ParserSpiritGrammar::base_type ( M_start ),
    M_start                               (),
    M_assignment                          (),
    //        M_command                             (),
    M_expression                          (),
    M_compare                             (),
    M_plusMinus                           (),
    M_multiplyDivide                      (),
    M_elevate                             (),
    M_element                             (),
    M_number                              (),
    M_function                            (),
    M_group                               (),
    M_variable                            ()
{
    M_start =
        (
            M_assignment
            //        |   M_command
            |  ( -qi::lit ('[') >> M_expression % ',' >> -qi::lit (']') )
        )
        ;

    M_assignment =
        (
            qi::raw[qi::lexeme[ (qi::alpha | '_') >> * (qi::alnum | '_')]]
            >>   qi::lit ('=')
            >>   M_expression
        )                                        [phoenix::bind (&ParserSpiritGrammar::assignVariable, this, qi::_1, qi::_2)]
        ;
    /*
    M_command =
        qi::lit("ShowMe")[phoenix::bind(&ParserSpiritGrammar::ShowMe, this)]
        ;
    */

    M_expression =
        *M_compare                                [qi::_val = qi::_1]
        ;

    M_compare =
        M_plusMinus                              [qi::_val = qi::_1]
        >> * (
            qi::lit (">=") >> M_plusMinus     [qi::_val = qi::_val >= qi::_1]
            |   qi::lit ("<=") >> M_plusMinus     [qi::_val = qi::_val <= qi::_1]
            |   qi::lit (">") >> M_plusMinus      [qi::_val = qi::_val > qi::_1]
            |   qi::lit ("<") >> M_plusMinus      [qi::_val = qi::_val < qi::_1]
        )
        ;

    M_plusMinus =
        M_multiplyDivide                         [qi::_val = qi::_1]
        >> * (
            qi::lit ('+') >> M_multiplyDivide [qi::_val += qi::_1]
            |   qi::lit ('-') >> M_multiplyDivide [qi::_val -= qi::_1]
        )
        ;

    M_multiplyDivide =
        M_elevate                                [qi::_val = qi::_1]
        >> * (
            qi::lit ('*') >> M_elevate        [qi::_val *= qi::_1]
            |   qi::lit ('/') >> M_elevate        [qi::_val /= qi::_1]
        )
        ;

    M_elevate =
        (
            qi::lit ('-') >> M_element            [qi::_val = qi::_1]
            >>  (
                qi::lit ('^') >> M_element        [qi::_val = -phoenix::bind (&ParserSpiritGrammar::pow,
                                                                              this, qi::_val, qi::_1)]
            )
            >> * (
                qi::lit ('^') >> M_element        [qi::_val = phoenix::bind (&ParserSpiritGrammar::pow,
                                                                             this, qi::_val, qi::_1)]
            )
        )
        |
        (
            M_element                            [qi::_val = qi::_1]
            >> * (
                qi::lit ('^') >> M_element        [qi::_val = phoenix::bind (&ParserSpiritGrammar::pow,
                                                                             this, qi::_val, qi::_1)]
            )
        )
        ;

    M_element =
        (
            qi::lit ('-') >> M_element            [qi::_val = -qi::_1]
        )
        |
        (
            M_number                             [qi::_val = qi::_1]
            |   M_function                           [qi::_val = qi::_1]
            |   M_variable                           [qi::_val = qi::_1]
            |   M_group                              [qi::_val = qi::_1]
        )
        ;

    M_number =
        (
            qi::double_
            //            ||   ('.' >> qi::double_)
            //            >>   -('.' >> qi::double_) | ('.' >> qi::double_)
        )
        ;

    M_function =
        (
            qi::lit ("sin")   >> M_group          [qi::_val = phoenix::bind (&ParserSpiritGrammar::sin, this,   qi::_1)]
            |   qi::lit ("cos")   >> M_group          [qi::_val = phoenix::bind (&ParserSpiritGrammar::cos, this,   qi::_1)]
            |   qi::lit ("tan")   >> M_group          [qi::_val = phoenix::bind (&ParserSpiritGrammar::tan, this,   qi::_1)]
            |   qi::lit ("sqrt")  >> M_group          [qi::_val = phoenix::bind (&ParserSpiritGrammar::sqrt, this,  qi::_1)]
            |   qi::lit ("exp")   >> M_group          [qi::_val = phoenix::bind (&ParserSpiritGrammar::exp, this,   qi::_1)]
            |   qi::lit ("log")   >> M_group          [qi::_val = phoenix::bind (&ParserSpiritGrammar::log, this,   qi::_1)]
            |   qi::lit ("log10") >> M_group          [qi::_val = phoenix::bind (&ParserSpiritGrammar::log10, this, qi::_1)]
        )
        ;

    M_group =
        (
            '('
            >>  M_expression                                 [qi::_val = qi::_1]
            >>  ')'
        )
        ;
}

template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >::ParserSpiritGrammar ( const ParserSpiritGrammar& spiritGrammar ) :
    ParserSpiritGrammar::base_type     ( spiritGrammar.M_start ),
    M_start                               ( spiritGrammar.M_start ),
    M_assignment                          ( spiritGrammar.M_assignment ),
    //        M_command                             ( spiritGrammar.M_command ),
    M_expression                          ( spiritGrammar.M_expression ),
    M_compare                             ( spiritGrammar.M_compare ),
    M_plusMinus                           ( spiritGrammar.M_plusMinus ),
    M_multiplyDivide                      ( spiritGrammar.M_multiplyDivide ),
    M_elevate                             ( spiritGrammar.M_elevate ),
    M_element                             ( spiritGrammar.M_element ),
    M_number                              ( spiritGrammar.M_number ),
    M_function                            ( spiritGrammar.M_function ),
    M_group                               ( spiritGrammar.M_group ),
    M_variable                            ( spiritGrammar.M_variable )
{
}

// ===================================================
// Operators
// ===================================================
template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >&
ParserSpiritGrammar< IteratorType, ResultsType >::operator= ( const ParserSpiritGrammar& spiritGrammar )
{
    if ( this != &spiritGrammar )
    {
        ParserSpiritGrammar::base_type::operator= ( spiritGrammar.M_start );
        M_start                               = spiritGrammar.M_start;
        M_assignment                          = spiritGrammar.M_assignment;
        //        M_command                             = spiritGrammar.M_command;
        M_expression                          = spiritGrammar.M_expression;
        M_compare                             = spiritGrammar.M_compare;
        M_plusMinus                           = spiritGrammar.M_plusMinus;
        M_multiplyDivide                      = spiritGrammar.M_multiplyDivide;
        M_elevate                             = spiritGrammar.M_elevate;
        M_element                             = spiritGrammar.M_element;
        M_number                              = spiritGrammar.M_number;
        M_function                            = spiritGrammar.M_function;
        M_group                               = spiritGrammar.M_group;
        M_variable                            = spiritGrammar.M_variable;
    }

    return *this;
}

// ===================================================
// Set Methods
// ===================================================
template < typename IteratorType, typename ResultsType >
inline void
ParserSpiritGrammar< IteratorType, ResultsType >::setDefaultVariables()
{
    M_variable.add ( "pi" , M_PI );
    M_variable.add ( "e", M_E );
}

template < typename IteratorType, typename ResultsType >
inline void
ParserSpiritGrammar< IteratorType, ResultsType >::setVariable ( const std::string& name, const Real& value )
{
    Real* p = M_variable.find ( name );
    if ( p != 0 )
    {
        *p = value;
    }
    else
    {
        M_variable.add ( name, value );
    }
}

#endif /* HAVE_BOOST_SPIRIT_QI || !ENABLE_SPIRIT_PARSER */

} // Namespace LifeV

#endif /* Parser_SpiritGrammar_H */
