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
 *  @brief File containing the Parser grammar
 *
 *  @date 05-02-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef Parser_SpiritGrammar_H
#define Parser_SpiritGrammar_H 1

#include <lifemc/lifecore/ParserDefinitions.hpp>

namespace LifeV
{

#ifdef HAVE_BOOST_SPIRIT_QI

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
 *  string = "a=5.12345 ; b=9.999999 ; (a*b*x, a/b*sqrt(x^2 + y^2), b*t)"
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

    typedef IteratorType                                        iterator_Type;
    typedef boost::iterator_range< iterator_Type >              iteratorRange_Type;
    typedef ResultsType                                         results_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit ParserSpiritGrammar();

    //! Copy constructor
    /*!
     * @param ParserSpiritGrammar ParserSpiritGrammar
     */
    explicit ParserSpiritGrammar( const ParserSpiritGrammar& spiritGrammar );

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
    ParserSpiritGrammar& operator=( const ParserSpiritGrammar& spiritGrammar );

    //@}


    //! @name Methods
    //@{

    /*! Assign a variable using a \c boost::iterator_range
     *
     * @param stringIteratorRange name of the parameter
     * @param value value of the parameter
     */
    void assignVariable( const iteratorRange_Type& stringIteratorRange, const Real& value ) { setVariable( std::string( stringIteratorRange.begin(), stringIteratorRange.end() ), value ); }

    //! Clear all the variables.
    void clearVariables() { Variable.clear(); }

    /*
        //! Show all the variables
        void ShowMe() { std::cout << "ShowMe called!" << std::endl; }
    */

    //@}


    //! @name Set Methods
    //@{

    //! Set default variables
    void setDefaultVariables();

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
    Real& variable( const std::string& name ) { return Variable.at( name ); }

    //@}

private:

    //! @name Private Methods
    //@

    /*! Phoenix wrapper for \c std::sin
     *
     * @param value input value
     * @return sin( value )
     */
    Real sin( const Real& value ) const { return std::sin( value ); }

    /*! Phoenix wrapper for \c std::cos
     *
     * @param value input value
     * @return cos( value )
     */
    Real cos( const Real& value ) const { return std::cos( value ); }

    /*! Phoenix wrapper for \c std::tan
     *
     * @param value input value
     * @return tan( value )
     */
    Real tan( const Real& value ) const { return std::tan( value ); }

    /*! Phoenix wrapper for \c std::pow
     *
     * @param base input base
     * @param exponent input exponent
     * @return pow( base, exponent )
     */
    Real pow( const Real& Base, const Real& Exponent ) const { return std::pow( Base, Exponent ); }

    /*! Phoenix wrapper for \c std::sqrt
     *
     * @param value input value
     * @return sqrt( value )
     */
    Real sqrt( const Real& value ) const { return std::sqrt( value ); }

    /*! Phoenix wrapper for \c std::exp
     *
     * @param value input value
     * @return exp( value )
     */
    Real exp( const Real& value ) const { return std::exp( value ); }

    /*! Phoenix wrapper for \c std::log
     *
     * @param value input value
     * @return log( value )
     */
    Real log( const Real& value ) const { return std::log( value ); }

    /*! Phoenix wrapper for \c std::log10
     *
     * @param value input value
     * @return log10( value )
     */
    Real log10( const Real& value ) const { return std::log10( value ); }

    //@}

    qi::rule< iterator_Type, results_Type(), ascii::space_type > Start;

    qi::rule< iterator_Type, void(), ascii::space_type >         Assignment;
//    qi::rule< iterator_Type, void(), ascii::space_type >         Command;

    qi::rule< iterator_Type, double(), ascii::space_type >       Expression;
    qi::rule< iterator_Type, double(), ascii::space_type >       Compare;
    qi::rule< iterator_Type, double(), ascii::space_type >       PlusMinus;
    qi::rule< iterator_Type, double(), ascii::space_type >       MultiplyDivide;
    qi::rule< iterator_Type, double(), ascii::space_type >       Elevate;
    qi::rule< iterator_Type, double(), ascii::space_type >       Element;
    qi::rule< iterator_Type, double(), ascii::space_type >       Number;
    qi::rule< iterator_Type, double(), ascii::space_type >       Function;
    qi::rule< iterator_Type, double(), ascii::space_type >       Group;

    qi::symbols<char, Real >                                     Variable;
};



// ===================================================
// Constructors & Destructor
// ===================================================
template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >::ParserSpiritGrammar() :
        ParserSpiritGrammar::base_type( Start ),
        Start                               (),
        Assignment                          (),
//        Command                             (),
        Expression                          (),
        Compare                             (),
        PlusMinus                           (),
        MultiplyDivide                      (),
        Elevate                             (),
        Element                             (),
        Number                              (),
        Function                            (),
        Group                               (),
        Variable                            ()
{
    Start
    =
        (
            Assignment
//        |  Command
        |  ( -qi::lit('[') >> Expression % ',' >> -qi::lit(']') )
        )
        ;

    Assignment
    =
        (    qi::raw[qi::lexeme[(qi::alpha | '_') >> *(qi::alnum | '_')]]
        >>   qi::lit('=')
        >>   Expression
        )              [phoenix::bind(&ParserSpiritGrammar::assignVariable,this, qi::_1, qi::_2)]
        ;
    /*
        Command
            =  qi::lit("ShowMe")[phoenix::bind(&ParserSpiritGrammar::ShowMe, this)]
        ;
    */

    Expression
    =  *Compare                                        [qi::_val = qi::_1]
       ;

    Compare
    =   PlusMinus                                      [qi::_val = qi::_1]
        >> *(
            qi::lit('>') >> PlusMinus                  [qi::_val = qi::_val > qi::_1]
        |   qi::lit('<') >> PlusMinus                  [qi::_val = qi::_val < qi::_1]
        )
        ;

    PlusMinus
    =   MultiplyDivide                                 [qi::_val = qi::_1]
        >> *(
            qi::lit('+') >> MultiplyDivide             [qi::_val += qi::_1]
        |   qi::lit('-') >> MultiplyDivide             [qi::_val -= qi::_1]
        )
        ;

    MultiplyDivide
    =   Elevate                                        [qi::_val = qi::_1]
        >> *(
            qi::lit('*') >> Elevate                    [qi::_val *= qi::_1]
        |   qi::lit('/') >> Elevate                    [qi::_val /= qi::_1]
        )
        ;

    Elevate
    =
        (
            qi::lit('-') >> Element                    [qi::_val = qi::_1]
            >>  (
                qi::lit('^') >> Element                [qi::_val = -phoenix::bind(&ParserSpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
            >> *(
                qi::lit('^') >> Element                [qi::_val = phoenix::bind(&ParserSpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
        )
        |
        (
            Element                                    [qi::_val = qi::_1]
            >> *(
                qi::lit('^') >> Element                [qi::_val = phoenix::bind(&ParserSpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
        )
        ;

    Element
    =
        (
            qi::lit('-') >> Element                    [qi::_val = -qi::_1]
        )
        |
        (
            Number                                     [qi::_val = qi::_1]
        |   Function                                   [qi::_val = qi::_1]
        |   Variable                                   [qi::_val = qi::_1]
        |   Group                                      [qi::_val = qi::_1]
        )
        ;

    Number
    =
        (
            qi::double_
//            ||   ('.' >> qi::double_)
//            >>   -('.' >> qi::double_) | ('.' >> qi::double_)
        )
        ;

    Function
    =
        (
            qi::lit("sin")   >> Group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::sin, this,   qi::_1)]
        |   qi::lit("cos")   >> Group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::cos, this,   qi::_1)]
        |   qi::lit("tan")   >> Group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::tan, this,   qi::_1)]
        |   qi::lit("sqrt")  >> Group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::sqrt, this,  qi::_1)]
        |   qi::lit("exp")   >> Group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::exp, this,   qi::_1)]
        |   qi::lit("log")   >> Group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::log, this,   qi::_1)]
        |   qi::lit("log10") >> Group                  [qi::_val = phoenix::bind(&ParserSpiritGrammar::log10, this, qi::_1)]
        )
        ;

    Group
    =
        (
            '('
        >>  Expression                                 [qi::_val = qi::_1]
        >>  ')'
        )
        ;
}

template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >::ParserSpiritGrammar( const ParserSpiritGrammar& spiritGrammar ) :
        ParserSpiritGrammar::base_type     ( spiritGrammar.Start ),
        Start                               ( spiritGrammar.Start ),
        Assignment                          ( spiritGrammar.Assignment ),
//    Command                             ( spiritGrammar.Command ),
        Expression                          ( spiritGrammar.Expression ),
        Compare                             ( spiritGrammar.Compare ),
        PlusMinus                           ( spiritGrammar.PlusMinus ),
        MultiplyDivide                      ( spiritGrammar.MultiplyDivide ),
        Elevate                             ( spiritGrammar.Elevate ),
        Element                             ( spiritGrammar.Element ),
        Number                              ( spiritGrammar.Number ),
        Function                            ( spiritGrammar.Function ),
        Group                               ( spiritGrammar.Group ),
        Variable                            ( spiritGrammar.Variable )
{
}

// ===================================================
// Operators
// ===================================================
template < typename IteratorType, typename ResultsType >
ParserSpiritGrammar< IteratorType, ResultsType >&
ParserSpiritGrammar< IteratorType, ResultsType >::operator=( const ParserSpiritGrammar& spiritGrammar )
{
    if ( this != &spiritGrammar )
    {
        ParserSpiritGrammar::base_type::operator=( spiritGrammar.Start );
        //ParserSpiritGrammar< IteratorType, ResultsType >::operator=( spiritGrammar );
        Start                               = spiritGrammar.Start;
        Assignment                          = spiritGrammar.Assignment;
//        Command                             = spiritGrammar.Command;
        Expression                          = spiritGrammar.Expression;
        Compare                             = spiritGrammar.Compare;
        PlusMinus                           = spiritGrammar.PlusMinus;
        MultiplyDivide                      = spiritGrammar.MultiplyDivide;
        Elevate                             = spiritGrammar.Elevate;
        Element                             = spiritGrammar.Element;
        Number                              = spiritGrammar.Number;
        Function                            = spiritGrammar.Function;
        Group                               = spiritGrammar.Group;
        Variable                            = spiritGrammar.Variable;
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
    Variable.add( "pi" , 3.141592653589793 );
    Variable.add( "e", 2.718281828459046 );
    //add( "pi", Pi ); //Using the tab of LifeV
    //add( "pi", 3.1415926535897932384626433832795 ); //Better only with long Real!
}

template < typename IteratorType, typename ResultsType >
inline void
ParserSpiritGrammar< IteratorType, ResultsType >::setVariable( const std::string& name, const Real& value )
{
    Real *p = Variable.find( name );
    if ( p != 0 )
        *p = value;
    else
        Variable.add( name, value );
}

//qi::rule< IteratorType, double(double), ascii::space_type > Power, MultiplyDivide, PlusMinus, Compare;
/*
        Expression
            =   ExpressionLevel_1                          [qi::_val = qi::_1]
            >> *PlusMinus(qi::_val)                        [qi::_val = qi::_1]
        ;

        ExpressionLevel_1
            =   Element                                    [qi::_val = qi::_1]
            >> *MultiplyDivide(qi::_val)                   [qi::_val = qi::_1]
        ;

        PlusMinus
            =   qi::eps(qi::_val = qi::_r1)
            >>  (
                qi::lit('+') >> ExpressionLevel_1          [qi::_val += qi::_1]
            |   qi::lit('-') >> ExpressionLevel_1          [qi::_val -= qi::_1]
                )
        ;

        MultiplyDivide
            =   qi::eps(qi::_val = qi::_r1)
            >>  (
                qi::lit('*') >> Element                    [qi::_val *= qi::_1]
            |   qi::lit('/') >> Element                    [qi::_val /= qi::_1]
                )
        ;
*/

#else

//! ParserSpiritGrammar - An empty implementation for boost version < 1.41
template < typename IteratorType = std::string::const_iterator, typename ResultsType = std::vector < Real > >
class ParserSpiritGrammar
{
public:

    typedef IteratorType                                        iterator_Type;
    typedef boost::iterator_range< iterator_Type >              iteratorRange_Type;
    typedef ResultsType                                         results_Type;

    ParserSpiritGrammar() : M_real(0.) {}
    ParserSpiritGrammar( const ParserSpiritGrammar& ) : M_real(0.) {}
    ~ParserSpiritGrammar() {}

    ParserSpiritGrammar& operator=( const ParserSpiritGrammar& ) { return *this; }

    void clearVariables() {}

    void setDefaultVariables() {}
    void setVariable( const std::string&, const Real& ) {}

    Real& variable( const std::string& ) { return M_real; }

    Real M_real;
};

#endif /* HAVE_BOOST_SPIRIT_QI */

} // Namespace LifeV

#endif /* Parser_SpiritGrammar_H */
