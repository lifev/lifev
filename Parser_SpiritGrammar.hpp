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
 * @brief Parser_SpiritGrammar
 *
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 * @date 05-02-2010
 */

#ifndef Parser_SpiritGrammar_H
#define Parser_SpiritGrammar_H 1

#include <lifemc/lifecore/Parser_Definitions.hpp>

namespace LifeV
{

#ifdef HAVE_BOOST_SPIRIT_QI

//! Parser_SpiritGrammar - A string parser grammar based on \c boost::spirit::qi
/*!
 *  @author(s) Cristiano Malossi
 *
 *  \c Parser_SpiritGrammar is a \c boost::spirit::qi based class to perform
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
 *  Currently Parser_SpiritGrammar works with the following operators:
 *  \verbatim
 *  +, -, *, /, ^, sqrt(), sin(), cos(), tan(), exp(), log(), log10(), >, <.
 *  \endverbatim
 *
 */
template <typename iterator>
class Parser_SpiritGrammar : public qi::grammar< iterator, Results_Type(), ascii::space_type >
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    Parser_SpiritGrammar();

    //! Copy constructor
    /*!
     * @param Parser_SpiritGrammar Parser_SpiritGrammar
     */
    Parser_SpiritGrammar( const Parser_SpiritGrammar& SpiritGrammar );

    //! Destructor
    ~Parser_SpiritGrammar() {}
    //@}


    //! @name Operators
    //@{

    //! Operator =
    /*!
     * @param SpiritGrammar Parser_SpiritGrammar
     * @return reference to a copy of the class
     */
    Parser_SpiritGrammar& operator=( const Parser_SpiritGrammar& SpiritGrammar );

    //@}


    //! @name Methods
    //@{

    /*! Assign a variable using a \c boost::iterator_range
     *
     * @param StringIteratorRange name of the parameter
     * @param Value value of the parameter
     */
    void AssignVariable( const boost::iterator_range< String_Iterator >& StringIteratorRange, const Real& Value );

    //! Clear all the variables.
    void ClearVariables();

    /*
        //! Show all the variables
        void ShowMe();
    */

    //@}


    //! @name Set Methods
    //@{

    //! Set default variables
    void SetDefaultVariables();

    /*! Set/replace a variable
     *
     * @param name name of the parameter
     * @param value value of the parameter
     */
    void SetVariable( const std::string& Name, const Real& Value );

    //@}

    //! @name Get Methods
    //@{

    /*! Get variable
     *
     * @param name name of the parameter
     * @return value of the variable
     */
    Real& GetVariable( const std::string& Name );

    //@}

private:

    //! @name Private Methods
    //@

    /*! Phoenix wrapper for \c std::sin
     *
     * @param Value input value
     * @return sin( Value )
     */
    Real sin( const Real& Value ) const;

    /*! Phoenix wrapper for \c std::cos
     *
     * @param Value input value
     * @return cos( Value )
     */
    Real cos( const Real& Value ) const;

    /*! Phoenix wrapper for \c std::tan
     *
     * @param Value input value
     * @return tan( Value )
     */
    Real tan( const Real& Value ) const;

    /*! Phoenix wrapper for \c std::pow
     *
     * @param base input base
     * @param exponent input exponent
     * @return pow( base, exponent )
     */
    Real pow( const Real& Base, const Real& Exponent ) const;

    /*! Phoenix wrapper for \c std::sqrt
     *
     * @param Value input value
     * @return sqrt( Value )
     */
    Real sqrt( const Real& Value ) const;

    /*! Phoenix wrapper for \c std::exp
     *
     * @param Value input value
     * @return exp( Value )
     */
    Real exp( const Real& Value ) const;

    /*! Phoenix wrapper for \c std::log
     *
     * @param Value input value
     * @return log( Value )
     */
    Real log( const Real& Value ) const;

    /*! Phoenix wrapper for \c std::log10
     *
     * @param Value input value
     * @return log10( Value )
     */
    Real log10( const Real& Value ) const;

    //@}

    qi::rule< iterator, Results_Type(), ascii::space_type > Start;

    qi::rule< iterator, void(), ascii::space_type >         Assignment;
//    qi::rule< iterator, void(), ascii::space_type >         Command;

    qi::rule< iterator, double(), ascii::space_type >       Expression;
    qi::rule< iterator, double(), ascii::space_type >       Compare;
    qi::rule< iterator, double(), ascii::space_type >       PlusMinus;
    qi::rule< iterator, double(), ascii::space_type >       MultiplyDivide;
    qi::rule< iterator, double(), ascii::space_type >       Elevate;
    qi::rule< iterator, double(), ascii::space_type >       Element;
    qi::rule< iterator, double(), ascii::space_type >       Number;
    qi::rule< iterator, double(), ascii::space_type >       Function;
    qi::rule< iterator, double(), ascii::space_type >       Group;

    qi::symbols<char, Real >                                Variable;
};



// ===================================================
// Constructors & Destructor
// ===================================================
template <typename iterator>
Parser_SpiritGrammar< iterator >::Parser_SpiritGrammar() :
        Parser_SpiritGrammar::base_type( Start ),
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
        )              [phoenix::bind(&Parser_SpiritGrammar::AssignVariable,this, qi::_1, qi::_2)]
        ;
    /*
        Command
            =  qi::lit("ShowMe")[phoenix::bind(&Parser_SpiritGrammar::ShowMe, this)]
        ;
    */

    Expression
    =  *Compare                                    [qi::_val = qi::_1]
       ;

    Compare
    =   PlusMinus                                  [qi::_val = qi::_1]
        >> *(
            qi::lit('>') >> PlusMinus                  [qi::_val = qi::_val > qi::_1]
            |   qi::lit('<') >> PlusMinus                  [qi::_val = qi::_val < qi::_1]
        )
        ;

    PlusMinus
    =   MultiplyDivide                             [qi::_val = qi::_1]
        >> *(
            qi::lit('+') >> MultiplyDivide             [qi::_val += qi::_1]
            |   qi::lit('-') >> MultiplyDivide             [qi::_val -= qi::_1]
        )
        ;

    MultiplyDivide
    =   Elevate                                    [qi::_val = qi::_1]
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
                qi::lit('^') >> Element                    [qi::_val = -phoenix::bind(&Parser_SpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
            >> *(
                qi::lit('^') >> Element                    [qi::_val = phoenix::bind(&Parser_SpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
        )
        |
        (
            Element                                    [qi::_val = qi::_1]
            >> *(
                qi::lit('^') >> Element                    [qi::_val = phoenix::bind(&Parser_SpiritGrammar::pow, this, qi::_val, qi::_1)]
            )
        )
        ;

    Element
    =
        (
            qi::lit('-') >> Element                        [qi::_val = -qi::_1]
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
            qi::lit("sin")   >> Group                  [qi::_val = phoenix::bind(&Parser_SpiritGrammar::sin, this,   qi::_1)]
            |   qi::lit("cos")   >> Group                  [qi::_val = phoenix::bind(&Parser_SpiritGrammar::cos, this,   qi::_1)]
            |   qi::lit("tan")   >> Group                  [qi::_val = phoenix::bind(&Parser_SpiritGrammar::tan, this,   qi::_1)]
            |   qi::lit("sqrt")  >> Group                  [qi::_val = phoenix::bind(&Parser_SpiritGrammar::sqrt, this,  qi::_1)]
            |   qi::lit("exp")   >> Group                  [qi::_val = phoenix::bind(&Parser_SpiritGrammar::exp, this,   qi::_1)]
            |   qi::lit("log")   >> Group                  [qi::_val = phoenix::bind(&Parser_SpiritGrammar::log, this,   qi::_1)]
            |   qi::lit("log10") >> Group                  [qi::_val = phoenix::bind(&Parser_SpiritGrammar::log10, this, qi::_1)]
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

template <typename iterator>
Parser_SpiritGrammar< iterator >::Parser_SpiritGrammar( const Parser_SpiritGrammar& SpiritGrammar ) :
        Parser_SpiritGrammar::base_type     ( SpiritGrammar.Start ),
        Start                               ( SpiritGrammar.Start ),
        Assignment                          ( SpiritGrammar.Assignment ),
//    Command                             ( SpiritGrammar.Command ),
        Expression                          ( SpiritGrammar.Expression ),
        Compare                             ( SpiritGrammar.Compare ),
        PlusMinus                           ( SpiritGrammar.PlusMinus ),
        MultiplyDivide                      ( SpiritGrammar.MultiplyDivide ),
        Elevate                             ( SpiritGrammar.Elevate ),
        Element                             ( SpiritGrammar.Element ),
        Number                              ( SpiritGrammar.Number ),
        Function                            ( SpiritGrammar.Function ),
        Group                               ( SpiritGrammar.Group ),
        Variable                            ( SpiritGrammar.Variable )
{
}

// ===================================================
// Operators
// ===================================================
template <typename iterator>
Parser_SpiritGrammar< iterator >&
Parser_SpiritGrammar< iterator >::operator=( const Parser_SpiritGrammar& SpiritGrammar )
{
    if ( this != &SpiritGrammar )
    {
        Parser_SpiritGrammar::base_type::operator=( SpiritGrammar.Start );
        //Parser_SpiritGrammar< iterator >::operator=( SpiritGrammar );
        Start                               = SpiritGrammar.Start;
        Assignment                          = SpiritGrammar.Assignment;
//        Command                             = SpiritGrammar.Command;
        Expression                          = SpiritGrammar.Expression;
        Compare                             = SpiritGrammar.Compare;
        PlusMinus                           = SpiritGrammar.PlusMinus;
        MultiplyDivide                      = SpiritGrammar.MultiplyDivide;
        Elevate                             = SpiritGrammar.Elevate;
        Element                             = SpiritGrammar.Element;
        Number                              = SpiritGrammar.Number;
        Function                            = SpiritGrammar.Function;
        Group                               = SpiritGrammar.Group;
        Variable                            = SpiritGrammar.Variable;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
template <typename iterator>
void
Parser_SpiritGrammar< iterator >::AssignVariable( const boost::iterator_range< String_Iterator >& StringIteratorRange, const Real& Value )
{
    SetVariable( std::string( StringIteratorRange.begin(), StringIteratorRange.end() ), Value );
}

template <typename iterator>
void
Parser_SpiritGrammar< iterator >::ClearVariables()
{
    Variable.clear();
}

/*
template <typename iterator>
void
Parser_SpiritGrammar< iterator >::ShowMe()
{
    std::cout << "ShowMe called!" << std::endl;
}
*/

// ===================================================
// Set Methods
// ===================================================
template <typename iterator>
void
Parser_SpiritGrammar< iterator >::SetDefaultVariables()
{
    Variable.add( "pi" , 3.141592653589793 );
    Variable.add( "e", 2.718281828459046 );
    //add( "pi", Pi ); //Using the tab of LifeV
    //add( "pi", 3.1415926535897932384626433832795 ); //Better only with long Real!
}

template <typename iterator>
void
Parser_SpiritGrammar< iterator >::SetVariable( const std::string& Name, const Real& Value )
{
    Real *p = Variable.find( Name );
    if ( p != 0 )
        *p = Value;
    else
        Variable.add( Name, Value );
}

// ===================================================
// Get Methods
// ===================================================
template <typename iterator>
Real&
Parser_SpiritGrammar< iterator >::GetVariable( const std::string& Name )
{
    return Variable.at( Name );
}

// ===================================================
// Private Methods
// ===================================================
template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::sin( const Real& Value ) const
{
    return std::sin( Value );
}

template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::cos( const Real& Value ) const
{
    return std::cos( Value );
}

template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::tan( const Real& Value ) const
{
    return std::tan( Value );
}

template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::pow( const Real& Base, const Real& Exponent ) const
{
    return std::pow( Base, Exponent );
}

template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::sqrt( const Real& Value ) const
{
    return std::sqrt( Value );
}

template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::exp( const Real& Value ) const
{
    return std::exp( Value );
}

template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::log( const Real& Value ) const
{
    return std::log( Value );
}

template <typename iterator>
Real
Parser_SpiritGrammar< iterator >::log10( const Real& Value ) const
{
    return std::log10( Value );
}

//qi::rule< iterator, double(double), ascii::space_type > Power, MultiplyDivide, PlusMinus, Compare;
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

//! Parser_SpiritGrammar - An empty implementation for boost version < 1.41
template <typename iterator>
class Parser_SpiritGrammar
{
public:

    Parser_SpiritGrammar() : M_real(0.) {}
    Parser_SpiritGrammar( const Parser_SpiritGrammar& ) : M_real(0.) {}
    ~Parser_SpiritGrammar() {}

    Parser_SpiritGrammar& operator=( const Parser_SpiritGrammar& ) { return *this; }

    void ClearVariables() {}

    void SetDefaultVariables() {}
    void SetVariable( const std::string&, const Real& ) {}

    Real& GetVariable( const std::string& ) { return M_real; }

    Real M_real;
};

#endif /* HAVE_BOOST_SPIRIT_QI */

} // Namespace LifeV

#endif /* Parser_SpiritGrammar_H */
