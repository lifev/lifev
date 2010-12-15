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
  @file
  @brief Complex assertion mechanism

  @date 16-01-2005
  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef SMART_ASSERT_H
#define SMART_ASSERT_H

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace LifeV
{
enum
{
    // default behavior - just loggs this assert
    // (a message is shown to the user to the console)
    lvl_warn = 100,

    // default behavior - asks the user what to do:
    // Ignore/ Retry/ etc.
    lvl_debug = 200,

    // default behavior - throws a SmartAssert_error
    lvl_error = 300,

    // default behavior - dumps all assert context to console,
    // and aborts
    lvl_fatal = 1000
};

/**
   \class AssertContext
   \brief contains details about a failed assertion
*/
class AssertContext
{
public:
    //! @name Public typedefs
    //@{
    typedef std::pair< std::string, std::string> val_and_str;
    typedef std::vector< val_and_str> vals_array;
    //@}

    //! @name Constructor
    //@{
    AssertContext() : M_level( lvl_debug) {}
    //@}

    //! @name Public methods
    //@{
    // adds one value and its corresponding string
    void add_val( const std::string & val, const std::string & str)
    {
        M_vals.push_back( val_and_str( val, str) );
    }
    //@}

    //! @name Set methods
    //@{
    void setFileLine( const char * file, int line)
    {
        M_file = file;
        M_line = line;
    }

    void setExpression( const std::string & str) { M_expression = str; }

    void setLevel( int nLevel) { M_level = nLevel; }

    void setLevelMsg( const char * strMsg)
    {
        if ( strMsg)
            M_msg = strMsg;
        else
            M_msg.erase();
    }
    //@}

    //! @name Get methods
    //@{
    const std::string & getContextFile() const __attribute__ ((deprecated)) { return M_file; }

    int getContextLine() const __attribute__ ((deprecated)) { return M_line; }

    const std::string & expression() const { return M_expression; }

    const vals_array & get_vals_array() const __attribute__ ((deprecated)) { return M_vals; }

    int get_level() const __attribute__ ((deprecated)) { return M_level; }

    const std::string & get_level_msg() const __attribute__ ((deprecated)) { return M_msg; }
    //@}

private:
    // where the assertion occured
    std::string M_file;
    int M_line;

    // expression and values
    std::string M_expression;
    vals_array M_vals;

    // level and message
    int M_level;
    std::string M_msg;
};


namespace SmartAssert
{

typedef void (*assert_function_type)( const AssertContext & context);

// helpers
std::std::string getTypeofLevel( int nLevel);
void dumpContextSummary( const AssertContext & context, std::ostream & out);
void dumpContextDetail( const AssertContext & context, std::ostream & out);

// defaults
void defaultWarnHandler( const AssertContext & context);
void defaultDebugHandler( const AssertContext & context);
void defaultErrorHandler( const AssertContext & context);
void defaultFatalHandler( const AssertContext & context);
void defaultLogger( const AssertContext & context);

} // namespace SmartAssert

namespace Private
{
void initAssert();
void setDefaultLogStream( std::ostream & out);
void setDefaultLogName( const char * str);

// allows finding if a value is of type 'const char *'
// and is null; if so, we cannot print it to an ostream
// directly!!!
template< class T>
struct isNullFinder
{
    bool is( const T &) const
    {
        return false;
    }
};

template<>
struct isNullFinder< char*>
{
    bool is( char * const & val)
    {
        return val == 0;
    }
};

template<>
struct isNullFinder< const char*>
{
    bool is( const char * const & val)
    {
        return val == 0;
    }
};

} // namespace Private

struct Assert
{
    //! @name Helpers and Typedefs
    //@{
    typedef SmartAssert::assert_function_type assert_function_type;

    // helpers, in order to be able to compile the code
    Assert & SMART_ASSERT_A;
    Assert & SMART_ASSERT_B;
    //@}

    //! @name Constructors and destructor
    //@{
    Assert( const char * expr):
        SMART_ASSERT_A( *this),
        SMART_ASSERT_B( *this),
        M_needs_handling( true)
    {
        M_context.setExpression( expr);

        if ( ( logger() == 0) || handlers().size() < 4)
        {
            // used before main!
            Private::initAssert();
        }
    }

    Assert( const Assert & other):
        SMART_ASSERT_A( *this),
        SMART_ASSERT_B( *this),
        M_context( other.M_context),
        M_needs_handling( true)
    {
        other.M_needs_handling = false;
    }

    ~Assert()
    {
        if ( M_needs_handling)
            handleAssert();
    }
    //@}

    //! @name Public methods
    //@{
    template< class type>
    Assert & printCurrentValue( const type & val, const char * msg)
    {
        std::ostringstream out;

        Private::isNullFinder< type> f;
        bool bIsNull = f.is( val);
        if ( !bIsNull)
            out << val;
        else
            // null string
            out << "null";
        M_context.add_val( out.str(), msg);
        return *this;
    }

    Assert & printContext( const char * file, int line)
    {
        M_context.setFileLine( file, line);
        return *this;
    }

    Assert & msg( const char * strMsg)
    {
        M_context.setLevelMsg( strMsg);
        return *this;
    }

    Assert & level( int nLevel, const char * strMsg = 0)
    {
        M_context.setLevel( nLevel);
        M_context.setLevelMsg( strMsg);
        return *this;
    }

    Assert & warn( const char * strMsg = 0)
    {
        return level( lvl_warn, strMsg);
    }

    Assert & debug( const char * strMsg = 0)
    {
        return level( lvl_debug, strMsg);
    }

    Assert & error( const char * strMsg = 0)
    {
        return level( lvl_error, strMsg);
    }

    Assert & error( const std::string strMsg )
    {
        return level( lvl_error, strMsg.c_str());
    }

    Assert & fatal( const char * strMsg = 0)
    {
        return  level( lvl_fatal, strMsg);
    }
    //@}

    //! @name Static methods
    //@{
    // in this case, we set the default logger, and make it
    // write everything to this file
    static void setLog( const char * strFileName)
    {
        Private::setDefaultLogName( strFileName);
        logger() = &SmartAssert::defaultLogger;
    }

    // in this case, we set the default logger, and make it
    // write everything to this log
    static void setLog( std::ostream & out)
    {
        Private::setDefaultLogStream( out);
        logger() = &SmartAssert::defaultLogger;
    }

    static void setLog( assert_function_type log)
    {
        logger() = log;
    }

    static void setHandler( int nLevel, assert_function_type handler)
    {
        handlers()[ nLevel] = handler;
    }
    //@}

private:
    //! @name Private methods
    //@{
    // handles the current assertion.
    void handleAssert()
    {
        logger()( M_context);
        get_handler( M_context.get_level() )( M_context);
    }

    /*
      IMPORTANT NOTE:
      The only reason logger & handlers are functions, are
      because you might use SMART_ASSERT before main().

      In this case, since they're statics, they might not
      be initialized. However, making them functions
      will make it work.
    */
    //@}

    //! @name Private static methods
    //@{
    // the log
    static assert_function_type & logger()
    {
        static assert_function_type inst;
        return inst;
    }

    // the handler
    typedef std::map< int, assert_function_type> handlers_collection;
    static handlers_collection & handlers()
    {
        static handlers_collection inst;
        return inst;
    }

    static assert_function_type get_handler( int nLevel)
    {
        handlers_collection::const_iterator found = handlers().find( nLevel);
        if ( found != handlers().end() )
            return found->second;
        else
            // we always assume the debug handler has been set
            return handlers().find( lvl_debug)->second;
    }
    //@}

private:
    AssertContext M_context;
    mutable bool M_needs_handling;
};

namespace SmartAssert
{
inline ::LifeV::Assert make_assert( const char * expr)
{
    return ::LifeV::Assert( expr);
}
} // namespace SmartAssert

} // LifeV Namespace

#ifdef LIFEV_SMART_ASSERT_DEBUG_MODE

#if LIFEV_SMART_ASSERT_DEBUG_MODE == 1
#define LIFEV_SMART_ASSERT_DEBUG
#else
#undef LIFEV_SMART_ASSERT_DEBUG
#endif

#else

// defaults
#ifndef NDEBUG
#define LIFEV_SMART_ASSERT_DEBUG
#else
#undef LIFEV_SMART_ASSERT_DEBUG
#endif

#endif


#ifdef LIFEV_SMART_ASSERT_DEBUG
// "debug" mode
#define LIFEV_SMART_ASSERT( expr) \
    if ( (expr) ) ; \
    else ::LifeV::SmartAssert::make_assert( #expr).printContext( __FILE__, __LINE__).SMART_ASSERT_A \
    /**/

#else
// "release" mode
#define LIFEV_SMART_ASSERT( expr) \
    if ( true ) ; \
    else ::LifeV::SmartAssert::make_assert( "").SMART_ASSERT_A \
    /**/

#endif // ifdef LIFEV_SMART_ASSERT_DEBUG

// LIFEV_ASSERT is a equivalent to LIFEV_SMART_ASSERT
#define LIFEV_ASSERT( expr) LIFEV_SMART_ASSERT(expr)


#define LIFEV_SMART_VERIFY( expr) \
    if ( (expr) ) ; \
    else ::LifeV::SmartAssert::make_assert( #expr).error().printContext( __FILE__, __LINE__).SMART_ASSERT_A \
    /**/
#define LIFEV_VERIFY( expr) LIFEV_SMART_VERIFY(expr)


#define SMART_ASSERT_A(x) LIFEV_SMART_ASSERT_OP(x, B)
#define SMART_ASSERT_B(x) LIFEV_SMART_ASSERT_OP(x, A)

#define LIFEV_SMART_ASSERT_OP(x, next) \
    SMART_ASSERT_A.printCurrentValue((x), #x).SMART_ASSERT_ ## next \
    /**/

#endif // SMART_ASSERT_H
