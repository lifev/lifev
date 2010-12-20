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

#include <cstdlib>

#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

#include <life/lifecore/SmartAssert.hpp>

// TODO: TO BE REMOVED: ARGUABLY USELESS AND BREAKS COMPILATION ON SOME
// SYSTEMS (BG/P etc.)
void breakIntoDebugger()
{
    // MSVC, BCB,
#if (defined _MSC_VER) || (defined __BORLANDC__)
    __asm { int 3 };
#elif defined(__GNUC__)
    // GCC
    __asm ("int $0x3");
#else
#  error Please supply instruction to break into code
#endif
}

namespace LifeV
{
namespace
{
// in case we're logging using the default logger...
struct StreamHolder
{
    StreamHolder() : M_outputStream( 0), M_ownership( false) {}
    ~StreamHolder()
    {
        if ( M_ownership)
            delete M_outputStream;
        M_outputStream = 0;
    }
    std::ostream * M_outputStream;
    bool M_ownership;
};
// information about the stream we write to, in case
// we're using the default logger
StreamHolder defaultLoggerInfo;

// intitializes the SMART_ASSERT library
struct assertInitializer
{
    assertInitializer()
    {
        Private::initAssert();
    }
}
init;
} // anonymous namespace

namespace SmartAssert
{

// returns a message corresponding to the type of level
std::string getTypeofLevel( int nLevel)
{
    switch ( nLevel)
    {
    case lvl_warn:
        return "Warning";
    case lvl_debug:
        return "Assertion failed";
    case lvl_error:
        return "Assertion failed (Error)";
    case lvl_fatal:
        return "Assertion failed (FATAL)";
    default:
    {
        std::ostringstream out;
        out << "Assertion failed (level=" << nLevel << ")";
        return out.str();
    }
    };
}

// helpers, for dumping the assertion context
void dumpContextSummary( const AssertContext & context, std::ostream & out)
{
    out << "\n" << getTypeofLevel( context.level() )
    << " in " << context.file() << ":" << context.contextLine() << '\n';
    if ( !context.message().empty())
        // we have a user-friendly message
        out << context.message();
    else
        out << "\nExpression: " << context.expression();
    out << std::endl;
}

void dumpContextDetail( const AssertContext & context, std::ostream & out)
{
    out << "\n" << getTypeofLevel( context.level() )
    << " in " << context.file() << ":" << context.contextLine() << '\n';
    if ( !context.message().empty())
        out << "User-friendly msg: '" << context.message() << "'\n";
    out << "\nExpression: '" << context.expression() << "'\n";

    typedef AssertContext::valueArray_Type valueArray_Type;
    const valueArray_Type & aVals = context.values();
    if ( !aVals.empty() )
    {
        bool bFirstTime = true;
        valueArray_Type::const_iterator first = aVals.begin(), last = aVals.end();
        while ( first != last)
        {
            if ( bFirstTime)
            {
                out << "Values: ";
                bFirstTime = false;
            }
            else
            {
                out << "        ";
            }
            out << first->second << "='" << first->first << "'\n";
            ++first;
        }
    }
    out << std::endl;
}

///////////////////////////////////////////////////////
// logger

void defaultLogger( const AssertContext & context)
{
    if ( defaultLoggerInfo.M_outputStream == 0)
        return;
    dumpContextDetail( context, *( defaultLoggerInfo.M_outputStream) );
}

///////////////////////////////////////////////////////
// handlers

// warn : just dump summary to console
void defaultWarnHandler( const AssertContext & context)
{
    dumpContextSummary( context, std::cout);
}


// debug: ask user what to do
void defaultDebugHandler( const AssertContext & context)
{
    static bool ignore_all = false;
    if ( ignore_all)
        // ignore All asserts
        return;
    typedef std::pair< std::string, int> fileAndLine;
    static std::set< fileAndLine> ignorer;
    if ( ignorer.find( fileAndLine( context.file(), context.contextLine())) != ignorer.end() )
        // this is Ignored Forever
        return;

    dumpContextSummary( context, std::cerr );
    std::cerr << "\nPress (I)gnore/ Igore (F)orever/ Ignore (A)ll/ (D)ebug/ A(b)ort: ";
    std::cerr.flush();
    char ch = 0;

    bool bContinue = true;
    while ( bContinue && std::cin.get( ch))
    {
        bContinue = false;
        switch ( ch)
        {
        case 'i':
        case 'I':
            // ignore
            break;

        case 'f':
        case 'F':
            // ignore forever
            ignorer.insert( fileAndLine( context.file(), context.contextLine()));
            break;

        case 'a':
        case 'A':
            // ignore all
            ignore_all = true;
            break;

        case 'd':
        case 'D':
            // break
            breakIntoDebugger();
            break;

        case 'b':
        case 'B':
            abort();
            break;

        default:
            bContinue = true;
            break;
        }
    }
}


// error : throw a runtime exception
void defaultErrorHandler( const AssertContext & context)
{
    std::ostringstream out;
    dumpContextSummary( context, out);
    throw std::runtime_error( out.str());
}


// fatal : dump error and abort
void defaultFatalHandler( const AssertContext & context)
{
    dumpContextDetail( context, std::cerr);
    abort();
}


} // namespace SmartAssert


namespace Private
{

void initAssert()
{
    ::LifeV::Assert::setLog( &::LifeV::SmartAssert::defaultLogger);
    ::LifeV::Assert::setHandler( lvl_warn, &::LifeV::SmartAssert::defaultWarnHandler);
    ::LifeV::Assert::setHandler( lvl_debug, &::LifeV::SmartAssert::defaultDebugHandler);
    ::LifeV::Assert::setHandler( lvl_error, &::LifeV::SmartAssert::defaultErrorHandler);
    ::LifeV::Assert::setHandler( lvl_fatal, &::LifeV::SmartAssert::defaultFatalHandler);
}

// sets the default logger to write to this stream
void setDefaultLogStream( std::ostream & out)
{
    defaultLoggerInfo.M_outputStream = &out;
    defaultLoggerInfo.M_ownership = false;
}

// sets the default logger to write to this file
void setDefaultLogName( const char * str)
{
    defaultLoggerInfo.M_ownership = false;
    defaultLoggerInfo.M_outputStream = new std::ofstream( str);
    defaultLoggerInfo.M_ownership = true;
}


} // namespace Private

}
