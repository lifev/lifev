/*
This file is part of the LifeV library.

Author: Christophe Prud'homme (christophe.prudhomme@epfl.ch)

Copyright (C) 2004 EPFL

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <cstring>
#include <cstdlib>
#include <errno.h>

#include <list>
#include <map>
#include <vector>

#include <iterator>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>


#if defined(HAVE_CONFIG_H)
# include <lifeconfig.h>
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_BACKTRACE
# include <execinfo.h>
#endif

#include <LDebug.hpp>

namespace LifeV
{
/*!
  \class LDebug
  \brief Area debugging tool
 
  \c LDebug() provides a debug stream to which you can pass a number, say 100, associated
  to an area of the code, say a class \c A.  In the implementation of the class \c A, you use
  debug statement like
 
  void A::f()
  {
    LDebug(100) << "A::f() is called.\n";
    // do something here
  }
 
  Now the debug message "A::f() is called." will be seen only if the area 100
  is defined in the environment(shell) variable \c DEBUG while executing a program
  \c A::f() is called \c runWithA that makes use of our class \c A.
 
  > runwithA
    --> no debug message related to A
  > export DEBUG="100"
  > runwithA
    A::f() is called.
 
   With this tool you can select the area you want to debug explicitly while keeping the
   others hidden.
 
  @author Christophe Prud'homme (christophe.prudhomme@epfl.ch)
*/
enum DebugLevels
{
    LDEBUG_INFO = 0,
    LDEBUG_WARN = 1,
    LDEBUG_ERROR = 2,
    LDEBUG_FATAL = 3
};
struct LDebugStream::Private
{
    Private()
            :
            debug( false ),
            __flush_function( 0 )
    {}
    bool debug;
    std::ostringstream _M_output;

    stprintf __flush_function;

    static bool _S_attached;
    static std::ofstream _S_logfile;
    //static std::map<int, std::ostream*> _S_logfile_per_area;

};
//
// getDescription
//
static std::map<uint, std::string>* DebugAreas = 0;
static std::string* StringNull = 0;
static std::list<int>* AREAS;
static std::string* DEBUG_AREA = 0;

namespace
{
// this Function makes sure that the static variables are initialized
// properly before being used
void
initDebugAreas ()
{
    static bool alloc = false;
    if ( alloc == false )
    {
        DEBUG_AREA = new std::string ( "" );
        AREAS = new std::list<int>;
        StringNull = new std::string ( "" );
        DebugAreas = new std::map<uint, std::string>;
        alloc = true;


        // read entries in sdebug.areas
        std::ostringstream __path;
        __path << LIFE_PREFIX << "/share/" << PACKAGE << "/ldebug.areas";

        std::ifstream fin ( __path.str().c_str() );
        if ( fin.fail () )
        {
            fin.open ( "../ldebug.areas" );
            if ( fin.fail() )
            {
                LDWarning() << "Debug areas not found.\n";
            }
        }
        while ( fin )
        {
            char __line[ 256 ];
            fin.getline ( __line, 256 );
            if ( __line[ 0 ] == '\0' ||
                    __line[ 0 ] == '\n' ||
                    __line[ 0 ] == '#' ||
                    isspace ( __line[ 0 ] ) )
                continue;
            std::istringstream __sentry ( __line );
            std::vector<std::string> __l;
            std::copy ( std::istream_iterator<std::string, char> ( __sentry ),
                        std::istream_iterator<std::string, char> (),
                        std::back_inserter ( __l ) );
            DebugAreas->insert ( std::make_pair ( std::atoi ( __l[ 0 ].c_str() ), __l[ 2 ] ) );
        }

        char * __env = getenv( "DEBUG" );
        if ( __env )
        {
            *DEBUG_AREA = __env;
        }
        std::istringstream __is ( *DEBUG_AREA );


        std::copy ( std::istream_iterator<int, char> ( __is ),
                    std::istream_iterator<int, char> (),
                    std::back_inserter ( *AREAS ) );

        DebugAreas->insert ( std::make_pair ( 0, "" ) );

    }
}
std::string const&
getDescription ( uint __area )
{
    if ( DebugAreas->empty() )
        return * StringNull;

    std::map<uint, std::string>::iterator entry_it = DebugAreas->find ( __area );

    if ( entry_it != DebugAreas->end() )
        return entry_it->second;
    else
        return DebugAreas->find ( 0 ) ->second;


}
}

//
// LDebugStream
//
LDebugStream::LDebugStream( int area, int level, bool print )
        :
        __p( new Private )
{
    initDebugAreas ();

    if ( DEBUG_AREA && ! DEBUG_AREA->empty() )
    {
        __p->debug = ( std::find ( AREAS->begin (), AREAS->end (), area ) != AREAS->end() &&
                       print ||
                       !area );
    }
    else
    {
        __p->debug = ( print && !area );
    }
    if ( __p->debug && getDescription( area ).empty() == false )
        __p->_M_output << getDescription ( area ) << ": ";

}
LDebugStream::LDebugStream( const char* initialString, int area, int level, bool print )
        :
        __p( new Private )
{
    initDebugAreas ();
    if ( DEBUG_AREA && ! DEBUG_AREA->empty() )
    {
        __p->debug = ( std::find ( AREAS->begin (), AREAS->end (), area ) != AREAS->end() &&
                       print ||
                       !area );
    }
    else
    {
        __p->debug = ( print && !area );
    }
    if ( __p->debug && getDescription( area ).empty() == false )
        __p->_M_output << getDescription ( area ) << ": "
        << initialString;
}
LDebugStream::LDebugStream( const LDebugStream& sd )
        :
        __p( new Private )
{
    __p->debug = sd.__p->debug;
    __p->__flush_function = sd.__p->__flush_function;
}
LDebugStream::~LDebugStream()
{
    delete __p;
}
LDebugStream&
LDebugStream::operator<<( char const* s )
{
    if ( __p->debug )
        __p->_M_output << s;
    if ( s[ strlen( s ) - 1 ] == '\n' )
        flush();
    return *this;
}
LDebugStream&
LDebugStream::operator<<( double s )
{
    if ( __p->debug )
    {
        __p->_M_output << s;
        flush();
    }
    return *this;
}
#if 0

LDebugStream&
LDebugStream::operator<<( ulong s )
{
    if ( __p->debug )
    {
        __p->_M_output << s;
        flush();
    }
    return *this;
}

LDebugStream&
LDebugStream::operator<<( long s )
{
    if ( __p->debug )
    {
        __p->_M_output << s;
        flush();
    }
    return *this;
}
LDebugStream&
LDebugStream::operator<<( int s )
{
    if ( __p->debug )
    {
        __p->_M_output << s;
        flush();
    }
    return *this;
}
#endif
LDebugStream&
LDebugStream::operator<<( std::string const& s )
{
    if ( __p->debug )
        __p->_M_output << s;
    if ( s[ s.size() - 1 ] == '\n' )
        flush();
    return *this;
}


LDebugStream&
LDebugStream::operator<<( LifeV::LManipFunction __f )
{
    if ( __p->debug )
    {
        ( *__f ) ( *this );
    }
    return *this;
}

void
LDebugStream::setFlush( stprintf func )
{
    __p->__flush_function = func;
}
void
LDebugStream::flush( )
{
    if ( !__p->_M_output.str().empty() )
    {
        __p->_M_output << std::ends;
        if ( Private::_S_attached )
        {
            Private::_S_logfile << __p->_M_output.str();
        }
        else if ( __p->__flush_function == 0 )
        {
            std::cerr << __p->_M_output.str();
        }
        else
        {
            __p->__flush_function( "%s", __p->_M_output.str().c_str() );
        }
        __p->_M_output.str( "" );
    }

}

bool LDebugStream::Private::_S_attached = false;
std::ofstream LDebugStream::Private::_S_logfile;
//std::map<int, std::ostream*> LDebugStream::_S_logfile_per_area;

void LDebugStream::attach( std::string const& __logfile )
{
    std::ostringstream __filename;
    __filename << __logfile << std::ends;

    if ( Private::_S_logfile.is_open() )
    {
        Private::_S_logfile.close();
    }

    Private::_S_logfile.open( __filename.str().c_str(), std::ios::out );

    if ( Private::_S_logfile.fail() )
    {
        LDWarning() << "LDebugStream::attach( " << __logfile.c_str() << " ) failed to open " << __filename.str() << "\n";
        LDWarning() << "Redirecting to default output\n";
        Private::_S_attached = false;
    }
    else if ( Private::_S_logfile.is_open() )
    {
        Private::_S_logfile << __filename.str() << " is opened for debug" << std::endl;
        Private::_S_attached = true;
    }
}
void
LDebugStream::attach( std::string const& __logfile, int __area )
{}
void
LDebugStream::detach( std::string const& __logfile, int __area )
{}

void
LDebugStream::detachAll()
{}



#ifndef NDEBUG_OLD
LDebugStream
LDebug( int area, LDebugStream::stprintf func )
{
    LDebugStream s( area, LDEBUG_INFO );
    s.setFlush( func );
    return s;
}

LDebugStream
LDebug( bool cond, int area, LDebugStream::stprintf func )
{
    if ( cond )
        return LDebugStream( area, LDEBUG_INFO );
    else
        return LDebugStream( 0, 0, false );
}
#endif

LDebugStream
LDWarning( int area )
{
    return LDebugStream( "WARNING: ", area, LDEBUG_WARN );
}

LDebugStream
LDWarning( bool cond, int area )
{
    if ( cond )
        return LDebugStream( "WARNING: ", area, LDEBUG_WARN );
    else
        return LDebugStream( 0, 0, false );

}

LDebugStream
LDError( int area )
{
    //LDebug () << LBacktrace() << "\n";
    return LDebugStream( "ERROR: ", area, LDEBUG_ERROR );
}

LDebugStream
LDError( bool cond, int area )
{
    //LDebug () << LBacktrace() << "\n";
    if ( cond )
        return LDebugStream( "ERROR: ", area, LDEBUG_ERROR );
    else
        return LDebugStream( 0, 0, false );

}

LDebugStream
LDFatal( int area )
{
    //LBacktrace();
    return LDebugStream( "FATAL: ", area, LDEBUG_FATAL );
}

LDebugStream
LDFatal( bool cond, int area )
{
    //LBacktrace();
    if ( cond )
        return LDebugStream( "FATAL: ", area, LDEBUG_FATAL );
    else
        return LDebugStream( 0, 0, false );
}

std::string
LBacktrace ()
{
    // show all backtrace
    return LBacktrace( -1 );
}
std::string
LBacktrace ( int __levels )
{
    std::ostringstream os;
#ifdef HAVE_BACKTRACE

    void* trace[ 256 ];
    int n = backtrace ( trace, 256 );
    char** strings = backtrace_symbols ( trace, n );

    if ( __levels != -1 )
        n = std::min ( n, __levels );
    os << "[\n";

    for ( int i = 0; i < n; ++i )
        os << i << ": " << strings[ i ] << "\n";
    os << "]\n";
    free ( strings );
#endif

    return os.str();
}
}

LifeV::LDebugStream&
perror( LifeV::LDebugStream& s )
{
    s << " " << strerror( errno );
    return s;
}
LifeV::LDebugStream&
endl( LifeV::LDebugStream& s )
{
    s << "\n";
    return s;
}
LifeV::LDebugStream&
flush( LifeV::LDebugStream& s )
{
    s.flush();
    return s;
}

