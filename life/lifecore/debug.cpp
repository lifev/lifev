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

#include <boost/lexical_cast.hpp>

#include <lifeconfig.h>

#ifdef HAVE_BACKTRACE
# include <execinfo.h>
#endif

#include <life/lifecore/debug.hpp>

namespace LifeV
{
/*!
  \class Debug
  \brief Area debugging tool

  \c Debug() provides a debug stream to which you can pass a number, say 100, associated
  to an area of the code, say a class \c A.  In the implementation of the class \c A, you use
  debug statement like

  void A::f()
  {
    Debug(100) << "A::f() is called.\n";
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
    DEBUG_INFO  = 0,
    DEBUG_WARN  = 1,
    DEBUG_ERROR = 2,
    DEBUG_FATAL = 3,
};
struct DebugStream::Private
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
static std::map<unsigned int, std::string>* DebugAreas = 0;
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
        DebugAreas = new std::map<unsigned int, std::string>;
        alloc = true;


        // read entries in sdebug.areas
        std::ostringstream __path;
        __path << LIFE_PREFIX << "/share/" << PACKAGE << "/debug.areas";

        std::ifstream fin ( __path.str().c_str() );
        if ( fin.fail () )
        {
            fin.open ( "../debug.areas" );
            if ( fin.fail() )
            {
                Warning() << "The file debug.areas was not found.\n"
                << "                 searched at ../debug.areas and\n"
                << "                 " << __path.str() << "\n";
            }
        }
        while ( fin )
        {
            char __line[256];
            fin.getline ( __line, 256 );
            if ( __line[ 0 ] == '\0' ||
                    __line[ 0 ] == '\n' ||
                    __line[ 0 ] == '#' ||
                    isspace ( __line[ 0 ] ) )
                continue;
            std::istringstream __sentry ( __line );
            std::vector<std::string> __l;
            std::copy ( std::istream_iterator<std::string,char> ( __sentry ),
                        std::istream_iterator<std::string,char> (),
                        std::back_inserter ( __l ) );
            DebugAreas->insert ( std::make_pair ( std::atoi ( __l[0].c_str() ), __l[2] ) );
        }

        char * __env = getenv("DEBUG");
        if ( __env )
        {
            *DEBUG_AREA = __env;
        }
        std::istringstream __is ( *DEBUG_AREA );


        std::copy ( std::istream_iterator<int,char> ( __is ),
                    std::istream_iterator<int,char> (),
                    std::back_inserter ( *AREAS ) );

        DebugAreas->insert ( std::make_pair ( 0, "" ) );

    }
}
std::string
getDescription ( unsigned int __area )
{
    if ( DebugAreas->empty() )
        return std::string( "Area " ) + boost::lexical_cast<std::string>(__area);

    std::map<unsigned int, std::string>::iterator entry_it = DebugAreas->find ( __area );

    if ( entry_it != DebugAreas->end() )
        return entry_it->second;
    else
        return std::string( "Area " ) + boost::lexical_cast<std::string>(__area);


}
}

//
// DebugStream
//
DebugStream::DebugStream( int area, int /*level*/, bool print )
        :
        __p( new Private )
{
    initDebugAreas ();

    if ( DEBUG_AREA && ! DEBUG_AREA->empty() )
    {
        __p->debug =  ( std::find ( AREAS->begin (), AREAS->end (), area ) != AREAS->end() &&
                        print ) ||
                      ( !area );
    }
    else
    {
        __p->debug =  ( print && !area );
    }
    if ( __p->debug )
        __p->_M_output << getDescription ( area ) << ": ";

}
DebugStream::DebugStream( const char* initialString, int area, int /*level*/, bool print )
        :
        __p( new Private )
{
    initDebugAreas ();
    if ( DEBUG_AREA && ! DEBUG_AREA->empty() )
    {
        __p->debug =  ( std::find ( AREAS->begin (), AREAS->end (), area ) != AREAS->end() &&
                        print ) ||
                      ( !area );
    }
    else
    {
        __p->debug =  ( print && !area );
    }
    if ( __p->debug )
        __p->_M_output << getDescription ( area ) << ": "
        << initialString;
}
DebugStream::DebugStream( const DebugStream& sd )
        :
        __p( new Private )
{
    __p->debug = sd.__p->debug;
    __p->__flush_function = sd.__p->__flush_function;
}
DebugStream::~DebugStream()
{
    delete __p;
}
DebugStream&
DebugStream::operator<<( char const* s)
{
    if ( __p->debug )
        __p->_M_output  << s;
    if ( s[strlen(s)-1] == '\n')
        flush();
    return *this;
}
DebugStream&
DebugStream::operator<<( double s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
#if 0

DebugStream&
DebugStream::operator<<( ulong s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}

DebugStream&
DebugStream::operator<<( long s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
DebugStream&
DebugStream::operator<<( int s)
{
    if ( __p->debug )
    {
        __p->_M_output  << s;
        flush();
    }
    return *this;
}
#endif
DebugStream&
DebugStream::operator<<( std::string const& s)
{
    if ( __p->debug )
        __p->_M_output  << s;
    if ( s[s.size() -1] == '\n')
        flush();
    return *this;
}


DebugStream&
DebugStream::operator<<( LifeV::LManipFunction __f )
{
    if ( __p->debug )
    {
        (*__f)( *this );
    }
    return *this;
}

void
DebugStream::setFlush( stprintf func )
{
    __p->__flush_function = func;
}
void
DebugStream::flush(  )
{
    if ( !__p->_M_output.str().empty() )
    {
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

bool DebugStream::Private::_S_attached = false;
std::ofstream DebugStream::Private::_S_logfile;
//std::map<int, std::ostream*> DebugStream::_S_logfile_per_area;

void DebugStream::attach( std::string const& __logfile )
{
    std::ostringstream __filename;
    __filename <<  __logfile;

    if ( Private::_S_logfile.is_open() )
    {
        Private::_S_logfile.close();
    }

    Private::_S_logfile.open( __filename.str().c_str(), std::ios::out );

    if ( Private::_S_logfile.fail() )
    {
        Warning() << "DebugStream::attach( " << __logfile.c_str() << " ) failed to open "  << __filename.str() << "\n";
        Warning() << "Redirecting to default output\n";
        Private::_S_attached = false;
    }
    else if ( Private::_S_logfile.is_open() )
    {
        Private::_S_logfile << __filename.str() << " is opened for debug" << std::endl;
        Private::_S_attached = true;
    }
}
void
DebugStream::attach( std::string const& /*__logfile*/, int /*__area*/ )
{

}
void
DebugStream::detach( std::string const& /*__logfile*/, int /*__area*/ )
{}

void
DebugStream::detachAll()
{}



#ifndef NDEBUG_OLD
DebugStream
Debug( int area, DebugStream::stprintf func )
{
    DebugStream s( area, DEBUG_INFO );
    s.setFlush( func );
    return s;
}

DebugStream
Debug( bool cond, int area, DebugStream::stprintf /*func*/ )
{
    if ( cond )
        return DebugStream( area, DEBUG_INFO );
    else
        return DebugStream( 0, 0, false );
}
#endif

DebugStream
Warning( int area )
{
    return DebugStream( "WARNING: ", area, DEBUG_WARN );
}

DebugStream
Warning( bool cond, int area )
{
    if ( cond )
        return DebugStream( "WARNING: ", area, DEBUG_WARN );
    else
        return DebugStream( 0, 0, false );

}

DebugStream
Error( int area )
{
    //Debug () << LBacktrace() << "\n";
    return DebugStream( "ERROR: ", area, DEBUG_ERROR );
}

DebugStream
Error( bool cond, int area )
{
    //Debug () << LBacktrace() << "\n";
    if ( cond )
        return DebugStream( "ERROR: ", area, DEBUG_ERROR );
    else
        return DebugStream( 0, 0, false );

}

DebugStream
Fatal( int area )
{
    //LBacktrace();
    return DebugStream( "FATAL: ", area, DEBUG_FATAL );
}

DebugStream
Fatal( bool cond, int area )
{
    //LBacktrace();
    if ( cond )
        return DebugStream( "FATAL: ", area, DEBUG_FATAL );
    else
        return DebugStream( 0, 0, false );
}

#ifdef HAVE_BACKTRACE
std::string
backtrace ()
{
    // show all backtrace
    return backtrace( -1 );
}

std::string
backtrace ( int __levels )
{
    std::ostringstream os;

    void* trace[256];
    int n = backtrace ( trace, 256 );
    char** strings = backtrace_symbols ( trace, n );

    if ( __levels != -1 )
        n = ( std::min ) ( n, __levels );
    os << "[\n";

    for (int i = 0; i < n; ++i)
        os << i << ": " << strings[i] << "\n";
    os << "]\n";
    free (strings);

    return os.str();
}
#endif
}

LifeV::DebugStream&
perror( LifeV::DebugStream& s )
{
    s << " " << strerror( errno );
    return s;
}
LifeV::DebugStream&
endl( LifeV::DebugStream& s )
{
    s << "\n";
    return s;
}
LifeV::DebugStream&
flush( LifeV::DebugStream& s )
{
    s.flush();
    return s;
}

