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
#ifndef __Debug_H
#define __Debug_H 1

#include <cstdio>
#include <iosfwd>

#include <string>
#include <sstream>

namespace LifeV
{
class DebugStream;
class NdebugStream;

typedef DebugStream & (*LManipFunction)( DebugStream &); // manipulator function
typedef NdebugStream & (*LNManipFunction)( NdebugStream&); // manipulator function

#ifdef __GNUC__
# define LIFEV_FUNCINFO "[" << __PRETTY_FUNCTION__ << "] "
#else
# define LIFEV_FUNCINFO "[" << __FILE__ << ":" << __LINE__ << "] "
#endif

#define LIFEV_LINEINFO "[" << __FILE__ << ":" << __LINE__ << "] "

class DebugStream
{
public:


    /** @name Internal Structures
     */
    //@{

    struct Private;

    typedef int (*stprintf)( const char* format, ... );
    //@}

    /** @name Constructors, destructor
     */
    //@{

    DebugStream(int area = 0, int level = 1, bool print = true);
    DebugStream(const char* initialString, int area = 0, int level = 1, bool print = true);
    DebugStream( DebugStream const& );
    ~DebugStream();

    //@}

    void setFlush( stprintf = 0 );

    /** @name  Methods
     */
    //@{

    static void attach( std::string const& __logfile );
    static void attach( std::string const& __logfile, int area );
    static void detach( std::string const& __logfile, int area );
    static void detachAll();
    void flush();


    DebugStream& operator<<( char const* );

    DebugStream& operator<<( double );

    DebugStream& operator<<( std::string const& );
    DebugStream& operator<<( LManipFunction f );
    //@}



protected:

private:
    Private* __p;

};

template<typename T>
DebugStream& operator<< ( DebugStream& __s, T const* __t )
{
    std::ostringstream __os;
    __os << __t;
    __s << __os.str();
    return __s;
}
#ifdef HAVE_BACKTRACE
std::string backtrace ();
std::string backtrace ( int );
#endif

class NdebugStream
{
public:
    /** @name Constructors, destructor
     */
    //@{
    typedef int (*stprintf)( const char* format, ... );

    NdebugStream() {}
    ~NdebugStream() {}

    //@}

    /** @name  Methods
     */
    //@{
    void flush( stprintf = 0 ) {}
    NdebugStream& operator<<( char const* ) { return *this; }
    NdebugStream& operator<<( std::string const& ) { return *this; }
    NdebugStream& operator<<( double ) { return *this; }
    NdebugStream& operator<<( LNManipFunction /*f*/ ) { return *this; }
    //@}
};

inline NdebugStream& perror( NdebugStream& s ) { return s; }
inline NdebugStream& endl( NdebugStream& s )   { return s; }
inline NdebugStream& flush( NdebugStream& s )    { return s; }

#ifndef NDEBUG_OLD
DebugStream Debug( int area = 0, DebugStream::stprintf = 0 );
DebugStream Debug( bool cond, int area = 0, DebugStream::stprintf = 0 );
#else
#define Debug Ndebug
inline NdebugStream Ndebug( int = 0, NdebugStream::stprintf = &printf ) { return NdebugStream(); }
#endif

DebugStream Warning( int area = 0 );
DebugStream Warning( bool cond, int area = 0 );

DebugStream Error( int area = 0 );
DebugStream Error( bool cond, int area = 0 );

DebugStream Fatal( int area = 0 );
DebugStream Fatal( bool cond, int area = 0 );

}




LifeV::DebugStream& perror( LifeV::DebugStream& s );
LifeV::DebugStream& endl( LifeV::DebugStream& s );
LifeV::DebugStream& flush( LifeV::DebugStream& );


#endif /* __Debug_H */
