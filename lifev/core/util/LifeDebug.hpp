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
  @brief Classes for debugging

  @date 13-12-2010
  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef LIFE_DEBUG_H
#define LIFE_DEBUG_H 1

#include <cstdio>
#include <iosfwd>
#include <sstream>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

// Forward declarations
class DebugStream;
class NdebugStream;

typedef DebugStream& (*LManipFunction) ( DebugStream&);  // manipulator function
typedef NdebugStream& (*LNManipFunction) ( NdebugStream&); // manipulator function

#ifdef __GNUC__
# define LIFEV_FUNCINFO "[" << __PRETTY_FUNCTION__ << "] "
#else
# define LIFEV_FUNCINFO "[" << __FILE__ << ":" << __LINE__ << "] "
#endif

#define LIFEV_LINEINFO "[" << __FILE__ << ":" << __LINE__ << "] "

class DebugStream
{
public:
    //! @name Public typedefs and structures
    //@{
    struct Private;
    typedef int (*stprintf) ( const char* format, ... );
    //@}

    /** @name Constructors, destructor
     */
    //@{
    DebugStream (int area = 0, int level = 1, bool print = true);
    DebugStream (const char* initialString, int area = 0, int level = 1, bool print = true);
    DebugStream ( DebugStream const& );
    ~DebugStream();
    //@}

    //! @name Operators
    //@{
    DebugStream& operator<< ( char const* c);
    DebugStream& operator<< ( double d);
    DebugStream& operator<< ( std::string const& str);
    DebugStream& operator<< ( LManipFunction f);
    //@}

    //! @name  Methods
    //@{
    void setFlush ( stprintf = 0 );
    void flush();

    static void attach ( std::string const& logfile );
    static void attach ( std::string const& logfile, int area );
    static void detach ( std::string const& logfile, int area );
    static void detachAll();
    //@}

private:
    Private* M_data;
};

// ===================================
// Debug Stream Implementation
// ===================================

template<typename T>
DebugStream& operator<< ( DebugStream& stream, T const* data )
{
    std::ostringstream out_stream;
    out_stream << data;
    stream << out_stream.str();
    return stream;
}

#ifdef HAVE_BACKTRACE
std::string backtrace ();
std::string backtrace ( int val);
#endif


class NdebugStream
{
public:
    //! @name Public typedefs
    //@{
    typedef int (*stprintf) ( const char* format, ... );
    //@}

    //! @name Constructors, destructor
    //@{
    NdebugStream() {}
    ~NdebugStream() {}
    //@}

    //! @name Operators
    //@{
    NdebugStream& operator<< ( char const* /*code*/ )
    {
        return *this;
    }
    NdebugStream& operator<< ( std::string const& /*str*/)
    {
        return *this;
    }
    NdebugStream& operator<< ( double /*code*/)
    {
        return *this;
    }
    NdebugStream& operator<< ( LNManipFunction /*f*/ )
    {
        return *this;
    }
    //@}

    //! @name  Methods
    //@{
    void flush ( stprintf = 0 ) {}
    //@}
};

inline NdebugStream& perror ( NdebugStream& s )
{
    return s;
}
inline NdebugStream& endl ( NdebugStream& s )
{
    return s;
}
inline NdebugStream& flush ( NdebugStream& s )
{
    return s;
}

#ifdef HAVE_LIFEV_DEBUG
DebugStream debugStream ( int area = 0, DebugStream::stprintf = 0 );
DebugStream debugStream ( bool cond, int area = 0, DebugStream::stprintf = 0 );
LIFEV_DEPRECATED ( DebugStream Debug ( int area = 0, DebugStream::stprintf func = 0 ) );
LIFEV_DEPRECATED ( DebugStream Debug ( bool cond, int area = 0, DebugStream::stprintf func = 0 ) );
#else
#define debugStream noDebugStream
inline NdebugStream noDebugStream ( int = 0, NdebugStream::stprintf = &printf )
{
    return NdebugStream();
}
#endif

DebugStream Warning ( int area = 0 );
DebugStream Warning ( bool cond, int area = 0 );

DebugStream Error ( int area = 0 );
DebugStream Error ( bool cond, int area = 0 );

DebugStream Fatal ( int area = 0 );
DebugStream Fatal ( bool cond, int area = 0 );

}

LifeV::DebugStream& perror ( LifeV::DebugStream& s );
LifeV::DebugStream& endl ( LifeV::DebugStream& s );
LifeV::DebugStream& flush ( LifeV::DebugStream& s );

#endif // LIFE_DEBUG_H
