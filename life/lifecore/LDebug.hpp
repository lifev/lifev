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
#ifndef __LDebug_H
#define __LDebug_H 1

#include <cstdio>
#include <iosfwd>

#include <string>
#include <sstream>

namespace LifeV
{
    class LDebugStream;
    class LNDebugStream;

    typedef LDebugStream & (*LManipFunction)( LDebugStream &); // manipulator function
    typedef LNDebugStream & (*LNManipFunction)( LNDebugStream&); // manipulator function

#ifdef __GNUC__
# define L_FUNCINFO "[" << __PRETTY_FUNCTION__ << "] "
#else
# define L_FUNCINFO "[" << __FILE__ << ":" << __LINE__ << "] "
#endif

#define L_LINEINFO "[" << __FILE__ << ":" << __LINE__ << "] "

    namespace detail
    {
#if 0
	class print
	{
	public:
	    explicit print ( std::ostream& __os ) : _M_os ( __os ) {}

	    template<typename T>
	    std::ostream& operator<< ( T const& __t )
		{
		    __print ( __t, St::SInt2Type<St::STypeTraits<T>::isFundamental>() );
		    return _M_os;
		}
	    template<typename T>
	    std::ostream& operator<< ( T const* __t )
		{
		    //__print ( __t, St::SInt2Type<St::STypeTraits<T>::isFundamental>() );
		    _M_os << __t;
		    return _M_os;
		}
	private:
	    template<typename T>
	    void __print ( T const& __t, St::SInt2Type<true> )
		{
		    _M_os << __t;
		}
	    template<typename T>
	    void __print ( T const& __t, St::SInt2Type<false> )
		{

		}
	private:
	    std::ostream& _M_os;
	};
#endif
    } // end namespace detail

    class LDebugStream
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

	LDebugStream(int area = 0, int level = 1, bool print = true);
	LDebugStream(const char* initialString, int area = 0, int level = 1, bool print = true);
	LDebugStream( LDebugStream const& );
	~LDebugStream();

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


	LDebugStream& operator<<( char const* );

	LDebugStream& operator<<( double );

	LDebugStream& operator<<( std::string const& );
	LDebugStream& operator<<( LManipFunction f );
	//@}



    protected:

    private:
	Private* __p;

    };

    template<typename T>
    LDebugStream& operator<< ( LDebugStream& __s, T const* __t )
    {
	std::ostringstream __os;
	__os << __t;
	__s << __os.str();
	return __s;
    }
    std::string LBacktrace ();
    std::string LBacktrace ( int );

    class LNDebugStream
    {
    public:
	/** @name Constructors, destructor
	 */
	//@{
	typedef int (*stprintf)( const char* format, ... );

	LNDebugStream(){}
	~LNDebugStream(){}

	//@}

	/** @name  Methods
	 */
	//@{
	void flush( stprintf = 0 ) {}
	LNDebugStream& operator<<( char const* ) { return *this; }
	LNDebugStream& operator<<( std::string const& ) { return *this; }
	LNDebugStream& operator<<( double ) { return *this; }
	LNDebugStream& operator<<( LNManipFunction f ) { return *this; }
	//@}
    };

    inline LNDebugStream& perror( LNDebugStream& s ) { return s; }
    inline LNDebugStream& endl( LNDebugStream& s )   { return s; }
    inline LNDebugStream& flush( LNDebugStream& s )    { return s; }

#ifndef NDEBUG_OLD
    LDebugStream LDebug( int area = 0, LDebugStream::stprintf = 0 );
    LDebugStream LDebug( bool cond, int area = 0, LDebugStream::stprintf = 0 );
#else
#define LDebug LNDebug
    inline LNDebugStream LNDebug( int = 0, LNDebugStream::stprintf = &printf ) { return LNDebugStream(); }
    inline LNDebugStream LNDebug( bool cond, int = 0, LNDebugStream::stprintf = &printf ) { return LNDebugStream(); }
#endif

    LDebugStream LDWarning( int area = 0 );
    LDebugStream LDWarning( bool cond, int area = 0 );

    LDebugStream LDError( int area = 0 );
    LDebugStream LDError( bool cond, int area = 0 );

    LDebugStream LDFatal( int area = 0 );
    LDebugStream LDFatal( bool cond, int area = 0 );

}




LifeV::LDebugStream& perror( LifeV::LDebugStream& s );
LifeV::LDebugStream& endl( LifeV::LDebugStream& s );
LifeV::LDebugStream& flush( LifeV::LDebugStream& );


#endif /* __LDebug_H */
