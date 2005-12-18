/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#ifndef _SIMPLEVECT_HH_
#define _SIMPLEVECT_HH_

#include <vector>

namespace LifeV
{
/******************************************************************
/#Version Luca Formaggia 30 Agu 1999 Experimental
/
/#Purpose: a stupid wrap up of stl vector Class to allow indexing from one.
/
/#Usage:
/                SimpleVect<float> a(10);
/                a(10)=90; // a[9] will contain 90.0
/                SimpleArray,int> b(3,5) // an arrray with 3 rows
/                                           and 5 columns
/                b(3,2)=5;
/                b.reshape(2,3) // now b is 2x3
/
******************************************************************/

template <typename T, int OFFSETVEC = 1>
class SimpleVect : public std::vector<T>
{
public:
    typedef std::vector<T> raw_container;
    typedef typename raw_container::size_type size_type;
    typedef typename raw_container::reference reference;
    typedef typename raw_container::const_reference const_reference;

    explicit SimpleVect( size_type i ) : raw_container( i )
    {}
    ;
    SimpleVect() : raw_container()
    {}
    ;
    SimpleVect( const SimpleVect<T, OFFSETVEC> & );
  
    ~SimpleVect() {}

    explicit SimpleVect( const raw_container & );
    SimpleVect<T, OFFSETVEC> & operator=( const SimpleVect<T, OFFSETVEC> & );
    T& fat( size_type i );
    const T& fat( size_type i ) const;
    inline reference operator() ( size_type const i )
    {
        return * ( this->begin() + ( i - OFFSETVEC ) );
    }
    inline const_reference operator() ( size_type const i ) const
    {
        return * ( this->begin() + ( i - OFFSETVEC ) );
    }
    //! Completely clear out the container, returning memory to the system
    void clean()
    {
        raw_container tmp;
        this->clear();
        this->swap( tmp );
    }
    inline bool bCheck( size_type const i ) const
    {
        return i >= OFFSETVEC && i < this->size() + OFFSETVEC ;
    }
};

template <typename T, int OFFSETVEC = 1>
class SimpleArray : public std::vector<T>
{
public:

    typedef std::vector<T> raw_container;
    typedef typename raw_container::size_type size_type;
    typedef typename raw_container::reference reference;
    typedef typename raw_container::const_reference const_reference;
    typedef typename raw_container::iterator iterator;
    typedef typename raw_container::const_iterator const_iterator;

    explicit SimpleArray( size_type ntot );
    explicit SimpleArray();
    explicit SimpleArray( size_type nrows, size_type ncols );

    ~SimpleArray() {}

    reference operator() ( size_type const i )
    {
        return * ( this->begin() + ( i - OFFSETVEC ) );
    } // Array is seen as vector (index from OFFSETVEC)

    const_reference operator() ( size_type const i ) const
    {
        return * ( this->begin() + ( i - OFFSETVEC ) );
    } // Array is seen as vector (index from OFFSETVEC)

    reference operator() ( size_type const i, size_type const j )
    {
        return * ( this->begin() + ( j - OFFSETVEC ) * _M_nrows + ( i - OFFSETVEC ) );
    } // from OFFSETVEC

    const_reference operator() ( size_type const i, size_type const j ) const
    {
        // from OFFSETVEC
        return * ( this->begin() + ( j - OFFSETVEC ) * _M_nrows + ( i - OFFSETVEC ) );
    }

    typename SimpleArray<T, OFFSETVEC>::iterator columnIterator( size_type const col )
    {
        if ( col > _M_ncols )
            return typename SimpleArray<T, OFFSETVEC>::iterator();
        else
            return this->begin() + ( col - OFFSETVEC ) * _M_nrows;
    }

    void reshape( SimpleArray<T, OFFSETVEC>::size_type const n, size_type const m );

    //! Completely clear out the container, returning memory to the system
    void clean()
    {
        raw_container tmp;
        this->clear();
        this->swap( tmp );
        _M_nrows = 0;
        _M_ncols = 0;
    }

    bool bCheck( size_type const i, size_type const j ) const
    {
        return i >= OFFSETVEC && i - OFFSETVEC + ( j - OFFSETVEC ) * _M_nrows < this->size() ;
    }

    inline size_type nrows() const
    {
        return _M_nrows;
    }
    inline size_type ncols() const
    {
        return _M_ncols;
    }

private:
    size_type _M_nrows;
    size_type _M_ncols;
};


template <typename T, int OFFSETVEC>
T& SimpleVect<T, OFFSETVEC>::fat( size_type i )
{
    if ( ! bCheck( i ) )
        abort();
    return *( this->begin() + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
const T& SimpleVect<T, OFFSETVEC>::fat( size_type i ) const
{
    if ( ! bCheck( i ) )
        abort();
    return *( this->begin() + ( i - OFFSETVEC ) );
}

template <typename T, int OFFSETVEC>
SimpleVect<T, OFFSETVEC>::SimpleVect( const SimpleVect<T, OFFSETVEC> & v ) : raw_container( v )
{}

template <typename T, int OFFSETVEC>
SimpleVect<T, OFFSETVEC> &
SimpleVect<T, OFFSETVEC>::operator=( const SimpleVect<T, OFFSETVEC> & v )
{
    raw_container::operator=( v );
    return *this;
}


// SIMPLE ARRAYS

template <typename T, int OFFSETVEC>
SimpleArray<T, OFFSETVEC>::SimpleArray()
        :
        std::vector<T>(),
        _M_nrows( 0 ),
        _M_ncols( 1 )
{}

template <typename T, int OFFSETVEC>
SimpleArray<T, OFFSETVEC>::SimpleArray( size_type ntot )
        :
        std::vector<T>( ntot ),
        _M_nrows( nrows ),
        _M_ncols( 1 )
{}

template <typename T, int OFFSETVEC>
SimpleArray<T, OFFSETVEC>::SimpleArray( size_type nrows, size_type ncols )
        :
        std::vector<T>( nrows*ncols ),
        _M_nrows( nrows ),
        _M_ncols( ncols )
{}

template <typename T, int OFFSETVEC>
void
SimpleArray<T, OFFSETVEC>::reshape( size_type nrows, size_type ncols )
{
    raw_container::resize( nrows * ncols ); // stl vector method
    _M_nrows = nrows;
    _M_ncols = ncols;
}
}
#endif

