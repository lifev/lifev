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
 *  @brief The file contains two classes implementing a wrap up
           of Standard Library vector class to allow indexing from one.
 *
 *  @date 30-08-1999
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *
 *  @contributor Laura Cattaneo
 *  @mantainer Laura Cattaneo
 */

#ifndef _SIMPLEVECT_HH_
#define _SIMPLEVECT_HH_

#include <cstdlib>
#include <vector>
#include <life/lifecore/life.hpp>

namespace LifeV
{
//! SimpleVect
/*!
    @author Luca Formaggia

    The class is a wrap up of Standard Library vector class to allow indexing from one.
    It implements one dimensional vector

    Example:

    SimpleVect<float> a(10);

    a(10)=90; // a[9] will contain 90.0

 */

template <typename DataType, int offsetVector = 1>
class SimpleVect : public std::vector<DataType>
{
public:

    //! @name Public Types
    //@{

    typedef DataType                                data_Type;
    typedef std::vector<DataType>                   vector_Type;
    typedef typename vector_Type::size_type         vectorSize_Type;
    typedef typename vector_Type::reference         vectorReference_Type;
    typedef typename vector_Type::const_reference   vectorConstReference_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    SimpleVect() : vector_Type() {};

    //! Constructor
    /*!
        @param size size of the vector
     */
    explicit SimpleVect( vectorSize_Type size ) : vector_Type( size ){};

    //! Copy constructor
    /*!
        @param vector SimpleVect vector to copy
     */
    SimpleVect( const SimpleVect<DataType, offsetVector> & vector );

    //! Constructor
    /*!
        Construct by copying a Standard Library vector
        @param vector Standard Library vector
     */
    explicit SimpleVect( const vector_Type & vector );

    //! Destructor
    ~SimpleVect() {}

    //@}


    //! @name Operators
    //@{

    //! Equivalence operator
    /*!
        Copies source SimpleVect vector into "this"
        @param vector SimpleVect vector
        @return Reference to a new SimpleVect vector with the same
                content of SimpleVect vector
     */
    SimpleVect<DataType, offsetVector> & operator=( const SimpleVect<DataType, offsetVector> & vector );

    //! Access operator
    /*!
        Example: a(10)=90; // a[9] will contain 90.0
        @param i index of the element of the SimpleVect vector
        @return a vector reference
     */
    vectorReference_Type operator() ( vectorSize_Type const i )
    {
        return ( this->operator[] ( i - offsetVector ) );
    }

    //! Const access operator
    /*!
        Example: a(10)=90; // a[9] will contain 90.0
        @param i index of the element of the SimpleVect vector
        @return a vector const reference
     */
    vectorConstReference_Type operator() ( vectorSize_Type const i ) const
    {
        return ( this->operator[] ( i - offsetVector ) );
    }

    //@}


    //! @name Methods
    //@{

    //! Completely clear out the container, returning memory to the system
    inline void clean();

    //! Check if the SimpleVect vector contains an element with index i
    /*!
        @param i index
        @return boolean
     */
    bool bCheck( vectorSize_Type const i ) const
    {
        return i >= offsetVector && i < this->size() + offsetVector ;
    }

    //! Return a reference to the element of type data_Type stored in the i-th position.
    //! The operator aborts if the SimpleVect vector doesn't contain an element with index i
    /*!
        @param i index
        @return reference to the data_Type element
     */
    data_Type& fat( vectorSize_Type i );

    //! Return a const reference to the element of type data_Type stored in the i-th position.
    //! The operator aborts if the SimpleVect vector doesn't contain an element with index i
    /*!
        @param i index
        @return const reference to the data_Type element
     */
    const data_Type& fat( vectorSize_Type i ) const;

    //@}

};

//============================================================================
// Constructors
//============================================================================
template <typename DataType, int offsetVector>
SimpleVect<DataType, offsetVector>::SimpleVect( const SimpleVect<DataType, offsetVector> & vector )
        :
        vector_Type( vector )
{}


//============================================================================
// Operators
//============================================================================
template <typename DataType, int offsetVector>
SimpleVect<DataType, offsetVector> &
SimpleVect<DataType, offsetVector>::operator=( const SimpleVect<DataType, offsetVector> & vector )
{
    vector_Type::operator=( vector );
    return *this;
}


//============================================================================
// Methods
//============================================================================
template <typename DataType, int offsetVector>
void
SimpleVect<DataType, offsetVector>::clean()
{
    vector_Type tmp;
    this->clear();
    this->swap( tmp );
}

template <typename DataType, int offsetVector>
typename SimpleVect<DataType, offsetVector>::data_Type&
SimpleVect<DataType, offsetVector>::fat( vectorSize_Type i )
{
    if ( ! bCheck( i ) )
        abort();
    return *( this->begin() + ( i - offsetVector ) );
}

template <typename DataType, int offsetVector>
const typename SimpleVect<DataType, offsetVector>::data_Type&
SimpleVect<DataType, offsetVector>::fat( vectorSize_Type i ) const
{
    if ( ! bCheck( i ) )
        abort();
    return *( this->begin() + ( i - offsetVector ) );
}


//! SimpleArray
/*!
    @author Luca Formaggia

    The class is a wrap up of Standard Library vector class to allow indexing from one.
    It defines vectors of n x m size.

    Example:

    SimpleArray<int> b(3,5) is an array with 3 rows and 5 columns

    b(3,2)=5 puts the number 5 in the third row, in the second column.

 */

template <typename DataType, int offsetVector = 1>
class SimpleArray : public std::vector<DataType>
{
public:

    //! @name Public Types
    //@{

    typedef DataType                                data_Type;
    typedef std::vector<DataType>                   vector_Type;
    typedef typename vector_Type::size_type         vectorSize_Type;
    typedef typename vector_Type::reference         vectorReference_Type;
    typedef typename vector_Type::const_reference   vectorConstReference_Type;
    typedef typename vector_Type::iterator          vectorIterator_Type;
    typedef typename vector_Type::const_iterator    vectorConstIterator_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    /*!
       Construct a SimpleArray vector of size (0,1)
     */
    explicit SimpleArray();

    //! Constructor
    /*!
        Construct a SimpleArray vector of size (size,1)
        @param size size of the SimpleArray vector
     */
    explicit SimpleArray( vectorSize_Type size );

    //! Constructor
    /*!
        Construct a SimpleArray object of size (numberOfRows,numberOfColumns)
        @param numberOfRows number of rows
        @param numberOfColumns number of columns
     */
    explicit SimpleArray( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns );

    //! Destructor
    ~SimpleArray() {}

    //@}


    //! @name Operators
    //@{

    //! Access operator
    /*!
        SimpleArray is seen as vector (index from one)
        @param i index of the element of the SimpleArray vector
        @return a vector reference
     */
    vectorReference_Type operator() ( vectorSize_Type const i )
    {
        return *( this->begin() + ( i - offsetVector ) );
    }

    //! Const access operator
    /*!
        SimpleArray is seen as vector (index from one)
        @param i index of the element of the SimpleArray vector
        @return a vector const reference
     */
    vectorConstReference_Type operator() ( vectorSize_Type const i ) const
    {
        return *( this->begin() + ( i - offsetVector ) );
    }

    //! Access operator
    /*!
        Indices start from one
        @param i row index
        @param j column index
        @return a vector reference
     */
    vectorReference_Type operator() ( vectorSize_Type const i, vectorSize_Type const j )
    {
        return *( this->begin() + ( j - offsetVector ) * M_numberOfRows + ( i - offsetVector ) );
    }

    //! Const access operator
    /*!
        Indices start from one
        @param i row index
        @param j column index
        @return a vector const reference
     */
    vectorConstReference_Type operator() ( vectorSize_Type const i, vectorSize_Type const j ) const
    {
        return *( this->begin() + ( j - offsetVector ) * M_numberOfRows + ( i - offsetVector ) );
    }

    //@}


    //! @name Methods
    //@{

    //!Return the column iterator
    /*!
        @param column column index
        @return SimpleArray iterator
     */
    inline typename SimpleArray<DataType, offsetVector>::iterator columnIterator( vectorSize_Type const column );

    //!Resize the SimpleArray vector
    /*!
        @param numberOfRows number of rows
        @param numberOfColumns number of columns
     */
    void reshape( vectorSize_Type const numberOfRows, vectorSize_Type const numberOfColumns );

    //! Completely clear out the container, returning memory to the system
    inline void clean();

    //! Check if the SimpleVect vector contains an element with row index i and column index j
    /*!
        @param i row index
        @param j column index
        @return boolean
     */
    bool bCheck( vectorSize_Type const i, vectorSize_Type const j ) const
    {
        return i >= offsetVector && i - offsetVector + ( j - offsetVector ) * M_numberOfRows < this->size();
    }

    //! Show the contents of the SimpleArray vector
    void showMe() const;

    //@}


    //! @name Get Methods
    //@{

    //!Return the number of rows
    /*!
        @return a vector size type variable holding the number of rows
     */
    vectorSize_Type nrows() const
    {
        return M_numberOfRows;
    }

    //!Return the number of columns
    /*!
        @return a vector size type variable holding the number of columns
     */
    vectorSize_Type ncols() const
    {
        return M_numberOfColumns;
    }

    //@}

private:

    //! Number of rows
    vectorSize_Type M_numberOfRows;

    //! Number of columns
    vectorSize_Type M_numberOfColumns;
};


//============================================================================
// Constructors
//============================================================================
template <typename DataType, int offsetVector>
SimpleArray<DataType, offsetVector>::SimpleArray()
        :
        std::vector<DataType>(),
        M_numberOfRows( 0 ),
        M_numberOfColumns( 1 )
{}

template <typename DataType, int offsetVector>
SimpleArray<DataType, offsetVector>::SimpleArray( vectorSize_Type size )
        :
        std::vector<DataType>( size ),
        M_numberOfRows( size ),
        M_numberOfColumns( 1 )
{}

template <typename DataType, int offsetVector>
SimpleArray<DataType, offsetVector>::SimpleArray( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns )
        :
        std::vector<DataType>( numberOfRows * numberOfColumns ),
        M_numberOfRows( numberOfRows ),
        M_numberOfColumns( numberOfColumns )
{}

//============================================================================
// Methods
//============================================================================
template <typename DataType, int offsetVector>
typename SimpleArray<DataType, offsetVector>::iterator
SimpleArray<DataType, offsetVector>::columnIterator( vectorSize_Type const column )
{
    if ( column > M_numberOfColumns )
        return typename SimpleArray<DataType, offsetVector>::iterator();
    else
        return this->begin() + ( column - offsetVector ) * M_numberOfRows;
}


template <typename DataType, int offsetVector>
void
SimpleArray<DataType, offsetVector>::reshape( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns )
{
    vector_Type::resize( numberOfRows * numberOfColumns ); // Standard Library vector method
    M_numberOfRows = numberOfRows;
    M_numberOfColumns = numberOfColumns;
}

template <typename DataType, int offsetVector>
void
SimpleArray<DataType, offsetVector>::clean()
{
    vector_Type tmp;
    this->clear();
    this->swap( tmp );
    M_numberOfRows = 0;
    M_numberOfColumns = 0;
}

template <typename DataType, int offsetVector>
void
SimpleArray<DataType, offsetVector>::showMe() const
{
    std::cout << " Number of rows: " << M_numberOfRows << std::endl;
    std::cout << " Number of columns: " << M_numberOfColumns << std::endl;
}

}// Namespace LifeV

#endif /* _SIMPLEVECT_HH_ */

