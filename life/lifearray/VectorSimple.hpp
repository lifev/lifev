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
#include <life/lifecore/LifeV.hpp>

namespace LifeV
{
//! VectorSimple
/*!
    @author Luca Formaggia

    The class is a wrap up of Standard Library vector class to allow indexing from one.
    It implements one dimensional vector

    Example:

    VectorSimple<float> a(10);

    a(10)=90; // a[9] will contain 90.0

 */

template <typename DataType, int offsetVector = 1>
class VectorSimple : public std::vector<DataType>
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
    VectorSimple() : vector_Type() {};

    //! Constructor
    /*!
        @param vectorSize size of the vector
     */
    explicit VectorSimple( vectorSize_Type vectorSize ) : vector_Type( vectorSize ){};

    //! Copy constructor
    /*!
        @param vector VectorSimple vector to copy
     */
    VectorSimple( const VectorSimple<DataType, offsetVector> & vector );

    //! Constructor
    /*!
        Construct by copying a Standard Library vector
        @param vector Standard Library vector 
     */
    explicit VectorSimple( const vector_Type & vector );

    //! Destructor
    ~VectorSimple() {}

    //@}


    //! @name Operators
    //@{

    //! Equivalence operator
    /*!
        Copies source VectorSimple vector into "this"
        @param vector VectorSimple vector
        @return Reference to a new VectorSimple vector with the same
                content of VectorSimple vector
     */
    VectorSimple<DataType, offsetVector> & operator=( const VectorSimple<DataType, offsetVector> & vector );

    //! Access operator
    /*!
        Example: a(10)=90; // a[9] will contain 90.0
        @param i index of the element of the VectorSimple vector
        @return a vector reference 
     */
    vectorReference_Type operator() ( vectorSize_Type const i )
    {
        return ( this->operator[] ( i - offsetVector ) );
    }

    //! Const access operator
    /*!
        Example: a(10)=90; // a[9] will contain 90.0
        @param i index of the element of the VectorSimple vector
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
    inline void clearVector();

    //! Check if the VectorSimple vector contains an element with index i
    /*!
        @param i index
        @return boolean 
     */
    bool checkIndex( vectorSize_Type const i ) const
    {
        return i >= offsetVector && i < this->size() + offsetVector ;
    }  

    //! Return a reference to the element of type data_Type stored in the i-th position.
    //! The operator aborts if the VectorSimple vector doesn't contain an element with index i
    /*!
        @param i index
        @return reference to the data_Type element
     */
    data_Type& returnElementAfterCheck( vectorSize_Type i );

    //! Return a const reference to the element of type data_Type stored in the i-th position.
    //! The operator aborts if the VectorSimple vector doesn't contain an element with index i
    /*!
        @param i index
        @return const reference to the data_Type element
     */
    const data_Type& returnElementAfterCheck( vectorSize_Type i ) const;

    //@}

};

//============================================================================
// Constructors
//============================================================================
template <typename DataType, int offsetVector>
VectorSimple<DataType, offsetVector>::VectorSimple( const VectorSimple<DataType, offsetVector> & vector ) 
        :
        vector_Type( vector )
{}


//============================================================================
// Operators
//============================================================================
template <typename DataType, int offsetVector>
VectorSimple<DataType, offsetVector> &
VectorSimple<DataType, offsetVector>::operator=( const VectorSimple<DataType, offsetVector> & vector )
{
    vector_Type::operator=( vector );
    return *this;
}


//============================================================================
// Methods
//============================================================================
template <typename DataType, int offsetVector> 
void
VectorSimple<DataType, offsetVector>::clearVector()
{
    vector_Type tmp;
    this->clear();
    this->swap( tmp );
}

template <typename DataType, int offsetVector>
typename VectorSimple<DataType, offsetVector>::data_Type& 
VectorSimple<DataType, offsetVector>::returnElementAfterCheck( vectorSize_Type i )
{
    if ( ! checkIndex( i ) )
        abort();
    return *( this->begin() + ( i - offsetVector ) );
}

template <typename DataType, int offsetVector>
const typename VectorSimple<DataType, offsetVector>::data_Type&  
VectorSimple<DataType, offsetVector>::returnElementAfterCheck( vectorSize_Type i ) const
{
    if ( ! checkIndex( i ) )
        abort();
    return *( this->begin() + ( i - offsetVector ) );
}


//! ArraySimple
/*!
    @author Luca Formaggia

    The class is a wrap up of Standard Library vector class to allow indexing from one.
    It defines vectors of n x m size.

    Example:

    ArraySimple<int> b(3,5) is an array with 3 rows and 5 columns

    b(3,2)=5 puts the number 5 in the third row, in the second column.
   
 */

template <typename DataType, int offsetVector = 1>
class ArraySimple : public std::vector<DataType>
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
       Construct a ArraySimple vector of size (0,1)
     */
    explicit ArraySimple();

    //! Constructor
    /*!
        Construct a ArraySimple vector of size (size,1)
        @param vectorSize size of the ArraySimple vector
     */
    explicit ArraySimple( vectorSize_Type vectorSize );

    //! Constructor
    /*!
        Construct a ArraySimple object of size (numberOfRows,numberOfColumns)
        @param numberOfRows number of rows
        @param numberOfColumns number of columns
     */
    explicit ArraySimple( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns );

    //! Destructor
    ~ArraySimple() {}

    //@}


    //! @name Operators
    //@{

    //! Access operator
    /*!
        ArraySimple is seen as vector (index from one)
        @param i index of the element of the ArraySimple vector
        @return a vector reference 
     */
    vectorReference_Type operator() ( vectorSize_Type const i )
    {
        return *( this->begin() + ( i - offsetVector ) );
    }

    //! Const access operator
    /*!
        ArraySimple is seen as vector (index from one)
        @param i index of the element of the ArraySimple vector
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

    //!Return the iterator of the column passed by argument
    /*!
        @param column column index 
        @return ArraySimple iterator
     */
    inline typename ArraySimple<DataType, offsetVector>::iterator columnIterator( vectorSize_Type const column );

    //!Resize the ArraySimple vector
    /*!
        @param numberOfRows number of rows
        @param numberOfColumns number of columns 
     */
    void reshape( vectorSize_Type const numberOfRows, vectorSize_Type const numberOfColumns );

    //! Completely clear out the container, returning memory to the system
    inline void clearArray();

    //! Check if the VectorSimple vector contains an element with row index i and column index j
    /*!
        @param i row index
        @param j column index
        @return boolean
     */
    bool checkIndex( vectorSize_Type const i, vectorSize_Type const j ) const
    {
        return i >= offsetVector && i - offsetVector + ( j - offsetVector ) * M_numberOfRows < this->size();
    }

    //! Show the contents of the ArraySimple vector
    void showMe() const;

    //@}


    //! @name Get Methods
    //@{

    //!Return the number of rows
    /*!
        @return a vector size type variable holding the number of rows
     */
    vectorSize_Type const numberOfRows() const
    {
        return M_numberOfRows;
    }

    //!Return the number of columns
    /*!
        @return a vector size type variable holding the number of columns
     */
    vectorSize_Type const numberOfColumns() const
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
ArraySimple<DataType, offsetVector>::ArraySimple()
        :
        std::vector<DataType>(),
        M_numberOfRows( 0 ),
        M_numberOfColumns( 1 )
{}

template <typename DataType, int offsetVector>
ArraySimple<DataType, offsetVector>::ArraySimple( vectorSize_Type vectorSize )
        :
        std::vector<DataType>( vectorSize ),
        M_numberOfRows( vectorSize ),
        M_numberOfColumns( 1 )
{}

template <typename DataType, int offsetVector>
ArraySimple<DataType, offsetVector>::ArraySimple( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns )
        :
        std::vector<DataType>( numberOfRows * numberOfColumns ),
        M_numberOfRows( numberOfRows ),
        M_numberOfColumns( numberOfColumns )
{}

//============================================================================
// Methods
//============================================================================
template <typename DataType, int offsetVector>
typename ArraySimple<DataType, offsetVector>::iterator 
ArraySimple<DataType, offsetVector>::columnIterator( vectorSize_Type const column )
{
    if ( column > M_numberOfColumns )
        return typename ArraySimple<DataType, offsetVector>::iterator();
    else
        return this->begin() + ( column - offsetVector ) * M_numberOfRows;
}


template <typename DataType, int offsetVector>
void
ArraySimple<DataType, offsetVector>::reshape( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns )
{
    vector_Type::resize( numberOfRows * numberOfColumns ); // Standard Library vector method
    M_numberOfRows = numberOfRows;
    M_numberOfColumns = numberOfColumns;
}

template <typename DataType, int offsetVector>
void 
ArraySimple<DataType, offsetVector>::clearArray()
{
    vector_Type tmp;
    this->clear();
    this->swap( tmp );
    M_numberOfRows = 0;
    M_numberOfColumns = 0;
}

template <typename DataType, int offsetVector>
void 
ArraySimple<DataType, offsetVector>::showMe() const
{
    std::cout << " Number of rows: " << M_numberOfRows << std::endl;
    std::cout << " Number of columns: " << M_numberOfColumns << std::endl;
}

}// Namespace LifeV

#endif /* _SIMPLEVECT_HH_ */

