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

#ifndef _SIMPLEARRAY_HH_
#define _SIMPLEARRAY_HH_

#include <lifev/core/LifeV.hpp>

namespace LifeV
{
//! ArraySimple
/*!
    @author Luca Formaggia

    The class is a wrap up of Standard Library vector class.
    It defines arrays of n x m size.
    This class is deprecated and it should be removed when the
    discussion is over. It is with COLUMNWISE ORDERING

    Example:

    ArraySimple<int> b(3,5) is an array with 3 rows and 5 columns

 */

template <typename DataType>
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
       Construct a ArraySimple vector of size (0,0)
     */
    explicit ArraySimple();

    //! Constructor
    /*!
        Construct a ArraySimple vector of size (size,1)
        @param vectorSize size of the ArraySimple vector
     */
    explicit ArraySimple ( vectorSize_Type vectorSize );

    //! Constructor
    /*!
        Construct a ArraySimple object of size (numberOfRows,numberOfColumns)
        @param numberOfRows number of rows
        @param numberOfColumns number of columns
     */
    explicit ArraySimple ( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns );

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
        return * ( this->begin() + ( i ) );
    }

    //! Const access operator
    /*!
        ArraySimple is seen as vector
        @param i index of the element of the ArraySimple vector
        @return a vector const reference
     */
    vectorConstReference_Type operator() ( vectorSize_Type const i ) const
    {
        return * ( this->begin() + ( i ) );
    }

    //! Access operator
    /*!
        @param i row index
        @param j column index
        @return a vector reference
     */
    vectorReference_Type operator() ( vectorSize_Type const i, vectorSize_Type const j )
    {
        return * ( this->begin() + ( j ) * M_numberOfRows + ( i ) );
    }

    //! Const access operator
    /*!
        @param i row index
        @param j column index
        @return a vector const reference
     */
    vectorConstReference_Type operator() ( vectorSize_Type const i, vectorSize_Type const j ) const
    {
        return * ( this->begin() + ( j ) * M_numberOfRows + ( i ) );
    }

    //@}


    //! @name Methods
    //@{

    //!Return the iterator of the column passed by argument
    /*!
        @param column column index
        @return ArraySimple iterator
     */
    inline typename ArraySimple<DataType>::iterator columnIterator ( vectorSize_Type const column );

    //!Resize the ArraySimple vector
    /*!
        @param numberOfRows number of rows
        @param numberOfColumns number of columns
     */
    void reshape ( vectorSize_Type const numberOfRows, vectorSize_Type const numberOfColumns );

    //! Completely clear out the container, returning memory to the system
    inline void clearArray();

    //! Check if the MeshEntityContainer vector contains an element with row index i and column index j
    /*!
        @param i row index
        @param j column index
        @return boolean
     */
    bool checkIndex ( vectorSize_Type const i, vectorSize_Type const j ) const
    {
        return i >= 0 && i + ( j ) * M_numberOfRows < this->size();
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
template <typename DataType>
ArraySimple<DataType>::ArraySimple()
    :
    std::vector<DataType>(),
    M_numberOfRows ( 0 ),
    M_numberOfColumns ( 0 )
{}

template <typename DataType>
ArraySimple<DataType>::ArraySimple ( vectorSize_Type vectorSize )
    :
    std::vector<DataType> ( vectorSize ),
    M_numberOfRows ( vectorSize ),
    M_numberOfColumns ( 1 )
{}

template <typename DataType>
ArraySimple<DataType>::ArraySimple ( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns )
    :
    std::vector<DataType> ( numberOfRows* numberOfColumns ),
    M_numberOfRows ( numberOfRows ),
    M_numberOfColumns ( numberOfColumns )
{}

//============================================================================
// Methods
//============================================================================
template <typename DataType>
typename ArraySimple<DataType>::iterator
ArraySimple<DataType>::columnIterator ( vectorSize_Type const column )
{
    if ( column > M_numberOfColumns )
    {
        return typename ArraySimple<DataType>::iterator();
    }
    else
    {
        return this->begin() + ( column ) * M_numberOfRows;
    }
}


template <typename DataType>
void
ArraySimple<DataType>::reshape ( vectorSize_Type numberOfRows, vectorSize_Type numberOfColumns )
{
    vector_Type::resize ( numberOfRows * numberOfColumns ); // Standard Library vector method
    M_numberOfRows = numberOfRows;
    M_numberOfColumns = numberOfColumns;
}

template <typename DataType>
void
ArraySimple<DataType>::clearArray()
{
    vector_Type tmp;
    this->clear();
    this->swap ( tmp );
    M_numberOfRows = 0;
    M_numberOfColumns = 0;
}

template <typename DataType>
void
ArraySimple<DataType>::showMe() const
{
    std::cout << " Number of rows: " << M_numberOfRows << std::endl;
    std::cout << " Number of columns: " << M_numberOfColumns << std::endl;
}

}// Namespace LifeV

#endif /* _SIMPLEARRAY_HH_ */
