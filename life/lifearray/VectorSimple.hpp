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

    The class is a wrap up of Standard Library vector class.
    This class is deprecated and it should be removed when the discussion is over.
    It implements one dimensional vector

    Example:

    VectorSimple<float> a(11);

    a(10)=90; // a[10] will contain 90.0

 */

template <typename DataType>
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
    VectorSimple( const VectorSimple<DataType> & vector );

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
    VectorSimple<DataType> & operator=( const VectorSimple<DataType> & vector );

    //! Access operator
    /*!
        Example: a(10)=90; // a[10] will contain 90.0
        @param i index of the element of the VectorSimple vector
        @return a vector reference 
     */
    vectorReference_Type operator() ( vectorSize_Type const i )
    {
        return ( this->operator[] ( i ) );
    }

    //! Const access operator
    /*!
        Example: a(10)=90; // a[10] will contain 90.0
        @param i index of the element of the VectorSimple vector
        @return a vector const reference  
     */
    vectorConstReference_Type operator() ( vectorSize_Type const i ) const
    {
        return ( this->operator[] ( i ) );
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
        return i >= 0 && i < this->size() ;
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
template <typename DataType>
VectorSimple<DataType>::VectorSimple( const VectorSimple<DataType> & vector ) 
        :
        vector_Type( vector )
{}


//============================================================================
// Operators
//============================================================================
template <typename DataType>
VectorSimple<DataType> &
VectorSimple<DataType>::operator=( const VectorSimple<DataType> & vector )
{
    vector_Type::operator=( vector );
    return *this;
}


//============================================================================
// Methods
//============================================================================
template <typename DataType> 
void
VectorSimple<DataType>::clearVector()
{
    vector_Type tmp;
    this->clear();
    this->swap( tmp );
}

template <typename DataType>
typename VectorSimple<DataType>::data_Type& 
VectorSimple<DataType>::returnElementAfterCheck( vectorSize_Type i )
{
    if ( ! checkIndex( i ) )
        abort();
    return *( this->begin() + ( i ) );
}

template <typename DataType>
const typename VectorSimple<DataType>::data_Type&  
VectorSimple<DataType>::returnElementAfterCheck( vectorSize_Type i ) const
{
    if ( ! checkIndex( i ) )
        abort();
    return *( this->begin() + ( i ) );
}

}// Namespace LifeV

#endif /* _SIMPLEVECT_HH_ */

