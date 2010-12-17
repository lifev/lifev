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
 *  @brief Containers Of Vectors
 *
 *  @date 29-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributor Simone Rossi <simone.rossi@epfl.ch>
 *  @mantainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */


#ifndef VectorContainer_H
#define VectorContainer_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{

//! VectorContainer - LifeV vector made of other vectors
/*!
 *  @author Cristiano Malossi
 *
 *  The VectorContainer class has the purpose to contains a certain number of
 *  different vectors as a standard container, but providing some general functions
 *  to perform common algebraic operations.
 *
 *
 *
 *  <b> Type of compatible vectors </b>
 *  This class has been developed in order to work with:
 *  <ol>
 *  	<li> LifeV::EpetraVector
 *  	<li> boost::numeric::ublas::vector
 *  </ol>
 *
 *
 *
 *  <b> Type of compatible containers </b>
 *  This class has been developed in order to work with:
 *  <ol>
 *  	<li> std::vector
 *  	<li> std::dequee
 *  	<li> std::list
 *  </ol>
 *
 *  <b> NOTES </b>
 *  <ol>
 *      <li> The VectorContainer contains shared_ptr to vectors: if you perform an operation (such as the
 *           element by element multiplication) among two containers which hold the same vectors (or at least one),
 *           the results can be different from what expected!
 *      <li> Up to now it has been strongly tested only with LifeV::EpetraVector. However the design should be
 *           more or less compatible with boost::numeric::ublas::vector.
 *  </ol>
 *
 */
template< class VectorType, class ContainerType = std::vector< boost::shared_ptr< VectorType > > >
class VectorContainer
{
public:

    //! @name Constructors & Destructor
    //@{

    typedef VectorType                               vector_Type;
    typedef boost::shared_ptr < vector_Type >        vectorPtr_Type;
    typedef ContainerType                            container_Type;
    typedef typename container_Type::iterator        iterator_Type;
    typedef typename container_Type::const_iterator  constIterator_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    VectorContainer();

    //! Copy constructor
    /*!
     * Makes a copy of all the objects contained in the shared_ptr(s).
     * @param VectorContainer - VectorContainer
     */
    VectorContainer( const VectorContainer& vectorContainer );

    //! Destructor
    virtual ~VectorContainer() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * Make a true copy of all the objects contained in the shared_ptr(s)
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& operator=( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Operator=
    /*!
     * Set all the components of the vector equal to a scalar quantity.
     * Cannot be done for an empty vector!
     * @param scalar - Scalar value
     */
    VectorContainer& operator=( const Real& scalar );

    //! Operator&=
    /*!
     * Make a copy of the shared_ptr(s), without duplicates the objects contained in them.
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& operator&=( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Operator[]
    /*!
     * The element access operator (by reference) has been implemented just for debugging or other special cases,
     * because it is not efficient with this type of implementation.
     *
     * @param id - Element id
     */
    Real& operator[]( const UInt& id ) const;

    //! Operator()
    /*!
     * Vector access operator (by shared_ptr).
     *
     * @param id - Vector id
     */
    vectorPtr_Type& operator()( const UInt& id );

    //! Operator+=
    /*!
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& operator+=( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Operator-=
    /*!
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& operator-=( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Operator*=
    /*!
     * Multiplication element by element of two vectors.
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& operator*=( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Operator/=
    /*!
     * Division element by element of two vectors.
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& operator/=( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Operator+
    /*!
     * @param vectorContainer - VectorContainer
     */
    const VectorContainer operator+( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const;

    //! Operator-
    /*!
     * @param vectorContainer - VectorContainer
     */
    const VectorContainer operator-( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const;

    //! Operator*
    /*!
     * Multiplication element by element of two vectors.
     * @param vectorContainer - VectorContainer
     */
    const VectorContainer operator*( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const;

    //! Operator/
    /*!
     * Division element by element of two vectors.
     * @param vectorContainer - VectorContainer
     */
    const VectorContainer operator/( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const;

    //! Operator+=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    VectorContainer& operator+=( const Real& scalar );

    //! Operator-=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    VectorContainer& operator-=( const Real& scalar );

    //! Operator*=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    VectorContainer& operator*=( const Real& scalar );

    //! Operator/=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    VectorContainer& operator/=( const Real& scalar );

    //! Operator*=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const VectorContainer operator+( const Real& scalar ) const;

    //! Operator/=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const VectorContainer operator-( const Real& scalar ) const;

    //! Operator*=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const VectorContainer operator*( const Real& scalar ) const;

    //! Operator/=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const VectorContainer operator/( const Real& scalar ) const;

    //! operator==
    /*!
     * Return a VectorContainer containing 1 where vector elements are == scalar;
     * @param scalar - Value for the comparison.
     */
    VectorContainer operator==( const Real& scalar );

    //! operator!=
    /*!
     * Return a VectorContainer containing 1 where vector elements are != scalar;
     * @param scalar - Value for the comparison.
     */
    VectorContainer operator!=( const Real& scalar );

    //! operator<
    /*!
     * Return a VectorContainer containing 1 where vector elements are < scalar;
     * @param scalar - Value for the comparison.
     */
    VectorContainer operator<( const Real& scalar );

    //! operator>
    /*!
     * Return a VectorContainer containing 1 where vector elements are > scalar;
     * @param scalar - Value for the comparison.
     */
    VectorContainer operator>( const Real& scalar );

    //! operator<=
    /*!
     * Return a VectorContainer containing 1 where vector elements are <= scalar;
     * @param scalar - Value for the comparison.
     */
    VectorContainer operator<=( const Real& scalar );

    //! operator>=
    /*!
     * Return a VectorContainer containing 1 where vector elements are >= scalar;
     * @param scalar - Value for the comparison.
     */
    VectorContainer operator>=( const Real& scalar );

    //! Logic operator&&
    /*!
     * Return a vector containing one where both the elements are different from zero;
     * @param vectorContainer - VectorContainer
     */
    VectorContainer operator&&( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Logic operator||
    /*!
     * Return a vector containing one if one of the two elements is different from zero;
     * @param vectorContainer - VectorContainer
     */
    VectorContainer operator||( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Logic operator!
    /*!
     * Return a vector containing one where the elements are equal to zero;
     */
    VectorContainer operator!();

    //@}


    //! @name Methods
    //@{

    //! Dot - Scalar product
    /*!
     * Scalar product of the vectors
     * @param vectorContainer - VectorContainer
     */
    Real dot( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const;

    //! Dot - Scalar product
    /*!
     * Scalar product of the vectors
     * @param vectorContainer - VectorContainer
     * @param scalarProduct - result
     */
    void dot( const VectorContainer< vector_Type, container_Type >& vectorContainer, Real& scalarProduct );

    //! Abs
    /*!
     * Replace all the elements in the vectorContainer with their abs.
     */
    void abs();

    /*!
     * Compute the abs of the VectorContainer and return it in a new container.
     * @param vectorContainer - The output vectorContainer containing the abs of the vector.
     */
    void abs( VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! Norm2
    /*!
     * Compute the weight Norm2 of the vector
     */
    Real weightNorm2();

    //! push_back
    /*!
     * Concatenate a VectorContainer to another
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& push_back( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! push_back
    /*!
     * Add a vector at the end of the VectorContainer
     * @param vector_ptr - Shared pointer to a vector
     */
    VectorContainer& push_back( const vectorPtr_Type& vector_ptr );

    //! push_front
    /*!
     * Concatenate a VectorContainer to another inserting it at the beginning
     * @param vectorContainer - VectorContainer
     */
    VectorContainer& push_front( const VectorContainer< vector_Type, container_Type >& vectorContainer );

    //! push_front
    /*!
     * Add a vector at the begin of the VectorContainer
     * @param vector_ptr - Shared pointer to a vector
     */
    VectorContainer& push_front( const vectorPtr_Type& vector_ptr );

    //! Replace
    /*!
     * @param vector_ptr - Shared pointer to the new vector
     * @param id - id of the vector that has to be replaced
     */
    void replace( const vectorPtr_Type& vector_ptr, const UInt& id ) { M_container[id] = vector_ptr; }

    //! resize
    /*!
     * @param size - New size of the container
     */
    void resize( const UInt& size ) { M_container.resize( size ); }

    //! Clear
    void clear() { M_container.clear(); }

    //! ShowMe
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Get Methods
    //@{

    //! Get the number of elements contained inside the VectorContainer object
    UInt size() const;

    //! Get the number of vectors contained inside the VectorContainer object
    UInt vectorsNumber() const
    {
        return static_cast< UInt > ( M_container.size() );
    }

    //@}

private:

    container_Type M_container;
};

// ===================================================
// Constructors
// ===================================================
template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >::VectorContainer():
        M_container()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::VectorContainer()" << "\n";
#endif

}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >::VectorContainer( const VectorContainer& vectorContainer ) :
        M_container()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::VectorContainer( vectorContainer )" << "\n";
#endif

    *this = vectorContainer;
}

// ===================================================
// Operators
// ===================================================
template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator=( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator=( vector_ptr )" << "\n";
#endif

    if ( this != &vectorContainer )
    {
        M_container.resize( vectorContainer.vectorsNumber() );

        vectorPtr_Type myVectorCopy;
        for ( UInt i( 0 ); i < vectorContainer.vectorsNumber(); ++i )
        {
            myVectorCopy.reset( new vector_Type( *( vectorContainer.M_container[i] ) ) );
            this->operator()( i ) = myVectorCopy;
        }
    }
    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator=( scalar )" << "\n";
#endif

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator=( scalar );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator&=( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator&=( vector_ptr )" << "\n";
#endif

    if ( this != &vectorContainer )
    {
        M_container = vectorContainer.M_container;
    }

    return *this;
}

template< class VectorType, class ContainerType >
Real&
VectorContainer< VectorType, ContainerType >::operator[]( const UInt& id ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator[]( id )" << "\n";
#endif

    UInt k( 0 );
    constIterator_Type i = M_container.begin();
    for ( ; i < M_container.end(); ++i )
    {
        if ( id <= k + ( *i )->size() - 1 )
            break;
        else
            k += static_cast< UInt > ( ( *i )->size() );
    }

    return ( *( *i ) )[id - k];
}

template< class VectorType, class ContainerType >
typename VectorContainer< VectorType, ContainerType >::vectorPtr_Type&
VectorContainer< VectorType, ContainerType >::operator()( const UInt& id )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator()( id )" << "\n";
#endif

    return M_container[id];
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator+=( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator+=( vectorContainer )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorContainer.vectorsNumber(); ++i )
        M_container[i]->operator+=( *( vectorContainer.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator-=( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator-=( vectorContainer )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorContainer.vectorsNumber(); ++i )
        M_container[i]->operator-=( *( vectorContainer.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator*=( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator*=( vectorContainer )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        M_container[i]->operator*=( *( vectorContainer.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator/=( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator/=( vectorContainer )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        M_container[i]->operator/=( *( vectorContainer.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator+( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator+( vectorContainer )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy += vectorContainer;

    return myVectorCopy;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator-( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator-( vectorContainer )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy -= vectorContainer;

    return myVectorCopy;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator*( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator*( vectorContainer )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy *= vectorContainer;

    return myVectorCopy;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator/( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator/( vectorContainer )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy /= vectorContainer;

    return myVectorCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator+=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator+=( scalar )" << "\n";
#endif

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator+=( scalar );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator-=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator-=( scalar )" << "\n";
#endif

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator+=( -scalar );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator*=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator*=( scalar )" << "\n";
#endif

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator*=( scalar );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::operator/=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator/=( scalar )" << "\n";
#endif

    this->operator*=( 1. / scalar );

    return *this;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator+( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator+( scalar )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy += scalar;

    return myVectorCopy;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator-( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator-( scalar )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy -= scalar;

    return myVectorCopy;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator*( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator*( scalar )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy *= scalar;

    return myVectorCopy;
}

template< class VectorType, class ContainerType >
const VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator/( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator/( scalar )" << "\n";
#endif

    VectorContainer myVectorCopy = *this;

    myVectorCopy /= scalar;

    return myVectorCopy;
}

// Multiplication by a scalar with the scalar on the left.
template< class ScalarType, class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
operator*( const ScalarType& scalar,
           const VectorContainer< VectorType, ContainerType >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator*( scalar, vectorContainer )" << "\n";
#endif

    VectorContainer< VectorType, ContainerType > vectorContainerCopy( vectorContainer );

    return vectorContainerCopy *= scalar;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator==( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator==( scalar )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy( *this );
    vectorPtr_Type myVectorCopy;

    for ( iterator_Type i = vectorContainerCopy.M_container.begin(); i < vectorContainerCopy.M_container.end(); ++i )
    {
        myVectorCopy.reset( new vector_Type( ( *i )->operator==( scalar ) ) );
        ( *i ) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator!=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator!=( scalar )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy( *this );
    vectorPtr_Type myVectorCopy;

    for ( iterator_Type i = vectorContainerCopy.M_container.begin(); i < vectorContainerCopy.M_container.end(); ++i )
    {
        myVectorCopy.reset( new vector_Type( ( *i )->operator!=( scalar ) ) );
        ( *i ) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator>( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator>( scalar )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy( *this );
    vectorPtr_Type myVectorCopy;

    for ( iterator_Type i = vectorContainerCopy.M_container.begin(); i < vectorContainerCopy.M_container.end(); ++i )
    {
        myVectorCopy.reset( new vector_Type( ( *i )->operator>( scalar ) ) );
        ( *i ) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator<( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator<( scalar )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy( *this );
    vectorPtr_Type myVectorCopy;

    for ( iterator_Type i = vectorContainerCopy.M_container.begin(); i < vectorContainerCopy.M_container.end(); ++i )
    {
        myVectorCopy.reset( new vector_Type( ( *i )->operator<( scalar ) ) );
        ( *i ) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator>=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator>=( scalar )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy( *this );
    vectorPtr_Type myVectorCopy;

    for ( iterator_Type i = vectorContainerCopy.M_container.begin(); i < vectorContainerCopy.M_container.end(); ++i )
    {
        myVectorCopy.reset( new vector_Type( ( *i )->operator>=( scalar ) ) );
        ( *i ) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator<=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator<=( scalar )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy( *this );
    vectorPtr_Type myVectorCopy;

    for ( iterator_Type i = vectorContainerCopy.M_container.begin(); i < vectorContainerCopy.M_container.end(); ++i )
    {
        myVectorCopy.reset( new vector_Type( ( *i )->operator<=( scalar ) ) );
        ( *i ) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator&&( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator&&( vectorContainer )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy;
    vectorPtr_Type myVectorCopy;

    vectorContainerCopy.resize( vectorContainer.vectorsNumber() );
    for ( UInt i( 0 ); i < vectorContainer.vectorsNumber(); ++i )
    {
        myVectorCopy.reset( new vector_Type( M_container[i]->operator&&( *( vectorContainer.M_container[i] ) ) ) );
        vectorContainerCopy.operator()(i) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator||( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator||( vectorContainer )" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy;
    vectorPtr_Type myVectorCopy;

    vectorContainerCopy.resize( vectorContainer.vectorsNumber() );
    for ( UInt i( 0 ); i < vectorContainer.vectorsNumber(); ++i )
    {
        myVectorCopy.reset( new vector_Type( M_container[i]->operator||( *( vectorContainer.M_container[i] ) ) ) );
        vectorContainerCopy.operator()(i) = myVectorCopy;
    }

    return vectorContainerCopy;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >
VectorContainer< VectorType, ContainerType >::operator!()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::operator!()" << "\n";
#endif

    VectorContainer< vector_Type, container_Type > vectorContainerCopy;
    vectorPtr_Type myVectorCopy;

    vectorContainerCopy.resize( vectorsNumber() );
    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
    {
        myVectorCopy.reset( new vector_Type( M_container[i]->operator!() ) );
        vectorContainerCopy.operator()(i) = myVectorCopy;
    }

    return vectorContainerCopy;
}

// ===================================================
// Methods
// ===================================================
template< class VectorType, class ContainerType >
Real
VectorContainer< VectorType, ContainerType >::dot( const VectorContainer< vector_Type, container_Type >& vectorContainer ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::dot( vectorContainer )" << "\n";
#endif

    Real scalarProduct = 0;

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        scalarProduct += M_container[i]->dot( *( vectorContainer.M_container[i] ) );

    return scalarProduct;
}

template< class VectorType, class ContainerType >
void
VectorContainer< VectorType, ContainerType >::dot( const VectorContainer< vector_Type, container_Type >& vectorContainer, Real& scalarProduct )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::dot( vectorContainer, scalarProduct )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        scalarProduct += M_container[i]->Dot( *( vectorContainer.M_container[i] ) );
}

template< class VectorType, class ContainerType >
void
VectorContainer< VectorType, ContainerType >::abs()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::abs()" << "\n";
#endif

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->abs();
}

template< class VectorType, class ContainerType >
void
VectorContainer< VectorType, ContainerType >::abs( VectorContainer< vector_Type,
                                                      container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::abs( vectorContainer )" << "\n";
#endif

    vectorContainer = *this;

    vectorContainer.abs();
}

template< class VectorType, class ContainerType >
Real
VectorContainer< VectorType, ContainerType >::weightNorm2()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::weightNorm2()" << "\n";
#endif

    Real PartialNorm, TotalNorm = 0;

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
    {
        ( *i )->norm2( PartialNorm );
        TotalNorm += PartialNorm * ( *i )->size();
    }

    return TotalNorm / this->size();
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::push_back( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::push_back( vector_ptr )" << "\n";
#endif

    for ( constIterator_Type i = vectorContainer.M_container.begin(); i < vectorContainer.M_container.end(); ++i )
        this->push_back( *i );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::push_back( const vectorPtr_Type& vector_ptr )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::push_back( vector_ptr )" << "\n";
#endif

    M_container.push_back( vector_ptr );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::push_front( const VectorContainer< vector_Type, container_Type >& vectorContainer )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::push_front( vector_ptr )" << "\n";
#endif

    for ( constIterator_Type i = vectorContainer.M_container.begin(); i < vectorContainer.M_container.end(); ++i )
        this->push_front( *i );

    return *this;
}

template< class VectorType, class ContainerType >
VectorContainer< VectorType, ContainerType >&
VectorContainer< VectorType, ContainerType >::push_front( const vectorPtr_Type& vector_ptr )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::push_front( vector_ptr )" << "\n";
#endif

    M_container.push_front( vector_ptr );

    return *this;
}

template< class VectorType, class ContainerType >
void
VectorContainer< VectorType, ContainerType >::showMe( std::ostream& output ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::showMe()" << "\n";
#endif

    output << "Number of vectors:  " << vectorsNumber() << std::endl;
    output << "Global vector size: " << size() << std::endl;

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->showMe( output );
//    for ( UInt i( 0 ); i < size(); ++i )
//        output << "V[" << i << "] = " << this->operator[]( i ) << std::endl;

    output << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
template< class VectorType, class ContainerType >
UInt
VectorContainer< VectorType, ContainerType >::size() const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "VectorContainer::size()" << "\n";
#endif

    UInt size = 0;

    for ( constIterator_Type i = M_container.begin(); i < M_container.end(); ++i )
        size += static_cast< UInt > ( ( *i )->size() );

    return size;
}

} // Namespace LifeV

#endif /* VectorContainer_H */
