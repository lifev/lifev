//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 * @file
 * @brief Containers Of Vectors
 *
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 * @date 29-09-2009
 */

#ifndef ContainerOfVectors_H
#define ContainerOfVectors_H 1

#include <life/lifecore/life.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV {

//! ContainerOfVectors - LifeV vector made of other vectors
/*!
 *  @author Cristiano Malossi
 *
 *  The ContainerOfVectors class has the purpose to contains a certain number of
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
 *      <li> The ContainerOfVectors contains shared_ptr to vectors: if you perform an operation (such as the
 *           element by element multiplication) among two containers which hold the same vectors (or at least one),
 *           the results can be different from what expected!
 *      <li> Up to now it has been strongly tested only with LifeV::EpetraVector. However the design should be
 *           more or less compatible with boost::numeric::ublas::vector.
 *  </ol>
 *
 */
template< class VectorType, class ContainerType = std::vector< boost::shared_ptr< VectorType > > >
class ContainerOfVectors
{
public:

    typedef typename ContainerType::iterator IT;
    typedef typename ContainerType::const_iterator constIT;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    ContainerOfVectors();

    //! Copy constructor
    /*!
     * Makes a copy of all the objects contained in the shared_ptr(s).
     * @param ContainerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors( const ContainerOfVectors& containerOfVectors );

    //! Destructor
    ~ContainerOfVectors() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * Make a true copy of all the objects contained in the shared_ptr(s)
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& operator=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Operator=
    /*!
     * Set all the components of the vector equal to a scalar quantity.
     * Cannot be done for an empty vector!
     * @param scalar - Scalar value
     */
    ContainerOfVectors& operator=( const Real& scalar );

    //! Operator&=
    /*!
     * Make a copy of the shared_ptr(s), without duplicates the objects contained in them.
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& operator&=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Operator[]
    /*!
     * The element access operator (by reference) has been implemented just for debugging or other special cases,
     * because it is not efficient with this type of implementation.
     *
     * @param ID - Element ID
     */
    Real& operator[]( const UInt& ID ) const;

    //! Operator()
    /*!
     * Vector access operator (by shared_ptr).
     *
     * @param ID - Vector ID
     */
    boost::shared_ptr< VectorType >& operator()( const UInt& ID );

    //! Operator+=
    /*!
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& operator+=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Operator-=
    /*!
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& operator-=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Operator*=
    /*!
     * Multiplication element by element of two vectors.
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& operator*=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Operator/=
    /*!
     * Division element by element of two vectors.
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& operator/=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Operator+
    /*!
     * @param containerOfVectors - ContainerOfVectors
     */
    const ContainerOfVectors operator+( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const;

    //! Operator-
    /*!
     * @param containerOfVectors - ContainerOfVectors
     */
    const ContainerOfVectors operator-( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const;

    //! Operator*
    /*!
     * Multiplication element by element of two vectors.
     * @param containerOfVectors - ContainerOfVectors
     */
    const ContainerOfVectors operator*( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const;

    //! Operator/
    /*!
     * Division element by element of two vectors.
     * @param containerOfVectors - ContainerOfVectors
     */
    const ContainerOfVectors operator/( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const;

    //! Operator+=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    ContainerOfVectors& operator+=( const Real& scalar );

    //! Operator-=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    ContainerOfVectors& operator-=( const Real& scalar );

    //! Operator*=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    ContainerOfVectors& operator*=( const Real& scalar );

    //! Operator/=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    ContainerOfVectors& operator/=( const Real& scalar );

    //! Operator*=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const ContainerOfVectors operator+( const Real& scalar ) const;

    //! Operator/=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const ContainerOfVectors operator-( const Real& scalar ) const;

    //! Operator*=
    /*!
     * Multiplication of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const ContainerOfVectors operator*( const Real& scalar ) const;

    //! Operator/=
    /*!
     * Division of the vector by a scalar.
     * @param scalar - Scalar value
     */
    const ContainerOfVectors operator/( const Real& scalar ) const;

    //! operator==
    /*!
     * Return a ContainerOfVectors containing 1 where vector elements are == scalar;
     * @param scalar - Value for the comparison.
     */
    ContainerOfVectors operator==( const Real& scalar );

    //! operator!=
    /*!
     * Return a ContainerOfVectors containing 1 where vector elements are != scalar;
     * @param scalar - Value for the comparison.
     */
    ContainerOfVectors operator!=( const Real& scalar );

    //! operator<
    /*!
     * Return a ContainerOfVectors containing 1 where vector elements are < scalar;
     * @param scalar - Value for the comparison.
     */
    ContainerOfVectors operator<( const Real& scalar );

    //! operator>
    /*!
     * Return a ContainerOfVectors containing 1 where vector elements are > scalar;
     * @param scalar - Value for the comparison.
     */
    ContainerOfVectors operator>( const Real& scalar );

    //! operator<=
    /*!
     * Return a ContainerOfVectors containing 1 where vector elements are <= scalar;
     * @param scalar - Value for the comparison.
     */
    ContainerOfVectors operator<=( const Real& scalar );

    //! operator>=
    /*!
     * Return a ContainerOfVectors containing 1 where vector elements are >= scalar;
     * @param scalar - Value for the comparison.
     */
    ContainerOfVectors operator>=( const Real& scalar );

    //! Logic operator&&
    /*!
     * Return a vector containing one where both the elements are different from zero;
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors operator&&( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Logic operator||
    /*!
     * Return a vector containing one if one of the two elements is different from zero;
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors operator||( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Logic operator!
    /*!
     * Return a vector containing one where the elements are equal to zero;
     */
    ContainerOfVectors operator!();

    //@}


    //! @name Methods
    //@{

    //! Dot - Scalar product
    /*!
     * Scalar product of the vectors
     * @param containerOfVectors - ContainerOfVectors
     */
    Real Dot( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const;

    //! Dot - Scalar product
    /*!
     * Scalar product of the vectors
     * @param containerOfVectors - ContainerOfVectors
     * @param scalarProduct - result
     */
    void Dot( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors, Real& scalarProduct );

    //! Abs
    /*!
     * Replace all the elements in the containerOfVectors with their abs.
     */
    void Abs();

    /*!
     * Compute the abs of the ContainerOfVectors and return it in a new container.
     * @param containerOfVectors - The output containerOfVectors containing the abs of the vector.
     */
    void Abs( ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! Norm2
    /*!
     * Compute the weight Norm2 of the vector
     */
    Real WeightNorm2();

    //! push_back
    /*!
     * Concatenate a ContainerOfVectors to another
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& push_back( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! push_back
    /*!
     * Add a vector at the end of the ContainerOfVectors
     * @param vector_ptr - Shared pointer to a vector
     */
    ContainerOfVectors& push_back( const boost::shared_ptr< VectorType >& vector_ptr );

    //! push_front
    /*!
     * Concatenate a ContainerOfVectors to another inserting it at the beginning
     * @param containerOfVectors - ContainerOfVectors
     */
    ContainerOfVectors& push_front( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors );

    //! push_front
    /*!
     * Add a vector at the begin of the ContainerOfVectors
     * @param vector_ptr - Shared pointer to a vector
     */
    ContainerOfVectors& push_front( const boost::shared_ptr< VectorType >& vector_ptr );

    //! Replace
    /*!
     * @param vector_ptr - Shared pointer to the new vector
     * @param ID - ID of the vector that has to be replaced
     */
    void Replace( const boost::shared_ptr< VectorType >& vector_ptr, const UInt& ID );

    //! resize
    /*!
     * @param size - New size of the container
     */
    void resize( const UInt& size );

    //! Clear
    void clear();

    //! ShowMe
    void ShowMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Get Methods
    //@{

    //! Get the number of elements contained inside the ContainerOfVectors object
    UInt size() const;

    //! Get the number of vectors contained inside the ContainerOfVectors object
    UInt vectorsNumber() const
    {
        return static_cast< UInt > ( M_container.size() );
    }

    //@}

private:

    ContainerType M_container;
};

// ===================================================
// Constructors
// ===================================================
template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >::ContainerOfVectors():
    M_container()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::ContainerOfVectors()" << "\n";
#endif

}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >::ContainerOfVectors( const ContainerOfVectors& containerOfVectors ) :
    M_container()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::ContainerOfVectors( containerOfVectors )" << "\n";
#endif

    *this = containerOfVectors;
}

// ===================================================
// Operators
// ===================================================
template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator=( vector_ptr )" << "\n";
#endif

    if ( this != &containerOfVectors )
    {
        M_container.resize( containerOfVectors.vectorsNumber() );

        boost::shared_ptr< VectorType > MyVectorCopy;
        for ( UInt i( 0 ); i < containerOfVectors.vectorsNumber(); ++i )
        {
            MyVectorCopy.reset( new VectorType( *( containerOfVectors.M_container[i] ) ) );
            this->operator()( i ) = MyVectorCopy;
        }
    }
    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator=( scalar )" << "\n";
#endif

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator=( scalar );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator&=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator&=( vector_ptr )" << "\n";
#endif

    if ( this != &containerOfVectors )
    {
        M_container = containerOfVectors.M_container;
    }

    return *this;
}

template< class VectorType, class ContainerType >
Real&
ContainerOfVectors< VectorType, ContainerType >::operator[]( const UInt& ID ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator[]( ID )" << "\n";
#endif

    UInt k( 0 );
    constIT i = M_container.begin();
    for ( ; i < M_container.end(); ++i )
    {
        if ( ID <= k + ( *i )->size() - 1 )
            break;
        else
            k += static_cast< UInt > ( ( *i )->size() );
    }

    return ( *( *i ) )[ID - k];
}

template< class VectorType, class ContainerType >
boost::shared_ptr< VectorType >&
ContainerOfVectors< VectorType, ContainerType >::operator()( const UInt& ID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator()( ID )" << "\n";
#endif

    return M_container[ID];
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator+=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator+=( containerOfVectors )" << "\n";
#endif

    for ( UInt i( 0 ); i < containerOfVectors.vectorsNumber(); ++i )
        M_container[i]->operator+=( *( containerOfVectors.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator-=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator-=( containerOfVectors )" << "\n";
#endif

    for ( UInt i( 0 ); i < containerOfVectors.vectorsNumber(); ++i )
        M_container[i]->operator-=( *( containerOfVectors.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator*=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator*=( containerOfVectors )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        M_container[i]->operator*=( *( containerOfVectors.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator/=( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator/=( containerOfVectors )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        M_container[i]->operator/=( *( containerOfVectors.M_container[i] ) );

    return *this;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator+( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator+( containerOfVectors )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy += containerOfVectors;

    return MyVectorCopy;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator-( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator-( containerOfVectors )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy -= containerOfVectors;

    return MyVectorCopy;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator*( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator*( containerOfVectors )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy *= containerOfVectors;

    return MyVectorCopy;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator/( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator/( containerOfVectors )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy /= containerOfVectors;

    return MyVectorCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator+=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator+=( scalar )" << "\n";
#endif

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator+=( scalar );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator-=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator-=( scalar )" << "\n";
#endif

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator+=( -scalar );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator*=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator*=( scalar )" << "\n";
#endif

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->operator*=( scalar );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::operator/=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator/=( scalar )" << "\n";
#endif

    this->operator*=( 1. / scalar );

    return *this;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator+( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator+( scalar )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy += scalar;

    return MyVectorCopy;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator-( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator-( scalar )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy -= scalar;

    return MyVectorCopy;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator*( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator*( scalar )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy *= scalar;

    return MyVectorCopy;
}

template< class VectorType, class ContainerType >
const ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator/( const Real& scalar ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator/( scalar )" << "\n";
#endif

    ContainerOfVectors MyVectorCopy = *this;

    MyVectorCopy /= scalar;

    return MyVectorCopy;
}

// Multiplication by a scalar with the scalar on the left.
template< class ScalarType, class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
operator*( const ScalarType& scalar,
           const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator*( scalar, containerOfVectors )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy( containerOfVectors );

    return containerOfVectorsCopy *= scalar;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator==( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator==( scalar )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy( *this );
    boost::shared_ptr< VectorType > MyVectorCopy;

    for ( IT i = containerOfVectorsCopy.M_container.begin(); i < containerOfVectorsCopy.M_container.end(); ++i )
    {
        MyVectorCopy.reset( new VectorType( ( *i )->operator==( scalar ) ) );
        ( *i ) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator!=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator!=( scalar )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy( *this );
    boost::shared_ptr< VectorType > MyVectorCopy;

    for ( IT i = containerOfVectorsCopy.M_container.begin(); i < containerOfVectorsCopy.M_container.end(); ++i )
    {
        MyVectorCopy.reset( new VectorType( ( *i )->operator!=( scalar ) ) );
        ( *i ) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator>( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator>( scalar )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy( *this );
    boost::shared_ptr< VectorType > MyVectorCopy;

    for ( IT i = containerOfVectorsCopy.M_container.begin(); i < containerOfVectorsCopy.M_container.end(); ++i )
    {
        MyVectorCopy.reset( new VectorType( ( *i )->operator>( scalar ) ) );
        ( *i ) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator<( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator<( scalar )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy( *this );
    boost::shared_ptr< VectorType > MyVectorCopy;

    for ( IT i = containerOfVectorsCopy.M_container.begin(); i < containerOfVectorsCopy.M_container.end(); ++i )
    {
        MyVectorCopy.reset( new VectorType( ( *i )->operator<( scalar ) ) );
        ( *i ) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator>=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator>=( scalar )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy( *this );
    boost::shared_ptr< VectorType > MyVectorCopy;

    for ( IT i = containerOfVectorsCopy.M_container.begin(); i < containerOfVectorsCopy.M_container.end(); ++i )
    {
        MyVectorCopy.reset( new VectorType( ( *i )->operator>=( scalar ) ) );
        ( *i ) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator<=( const Real& scalar )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator<=( scalar )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy( *this );
    boost::shared_ptr< VectorType > MyVectorCopy;

    for ( IT i = containerOfVectorsCopy.M_container.begin(); i < containerOfVectorsCopy.M_container.end(); ++i )
    {
        MyVectorCopy.reset( new VectorType( ( *i )->operator<=( scalar ) ) );
        ( *i ) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator&&( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator&&( containerOfVectors )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy;
    boost::shared_ptr< VectorType > MyVectorCopy;

    containerOfVectorsCopy.resize( containerOfVectors.vectorsNumber() );
    for ( UInt i( 0 ); i < containerOfVectors.vectorsNumber(); ++i )
    {
        MyVectorCopy.reset( new VectorType( M_container[i]->operator&&( *( containerOfVectors.M_container[i] ) ) ) );
        containerOfVectorsCopy.operator()(i) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator||( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator||( containerOfVectors )" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy;
    boost::shared_ptr< VectorType > MyVectorCopy;

    containerOfVectorsCopy.resize( containerOfVectors.vectorsNumber() );
    for ( UInt i( 0 ); i < containerOfVectors.vectorsNumber(); ++i )
    {
        MyVectorCopy.reset( new VectorType( M_container[i]->operator||( *( containerOfVectors.M_container[i] ) ) ) );
        containerOfVectorsCopy.operator()(i) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >
ContainerOfVectors< VectorType, ContainerType >::operator!()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::operator!()" << "\n";
#endif

    ContainerOfVectors< VectorType, ContainerType > containerOfVectorsCopy;
    boost::shared_ptr< VectorType > MyVectorCopy;

    containerOfVectorsCopy.resize( vectorsNumber() );
    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
    {
        MyVectorCopy.reset( new VectorType( M_container[i]->operator!() ) );
        containerOfVectorsCopy.operator()(i) = MyVectorCopy;
    }

    return containerOfVectorsCopy;
}

// ===================================================
// Methods
// ===================================================
template< class VectorType, class ContainerType >
Real
ContainerOfVectors< VectorType, ContainerType >::Dot( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::Dot( containerOfVectors )" << "\n";
#endif

    Real scalarProduct = 0;

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        scalarProduct += M_container[i]->Dot( *( containerOfVectors.M_container[i] ) );

    return scalarProduct;
}

template< class VectorType, class ContainerType >
void
ContainerOfVectors< VectorType, ContainerType >::Dot( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors,
                                                            Real& scalarProduct )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::Dot( containerOfVectors, scalarProduct )" << "\n";
#endif

    for ( UInt i( 0 ); i < vectorsNumber(); ++i )
        scalarProduct += M_container[i]->Dot( *( containerOfVectors.M_container[i] ) );
}

template< class VectorType, class ContainerType >
void
ContainerOfVectors< VectorType, ContainerType >::Abs()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::Abs()" << "\n";
#endif

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->Abs();
}

template< class VectorType, class ContainerType >
void
ContainerOfVectors< VectorType, ContainerType >::Abs( ContainerOfVectors< VectorType,
        ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::Abs( containerOfVectors )" << "\n";
#endif

    containerOfVectors = *this;

    containerOfVectors.Abs();
}

template< class VectorType, class ContainerType >
Real
ContainerOfVectors< VectorType, ContainerType >::WeightNorm2()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::WeightNorm2()" << "\n";
#endif

    Real PartialNorm, TotalNorm = 0;

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
    {
        ( *i )->Norm2( PartialNorm );
        TotalNorm += PartialNorm * ( *i )->size();
    }

    return TotalNorm / this->size();
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::push_back( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::push_back( vector_ptr )" << "\n";
#endif

    for ( constIT i = containerOfVectors.M_container.begin(); i < containerOfVectors.M_container.end(); ++i )
        this->push_back( *i );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::push_back( const boost::shared_ptr< VectorType >& vector_ptr )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::push_back( vector_ptr )" << "\n";
#endif

    M_container.push_back( vector_ptr );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::push_front( const ContainerOfVectors< VectorType, ContainerType >& containerOfVectors )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::push_front( vector_ptr )" << "\n";
#endif

    for ( constIT i = containerOfVectors.M_container.begin(); i < containerOfVectors.M_container.end(); ++i )
        this->push_front( *i );

    return *this;
}

template< class VectorType, class ContainerType >
ContainerOfVectors< VectorType, ContainerType >&
ContainerOfVectors< VectorType, ContainerType >::push_front( const boost::shared_ptr< VectorType >& vector_ptr )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::push_front( vector_ptr )" << "\n";
#endif

    M_container.push_front( vector_ptr );

    return *this;
}

template< class VectorType, class ContainerType >
void
ContainerOfVectors< VectorType, ContainerType >::Replace( const boost::shared_ptr< VectorType >& vector_ptr,
                                                          const UInt& ID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::Replace( vector_ptr, ID )" << "\n";
#endif

    M_container[ID] = vector_ptr;
}

template< class VectorType, class ContainerType >
void
ContainerOfVectors< VectorType, ContainerType >::resize( const UInt& size )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::resize( size )" << "\n";
#endif

    M_container.resize( size );
}

template< class VectorType, class ContainerType >
void
ContainerOfVectors< VectorType, ContainerType >::clear()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::clear()" << "\n";
#endif

    M_container.clear();
}

template< class VectorType, class ContainerType >
void
ContainerOfVectors< VectorType, ContainerType >::ShowMe( std::ostream& output ) const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::ShowMe()" << "\n";
#endif

    output << "Number of vectors:  " << vectorsNumber() << std::endl;
    output << "Global vector size: " << size() << std::endl;

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
        ( *i )->ShowMe( output );
//    for ( UInt i( 0 ); i < size(); ++i )
//        output << "V[" << i << "] = " << this->operator[]( i ) << std::endl;

    output << std::endl;
}

// ===================================================
// Get Methods
// ===================================================
template< class VectorType, class ContainerType >
UInt
ContainerOfVectors< VectorType, ContainerType >::size() const
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 3100 ) << "ContainerOfVectors::size()" << "\n";
#endif

    UInt size = 0;

    for ( constIT i = M_container.begin(); i < M_container.end(); ++i )
        size += static_cast< UInt > ( ( *i )->size() );

    return size;
}

} // Namespace LifeV

#endif /* ContainerOfVectors_H */
