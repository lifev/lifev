//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains a simple \f$ R^3 \f$ point

     @date 06/2011
     @author A. Cervone <ant.cervone@gmail.com>
 */

#ifndef _VECTORSMALL_H_
#define _VECTORSMALL_H_ 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/RNM.hpp>

// LifeV namespace.
namespace LifeV
{
//! class VectorSmall   This class implements a simple \f$ R^n \f$ vector

/*!
  @author A. Cervone <ant.cervone@gmail.com>

  This class implements a simple \f$ R^n \f$ vector.
  <br>
  It allows all kind of geometric operations on the node,
  such as summation, multiplication by scalar, scalar product,
  cross product, norm, etc.
  <br>
  The implementation is oriented to best perform with small (less than 30)
  values of \f$ n \f$.
  <br>
  The interface of the class is designed to stay compatible with the Eigen
  library Matrix class.

*/

template <UInt Dim>
class VectorSmall
{

public:

    //! @name Constructors and destructors
    //@{

    //! Empty constructor (all components are set to zero)
    VectorSmall()
    {
        for ( UInt i = 0; i < Dim; i++ )
        {
            M_coords[ i ] = 0.;
        }
    }

    //! Assignment operator
    VectorSmall<Dim>& operator= ( VectorSmall<Dim> const& vector )
    {
        for ( UInt i = 0; i < Dim; i++ )
        {
            M_coords[ i ] = vector.M_coords[ i ];
        }
        return *this;
    }

    //! Copy constructor
    VectorSmall ( VectorSmall<Dim> const& vector )
    {
        *this = vector;
    }

    //@}

    //! @name Static initializations
    //@{

    //! Constant initialization
    static VectorSmall<Dim> Constant ( Real const& value )
    {
        VectorSmall<Dim> v;
        for ( UInt i = 0; i < Dim; i++ )
        {
            v[ i ] = value;
        }
        return v;
    }

    //! Zero initialization
    static VectorSmall<Dim> Zero ()
    {
        return VectorSmall<Dim>::Constant ( 0. );
    }

    //@}

    //! @name Overloaded operators
    //@{

    //! Operator +=
    VectorSmall<Dim>& operator+= ( VectorSmall<Dim> const& vector )
    {
        for ( UInt i = 0; i < Dim; i++ )
        {
            M_coords[ i ] += vector.M_coords[ i ];
        }
        return *this;
    }

    //! Operator +
    VectorSmall<Dim> operator+ ( VectorSmall<Dim> const& vector ) const
    {
        VectorSmall<Dim> tmp ( *this );
        return tmp += vector;
    }

    //! Operator -=
    VectorSmall<Dim>& operator-= ( VectorSmall<Dim> const& vector )
    {
        for ( UInt i = 0; i < Dim; i++ )
        {
            M_coords[ i ] -= vector.M_coords[ i ];
        }
        return *this;
    }

    //! Operator -
    VectorSmall<Dim> operator- ( VectorSmall<Dim> const& vector ) const
    {
        VectorSmall tmp ( *this );
        return tmp -= vector;
    }

    //! Operator *= (multiplication by scalar)
    VectorSmall<Dim>&   operator*= ( Real const& factor )
    {
        for ( UInt i = 0; i < Dim; i++ )
        {
            M_coords[ i ] *= factor;
        }
        return *this;
    }

    //! Operator /= (division by scalar)
    VectorSmall<Dim>& operator/= ( Real const& factor )
    {
        ASSERT ( factor != 0. , "Division by zero!" );
        *this *= 1. / factor;
        return *this;
    }

    //! Operator / (division by scalar)
    VectorSmall<Dim> operator/ ( Real const& factor ) const
    {
        VectorSmall<Dim> tmp ( *this );
        return tmp /= factor;
    }

    //! Operator []
    Real const& operator[] ( UInt const& i ) const
    {
        ASSERT ( i < Dim, "trying to access an index that exceeds the dimension of the array" );
        return M_coords [ i ];
    }

    //! Operator []
    Real& operator[] ( UInt const& i )
    {
        ASSERT ( i < Dim, "trying to set an index that exceeds the dimension of the array" );
        return M_coords [ i ];
    }

    //! Operator ()
    Real const& operator() ( UInt const& i ) const
    {
        ASSERT ( i < Dim, "trying to access an index that exceeds the dimension of the array" );
        return M_coords [ i ];
    }

    //! Operator ()
    Real& operator() ( UInt const& i )
    {
        ASSERT ( i < Dim, "trying to set an index that exceeds the dimension of the array" );
        return M_coords [ i ];
    }

    //@}

    //! @name Geometric Methods
    //@{

    //! Scalar product
    /*!
    @param vector second operand
    @return scalar product value
    */
    Real dot ( VectorSmall<Dim> const& vector ) const
    {
        Real scalarProduct = 0.;
        for ( UInt i = 0; i < Dim; i++ )
        {
            scalarProduct += M_coords[ i ] * vector.M_coords[ i ];
        }
        return scalarProduct;
    }

    //! \f$ L^2 \f$ norm
    /*!
    @return norm value
    */
    Real norm () const
    {
        return std::sqrt ( this->dot ( *this ) );
    }

    //! Normalize vector
    void normalize ()
    {
        *this /= norm ();
    }

    //! Create the versor associated to this VectorSmall
    /*!
    @return the versor associated to this VectorSmall
    */
    VectorSmall<Dim> normalized ()
    {
        return VectorSmall<Dim> ( ( *this ) / norm () );
    }

    //@}

    //! @name Tools
    //@{

    //! function to get the size of the VectorSmall ( for compatibility with Eigen)
    /*!
    @return the fixed size of the VectorSmall
    */
    static UInt size()
    {
        return Dim;
    }

    //@}

private:

    //! @name Data
    //@{

    //! Data storage
    Real M_coords[ Dim ];

    //@}
};

//! @name External overloaded operators
//@{

//! Operator * (multiplication by scalar on the right)
template <UInt Dim>
inline VectorSmall<Dim> operator* ( VectorSmall<Dim> const& vector, Real const& factor )
{
    VectorSmall<Dim> tmp ( vector );
    return tmp *= factor;
}

//! Operator * (multiplication by scalar on the left)
template <UInt Dim>
inline VectorSmall<Dim> operator* ( Real const& factor, VectorSmall<Dim> const& vector )
{
    VectorSmall<Dim> tmp ( vector );
    return tmp *= factor;
}

//! Operator <<
template <UInt Dim>
inline std::ostream& operator<< ( std::ostream& out , VectorSmall<Dim> const& point )
{
    out << "( ";
    for ( UInt i = 0; i < Dim; i++ )
    {
        out << point[ i ] << " ";
    }
    out << ")";
    return out;
}

//@}

//! @name Conversion free-functions
//@{

//! Conversion of an array (std::vector, KN, ecc.) to a VectorSmall
/*!
@param coords vector of point coordinates with operator[] available
@return the VectorSmall that corresponds to the input
*/
template <UInt Dim, typename Vector>
inline VectorSmall<Dim> castToVectorSmall ( Vector const& coords )
{
    ASSERT ( coords.size() == Dim , "the input vector has the wrong dimension" );
    VectorSmall<Dim> tmp;
    for ( UInt i = 0; i < Dim; i++ )
    {
        tmp[ i ] = coords[ i ];
    }
    return tmp;
}

//@}

//! class VectorSmall<3>   Partial specialization for the 3D case
template <>
class VectorSmall<3>
{

public:

    //! @name Constructors and destructors
    //@{
    //! Empty constructor (all components are set to zero)
    VectorSmall()
    {
        M_coords[ 0 ] = M_coords[ 1 ] = M_coords[ 2 ] = 0.;
    }

    //! Full constructor with all components explicitly initialized
    /*!
    @param x x-component of the point
    @param y y-component of the point
    @param z z-component of the point
    */
    VectorSmall ( Real const& x, Real const& y, Real const& z )
    {
        M_coords[ 0 ] = x;
        M_coords[ 1 ] = y;
        M_coords[ 2 ] = z;
    }

    //! Assignment operator
    VectorSmall<3>& operator= ( VectorSmall<3> const& vector )
    {
        M_coords[ 0 ] = vector.M_coords[ 0 ];
        M_coords[ 1 ] = vector.M_coords[ 1 ];
        M_coords[ 2 ] = vector.M_coords[ 2 ];
        return *this;
    }

    //! Copy constructor
    VectorSmall ( VectorSmall<3> const& vector )
    {
        *this = vector;
    }

    //@}

    //! @name Overloaded operators
    //@{

    //! Operator +=
    VectorSmall<3>& operator+= ( VectorSmall<3> const& vector )
    {
        M_coords[ 0 ] += vector.M_coords[ 0 ];
        M_coords[ 1 ] += vector.M_coords[ 1 ];
        M_coords[ 2 ] += vector.M_coords[ 2 ];
        return *this;
    }

    //! Operator +
    VectorSmall<3> operator+ ( VectorSmall<3> const& vector ) const
    {
        VectorSmall<3> tmp ( *this );
        return tmp += vector;
    }

    //! Operator -=
    VectorSmall<3>& operator-= ( VectorSmall<3> const& vector )
    {
        M_coords[ 0 ] -= vector.M_coords[ 0 ];
        M_coords[ 1 ] -= vector.M_coords[ 1 ];
        M_coords[ 2 ] -= vector.M_coords[ 2 ];
        return *this;
    }

    //! Operator -
    VectorSmall<3> operator- ( VectorSmall<3> const& vector ) const
    {
        VectorSmall<3> tmp ( *this );
        return tmp -= vector;
    }

    //! Operator *= (multiplication by scalar)
    VectorSmall<3>&   operator*= ( Real const& factor )
    {
        M_coords[ 0 ] *= factor;
        M_coords[ 1 ] *= factor;
        M_coords[ 2 ] *= factor;
        return *this;
    }

    //! Operator * (multiplication by scalar on the right)
    VectorSmall<3> operator* ( Real const& factor ) const
    {
        VectorSmall<3> tmp ( *this );
        return tmp *= factor;
    }

    //! Operator /= (division by scalar)
    VectorSmall<3>& operator/= ( Real const& factor )
    {
        ASSERT ( factor != 0. , "Division by zero!" );
        *this *= 1. / factor;
        return *this;
    }

    //! Operator / (division by scalar)
    VectorSmall<3> operator/ ( Real const& factor ) const
    {
        VectorSmall<3> tmp ( *this );
        return tmp /= factor;
    }

    //! Operator []
    Real const& operator[] ( UInt const& i ) const
    {
        ASSERT ( i < 3 , "trying to access an index different from 0,1,2" );
        return M_coords [ i ];
    }

    //! Operator []
    Real& operator[] ( UInt const& i )
    {
        ASSERT ( i < 3 , "trying to set an index different from 0,1,2" );
        return M_coords [ i ];
    }

    //! Operator ()
    Real const& operator() ( UInt const& i ) const
    {
        ASSERT ( i < 3 , "trying to access an index different from 0,1,2" );
        return M_coords [ i ];
    }

    //! Operator ()
    Real& operator() ( UInt const& i )
    {
        ASSERT ( i < 3 , "trying to set an index different from 0,1,2" );
        return M_coords [ i ];
    }

    //@}

    //! @name Geometric Methods
    //@{

    //! Scalar product
    /*!
    @param vector second operand
    @return scalar product value
    */
    Real dot ( VectorSmall<3> const& vector ) const
    {
        return ( M_coords[ 0 ] * vector.M_coords[ 0 ]
                 + M_coords[ 1 ] * vector.M_coords[ 1 ]
                 + M_coords[ 2 ] * vector.M_coords[ 2 ] );
    }

    //! Cross product
    /*!
    @param vector second operand
    */
    VectorSmall<3> cross ( VectorSmall<3> const& vector ) const
    {
        return VectorSmall ( M_coords[ 1 ] * vector.M_coords[ 2 ]
                             - M_coords[ 2 ] * vector.M_coords[ 1 ],
                             M_coords[ 2 ] * vector.M_coords[ 0 ]
                             - M_coords[ 0 ] * vector.M_coords[ 2 ],
                             M_coords[ 0 ] * vector.M_coords[ 1 ]
                             - M_coords[ 1 ] * vector.M_coords[ 0 ] );
    }

    //! \f$ L^2 \f$ norm
    /*!
    @return norm value
    */
    Real norm () const
    {
        return std::sqrt ( this->dot ( *this ) );
    }

    //! Normalize vector
    void normalize ()
    {
        *this /= norm ();
    }

    //! Create the versor associated to this VectorSmall
    /*!
    @return the versor associated to this VectorSmall
    */
    VectorSmall<3> normalized ()
    {
        return VectorSmall<3> ( ( *this ) / norm () );
    }

    //@}

    //! @name Tools
    //@{

    //! function to get the size of the VectorSmall ( for compatibility with Eigen)
    /*!
    @return the fixed size of the VectorSmall
    */
    UInt size() const
    {
        return 3;
    }

    //@}

private:

    //! @name Data
    //@{

    //! Data storage
    Real M_coords[3];

    //@}
};

//! @name External overloaded operators
//@{

//! Operator * (multiplication by scalar on the left)
inline VectorSmall<3> operator* ( Real const& factor, VectorSmall<3> const& vector )
{
    VectorSmall<3> tmp ( vector );
    tmp *= factor;
    return tmp;
}

//@}

//! @name Conversion free-functions
//@{

//! Conversion of an array (std::vector, KNM, ecc.) to a VectorSmall
/*!
@param coords vector of point coordinates with operator[] available
@return the VectorSmall that corresponds to the input
*/
template <typename Vector>
inline VectorSmall<3> castToVector3D ( Vector const& coords )
{
    ASSERT ( coords.size() == 3 , "the input vector has the wrong dimension" );
    return VectorSmall<3> ( coords[ 0 ], coords[ 1 ], coords[ 2 ] );
}

//@}

typedef VectorSmall<3> Vector3D;

} // namespace LifeV


#endif //_VECTORSMALL_H_

// -*- mode: c++ -*-
