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

#ifndef _BAREPOINT_H_
#define _BAREPOINT_H_ 1

#include "life/lifecore/LifeV.hpp"
#include "life/lifemesh/MeshVertex.hpp"
#include "life/lifearray/RNM.hpp"

#include <vector>
#include <cmath>

// LifeV namespace.
namespace LifeV
{
//! class BarePoint   This class implements a simple \f$ R^3 \f$ vector

/*!
  @author A. Cervone <ant.cervone@gmail.com>

  This class implements a simple R^3 vector.
  <br>
  It allows all kind of geometric operations on the node,
  such as summation, multiplication by scalar, scalar product,
  cross product, norm, etc.

*/

class BarePoint
{

public:

    //! @name Constructors and destructors
    //@{

    //! Empty constructor (all components are set to zero)
    BarePoint()
    {
        M_coords[ 0 ] = M_coords[ 1 ] = M_coords[ 2 ] = 0.;
    }

    //! Full contructor with all components explicitly initialized
    /*!
    @param x x-component of the point
    @param y y-component of the point
    @param z z-component of the point
    */
    BarePoint( Real const & x, Real const & y, Real const & z )
    {
        M_coords[ 0 ] = x;
        M_coords[ 1 ] = y;
        M_coords[ 2 ] = z;
    }

    //! Assignment operator
    BarePoint & operator= ( BarePoint const & vector )
    {
        this->M_coords[ 0 ] = vector.M_coords[ 0 ];
        this->M_coords[ 1 ] = vector.M_coords[ 1 ];
        this->M_coords[ 2 ] = vector.M_coords[ 2 ];
        return *this;
    }

    //! Copy constructor
    BarePoint( BarePoint const & vector )
    {
        *this = vector;
    }

    //@}

    //! @name Overloaded operators 
    //@{

    //! Operator +=
    BarePoint & operator+= ( BarePoint const & vector )
    {
        this->M_coords[ 0 ] += vector.M_coords[ 0 ];
        this->M_coords[ 1 ] += vector.M_coords[ 1 ];
        this->M_coords[ 2 ] += vector.M_coords[ 2 ];
        return *this;
    }

    //! Operator +
    BarePoint operator+ ( BarePoint const & vector )
    {
        BarePoint tmp ( *this ); return tmp += vector;
    }

    //! Operator -=
    BarePoint & operator-= ( BarePoint const & vector )
    {
        this->M_coords[ 0 ] -= vector.M_coords[ 0 ];
        this->M_coords[ 1 ] -= vector.M_coords[ 1 ];
        this->M_coords[ 2 ] -= vector.M_coords[ 2 ];
        return *this;
    }

    //! Operator -
    BarePoint operator- ( BarePoint const & vector )
    {
        BarePoint tmp ( *this ); return tmp -= vector;
    }

    //! Operator *= (multiplication by scalar)
    BarePoint &  operator*= ( Real const & factor )
    {
        this->M_coords[ 0 ] *= factor;
        this->M_coords[ 1 ] *= factor;
        this->M_coords[ 2 ] *= factor;
        return *this;
    }

    //! Operator * (multiplication by scalar on the right)
    BarePoint operator* ( Real const & factor )
    {
        BarePoint tmp ( *this ); return tmp *= factor;
    }

    //! Operator /= (division by scalar)
    BarePoint & operator/= ( Real const & factor )
    {
        ASSERT ( factor != 0. , "Division by zero!" );
        *this *= 1. / factor;
        return *this;
    }

    //! Operator / (division by scalar)
    BarePoint operator/ ( Real const & factor )
    {
        BarePoint tmp ( *this ); return tmp /= factor;
    }

    //! Operator []
    Real const & operator[] ( UInt const & i ) const
    {
        ASSERT ( i < 3 , "trying to set an index different from 0,1,2" );
        return M_coords [ i ];
    }

    //! Operator []
    Real & operator[] ( UInt const & i )
    {
        ASSERT ( i < 3 , "trying to set an index different from 0,1,2" );
        return M_coords [ i ];
    }

    //! Operator ()
    Real const & operator() ( UInt const & i ) const
    {
        ASSERT ( i < 3 , "trying to set an index different from 0,1,2" );
        return M_coords [ i ];
    }

    //! Operator ()
    Real & operator() ( UInt const & i )
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
    Real dot ( BarePoint const & vector ) const
    {
        return ( this->M_coords[ 0 ] * vector.M_coords[ 0 ]
               + this->M_coords[ 1 ] * vector.M_coords[ 1 ]
               + this->M_coords[ 2 ] * vector.M_coords[ 2 ] );
    }

    //! Cross product
    /*!
    @param vector second operand
    @return a BarePoint with the cross product result
    */
    BarePoint cross ( BarePoint const & vector ) const
    {
        return BarePoint ( this->M_coords[ 1 ] * vector.M_coords[ 2 ]
                         - this->M_coords[ 2 ] * vector.M_coords[ 1 ],
                           this->M_coords[ 2 ] * vector.M_coords[ 0 ]
                         - this->M_coords[ 0 ] * vector.M_coords[ 2 ],
                           this->M_coords[ 0 ] * vector.M_coords[ 1 ]
                         - this->M_coords[ 1 ] * vector.M_coords[ 0 ] );
    }

    //! \f$ L^2 \f$ norm
    /*!
    @return norm value
    */
    Real norm () const
    {
        return std::sqrt( M_coords[ 0 ] * M_coords[ 0 ] +
                          M_coords[ 1 ] * M_coords[ 1 ] +
                          M_coords[ 2 ] * M_coords[ 2 ] );
    }

    //! Normalize vector
    void normalize ()
    {
        *this /= norm ();
    }

    //! Create the versor associated to this BarePoint
    /*!
    @return the versor associated to this BarePoint
    */
    BarePoint normalized ()
    {
        return BarePoint ( ( *this ) / norm () );
    }

    //@}

    //! @name Tools
    //@{

    //! function to get the size of the barePoint ( for compatibility with Eigen)
    /*!
    @return the fixed size of the BarePoint
    */
    const UInt size() const { return 3;}

    //! Operator <<
    friend std::ostream & operator<< ( std::ostream & out , BarePoint const & point );

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
inline BarePoint operator* ( Real const & factor, BarePoint const & vector )
{
    BarePoint tmp ( vector ); return tmp *= factor;
}

//! Operator <<
inline std::ostream & operator<< ( std::ostream & out , BarePoint const & point )
{
    out << "(" << point.M_coords[ 0 ] << ", " 
               << point.M_coords[ 1 ] << ", " 
               << point.M_coords[ 2 ] << ")";
    return out;
}

//@}

//! @name Conversion free-functions
//@{

//! Conversion of an STL-vector of coordinates to a BarePoint

/*!
@param coords vector of point coordinates
@return the BarePoint that corresponds to the input
*/
inline BarePoint castToBarePoint ( std::vector<Real> const & coords )
{
    ASSERT ( coords.size() == 3 , "the inpunt vector is of the wrong dimension" );
    return BarePoint( coords[ 0 ], coords[ 1 ], coords[ 2 ] );
}

//! Conversion of a MeshVertex to a BarePoint

/*!
@param vertex MeshVertex original object to be copied
@return the BarePoint that corresponds to the input
*/
inline BarePoint castToBarePoint ( MeshVertex const & vertex )
{
    return BarePoint( vertex.x(), vertex.y(), vertex.z() );
}

//! Conversion of a KN<Real> vector to a BarePoint

/*!
@param vector vector of point coordinates
@return the BarePoint that corresponds to the input
*/
inline BarePoint castToBarePoint ( KN<Real> const & vector )
{
    ASSERT ( vector.size() == 3 , "the inpunt vector is of the wrong dimension" );
    return BarePoint( vector[ 0 ], vector[ 1 ], vector[ 2 ]);
}

} // namespace LifeV


#endif //_BAREPOINT_H_

// -*- mode: c++ -*-
