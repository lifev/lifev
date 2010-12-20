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
    @file
    @brief Definition of the quadPoint class, usefull for the quadrature rules.

    @author Jean-Frederic Gerbeau
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 18-05-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef QUADPOINT_H
#define QUADPOINT_H 1

#include <life/lifecore/life.hpp>

#include <boost/array.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <iostream>

namespace LifeV
{
//! @name Public typedefs
//@{
typedef boost::numeric::ublas::vector<Real> GeoVector;
//@}

//! QuadPoint - Simple container for a point of a quadrature rule.
/*!
    @author Samuel Quinodoz
    @date 05/2010
    @version 2.0
    @note This class is based on a previous implementation, due to J.-F. Gerbeau (04/2002).

    <b> Definition </b>

    The QuadPoint class consists basically in a real number (the weight of the point) and a vector of real numbers (the coordinates of the point). To enable fast computations if needed, blas vectors (look in the file /lifearray/tab.hpp to know precisely the type) are used internally.

    <b> Create a QuadPoint </b>

    A QuadPoint can be defined using directly 1,2 or 3 coordinates and the weight.

    \code
    QuadPoint myPoint(1,0,1,0.5); // Point(1,0,1) with weight 0.5
    \endcode

    However, this use should be avoided if possible, because it assumes (for backward compatibility purpose) that the dimension of the QuadPoint is 3, even if only 1 or 2 coordinates are passed!

    \code
    QuadPoint myPoint(1,0.3); // Point(1,0,0) with weight 0.3
    \endcode

    To create properly a QuadPoint, one should use the native vector (GeoVector):

    \code
    GeoVector myCoordinates(1); // Create a vector with 1 component
    myCoordinates[0]=1;         // Put 1 in the first component
    QuadPoint myPoint(myCoordinates,0.3);   // Point(1) with weight 0.3
    \endcode

    This makes the code "surprise-free" and also easier to generalize to any dimension.

    <b> Dimension of the QuadPoint </b>

    The QuadPoint has naturally a dimension, the number of coordinates of the point. To change the dimension of the QuadPoint, use the copy constructor where you can specify the new dimension.

    If you intend to use your code for different dimensions, take care of the methods that you call: avoid using x(),y(),z() and replace them by the coor() method.


 */
class QuadPoint
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty constructor (all zero data).
    /*!
      This constructor builds a QuadPoint with
      3 components (for the location) and a weight
      that are all set to zero
     */
    QuadPoint();

    //! Full constructor for 3D
    /*!
      This builds a quadrature with 3D coordinates.
      @param x First coordinate of the point
      @param y Second coordinate of the point
      @param z Third coordinate of the point
      @param weight Weight of the point
     */
    QuadPoint( Real x, Real y, Real z, Real weight );

    //! Full constructor for 2D
    /*!
      @param x First coordinate of the point
      @param y Second coordinate of the point
      @param weight Weight of the point
     */
    QuadPoint( Real x, Real y, Real weight );

    //! Full constructor for 1D
    /*!
      @param x First coordinate of the point
      @param weight Weight of the point
     */
    QuadPoint( Real x, Real weight );

    //! Full multidimension constructor
    /*!
      @param coor Coordinates of the point
      @param weight Weight of the point
     */
    QuadPoint(const GeoVector& coor, const Real& weight);

    //! Multidimension constructor with specified dimension
    /*!
      With this constructor, one can specify the dimension
      of the QuadPoint (that can therefore be different from
      coor.size()).

      @param coor Coordinates of the point
      @param weight Weight of the point
      @param spaceDim The new dimension of the point
     */
    QuadPoint(const GeoVector& coor, const Real& weight, const UInt& spaceDim);


    //! Simple copy constructor
    /*!
      @param qp The quadrature point to copy
     */
    QuadPoint(const QuadPoint& qp);

    //! Import from another dimension
    /*!
      With this constructor, one can change the dimension of the
      space where the quadrature point is living. For example, we can see a
      2D quadrature point as a 3D quadrature point with 0 for the third component.
      @param qp The quadrature point to import
      @param spaceDim The dimension of the space where the quadrature point is defined
     */
    QuadPoint(const QuadPoint& qp, const UInt spaceDim);

    //! Destructor
    virtual ~QuadPoint() {};

    //@}


    //! @name Methods
    //@{

    //! Returns the dimension of the quadPoint
    Real dimension() const
    {
        return M_coor.size();
    }

    //@}


    //! @name Operator
    //@{

    //! Output operator
    friend std::ostream & operator << ( std::ostream& out, const QuadPoint& qp )
    {
        out << " Coordinates = ";
        for (UInt i(0); i<qp.M_coor.size(); ++i)
        {
            out << qp.M_coor[i] << ", ";
        };
        out << " Weight = " << qp.M_weight;
        return out;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the weight
    const Real& weight() const
    {
        return M_weight;
    }

    //! Getter for the first coordinate
    const Real& x() const
    {
        ASSERT(0<M_coor.size()," No x coordinate for this quadrature point");
        return M_coor[ 0 ];
    }

    //! Getter for the second coordinate
    const Real& y() const
    {
        ASSERT(1<M_coor.size()," No y coordinate for this quadrature point");
        return M_coor[ 1 ];
    }

    //! Getter for the third coordinate
    const Real& z() const
    {
        ASSERT(2<M_coor.size()," No z coordinate for this quadrature point");
        return M_coor[ 2 ];
    }

    //! Getter for the coordinate (0<=i)
    const Real& coor(const UInt& i ) const
    {
        ASSERT(i<M_coor.size()," Error in the coordinate for this quadrature point");
        return M_coor[ i ];
    }

    //! Getter for the full vector of coordinates
    const GeoVector& coor() const
    {
        return M_coor;
    }

    //@}

private:
    Real M_weight;
    GeoVector M_coor;
};


} // Namespace LifeV

#endif /* QUADPOINT_H */
