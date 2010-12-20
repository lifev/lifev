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
    @brief Zero dimensional entity

    @author Luca Formaggia <luca.formaggia@polimi.it>
    @contributor Marta D'Elia <mdelia2@mathcs.emory.edu>
    @maintainer Marta D'Elia <mdelia2@mathcs.emory.edu>

    @date 00-00-0000

*/


#ifndef _GEO0D_HH_
#define _GEO0D_HH_

#include <boost/array.hpp>
#include <life/lifemesh/meshEntity.hpp>
#include <life/lifemesh/basisElSh.hpp>

namespace LifeV
{
//! Geo0D -  Zero dimensional entity.
/*!
    @author
	Intermediate class used to build the actual Geometry classes; it stores boundary information.

	@warning Geo1D/2D/3D are template classes; in fact, information might not be known a priori.
	All vector dimensions are determined at compile time to enhance memory access time.
	A coherent GeoShape has to be provided by the user.

 */
class Geo0D : public MeshEntityWithBoundary
{
public:

    //! @name Public Types
    //@{

    typedef GeoPoint geoShape_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    Geo0D();

    //! Declares item identity and states if it is on boundary
    /*!
    	@param identity Element identity
        @param boundary True if the element is on boundary
     */
    explicit Geo0D( ID identity, bool boundary = false );

    //! Declares item identity, provides coordinate and states if it is on boundary
    /*!
    	@param identity Element identity
    	@param x Element x coordinate
    	@param y Element y coordinate
    	@param z Element z coordinate
    	@param boundary True if the element is on boundary
     */
    Geo0D( ID identity, Real x, Real y, Real z, bool boundary = false );

    //! Copy constructor
    /*!
        @param Element Geo0D to be copied
     */
    Geo0D( Geo0D const & Element );

    //! Destructor
    virtual ~Geo0D()
    {
        // nothing to be done
    }

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        @param Element Equivalent GeoElement0D
        @return Reference to a new GeoElement0D with the same content of GeoElement0D Element
     */
    Geo0D & operator=( Geo0D const & Element );

    //@}

    //! @name Methods
    //@{

    //! Display general information about the content of the class
    /*!
        List of things displayed in the class
        @param output specify the output format (std::cout by default)
     */
    std::ostream & showMe( bool Verbose = false, std::ostream & coordinateVector = std::cout ) const;

    //! Returns the pointer to the coordinates vector
    /*!
    	@return Pointer to coordinate vector
     */
    Real const * coordinatesArray() const
    {
        return M_coordinates.data();
    };
    //! Returns the pointer to the coordinates vector
    /*!
    	@return Pointer to coordinate vector
     */
//    Real const * __attribute__ ((__deprecated__)) coor() const
//    {
//        return coordinatesArray();
//    };

    //! Returns the reference to the x-coordinate
    /*!
    	Used to provide coordinates to object created using a constructor with no coordinates given, or to modify existing coordinates
    	@return Reference to element x-coordinate
     */
    Real & x()
    {
        return M_coordinates[ 0 ];
    }
    //! Returns the reference to the y-coordinate
    /*!
      	Used to provide coordinates to object created using a constructor with no coordinates given, or to modify existing coordinates
    	@return Reference to element y-coordinate
     */
    Real & y()
    {
        return M_coordinates[ 1 ];
    }
    //! Returns the reference to the z-coordinate and checks if working in two dimensions
    /*!
    	Used to provide coordinates to object created using a constructor with no coordinates given, or to modify existing coordinates
    	@return Reference to element z-coordinate
     */
    Real & z()
    {
#ifdef TWODIM
        ERROR_MSG( "z coordinate may be modified only in a 3D problem" );
#else
        return M_coordinates[ 2 ];
#endif

    }
    //! Returns the x-coordinate
    /*!
    	@return Element x-coordinate
     */
    Real x() const
    {
        return M_coordinates[ 0 ];
    }
    //! Returns the y-coordinate
    /*!
    	@return Element y-coordinate
     */
    Real y() const
    {
        return M_coordinates[ 1 ];
    };
    //! Returns the z-coordinate and checks if working in two dimensions
    /*!
    	@return Element z-coordinate
     */
    Real z() const
    {
#ifdef TWODIM
        return 0;
#else
        return M_coordinates[ 2 ];
#endif

    }

    //! Returns the coordinate specified in the argument
    /*!
        The method allows to access the coordinate specified in the argument
    	@param coordinate x, y, or z coordinate to be returned
    	@return Coordinate specified in the argument
     */
    Real coordinate ( ID const coordinate ) const
    {
        ASSERT_BD( coordinate > 0 && coordinate <= NDIM ) ;
        return M_coordinates[ coordinate -1 ]; // indexing from 1
    }
    //! Returns the reference to the coordinate specified in the argument
    /*!
        The method allows to modify the coordinate specified in the argument
    	@param coordinate x, y, or z coordinate to be returned
    	@return Reference to the coordinate specified in the argument
     */
    Real & coordinate ( ID const coordinate )
    {
        ASSERT_BD( coordinate > 0 && coordinate <= NDIM ) ;
        return M_coordinates[ coordinate -1 ];
    }

    //@}

    //! @name Get Methods
    //@{

    //! Returns the coordinates vector
    /*!
        The method allows to access coordinates and modify them
    	@return Coordinates array
     */
    boost::array<Real,NDIM>& coordinates ( void )
    {
        return M_coordinates;
    }
    //! Returns the coordinates vector
    /*!
        The method allows to access coordinates and modify them
    	@return Coordinates array
     */
//    boost::array<Real,NDIM>& __attribute__ ((__deprecated__)) coordinate ( void )
//    {
//        return coordinates();
//    }

    //@}

private:
    boost::array<Real,NDIM> M_coordinates;
};

}
#endif
