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


#ifndef MESHVERTEX_H
#define MESHVERTEX_H

#include <lifev/core/array/VectorSmall.hpp>
#include <lifev/core/mesh/MeshEntity.hpp>
#include <lifev/core/mesh/ElementShapes.hpp>

namespace LifeV
{
//! MeshVertex -  Zero dimensional entity.
/*!
    @author
    Intermediate class used to build the actual Geometry classes; it stores boundary information.

    @warning MeshVertex is a template class; in fact, information might not be known a priori.
    All vector dimensions are determined at compile time to enhance memory access time.
    A coherent GeoShape has to be provided by the user.

 */
class MeshVertex : public MeshEntity
{
public:

    //! @name Public Types
    //@{

    typedef GeoPoint geoShape_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MeshVertex();

    //! Declares item identity and states if it is on boundary
    /*!
        Local and global id are set equal. To change global id
        use setId() method (inherited from MeshEntity)
        @param identity Element identity (local and global)
        @param boundary True if the element is on boundary
     */
    explicit MeshVertex ( ID identity, bool boundary = false );

    //! Declares item identity, provides coordinate and states if it is on boundary
    /*!
        @param identity Element identity
        @param x Element x coordinate
        @param y Element y coordinate
        @param z Element z coordinate
        @param boundary True if the element is on boundary
     */
    MeshVertex ( ID identity, Real x, Real y, Real z, bool boundary = false );

    //! Destructor
    virtual ~MeshVertex()
    {
        // nothing to be done
    }

    //@}

    //! @name Operators
    //@{
    //@}

    //! @name Methods
    //@{

    //! Display general information about the content of the class
    /*!
        List of things displayed in the class
        @param output specify the output format (std::cout by default)
     */
    std::ostream& showMe ( bool Verbose = false, std::ostream& coordinateVector = std::cout ) const;

    //! Returns the pointer to the coordinates vector
    /*!
        @return Pointer to coordinate vector
     */
    Real const* coordinatesArray() const
    {
        return &M_coordinates[0];
    };

    //! Returns the reference to the x-coordinate
    /*!
        Used to provide coordinates to object created using a constructor with no coordinates given, or to modify existing coordinates
        @return Reference to element x-coordinate
     */
    Real& x()
    {
        return M_coordinates[ 0 ];
    }
    //! Returns the reference to the y-coordinate
    /*!
        Used to provide coordinates to object created using a constructor with no coordinates given, or to modify existing coordinates
        @return Reference to element y-coordinate
     */
    Real& y()
    {
        return M_coordinates[ 1 ];
    }
    //! Returns the reference to the z-coordinate and checks if working in two dimensions
    /*!
        Used to provide coordinates to object created using a constructor with no coordinates given, or to modify existing coordinates
        @return Reference to element z-coordinate
     */
    Real& z()
    {
        return M_coordinates[ 2 ];
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
        return M_coordinates[ 2 ];
    }

    //! Returns the coordinate specified in the argument
    /*!
        The method allows to access the coordinate specified in the argument
        @param coordinate x, y, or z coordinate to be returned
        @return Coordinate specified in the argument
     */
    Real coordinate ( ID const coordinate ) const
    {
        ASSERT_BD ( coordinate < NDIM ) ;
        return M_coordinates[ coordinate ];
    }
    //! Returns the reference to the coordinate specified in the argument
    /*!
        The method allows to modify the coordinate specified in the argument
        @param coordinate x, y, or z coordinate to be returned
        @return Reference to the coordinate specified in the argument
     */
    Real& coordinate ( ID const coordinate )
    {
        ASSERT_BD ( coordinate < NDIM ) ;
        return M_coordinates[ coordinate ];
    }

    //@}

    //! @name Get Methods
    //@{

    //! Returns the coordinates vector
    /*!
        The method allows to access coordinates and modify them
        @return Coordinates array
    */
    Vector3D const& coordinates () const
    {
        return M_coordinates;
    }

    //@}

private:
    Vector3D M_coordinates;
};

}
#endif
