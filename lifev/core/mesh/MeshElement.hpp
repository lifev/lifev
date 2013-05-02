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
#ifndef MESHELEMENT_H
#define MESHELEMENT_H

#include <lifev/core/mesh/MeshVertex.hpp>
#include <lifev/core/mesh/MeshEntity.hpp>
#include <lifev/core/array/VectorSmall.hpp>
#include <algorithm>

namespace LifeV
{

//! MeshVertex -  Zero dimensional entity.
/*!
    @author Luca Formaggia

    Base class for Multidimensional basis Geometric Entities.

    @warning It has no boundary information, in fact GeoXD boundary items are stored in the corresponding RegionMesh List.

 */
template <typename GeoShape, typename PointType = MeshVertex>
class MeshElement :
    public MeshEntity,
    public GeoShape
{
public:

    //! @name Public Types
    //@{

    //! Number of points associated to the entity
    static const UInt S_numLocalPoints = GeoShape::S_numPoints;
    //! Number of Vertices associated to the entity
    static const UInt S_numLocalVertices = GeoShape::S_numVertices;

    typedef GeoShape geoShape_Type;
    typedef PointType point_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MeshElement();
    //! Declares item identity
    /*!
        @param Identity Element identity
     */
    explicit MeshElement ( ID identity );

    //! Copy constructor
    /*!
        @param Element MeshElement to be copied
     */
    MeshElement ( const MeshElement<GeoShape, PointType>& element);

    //! Destructor
    ~MeshElement();

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
        @param Element Equivalent MeshElementMarked0D
        @return Reference to a new MeshElementMarked0D with the same content of MeshElementMarked0D Element
     */
    MeshElement& operator= ( MeshElement const& element );

    //@}

    //! @name Methods
    //@{

    //! Display general information about the content of the class
    /*!
        List of things displayed in the class
        @param verbose If true more information is displayed
        @param c Output
     */
    std::ostream& showMe ( bool verbose = false, std::ostream& c = std::cout ) const;

    //! Swap Points
    /*!
        @param firstIdentity Identity of the first point to be swapped
        @param secondIdentity Identity of the second point to be swapped
        @warning Function to be used only by routines for checking or amending meshes
     */
    void swapPoints ( const ID& firstIdentity, const ID& secondIdentity );

    //! Reverses the points of the element.
    /*!
        Useful in the phase of building a mesh if you want to change the orientation of
        the element. First point becomes last etc.etc

        \note Beware that the some adjacency list may be invalidated. If you revert the point
        of a 3D element, for instance, the local faces change and a elementToFace list may be
        now incorrect. So this method must be used with care and only when building the basic entities
        of a mesh.
     */
    void reversePoints();
    //! Exchange points
    /*!
         Exchanges points according to a list of old-to-new local identity numbering
        @param oldToNew New local identity of a point
        @warning Function to be used only by routines for checking or amending meshes
     */
    void exchangePoints ( const ID oldToNew[ GeoShape::S_numPoints ] );

    //! Returns the point of identity indicated in the argument
    /*!
        @param identity Identity of the point
        @return reference to a point object, possibly derived from point_Type
    */
    point_Type const& point ( ID const identity ) const;
    //! Returns the point of identity indicated in the argument
    /*!
        It starts from the last point and it follows the rule: vertices first.
        It may be used to access the points of a Geometry Element in a reverse way
        (i.e. with the opposite MeshElementMarked orientation)
            @param identity Identity of the point
            @return reference to a point object, possibly derived from point_Type
    */
    point_Type const& reversepoint ( ID const identity ) const;

    //@}

    //! @name Set Methods
    //@{

    //! Inserts a point using point references
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
     */
    void setPoint ( ID const identity, point_Type const& point );
    //!Inserts a point using point references with forced boundary check
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
        @return TRUE if the point is set
     */
    bool setPointWithBoundaryCheck ( ID const identity, point_Type const& point );

    //!Inserts a point using pointers
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
     */
    void setPoint ( ID const identity, point_Type const* point );
    //!Inserts a point using pointers with forced boundary check
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
        @return %%%%
     */
    bool setPointWithBoundaryCheck ( ID const identity, point_Type const* point );

    //! Sets the Marker ID of a point
    /*!
        Sets the Marker ID the stronger between the stored one and the one provided by the argument
        @param identity Elemental numbering of the point to be inserted
        @param point Point to be inserted
        @return TRUE if the point is set
        @warning A const_cast to M_points is done
    */
    markerID_Type setStrongerMarkerIDAtPoint ( const ID& identity, markerID_Type const& flag ) const;

    //@}

    //! @name Get Methods
    //@{

    //! Returns the points vector
    /*!
        The method allows to access coordinates but not modify them
        @return Points vector
     */
    point_Type const* points () const
    {
        return M_points;
    }

    //@}

private:
    point_Type const* M_points[ GeoShape::S_numPoints ];
};


/*--------------------------------------------------------------
                 MeshElement
---------------------------------------------------------------*/
template <typename GeoShape, typename PointType>
const UInt MeshElement<GeoShape, PointType>::S_numLocalPoints;

template <typename GeoShape, typename PointType>
const UInt MeshElement<GeoShape, PointType>::S_numLocalVertices;

template <typename GeoShape, typename PointType>
MeshElement<GeoShape, PointType>::MeshElement() :
    MeshEntity ( NotAnId )
{

}

template <typename GeoShape, typename PointType>
MeshElement<GeoShape, PointType>::MeshElement ( ID identity ) :
    MeshEntity ( identity )
{

}

template <typename GeoShape, typename PointType>
MeshElement<GeoShape, PointType>::MeshElement ( MeshElement<GeoShape, PointType> const& element ) :
    MeshEntity ( element )
{
    setLocalId ( element.localId() );
    for ( UInt i = 0; i < MeshElement<GeoShape, PointType>::S_numLocalPoints; ++i )
    {
        M_points[ i ] = element.M_points[ i ];
    }
}

template <typename GeoShape, typename PointType>
MeshElement<GeoShape, PointType>::~MeshElement()
{

}

template <typename GeoShape, typename PointType>
MeshElement<GeoShape, PointType>&
MeshElement<GeoShape, PointType>::operator= ( MeshElement<GeoShape, PointType> const& element )
{
    if ( this != &element )
    {
        MeshEntity::operator= ( element );
        for ( UInt i = 0; i < MeshElement<GeoShape, PointType>::S_numLocalPoints; ++i )
        {
            M_points[ i ] = element.M_points[ i ];
        }
    }
    return *this;
}

template <typename GeoShape, typename PointType>
std::ostream& MeshElement<GeoShape, PointType>::
showMe ( bool verbose, std::ostream& out ) const
{
    out << "----- MeshElement object -----" << std::endl;
    out << " Number of Vertices = " << GeoShape::S_numVertices << std::endl;
    out << " Number of Points   = " << GeoShape::S_numPoints << std::endl;
    out << " ID                 = " << id() << std::endl;
    out << " local ID           = " << localId() << std::endl;
    if ( verbose )
    {
        out << " POINTS INFORMATION" << std::endl << std::endl;
        for ( UInt i = 0 ; i < GeoShape::S_numVertices; i++ )
        {
            out << "POINT ID. " << i << std::endl;
            out << point ( i ).showMe ( verbose, out );
        }
    }
    out << "----- END OF MeshElement data ---" << std::endl << std::endl;
    return out;
}

template <typename GeoShape, typename PointType>
void MeshElement<GeoShape, PointType>::swapPoints ( const ID& firstIdentity, const ID& secondIdentity )
{
    std::swap (M_points[ firstIdentity ], M_points[ secondIdentity ]);
}

template <typename GeoShape, typename PointType>
void MeshElement<GeoShape, PointType>::reversePoints()
{
    static bool first (true);
    static ID oldToNew[ GeoShape::S_numPoints ];
    if (first)
    {
        for (ID i = 0; i < GeoShape::S_numPoints; ++i)
        {
            oldToNew[i] = reversePoint<GeoShape> (i);
        }
        first = false;
    }
    this->exchangePoints (oldToNew);
}


template <typename GeoShape, typename PointType>
void MeshElement<GeoShape, PointType>::exchangePoints ( const ID oldToNew[ GeoShape::S_numPoints ] )
{
    point_Type const* tmp[ GeoShape::S_numPoints ];
    for ( UInt i = 0; i < GeoShape::S_numPoints; ++i )
    {
        tmp[ i ] = M_points[ i ];
    }
    for ( UInt i = 0; i < GeoShape::S_numPoints; ++i )
    {
        M_points[ i ] = tmp[ oldToNew[ i ] ];
    }
}

template <typename GeoShape, typename PointType>
inline
PointType const& MeshElement<GeoShape, PointType>::point ( ID const identity ) const
{
    ASSERT_BD ( ( identity < MeshElement<GeoShape, PointType>::S_numLocalPoints ) );
    return * ( static_cast<PointType const*> ( M_points[ identity ] ) );
}

template <typename GeoShape, typename PointType>
inline
PointType const& MeshElement<GeoShape, PointType>::reversepoint ( ID const identity ) const
{
    ASSERT_BD ( ( identity < MeshElement<GeoShape, PointType>::S_numLocalPoints ) );
    return * ( static_cast<PointType const*> ( M_points[ reversePoint<GeoShape> ( identity ) ] ) );
}

template <typename GeoShape, typename PointType>
inline
void MeshElement<GeoShape, PointType>::setPoint ( ID const identity, point_Type const& point )
{
    ASSERT_BD ( ( identity < MeshElement<GeoShape, PointType>::S_numLocalPoints ) ) ;
    M_points[ identity ] = ( &point );
}

template <typename GeoShape, typename PointType>
bool MeshElement<GeoShape, PointType>::setPointWithBoundaryCheck ( ID const identity, point_Type const& point )
{
    ASSERT_BD0 ( ( identity < MeshElement<GeoShape, PointType>::S_numLocalPoints ) ) ;
    if ( identity >= MeshElement<GeoShape, PointType>::S_numLocalVertices )
    {
        return false;
    }
    M_points[ identity ] = ( &point );
    return true;
}

template <typename GeoShape, typename PointType>
inline
void MeshElement<GeoShape, PointType>::setPoint ( ID const identity, point_Type const* point )
{
    ASSERT_BD ( ( identity < MeshElement<GeoShape, PointType>::S_numLocalPoints ) ) ;
    M_points[ identity ] = ( point );
}

template <typename GeoShape, typename PointType>
bool MeshElement<GeoShape, PointType>::setPointWithBoundaryCheck ( ID const identity, point_Type const* point )
{
    ASSERT_BD0 ( ( identity < MeshElement<GeoShape, PointType>::S_numLocalPoints ) ) ;
    if ( identity >= MeshElement<GeoShape, PointType>::S_numLocalVertices )
    {
        return false;
    }
    M_points[ identity ] = ( point );
    return true;
}

template <typename GeoShape, typename PointType>
markerID_Type MeshElement<GeoShape, PointType>::setStrongerMarkerIDAtPoint ( const ID& identity, markerID_Type const& flag )
const
{
    return (const_cast<PointType*> ( M_points[identity]) ) -> setStrongerMarkerID (flag);
}

template <typename PointType>
Real edgeLength ( const MeshElement<LinearLine, PointType>& edge )
{
    //    Real deltaX, deltaY, deltaZ;
    //
    //    deltaX = ( edge.point( 1 ) ).x() - ( edge.point( 0 ) ).x();
    //    deltaY = ( edge.point( 1 ) ).y() - ( edge.point( 0 ) ).y();
    //    deltaZ = ( edge.point( 1 ) ).z() - ( edge.point( 0 ) ).z();
    //
    //    return std::sqrt( deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ );

    Vector3D dVec = castToVector3D (edge.point ( 1 ).coordinates() );
    dVec -= castToVector3D (edge.point ( 0 ).coordinates() );

    return dVec.norm();
}

}
#endif

