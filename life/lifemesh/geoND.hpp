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
#ifndef _GEOND_HH_
#define _GEOND_HH_

#include <life/lifemesh/geo0D.hpp>
#include <life/lifemesh/meshEntity.hpp>

namespace LifeV
{

//! Geo0D -  Zero dimensional entity.
/*!
    @author Luca Formaggia

	Base class for Multidimensional basis Geometric Entities.

	@warning It has no boundary information, in fact GeoXD boundary items are stored in the corresponding RegionMesh List.

 */
template <typename GEOSHAPE, typename POINTTYPE = Geo0D>
class GeoND :
        public MeshEntity,
        public GEOSHAPE
{
public:

    //! @name Public Types
    //@{

    //! Number of points associated to the entity
    static const UInt numLocalPoints = GEOSHAPE::numPoints;
    //! Number of Vertices associated to the entity
    static const UInt numLocalVertices = GEOSHAPE::numVertices;

    typedef GEOSHAPE geoShape_Type;
    typedef POINTTYPE point_Type;
    typedef GEOSHAPE GeoShape; //to be removed
    typedef POINTTYPE PointType; //to be removed
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    GeoND();
    //! Declares item identity
    /*!
    	@param Identity Element identity
     */
    explicit GeoND( ID identity );

    //! Declares item identity and item local identiy
    /*!
    	@param identity Element identity
    	@param localIdentity Element local identity
     */
    explicit GeoND( ID identity, ID localIdentity );

    //! Copy constructor
    /*!
        @param Element GeoND to be copied
     */
    GeoND( const GeoND<GEOSHAPE, POINTTYPE> & element);

    //! Destructor
    ~GeoND();

    //@}

    //! @name Operators
    //@{

    //! The equivalence operator
    /*!
    	@param Element Equivalent GeoElement0D
        @return Reference to a new GeoElement0D with the same content of GeoElement0D Element
     */
    GeoND & operator=( GeoND const & element );

    //@}

    //! @name Methods
    //@{

    //! Display general information about the content of the class
    /*!
        List of things displayed in the class
        @param verbose If true more information is displayed
        @param c %%%%
     */
    std::ostream & showMe( bool verbose = false, std::ostream & c = std::cout ) const;

    //! Swap Points
    /*!
        @param firstIdentity Identity of the first point to be swapped (numbering starts from 1)
        @param secondIdentity Identity of the second point to be swapped (numbering starts from 1)
        @warning Function to be used only by routines for checking or amending meshes
     */
    void swapPoints( const ID & firstIdentity, const ID & secondIdentity );

    //! Exchange points
    /*!
     	 Exchanges points according to a list of old-to-new local identity numbering
        @param oldToNew New local identity of a point (numbering starts from 1)
        @warning Function to be used only by routines for checking or amending meshes
     */
    void exchangePoints( const ID oldToNew[ GEOSHAPE::numPoints ] );

    //@}

    //! @name Set Methods
    //@{

    //! Inserts a point using point references
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
     */
    void setPoint( ID const identity, PointType const & point );
    //!Inserts a point using point references with forced bound check
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
        @return %%%%
     */
    bool setPointBD( ID const identity, PointType const & point );
    //!Inserts a point using pointers
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
     */
    void setPoint( ID const identity, PointType const * point );
    //!Inserts a point using pointers with forced bound check
    /*!
        @param identity Identity of the point to be inserted
        @param point Point to be inserted
        @return %%%%
     */
    bool setPointBD( ID const identity, PointType const * point );

    //! Sets the flag of a point
    /*!
    	Sets the flag to the stronger between the stored one and the one provided by the argument
    	@param identity Identity of the point to be inserted
        @param point Point to be inserted
        @return %%%%%
    	@warning A const_cast to M_points is done in order to change the flag
    */
    EntityFlag setStrongerMarkerAtPoint( const ID& identity, EntityFlag const & flag );

    //@}

    //! @name Get Methods
    //@{

    //! Returns the point of identity indicated in the argument
    /*!
    	@param identity Identity of the point (numbering starts from 1)
        @return reference to a point object, possibly derived from PointType
    */
    PointType const & point ( ID const identity ) const;
    //! Returns the point of identity indicated in the argument
    /*!
     	It starts from the last point and it follows the rule: vertices first.
     	It may be used to access the points of a Geometry Element in a reverse way
     	(i.e. with the opposite GeoElement orientation)
    		@param identity Identity of the point (numbering starts from 1)
        	@return reference to a point object, possibly derived from PointType
    */
    PointType const & reversepoint ( ID const identity ) const;

    //@}

private:
    PointType const* M_points[ GEOSHAPE::numPoints ];
};


/*--------------------------------------------------------------
                 GeoND
---------------------------------------------------------------*/
template <typename GEOSHAPE, typename POINTTYPE>
const UInt GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints;

template <typename GEOSHAPE, typename POINTTYPE>
const UInt GeoND<GEOSHAPE, POINTTYPE>::numLocalVertices;

template <typename GEOSHAPE, typename POINTTYPE>
GeoND<GEOSHAPE, POINTTYPE>::GeoND() :
        MeshEntity( 0 )
{

}

template <typename GEOSHAPE, typename POINTTYPE>
GeoND<GEOSHAPE, POINTTYPE>::GeoND( ID identity ) :
        MeshEntity( identity, identity )
{

}

template <typename GEOSHAPE, typename POINTTYPE>
GeoND<GEOSHAPE, POINTTYPE>::GeoND( ID identity, ID localIdentity ) :
        MeshEntity( identity, localIdentity )
{

}

template <typename GEOSHAPE, typename POINTTYPE>
GeoND<GEOSHAPE, POINTTYPE>::GeoND( GeoND<GEOSHAPE, POINTTYPE> const & element ) :
        MeshEntity( element.id(), element.localId() )
{
    for ( UInt i = 0; i < GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints; ++i )
    {
        M_points[ i ] = element.M_points[ i ];
    }
}

template <typename GEOSHAPE, typename POINTTYPE>
GeoND<GEOSHAPE, POINTTYPE>::~GeoND()
{

}

template <typename GEOSHAPE, typename POINTTYPE>
GeoND<GEOSHAPE, POINTTYPE> &
GeoND<GEOSHAPE, POINTTYPE>::operator=( GeoND<GEOSHAPE, POINTTYPE> const & element )
{
    if ( this != &element )
    {
        this->setId     (element.id());
        this->setLocalId(element.localId());
        for ( UInt i = 0; i < GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints; ++i )
        {
            M_points[ i ] = element.M_points[ i ];
        }
    }
    return *this;
}

template <typename GEOSHAPE, typename POINTTYPE>
std::ostream & GeoND<GEOSHAPE, POINTTYPE>::
showMe( bool verbose, std::ostream & out ) const
{
    out << "----- GeoND object -----" << std::endl;
    out << " Number of Vertices = " << GEOSHAPE::numVertices << std::endl;
    out << " Number of Points   = " << GEOSHAPE::numPoints << std::endl;
    out << " ID                 = " << id() << std::endl;
    out << " local ID           = " << localId() << std::endl;
    if ( verbose )
    {
        out << " POINTS INFORMATION" << std::endl << std::endl;
        for ( unsigned i = 1 ; i <= GEOSHAPE::numVertices; i++ )
        {
            out << "POINT ID. " << i << std::endl;
            out << point( i ).showMe( verbose, out );
        }
    }
    out << "----- END OF GeoND data ---" << std::endl << std::endl;
    return out;
}

template <typename GEOSHAPE, typename POINTTYPE>
void GeoND<GEOSHAPE, POINTTYPE>::swapPoints( const ID & firstIdentity, const ID & secondIdentity )
{
    PointType const* tmp( M_points[ firstIdentity - 1 ] );
    M_points[ firstIdentity - 1 ] = M_points[ secondIdentity - 1 ];
    M_points[ secondIdentity - 1 ] = tmp;
}

template <typename GEOSHAPE, typename POINTTYPE>
void GeoND<GEOSHAPE, POINTTYPE>::exchangePoints( const ID oldToNew[ GEOSHAPE::numPoints ] )
{
    PointType const* tmp[ GEOSHAPE::numPoints ];
    for ( UInt i = 0; i < GEOSHAPE::numPoints; ++i )
    {
        tmp[ i ] = M_points[ i ];
    }
    for ( UInt i = 0; i < GEOSHAPE::numPoints; ++i )
    {
        M_points[ i ] = tmp[ oldToNew[ i ] - 1 ];
    }
}

template <typename GEOSHAPE, typename POINTTYPE>
INLINE
void GeoND<GEOSHAPE, POINTTYPE>::setPoint( ID const identity, PointType const & point )
{
    ASSERT_BD( ( identity > 0 && identity <= GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints ) ) ;
    M_points[ identity - 1 ] = ( &point );
}

template <typename GEOSHAPE, typename POINTTYPE>
bool GeoND<GEOSHAPE, POINTTYPE>::setPointBD( ID const identity, PointType const & point )
{
    ASSERT_BD0( ( identity > 0 && identity <= GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints ) ) ;
    if ( identity <= 0 || identity > GeoND<GEOSHAPE, POINTTYPE>::numLocalVertices )
        return false;
    M_points[ identity -1 ] = ( &point );
    return true;
}

template <typename GEOSHAPE, typename POINTTYPE>
INLINE
void GeoND<GEOSHAPE, POINTTYPE>::setPoint( ID const identity, PointType const * point )
{
    ASSERT_BD( ( identity > 0 && identity <= GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints ) ) ;
    M_points[ identity - 1 ] = ( point );
}

template <typename GEOSHAPE, typename POINTTYPE>
bool GeoND<GEOSHAPE, POINTTYPE>::setPointBD( ID const identity, PointType const * point )
{
    ASSERT_BD0( ( identity > 0 && identity <= GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints ) ) ;
    if ( identity <= 0 || identity > GeoND<GEOSHAPE, POINTTYPE>::numLocalVertices )
        return false;
    M_points[ identity -1 ] = ( point );
    return true;
}

template <typename GEOSHAPE, typename POINTTYPE>
EntityFlag GeoND<GEOSHAPE, POINTTYPE>::setStrongerMarkerAtPoint( const ID& identity, EntityFlag const & flag )
{
    return (const_cast<POINTTYPE *> ( M_points[identity -1]) ) -> setStrongerMarker(flag);
}

template <typename GEOSHAPE, typename POINTTYPE>
INLINE
POINTTYPE const & GeoND<GEOSHAPE, POINTTYPE>::point( ID const identity ) const
{
    ASSERT_BD( ( identity > 0 && identity <= GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints ) );
    return *( static_cast<POINTTYPE const*>( M_points[ identity - 1 ] ) );
}

template <typename GEOSHAPE, typename POINTTYPE>
INLINE
POINTTYPE const & GeoND<GEOSHAPE, POINTTYPE>::reversepoint( ID const identity ) const
{
    ASSERT_BD( ( identity > 0 && identity <= GeoND<GEOSHAPE, POINTTYPE>::numLocalPoints ) );
    return *( static_cast<POINTTYPE const*>( M_points[ reversePoint<GEOSHAPE>::operate( identity ) - 1 ] ) );
}

}
#endif
