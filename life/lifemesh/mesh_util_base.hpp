/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __MESH_UTIL_BASE__
#define __MESH_UTIL_BASE__
#include <vector>
#include <algorithm>
#include <set>

#include "lifeV.hpp"
#include "regionMesh3D.hpp"

namespace LifeV
{
/*!
  \brief Base tilities operating on meshes
 
 
  This file contains a set of base utilities used to test mesh entities or
  operate on them
 
 */

//! A locally used structure, not meant for general use
typedef std::map<BareFace, std::pair<ID, ID >, cmpBareItem<BareFace> > TempFaceContainer;

//! A locally used structure, not meant for general use
typedef std::map<BareEdge, std::pair<ID, ID>, cmpBareItem<BareEdge> > TempEdgeContainer;

/*
*************************************************************************************
                            FUNCTORS
*************************************************************************************
*/ 
//! \defgroup Test_Functors Some useful functors to be used for test mesh entities

/*! \ingroup Test_Functors
  \briefFunctor to check if a Point, Face or Edge is on the boundary.
 
  \precond It assumes that boundary points in RegionMesh are correctly set.
  \precond   the RegionMesh must export the typenames
  PointType, FaceType and EdgeType.
*/
template <typename RegionMesh>
class EnquireBEntity
{
public:
    EnquireBEntity( RegionMesh const & mesh ) : pmesh( &mesh )
    {}
    ;
    typedef typename RegionMesh::FaceType FaceType;
    typedef typename RegionMesh::EdgeType EdgeType;
    typedef typename RegionMesh::PointType PointType;

    bool operator() ( FaceType & face )
    {
        bool isboundary = true;
        for ( UInt k = 1;k <= FaceType::numVertices;++k )
        {
            isboundary = isboundary & face.point( k ).boundary();
        }
        return isboundary;
    }

    bool operator() ( EdgeType & edge )
    {
        bool isboundary = true;
        for ( UInt k = 1;k <= EdgeType::numVertices;++k )
        {
            isboundary = isboundary & edge.point( k ).boundary();
        }
        return isboundary;
    }

    INLINE bool operator() ( PointType & point )
    {
        return point.boundary();
    }

private:
    EnquireBEntity()
    {}
    RegionMesh const * pmesh;
};

//!\ingroup Test_Functors
/*! Functor to check if a Face is on the boundary, by using the information
 contained in a TempFaceContainer produced by findBoundaryFaces(). It does
 not use the information contained in the mesh PointList, so it differs
 from EnquireBEntity.
 \precond bfaces have been previously set by a call to FindBoundaryFaces
*/
template <typename RegionMesh>
class EnquireBFace
{
public:
    typedef typename RegionMesh::FaceType FaceType;
    typedef typename RegionMesh::FaceShape FaceShape;

    EnquireBFace( RegionMesh const & mesh, TempFaceContainer const & bfaces ) :
            pmesh( &mesh ), pbfaces( &bfaces )
    {}

    bool operator() ( FaceType & f )
    {
        ID i1, i2, i3, i4;
        BareFace bface;

        i1 = f.point( 1 ).id();
        i2 = f.point( 2 ).id();
        i3 = f.point( 3 ).id();
        if ( FaceShape::numVertices == 4 )
        {
            i4 = f.point( 4 ).id();
            bface = ( makeBareFace( i1, i2, i3, i4 ) ).first;
        }
        else
        {
            bface = ( makeBareFace( i1, i2, i3 ) ).first;
        }
        return pbfaces->find( bface ) != pbfaces->end();
    }
private:
    EnquireBFace()
    {}
    ;
    RegionMesh const * pmesh;
    TempFaceContainer const * pbfaces;
};

//!\ingroup Test_Functors
/*! Functor to check if an edge is on the boundary, by using the information
 contained in a TempFaceContainer produced by findBoundaryEdges(). It does
 not use the information contained in the mesh PointList, so it differs
 from EnquireBEntity.
 
 \precond bedges have been previously set by a call to FindBoundaryEdges()
 
*/
template <typename RegionMesh>
class EnquireBEdge
{
public:
    typedef typename RegionMesh::EdgeType EdgeType;
    typedef typename RegionMesh::EdgeShape EdgeShape;

    EnquireBEdge( RegionMesh const & mesh, TempEdgeContainer const & bedges ) :
            pmesh( &mesh ), pbedges( &bedges )
    {}

    bool operator() ( EdgeType & f )
    {
        ID i1, i2;
        BareEdge bedge;

        i1 = f.point( 1 ).id();
        i2 = f.point( 2 ).id();
        bedge = ( makeBareEdge( i1, i2 ) ).first;
        return pbedges->find( bedge ) != pbedges->end();
    }

private:
    EnquireBEdge()
    {}
    ;
    RegionMesh const * pmesh;
    TempEdgeContainer const * pbedges;
};

//! \ingroup Test_Functors
/*! Functor to check if a mesh entity with boundary indicator (for instance a GeoPoint)
  is on the boundary, by enquiring its boundary flag.
  \warning It assumes that boundary points are correctly set.
*/
template <typename RegionMesh>
class EnquireBPoint
{
public:
    EnquireBPoint( RegionMesh & mesh ) : pmesh( &mesh )
    {}
    ;
    bool operator() ( MeshEntityWithBoundary & e )
    {
        return e.boundary();
    }
private:
    EnquireBPoint()
    {}
    ;
    RegionMesh * pmesh;
};


//! \ingroup Test_Functors
/*! This functor is used to do some geometry checks It returns a coordinate
 */
class GetCoordComponent
{
public:
    GetCoordComponent();
    GetCoordComponent( int i );
    void operator() ( Real const x, Real const y, Real const z, Real ret[ 3 ] ) const;
private:
    int comp;
};

//! \ingroup Test_Functors
/*! This functor is used to do some geometry checks It returns a vector of ones
 */
class GetOnes
{
public:
    void operator() ( Real const x, Real const y, Real const z, Real ret[ 3 ] ) const;
};
/*
*************************************************************************************
                            EDGES/FACES FINDERS
*************************************************************************************
*/

//! Finds boundary faces.
/*!  A low level routine, not meant to be called directly. It creates a
container with all the information needed to set up properly the boundary
faces connectivities.
 
\param mesh A 3D mesh.
 
\param NumInternalFaces. A reference to an integer returning the number of internal faces found.
 
\param bfaces This container will eventually contain a map whose key are
the BareFace corresponding to the boundary faces and the data a pair of
IDs: the ID of the adjacent element and the relative position of the face
in the element.
 
\param allFaces When this bool is set true the function will also construct the set of internale faces, stored in intfaces.
 
\param intfaces A container that will possibly contain a map whose keys are
the BareFace corresponding to an internal faces and the data a pair of IDs:
the ID of the two elements adjacent to the face.
 
\return Number of boundary faces found
*/
template <typename RegionMesh3D>
UInt findFaces( const RegionMesh3D & mesh, TempFaceContainer & bfaces, UInt & numInternalFaces, TempFaceContainer & intfaces, bool allFaces = false )
{
    UInt i1, i2, i3, i4;
    BareFace bface;
    typename RegionMesh3D::VolumeShape ele;
    typedef typename RegionMesh3D::Volumes Volumes;
    TempFaceContainer::iterator fi;

    // clean first in case it has been alredy used

    bfaces.clear();
    if ( allFaces )
        intfaces.clear();
    numInternalFaces = 0;

    for ( typename Volumes::const_iterator iv = mesh.volumeList.begin();
            iv != mesh.volumeList.end(); ++iv )
    {
        for ( ID j = 1;j <= mesh.numLocalFaces();++j )
        {
            i1 = ele.fToP( j, 1 );
            i2 = ele.fToP( j, 2 );
            i3 = ele.fToP( j, 3 );
            // go to global
            i1 = ( iv->point( i1 ) ).id();
            i2 = ( iv->point( i2 ) ).id();
            i3 = ( iv->point( i3 ) ).id();
            if ( RegionMesh3D::FaceShape::numVertices == 4 )
            {
                i4 = ele.fToP( j, 4 );
                i4 = ( iv->point( i4 ) ).id();
                bface = ( makeBareFace( i1, i2, i3, i4 ) ).first;
            }
            else
            {
                bface = ( makeBareFace( i1, i2, i3 ) ).first;
            }

            if ( ( fi = bfaces.find( bface ) ) == bfaces.end() )
            {
                bfaces.insert( make_pair( bface, make_pair( iv->id(), j ) ) );
            }
            else
            {
                if ( allFaces && i1 > i2 )
                    intfaces.insert( ( make_pair( bface, make_pair( iv->id(), j ) ) ) );
                bfaces.erase( fi ); // counted twice: internal face
                ++numInternalFaces;
            }
        }
    }
    return bfaces.size();
}

template <typename RegionMesh3D>
UInt findBoundaryFaces( const RegionMesh3D & mesh, TempFaceContainer & bfaces, UInt & numInternalFaces )
{
    TempFaceContainer dummy;
    return findFaces( mesh, bfaces, numInternalFaces, dummy, false );
}



//! Finds boundary edges.
/*!  A low level routine, not meant to be called directly. It creates a
container with all the information needed to set up properly the boundary
edges connectivities.
 
\param mesh A 3D mesh.
 
\param bedges This container contains a set with the BareEdge of the
boundary edges.
 
\return Number of boundary edges found.
 
\pre The list of boundary faces must be correctly set.
*/
template <typename RegionMesh3D>
UInt findBoundaryEdges( const RegionMesh3D & mesh, TempEdgeContainer & bedges )
{
    UInt i1, i2;
    BareEdge bedge;
    typedef typename RegionMesh3D::FaceShape FaceShape;
    typedef typename RegionMesh3D::Faces Faces;


    if ( ! mesh.hasFaces() )
        return 0;

    // clean first in case it has been alredy used
    bedges.clear();

    for ( typename Faces::const_iterator ifa = mesh.faceList.begin();
            ifa != mesh.faceList.begin() + mesh.numBFaces(); ++ifa )
    {
        for ( ID j = 1;j <= mesh.numLocalEdgesOfFace();++j )
        {
            i1 = FaceShape::eToP( j, 1 );
            i2 = FaceShape::eToP( j, 2 );
            // go to global
            i1 = ( ifa->point( i1 ) ).id();
            i2 = ( ifa->point( i2 ) ).id();
            bedge = ( makeBareEdge( i1, i2 ) ).first;
            bedges.insert( make_pair( bedge, make_pair( ifa->id(), j ) ) );
        }
    }
    return bedges.size();
}

//! Finds all  edges.
/*!  A low level routine, not meant to be called directly. It creates a
container with all the information needed to set up properly the edge connectivities.
 
\param mesh A 3D mesh.
 
\param bedges This container contains a set of  BareEdges for all mesh edges.
 
\return Number of edges found.
 
*/

template <typename RegionMesh3D>
UInt findInternalEdges( const RegionMesh3D & mesh, const TempEdgeContainer & boundary_edges, TempEdgeContainer & internal_edges )
{
    UInt i1, i2;
    BareEdge bedge;
    typedef typename RegionMesh3D::ElementShape VolumeShape;
    typedef typename RegionMesh3D::Volumes Volumes;


    ASSERT0( mesh.numVolumes() > 0, "We must have some 3D elements stored n the mesh to use this function!" );

    internal_edges.clear();


    for ( typename Volumes::const_iterator ifa = mesh.volumeList.begin();
            ifa != mesh.volumeList.end(); ++ifa )
    {
        for ( ID j = 1;j <= mesh.numLocalEdges();++j )
        {
            i1 = VolumeShape::eToP( j, 1 );
            i2 = VolumeShape::eToP( j, 2 );
            // go to global
            i1 = ( ifa->point( i1 ) ).id();
            i2 = ( ifa->point( i2 ) ).id();
            bedge = ( makeBareEdge( i1, i2 ) ).first;
            if ( boundary_edges.find( bedge ) == boundary_edges.end() )
                internal_edges.insert( make_pair( bedge, make_pair( ifa->id(), j ) ) );
        }
    }
    return internal_edges.size();
}
/*
*************************************************************************************
                            MARKERS HANDLERS
*************************************************************************************
*/ 
//! \defgroup marker_handlers Used to manage missing handlers

/*! \ingroup marker_handlers
 
//! \brief Sets the marker flag of a GeoElement of dimension greater one
 
 It gets the stronger marker of the GeoElement points. The marker
hierarchy is defined in the marker.h file.  It returns a bool indicating if
the flag has changed. If any of the vertices has an unset marker the result
is an unset flag for the GeoElement.
 
\warning It overrides the original marker flag.
*/
template <typename GeoElement>
EntityFlag inheritStrongerMarker( GeoElement & fp )
{
    ASSERT_PRE( GeoElement::nDim > 0, "A GeoElement with ndim<1 cannot inherit marker flags" );

    fp.setMarker( fp.point( 1 ).marker() );
    for ( ID j = 2;j <= GeoElement::numVertices;++j )
        fp.setStrongerMarker( fp.point( j ).marker() );
    return fp.marker();

}


/*! \ingroup marker_handlers
 
//! \brief Sets the marker flag of a GeoElement of dimension greater one
 
 It gets the weaker marker of the GeoElement points. The marker
hierarchy is defined in the marker.h file.  It returns a bool indicating if
the flag has changed. If any of the vertices has an unset marker the result
is an unset flag for the GeoElement.
 
  \warning It overrides the original marker flag.*/
template <typename GeoElement>
EntityFlag inheritWeakerMarker( GeoElement & fp )
{
    ASSERT_PRE( GeoElement::nDim > 0, "A GeoElement with ndim<1 cannot inherit marker flags" );
    fp.setMarker( fp.point( 1 ).marker() );
    for ( ID j = 2;j <= GeoElement::numVertices;++j )
        fp.setWeakerMarker( fp.point( j ).marker() );
    return fp.marker();

}
/*!
  This routine tests if the topological descrption of boundary face is sane. In particular
  all boundary edges must be adjacent to only 2 surface elements and the orientation must be correct.
 
  \param mesh a mesh
  \param numBedges  The function also returns the number of boundary edges in numBedges.
  \return It it returns 0 the test has been passed. If  not it returns the number of of wrong boundary edges.
  \warning numBEdges is properly set only if the test has been passed.
*/
template <typename RegionMesh3D>
UInt testClosedDomain_Top( RegionMesh3D const & mesh, UInt & numBEdges )
{

    typedef std::set
        <BareEdge, cmpBareItem<BareEdge> > TempEdgeContainer2;
    TempEdgeContainer2 bedges;
    UInt i1, i2;
    BareEdge bedge;
    typename RegionMesh3D::BElementShape ele;
    typedef typename RegionMesh3D::Faces Faces;
    typedef typename RegionMesh3D::FaceType FaceType;
    TempEdgeContainer2::iterator ed;


    // clean first in case it has been alredy used

    typename Faces::const_iterator iv = mesh.faceList.begin();

    for ( UInt k = 0;k < mesh.numBFaces();++k )
    {
        ASSERT( iv != mesh.faceList.end(), " Trying to get not existing face" << k << " " << mesh.numBFaces() );

        for ( ID j = 1;j <= FaceType::numEdges;++j )
        {
            i1 = ele.eToP( j, 1 );
            i2 = ele.eToP( j, 2 );
            // go to global
            i1 = ( iv->point( i1 ) ).id();
            i2 = ( iv->point( i2 ) ).id();
            bedge = ( makeBareEdge( i1, i2 ) ).first;

            if ( ( ed = bedges.find( bedge ) ) == bedges.end() )
            {
                bedges.insert( bedge );
                ++numBEdges;
            }
            else
            {
                bedges.erase( ed );
            }
        }
        ++iv;
    }
    return bedges.size();
}
/*
*****************************************************************************
                                MARKERS FIXING
*****************************************************************************
*/

//! Check wether all markers of a the goemetry entities stored in a list are set
template <typename MeshEntityList>
bool checkMarkerSet( const MeshEntityList & list )
{
    typedef typename MeshEntityList::const_iterator C_Iter;
    bool ok( true );
    for ( C_Iter l = list.begin();l != list.end();++l )
        ok = ( ok & l->isMarkerSet() );
    return ok;
}

//! Sets the marker flag for all boundary edges by inheriting them from boundary points.
/*! The paradigm is that an edge <B>WHOSE MARKER HAS NOT ALREADY BEEN
  SET</B> will get the WEAKER marker flag among its VERTICES. For instance
  is a vertex is assigned to an Essential B.C and the other to a Natural
  B.C. the edge will get the flag related to the Natural B.C.
 
  /param mesh A mesh
  /param clog ostream to which the logging of the map of the newly assigned marked will be output
  /param err ostream to which error messages will be sent
 
  /todo better handling of flags: all function handling flags should be
  wrapped into a class
*/

template <typename RegionMesh>
void
setBEdgesMarker( RegionMesh & mesh, std::ostream & clog = std::cout, std::ostream & err = std::cerr, bool verbose = true )
{
    typename RegionMesh::EdgeType * fp = 0;
    unsigned int count( 0 );

    if ( verbose )
        clog << "NEW EDGE MARKER MAP" << std::endl << " ID->New Marker" << std::endl;

    for ( ID k = 1; k <= mesh.numBEdges(); ++k )
    {
        fp = &( mesh.edge( k ) );
        if ( fp->isMarkerUnset() )
        {
            inheritWeakerMarker( *fp );
            if ( verbose )
            {
                clog << fp->id() << " -> ";
                fp->printFlag( clog );
                clog << " ";
                if ( ++count % 3 == 0 )
                    clog << std::endl;
            }
        }
    }
    if ( verbose )
        clog << std::endl;
}


//! Sets the marker flag for all boundary faces by inheriting them from boundary points.
/*! The paradigm is that a face WHOSE MARKER HAS NOT ALREADY BEEN SET will
  get the WEAKER marker flag among its VERTICES. For instance if a vertex
  is assigned to a Natural B.C and the others to a Natural B.C. the face
  will get the flag related to the Natural B.C.
 
  /todo better handling of flags: all function handling flags should be
  wrapped into a class
*/
template <typename RegionMesh>
void
setBFacesMarker( RegionMesh & mesh, std::ostream & clog = std::cout, std::ostream & err = std::cerr, bool verbose = true )
{
    typename RegionMesh::FaceType * fp = 0;
    unsigned int count( 0 );

    if ( verbose )
        clog << "NEW FACE MARKER MAP" << std::endl << " ID->New Marker" << std::endl;

    for ( UInt k = 1;k <= mesh.numBFaces();++k )
    {
        fp = &( mesh.face( k ) );
        if ( fp->isMarkerUnset() )
        {
            inheritWeakerMarker( *fp );
            if ( verbose )
            {
                clog << fp->id() << " -> ";
                fp->printFlag( clog );
                clog << " ";
                if ( ++count % 3 == 0 )
                    clog << std::endl;
            }
        }
    }
    if ( verbose )
        clog << std::endl;
}

//! It sets the marker flag of boundary points, by inheriting it from boundary elements.
/*! The paradigm is that a point whose marker flag is unset will inherhit
  the strongest marker flag of the surrounding Boundary Elements, with the
  convention that if the marker flag of one of the surrounding boundary
  elements is null is ignored.
*/
template <typename RegionMesh>
void
setBPointsMarker( RegionMesh & mesh, std::ostream & clog = std::cout, std::ostream& err = std::cerr, bool verbose = false )
{
    // First looks at points whose marker has already been set
    std::vector<bool> markset( mesh.storedPoints(), false );

    typedef typename RegionMesh::Points::iterator PointIterator;
    typedef typename RegionMesh::BElementShape BElementShape;

    std::vector<bool>::iterator pm = markset.begin();

    for ( PointIterator p = mesh.pointList.begin();p != mesh.pointList.end();++p )
        *( pm++ ) = p->isMarkerSet();

    typename RegionMesh::BElementType * fp = 0;
    for ( UInt k = 1;k <= mesh.numBElements();++k )
    {
        fp = &( mesh.bElement( k ) );
        if ( fp->isMarkerSet() )
        {
            for ( UInt j = 1;j <= BElementShape::numPoints;++j )
            {
                if ( !markset[ ( fp->point( j ).id() ) - 1 ] )
                    fp->point( j ).setStrongerMarker( fp->marker() );
            }
        }
    }
    unsigned int count( 0 );
    if ( verbose )
    {
        clog << "**** NEW POINTS MARKERS **************" << std::endl;
        clog << "id->marker    id->marker     id->marker" << std::endl;
        pm = markset.begin();
        for ( PointIterator p = mesh.pointList.begin();p != mesh.pointList.end();++p )
        {
            if ( *pm++ )
            {
                clog << p->id() << " -> ";
                p->printFlag( clog );
                clog << " ";
                if ( ++count % 3 )
                    clog << std::endl;
            }
        }
    }
}
/*
*****************************************************************************
                                FIXING ID AND COUNTERS
*****************************************************************************
*/ 
//! \brief Verifies if a list of mesh entities hav ethe ID properly set.
/* More precisely, the id() must correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.
 
   \pre The template argument MeshEntityList must be a stl
   compliant container and its elements must have the method id().
*/
template <typename MeshEntityList>
bool checkIdnumber( const MeshEntityList & list )
{
    typedef typename MeshEntityList::const_iterator C_Iter;
    bool ok( true );
    unsigned int count( 1 );
    for ( C_Iter l = list.begin();l != list.end() && ok;++l, ++count )
        ok = ( l->id() == count );
    return ok;
}

//! \brief Fixes a a list of mesh entities so that the ID is properly set.
/* \post  The id will correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.
 
   \pre The template argument MeshEntityList must be a stl
   compliant container and its elements must have the method UInt &id().
*/
template
<typename MeshEntityList>
void fixIdnumber( MeshEntityList & list )
{
    unsigned int count( 0 );
    typedef typename MeshEntityList::iterator Iter;
    for ( Iter l = list.begin() ;l != list.end(); ++l )
        l->id() = ++count;
}

/*! \brief Fixes boundary points counter
  It fix the boundary points counter by counting
  how many points have te boundary flag set.
  It also reset the Bpoints list.
 
  \pre It assumes that the points have the boundary flag corretly set
*/

template <typename RegionMesh>
void
setBPointsCounters( RegionMesh & mesh )
{

    unsigned int countBP( 0 );
    unsigned int countBV( 0 );

    mesh._bPoints.clear();

    for ( UInt k = 1;k <= mesh.numVertices();++k )
    {
        if ( mesh.isBoundaryPoint( k ) )
        {
            ++countBP;
            ++countBV;
        }
    }

    for ( UInt k = mesh.numVertices() + 1;k <= mesh.numPoints();++k )
    {
        if ( mesh.isBoundaryPoint( k ) )
        {
            ++countBP;
        }
    }

    mesh.numBVertices() = countBV;
    mesh.setNumBPoints( countBP );
    mesh._bPoints.reserve( countBP );

    for ( UInt k = 1;k <= mesh.numPoints();++k )
    {
        if ( mesh.isBoundaryPoint( k ) )
            mesh._bPoints.push_back( &mesh.point( k ) );
    }
}

/*
*****************************************************************************
                                BOUNDARY INDICATOR FIXING
*****************************************************************************
*/ 
//! It fixes boundary flag on points laying on boundary faces.
/*!
  \param mesh a mesh
  \param clog logging stream
  \param err error stream
  \param verbose If true you have a verbose output
 
  \pre mesh point list must exists and boundary face lsist  must have been set properly.
*/
template <typename RegionMesh>
void
fixBPoints( RegionMesh & mesh, std::ostream & clog = std::cout, std::ostream & err = std::cerr, bool verbose = true )
{
    ASSERT_PRE( mesh.numPoints() > 0, "The point list should not be empty" );
    ASSERT_PRE( mesh.numBElements() > 0, "The BElements list should not be empty" );

    typedef typename RegionMesh::BElements BElements;
    typedef typename RegionMesh::BElementShape BElementShape;
    typename RegionMesh::BElementType * fp;

    if ( verbose )
        clog << "New BPoints Found " << std::endl;
    for ( UInt k = 1;k <= mesh.numBElements();++k )
    {
        fp = &( mesh.bElement( k ) );
        for ( UInt j = 1;j <= BElementShape::numPoints;++j )
        {
            if ( verbose && !fp->point( j ).boundary() )
                clog << "ID: " << fp->point( j ).id() << std::endl;
            fp->point( j ).boundary() = true;
        }
    }
    // Fix now the number of vertices/points
    setBPointsCounters( mesh );
}

//!It makes sure that boundary edges are stored first
/*!
\pre It assumes that boundary points are properly stored in the mesh
*/
template <typename RegionMesh>
void
setBoundaryEdgesFirst( RegionMesh & mesh )
{

    typedef typename RegionMesh::Edges Edges;
    // set the functor
    EnquireBEntity<RegionMesh > enquireBEdge( mesh );

    std::partition( mesh.edgeList.begin(), mesh.edgeList.end(), enquireBEdge );
    fixIdnumber( mesh.edgeList );
}

//!It makes sure that boundary faces are stored first
/*!
\pre It assumes that boundary points are properly stored in the mesh
*/
template <typename RegionMesh>
void
setBoundaryFacesFirst( RegionMesh & mesh )
{

    typedef typename RegionMesh::Faces Faces;
    // set the functor
    EnquireBEntity<RegionMesh> enquireBFace( mesh );

    std::partition( mesh.faceList.begin(), mesh.faceList.end(), enquireBFace );
    fixIdnumber( mesh.faceList );
}

//! Tests if boundary faces are stored first
/*! \return true if boundary faces are indeed stored first
  \pre It assumes that boundary points are set */
template <typename RegionMesh>
bool checkBoundaryFacesFirst( const RegionMesh & mesh )
{

    typedef typename RegionMesh::Faces Faces;

    // set the functor
    EnquireBEntity<RegionMesh> enquireBFace( mesh );
    typename RegionMesh::FaceType * fp;
    bool ok( true );

    for ( UInt k = 1;k <= mesh.numBElements();++k )
        ok = ok && enquireBFace( mesh.boundaryFace( k ) );
    for ( UInt k = mesh.numBElements() + 1;k <= mesh.storedFaces();++k )
        ok = ok && ! enquireBFace( mesh.face( k ) );

    return ok;
}

//! Tests if boundary edges are stored first
/*! \return true if boundary edges are indeed stored first
  \pre It assumes that boundary points are set */
template <typename RegionMesh>
bool checkBoundaryEdgesFirst( const RegionMesh & mesh )
{

    typedef typename RegionMesh::Edges Edges;

    // set the functor
    EnquireBEntity<RegionMesh> enquireBEdge( mesh );
    typename RegionMesh::EdgeType * fp;
    bool ok( true );

    for ( UInt k = 1;k <= mesh.numBEdges();++k )
        ok = ok && enquireBEdge( mesh.boundaryEdge( k ) );
    for ( UInt k = mesh.numBEdges() + 1;k <= mesh.storedEdges();++k )
        ok = ok && ! enquireBEdge( mesh.edge( k ) );
    return ok;
}

/*
*****************************************************************************
 UTILITIES TO VERIFY/CREATE FACES/EDGES
*****************************************************************************
*/ 
//! It fixes boundary faces so that they are consistently numbered with volumes.

/*! An important step for building degrees of freedom on faces.  It also
  fixes other face related data.
\param mesh a mesh
\param err  ostream for error messages
\param sw A switch that will contain information on what has been done
Possible values are
<ol>
<li>NUM_FACES_MISMATCH</li>
<li>FIXED_FACE_COUNTER</li>
<li>BFACE_MISSING</li>
<li>BFACE_STORED_MISMATCH</li>
<li>BELEMENT_COUNTER_UNSET</li>
<li>BFACE_STORED_MISMATCH</li>
<li>FIXED_MAX_NUM_FACES</li>
</ol>
 
\param fixMarker If set to the true value all faces without a markerFlag set will inherit it from the points.
 
\param clog ostream that will all information regarding the markers
 
\param verbose if falso nothng is written to clog
 
\param numFaces It returns the number of faces found by the function
 
\param bfaces_found It returns the number of boundary faces found by the function
 
\param ext_container. If not NULL it is a pointer to an external map of bondary faces, already
  produced by a call to findBoundaryFaces(). This parameter may be used to save al lot of computational work, since
  findBoundaryFaces() is rather expensive.
 
\pre Boundary faces list must be properly set.
*/

template <class RegionMesh3D>
bool fixBoundaryFaces( RegionMesh3D & mesh,
                       std::ostream & clog, std::ostream &err, Switch & sw,
                       UInt & numFaces, UInt & bfaces_found,
                       bool fixMarker = false, bool verbose = false,
                       TempFaceContainer * ext_container )
{

    typedef typename RegionMesh3D::Volumes Volumes;
    typedef typename RegionMesh3D::VolumeType VolumeType;
    typedef typename RegionMesh3D::Faces Faces;
    typedef typename RegionMesh3D::FaceType FaceType;

    UInt i1, i2, i3, i4;
    BareFace bface;
    VolumeType * pv;
    typename Faces::iterator fit;
    typename RegionMesh3D::VolumeShape ele;
    TempFaceContainer * bfaces;
    TempFaceContainer::iterator fi;
    std::pair<ID, ID>info;
    ID j;
    ID vol;
    UInt numInternalFaces;
    bool notfound( false );
    bool extcont( false );

    if ( extcont = ( ext_container != 0 ) )
    {
        bfaces = ext_container;
        bfaces_found = bfaces->size();
    }
    else
    {
        bfaces = new TempFaceContainer;
        bfaces_found = findBoundaryFaces( mesh, *bfaces, numInternalFaces );
        numFaces = bfaces_found + numInternalFaces;
    }


    bool notEnough = mesh.storedFaces() < bfaces_found;



    if ( notEnough )
    {
        err << "WARNING: number of B. Faces stored smaller" << std::endl;
        err << "than the number of bfaces found  and build is not set" << std::endl;
        err << "POSSIBLE ERROR" << std::endl;
        sw.create( "BFACE_STORED_MISMATCH", true );
    }

    if ( mesh.numBElements() == 0 )
    {
        err << "ERROR: Boundary Element counter was not set" << std::endl;
        err << "I Cannot proceed because the situation is ambiguous" << std::endl;
        err << "Please check and eventually either: (a) call buildBoundaryFaces()" << std::endl;
        err << "or (b) set the correct number of bfaces in the mesh using mesh.numBElements()" << std::endl;
        err << "ABORT";
        sw.create( "BELEMENT_COUNTER_UNSET", true );
    }

    if ( mesh.numBFaces() != bfaces_found )
    {
        err << "WARNING: B Face counter in mesh is set to " << mesh.numBFaces();
        err << "While I have found " << bfaces_found << " B. Elements in mesh" << std::endl;
        err << "Plese check... I continue anyway" << std::endl;
        sw.create( "BFACE_COUNTER_MISMATCH", true );
    }

    if ( verbose )
    {
        clog << "**** Marker Flags for Fixed Boundary Faces ***" << std::endl;
        clog << " (it only contains those that were fixed because unset !" << std::endl;
        clog << "id->marker   id->marker  id->marker" << std::endl;
    }

    UInt count( 0 );

    fit = mesh.faceList.begin();
    for ( UInt facid = 0;facid < mesh.numBElements();++facid )
    {
        i1 = ( fit->point( 1 ) ).id();
        i2 = ( fit->point( 2 ) ).id();
        i3 = ( fit->point( 3 ) ).id();
        if ( RegionMesh3D::FaceShape::numVertices == 4 )
        {
            i4 = ( fit->point( 4 ) ).id();
            bface = ( makeBareFace( i1, i2, i3, i4 ) ).first;
        }
        else
        {
            bface = ( makeBareFace( i1, i2, i3 ) ).first;
        }
        fi = bfaces->find( bface );
        if ( fi == bfaces->end() )
        {
            notfound = true;
        }
        else
        {
            info = fi->second;
            vol = info.first; // Element ID
            pv = &mesh.volume( vol ); // Element
            j = info.second;       // The local ID of face on element
            // Reset face point definition to be consistent with face.
            for ( UInt k = 1;k <= FaceType::numPoints;++k )
            {
                fit->setPoint( k, pv->point( ele.fToP( j, k ) ) );
            }
            // Correct extra info
            fit->ad_first() = vol;
            fit->pos_first() = j;
            if ( fit->markerUnset() )
            {
                inheritWeakerMarker( *fit );
                if ( verbose )
                {
                    clog << fit->id() << " -> ";
                    fit->printFlag( clog );
                    clog << " ";
                    if ( ++count % 3 == 0 )
                        clog << std::endl;
                }
            }
            // Take out face from temporary container
            bfaces->erase( fi );
        }
        ++fit;
    }

    if ( !extcont )
        delete bfaces;

    if ( notfound )
    {
        err << "WARNING: At least one boundary face has not been found on the list stored in RegionMesh3D\n";
        sw.create( "BFACE_MISSING", true );
    }

    if ( verbose )
    {
        clog << std::endl << "  *****  END OF LIST ****" << std::endl;
    }

    // Here I calculate the number of faces,

    if ( mesh.numFaces() != numFaces )
    {
        err << "WARNING: faces counter in mesh  should be " << numFaces << std::endl;
        err << "(bfaces->size()+numInternalFaces)" << std::endl;
        err << "it is instead " << mesh.numFaces();
        sw.create( "NUM_FACES_MISMATCH", true );
    }
    mesh.setLinkSwitch( std::string( "HAS_BOUNDARY_FACES" ) );

    return true;
}

//! Builds faces
/*! This function may alternatively be used to build the compulsory boundary faces, all the mesh faces, or just add to an
  existing list of just boundary faces the internal ones.
 
  \param mesh A mesh
 
  \param clog Log file for information on the newly created markers
 
  \param err  Error stream
 
  \param buildbounary if true the function builds boundary faces
 
  \param buildinternal if true the function builds internal faces
 
  \param verbose. If true markerFrlags info is written on clog.
 
  \param numInternalFaces It returns the number of internal faces (only if ext_container is not provided!)
 
  \param bfaces_found It returns the number of boundary faces
 
  \param ext_container. If not NULL it is a pointer to an external map of bondary faces, already
  produced by a call to findBoundaryFaces(). This parameter may be used to save al lot of computational work, since
  findBoundaryFaces() is rather expensive.
 
  \pre If buildinternal=true and buildboundary=false the mesh must contain a proper list
  of boundary faces
 
  \note By setting buildinternal=true and buildboundary=true the function just fixes the counters
  with the number of faces in the mesh
 */
template <class RegionMesh3D>
bool buildFaces( RegionMesh3D & mesh,
                 std::ostream & clog, std::ostream &err, UInt & bfaces_found,
                 UInt & numInternalFaces,
                 bool buildboundary = true,
                 bool buildinternal = false,
                 bool verbose = false,
                 TempFaceContainer * ext_container = 0 )
{
    UInt i1, i2, i3, i4;
    typename RegionMesh3D::VolumeShape ele;
    typedef typename RegionMesh3D::Volumes Volumes;
    typedef typename RegionMesh3D::VolumeType VolumeType;
    typedef typename RegionMesh3D::Faces Faces;
    typedef typename RegionMesh3D::FaceType FaceType;
    VolumeType * pv;
    TempFaceContainer* bfaces;
    TempFaceContainer::iterator fi;
    bool extcont( false );

    std::pair<ID, ID>info;
    ID j, id;
    ID vol;

    if ( extcont = ( ext_container != 0 ) )
    {
        bfaces = ext_container;
        bfaces_found = bfaces->size();
    }
    else
    {
        bfaces = new TempFaceContainer;
        bfaces_found = findBoundaryFaces( mesh, *bfaces, numInternalFaces );
    }

    if ( buildboundary )
        mesh.faceList.clear();
    mesh.setNumBFaces( bfaces_found );
    if ( !buildinternal )
    {
        mesh.setMaxNumFaces( bfaces_found, false );
        mesh.numFaces() = numInternalFaces + bfaces_found;
    }
    else
    {
        mesh.setMaxNumFaces( numInternalFaces + bfaces_found, true );
    }

    FaceType face;

    if ( buildboundary )
    {

        if ( verbose )
        {
            clog << "**** Marker Flags for Newly Created Boundary Faces ***" << std::endl;
            clog << "id->marker   id->marker  id->marker" << std::endl;
        }

        for ( fi = bfaces->begin();fi != bfaces->end();++fi )
        {
            info = fi->second;
            vol = info.first; // Element ID
            pv = &mesh.volume( vol ); // Element
            j = info.second;       // The local ID of face on element

            for ( UInt k = 1;k <= FaceType::numPoints;++k )
                face.setPoint( k, pv->point( ele.fToP( j, k ) ) );
            // Add extra info
            face.ad_first() = vol;
            face.pos_first() = j;
            // Get marker value
            inheritWeakerMarker( face );
            id = mesh.addFace( face, true ).id();
            if ( verbose )
            {
                if ( id % 3 == 0 )
                    clog << std::endl;
                clog << id << " -> ";
                face.printFlag( clog );
                clog << " ";
            }
        }
        mesh.setLinkSwitch( std::string( "HAS_BOUNDARY_FACES" ) );
        if ( ! buildinternal )
            mesh.unsetLinkSwitch( std::string( "HAS_ALL_FACES" ) );
        mesh.setLinkSwitch( std::string( "FACES_HAVE_ADIACENCY" ) );
    }

    if ( !extcont )
        delete bfaces;

    if ( ! buildinternal )
        return true;


    if ( !buildboundary )
    {
        if ( mesh.storedFaces() < mesh.numBFaces() )
        {
            err << "ERROR: mesh has not boundary faces, cannot just create internal ones!!!" << std::endl;
            err << "ABORT CONDITION" << std::endl;
            return false;
        }
        else if ( mesh.storedFaces() > mesh.numBFaces() )
        {
            mesh.faceList.resize( mesh.numBFaces() );
        }
    }


    // I may get rid of the bfaces container. Unfortunately now I need a more complex structure, a BareItemsHandel,
    // in order to generate the internal faces id. An alternative would be to use the point data to identify boundary faces
    // as the ones with all point on the boundary. Yet in this function we do not want to use a priori infromation, so that
    // it might work even if the points boundary flag is not properly set.


    BareItemsHandler<BareFace> _be;
    std::pair<UInt, bool> e;
    std::pair<BareFace, bool> _face;

    for ( UInt j = 0; j < mesh.faceList.size();++j )
    {
        i1 = ( mesh.faceList[ j ].point( 1 ) ).id();
        i2 = ( mesh.faceList[ j ].point( 2 ) ).id();
        i3 = ( mesh.faceList[ j ].point( 3 ) ).id();
        if ( RegionMesh3D::FaceShape::numVertices == 4 )
        {
            i4 = ( mesh.faceList[ j ].point( 4 ) ).id();
            _face = makeBareFace( i1, i2, i3, i4 );
        }
        else
        {
            _face = makeBareFace( i1, i2, i3 );
        }
        _be.addIfNotThere( _face.first );
    }

    EntityFlag mm( mesh.marker() );
    ID vid;

    for ( typename Volumes::iterator iv = mesh.volumeList.begin();
            iv != mesh.volumeList.end(); ++iv )
    {
        vid = iv->id();
        // REMEMBER: numbering from 1
        for ( UInt j = 1;j <= mesh.numLocalFaces();j++ )
        {
            i1 = ele.fToP( j, 1 );
            i2 = ele.fToP( j, 2 );
            i3 = ele.fToP( j, 3 );
            // go to global
            i1 = ( iv->point( i1 ) ).id();
            i2 = ( iv->point( i2 ) ).id();
            i3 = ( iv->point( i3 ) ).id();
            if ( RegionMesh3D::FaceShape::numVertices == 4 )
            {
                i4 = ele.fToP( j, 4 );
                i4 = ( iv->point( i4 ) ).id();
                _face = makeBareFace( i1, i2, i3, i4 );
            }
            else
            {
                _face = makeBareFace( i1, i2, i3 );
            }
            e = _be.addIfNotThere( _face.first );
            if ( e.second )
            {
                // a new face It must be internal.
                for ( UInt k = 1;k <= FaceType::numPoints;++k )
                    face.setPoint( k, iv->point( ele.fToP( j, k ) ) );
                face.ad_first() = vid;
                face.pos_first() = j;
                // gets the marker from the RegionMesh
                face.setMarker( mm );
                mesh.addFace( face, false ); //The id should be correct
            }
            else
            {
                if ( e.first > bfaces_found )  // internal
                {
                    mesh.faceList( e.first ).ad_second() = vid;
                    mesh.faceList( e.first ).pos_second() = j;
                }
            }
        }
    }
    mesh.setLinkSwitch( std::string( "HAS_ALL_FACES" ) );
    return true;
}

//! It builds edges.

/*! This function may alternatively be used to build the boundary edges, all the mesh faces, or just add the internal edges
  to an existing list of just boundary edges.
 
  \param mesh A mesh
 
  \param clog Log file for information on the newly created markers for boundary edges
 
  \param err  Error stream
 
  \param bedges_found Returns the number of boundary edges
 
  \param iedges_found Returns the number of internal edges
 
  \param buildbounary if true the function builds boundary edges
 
  \param buildinternal if true the function builds internal edges
 
  \param verbose. If true markerFlags info is written on clog.
 
  \param ext_container. If not NULL it is a pointer to an external map of bondary edges, already
  produced by a call to findBoundaryEdges(). This parameter may be used to save al lot of computational work, since
  findBoundaryEdges() is rather expensive.
 
  \pre If buildinternal=true and buildboundary=false the mesh must contain a proper list
  of boundary edges
  \pre The mesh must copntain a proper list of boundary faces
 
  \note By setting buildinternal=true and buildboundary=true the function just fixes the counters
  with the number of edges in the mesh
 */

template <typename RegionMesh3D>
bool buildEdges( RegionMesh3D & mesh, std::ostream & clog, std::ostream &err, UInt & bedges_found,
                 UInt & iedges_found, bool buildboundary = true, bool buildinternal = false, bool verbose = false,
                 TempEdgeContainer * ext_container = 0 )
{
    typedef typename RegionMesh3D::Volumes Volumes;
    typedef typename RegionMesh3D::Faces Faces;
    typedef typename RegionMesh3D::VolumeType VolumeType;
    typedef typename RegionMesh3D::VolumeShape VolumeShape;
    typedef typename RegionMesh3D::Edges Edges;
    typedef typename RegionMesh3D::EdgeType EdgeType;
    typedef typename RegionMesh3D::FaceType FaceType;
    typedef typename RegionMesh3D::FaceShape FaceShape;
    typename RegionMesh3D::FaceType * pf;
    typename RegionMesh3D::VolumeType * pv;

    TempEdgeContainer * bedges;
    TempEdgeContainer iedges;
    std::pair<ID, ID>info;
    ID j, id;
    ID facID;


    bool extcont( false );


    if ( extcont = ( ext_container != 0 ) )
    {
        bedges = ext_container;
        bedges_found = bedges->size();
    }
    else
    {
        bedges = new TempEdgeContainer;
        bedges_found = findBoundaryEdges( mesh, *bedges );
    }

    iedges_found = findInternalEdges( mesh, *bedges, iedges );
    // free some memory if not needed!
    if ( !buildinternal )
        iedges.clear();
    if ( !buildboundary && buildinternal )
    {
        if ( mesh.storedEdges() < bedges_found )
        {
            err << "ERROR in buildedges(): mesh does not contain boundary edges" << std::endl;
            err << "I need to set buildboundary=true" << std::endl;
            err << "ABORT CONDITION" << std::endl;
            return false;
        }
        else if ( mesh.storedEdges() > bedges_found )
        {
            mesh.edgeList.resize( bedges_found );
        }
    }
    mesh.setNumBEdges( bedges_found );
    mesh.numEdges() = ( bedges_found + iedges_found );
    if ( buildboundary )
        mesh.edgeList.clear();
    if ( buildboundary && ! buildinternal )
        mesh.setMaxNumEdges( bedges_found, false );
    if ( buildinternal )
        mesh.setMaxNumEdges( bedges_found + iedges_found, true );
    err << "Building edges from scratch" << std::endl;

    EdgeType edge;

    if ( buildboundary )
    {

        if ( verbose )
        {
            clog << "**** Marker Flags for Newly Created Boundary Edges ***" << std::endl;
            clog << "id->marker   id->marker   id->marker" << std::endl;
        }

        // First boundary.
        for ( TempEdgeContainer::iterator ei = bedges->begin();
                ei != bedges->end();++ei )
        {
            info = ei->second;
            facID = info.first; // Face ID
            pf = &mesh.face( facID ); // Face
            j = info.second;       // The local ID of edge on face
            for ( UInt k = 1;k <= EdgeType::numPoints;++k )
            {
                edge.setPoint( k, pf->point( FaceShape::eToP( j, k ) ) );
            }

            inheritWeakerMarker( edge ); // Get marker value inheriting from points

            id = mesh.addEdge( edge, true ).id();
            if ( verbose )
            {
                if ( id % 3 == 0 )
                    clog << std::endl;
                clog << id << " -> ";
                edge.printFlag( clog );
                clog << " ";
            }
        }

        if ( verbose )
            clog << std::endl << "  *****  END OF LIST OF BOUNDARY EDGES ****" << std::endl;

        mesh.setLinkSwitch( std::string( "HAS_BOUNDARY_EDGES" ) );
    }

    if ( !extcont )
        delete bedges;

    if ( !buildinternal )
    {
        mesh.unsetLinkSwitch( std::string( "HAS_ALL_EDGES" ) );
        return true;
    }



    // Now internal edges
    // free some memory

    for ( TempEdgeContainer::iterator ei = iedges.begin();
            ei != iedges.end();++ei )
    {
        info = ei->second;
        facID = info.first; // Volume ID
        pv = &mesh.volume( facID ); // Volume that generated the edge
        j = info.second;       // The local ID of edge on volume
        for ( UInt k = 1;k <= EdgeType::numPoints;++k )
            edge.setPoint( k, pv->point( VolumeShape::eToP( j, k ) ) );
        edge.setMarker( mesh.marker() ); // Get marker value: that of the mesh
        mesh.addEdge( edge, false );
    }

    mesh.setLinkSwitch( std::string( "HAS_ALL_EDGES" ) );

    return true;
}


/*
*****************************************************************************
 UTILITIES TO TRANSFORM A MESH
*****************************************************************************
 
*/ 
//! It builds a P2 mesh from P1 data.
/*! \author L.Formaggia.
  \version Version 1.0
  \pre All compulsory structures in mesh must have been already set: volumes and boundary faces.
  \pre Points list MUST have been dimensioned correctly!!!
  \note the function takes advantage of the fact that
*/
template <typename RegionMesh>
void
p1top2( RegionMesh & mesh, std::ostream & out = std::cout )
{

    typedef typename RegionMesh::ElementShape GeoShape;
    typedef typename RegionMesh::BElementShape GeoBShape;
    ASSERT_PRE( GeoShape::numPoints > 4, "p1top2 ERROR: we need a P2 mesh" );

    out << "Building P2 mesh points and connectivities from P1 data" << std::endl;


    typename RegionMesh::PointType * pp = 0;
    typename RegionMesh::EdgeType * pe = 0;
    typename RegionMesh::ElementType * pv = 0;
    typename RegionMesh::BElementType * pbe = 0;
    typedef typename RegionMesh::Elements Elements;
    typedef typename RegionMesh::BElements BElements;

    BareItemsHandler<BareEdge> _be;
    std::pair<UInt, bool> _edgeid;
    UInt i1, i2, e_id;
    std::pair<BareEdge, bool> _edge;
    typename RegionMesh::ElementShape ele;
    out << "Processing " << mesh.storedEdges() << " P1 Edges" << std::endl;
    UInt nbe = mesh.numBEdges();
    for ( UInt j = 1; j <= mesh.storedEdges();++j )
    {
        pe = & mesh.edge( j );
        i1 = ( pe->point( 1 ) ).id();
        i2 = ( pe->point( 2 ) ).id();
        pp = & mesh.addPoint( j <= nbe ); // true for boundary points
        pp->x() = ( ( pe->point( 1 ) ).x() +
                    ( pe->point( 2 ) ).x() ) * .5;
        pp->y() = ( ( pe->point( 1 ) ).y() +
                    ( pe->point( 2 ) ).y() ) * .5;
        pp->z() = ( ( pe->point( 1 ) ).z() +
                    ( pe->point( 2 ) ).z() ) * .5;

        // If we have set a marker for the boundary edge, that marker is inherited by the new created point
        // Otherwise the edge (and the new created point) gets the WORST marker among the two end Vertices
        /*
          JFG 07/2002:
          if the mesh file do not contain the edges (inria files), they are built in
          fixBoundaryEdges(...), but I suspect that this function does not attribute the right
          marker to the edges (maybe a problem in setWorseMarkerOfEntity, or something like that...)
          If you do #undef JFG : no change, if you do #define JFG, we do not consider the (wrong) marker
          of the edge to define the marker of the added node (I arbitrarily take the marker of the first
          node)
         */ 
        //#define JFG
        //#ifndef JFG
        // original version: DOES NOT work when the edges are not give in the mesh file
        if ( pe->markerUnset() )
            inheritWeakerMarker( *pe );
        pp->setMarker( pe->marker() );
        //#else
        // temporary version that works when the edges are not given in the mesh file
        // pe->setMarker(mesh.point(i1).marker());
        // pp->setMarker(pe->marker());
        //#endif
        // pp->id()=++i; // Indexing from 1
        // if I storealso non B. edges:
        // pointList[i].boundary()=
        // (edgeList[j].point(1)).boundary() &&(edgeList[j].point(2)).boundary();
        pe->setPoint( 3, pp ); //use overloaded version that takes a pointer
        _edge = makeBareEdge( i1, i2 );
        _edgeid = _be.addIfNotThere( _edge.first, pp->id() );
    }
    // Now the other edges, of which I do NOT build the global stuff
    // (I would need to check the switch but I will do that part later on)
    if ( GeoShape::nDim == 3 )
    {
        UInt nbf = mesh.numBFaces();
        UInt nbv = GeoBShape::numVertices;
        out << "Processing " << mesh.storedFaces() << " Face Edges" << std::endl;
        for ( UInt k = 1; k <= mesh.storedFaces(); ++k )
        {
            pbe = &mesh.face( k );
            for ( UInt j = 1;j <= mesh.numLocalEdgesOfFace();j++ )
            {
                i1 = GeoBShape::eToP( j, 1 );
                i2 = GeoBShape::eToP( j, 2 );
                i1 = ( pbe->point( i1 ) ).id();
                i2 = ( pbe->point( i2 ) ).id();
                _edge = makeBareEdge( i1, i2 );
                e_id = _be.id( _edge.first );
                if ( e_id != 0 )
                {
                    pp = &mesh.point( e_id );
                }
                else
                {
                    // new edge -> new Point
                    pp = &mesh.addPoint( k <= nbf ); // true for boundary points
                    _edgeid = _be.addIfNotThere( _edge.first, pp->id() );
                    pp->x() = ( mesh.point( i1 ).x() +
                                mesh.point( i2 ).x() ) * .5;
                    pp->y() = ( mesh.point( i1 ).y() +
                                mesh.point( i2 ).y() ) * .5;
                    pp->z() = ( mesh.point( i1 ).z() +
                                mesh.point( i2 ).z() ) * .5;
                    // If we have set a marker for the face, that marker is inherited by the new created point
                    pp->setMarker( pbe->marker() );
                }
                pbe->setPoint( nbv + j, pp );
            }
        }
    }

    out << "Processing " << mesh.numElements() << " Mesh Elements" << std::endl;
    UInt nev = GeoShape::numVertices;
    for ( UInt k = 1; k <= mesh.numElements(); ++k )
    {
        pv = &mesh.element( k );
        for ( UInt j = 1;j <= mesh.numLocalEdges();j++ )
        {
            i1 = ele.eToP( j, 1 );
            i2 = ele.eToP( j, 2 );
            i1 = ( pv->point( i1 ) ).id();
            i2 = ( pv->point( i2 ) ).id();
            _edge = makeBareEdge( i1, i2 );
            e_id = _be.id( _edge.first );
            if ( e_id != 0 )
            {
                pp = &mesh.point( e_id );
            }
            else
            {
                pp = &mesh.addPoint( false ); // cannot be on boundary is the mesh is proper!
                _edgeid = _be.addIfNotThere( _edge.first, pp->id() );
                pp->x() = ( mesh.point( i1 ).x() +
                            mesh.point( i2 ).x() ) * .5;
                pp->y() = ( mesh.point( i1 ).y() +
                            mesh.point( i2 ).y() ) * .5;
                pp->z() = ( mesh.point( i1 ).z() +
                            mesh.point( i2 ).z() ) * .5;
                pp->setMarker( pe->marker() );
            }
            pv->setPoint( nev + j, pp );
        }
    }
    /*=============================*/
    out << " ******* Done Construction of P2 Mmesh *******" << std::endl << std::endl;
}

//! Fix mesh switches
/*!
  Using some heuristics it tries to fix mesh switches
 */

// template<typename RegionMesh>
// void
// fixSwitches(RegionMesh ^ mesh, std::ostream & clog=std::cout, bool verbose=false)
// {

//   clog<<" ************** FIXING MESH SWITCHES **********************"<<std::endl;
//   clog<<"            Mesh switches Status before fixing"<<std::endl;
//   mesh.showLinkSwitch(verbose,clog);
//   if (mesh.storedFaces()> mesh.numBFaces()){
//     mesh.setLinkSwitch("HAS_ALL_FACES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_ALL_FACES");
//   }
//   if (mesh.numBFaces>0){
//     mesh.setLinkSwitch("HAS_BOUNDARY_FACES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDARY_FACES");
//   }
//   if (mesh.storedEdges()> mesh.numBEdges()){
//     mesh.setLinkSwitch("HAS_ALL_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_ALL_EDGES");
//   }
//   if (mesh.numBEdges()> 0){
//     mesh.setLinkSwitch("HAS_BOUNDARY_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDAY_EDGES");
//   }
//   if (mesh.storedFaces()>0){
//     if(mesh.face(1).ad_first(
//     mesh.setLinkSwitch("HAS_BOUNDARY_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDAY_EDGES");
//   }


}


#endif
