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
    @brief class to manage ghost data across procs

    @author Antonio Cervone <ant.cervone@gmail.com>

    @date 27-10-2011
*/

#ifndef _GHOSTHANDLER_HPP_
#define _GHOSTHANDLER_HPP_

#include <bitset>

#ifdef HAVE_HDF5
#include <EpetraExt_HDF5.h>
#endif

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChronoManager.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/NeighborMarker.hpp>

#ifdef HAVE_LIFEV_DEBUG
//#define LIFEV_GHOSTHANDLER_DEBUG 1
#endif

namespace LifeV
{

typedef std::bitset<4> NeighborType;

NeighborType const   POINT_NEIGHBORS = 0x1;
NeighborType const   RIDGE_NEIGHBORS = 0x2;
NeighborType const   FACET_NEIGHBORS = 0x4;
NeighborType const ELEMENT_NEIGHBORS = 0x8;
NeighborType const     ALL_NEIGHBORS = POINT_NEIGHBORS | RIDGE_NEIGHBORS | FACET_NEIGHBORS | ELEMENT_NEIGHBORS;

//! GhostHandler
/*!
  This class manages neighborhood information across processes.
  The aim is to have the possibility to build overlapping maps
  in order to ease the retrieving of ghosted values. The class
  offers also the possibility to build phisically overlapped
  meshes that do not require communication to retrieve mesh
  information from adjacent elements.
 */
template <typename MeshType>
class GhostHandler
{
public:

    //! @name Public Types
    //@{

    typedef MeshType mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> commPtr_Type;
    typedef MapEpetra map_Type;
    typedef boost::shared_ptr<map_Type> mapPtr_Type;
    typedef std::vector<Int> idList_Type;
    typedef boost::shared_ptr<idList_Type> idListPtr_Type;
    typedef std::vector<idList_Type> graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;
    typedef std::vector<idListPtr_Type> vertexPartition_Type;
    typedef boost::shared_ptr<vertexPartition_Type> vertexPartitionPtr_Type;
    typedef std::vector<markerID_Type> markerIDList_Type;
    typedef std::vector<int> markerIDListSigned_Type;

    //@}

    //! @name Constructors & Destructors
    //@{
    //! Constructor
    /*!
     * @param comm. Communicator
     */
    explicit GhostHandler ( commPtr_Type const& comm );

    //! Constructor
    /*!
     * @param fullMesh. Original mesh, before partitioning
     * @param comm. Communicator
     */
    GhostHandler ( meshPtr_Type fullMesh, commPtr_Type const& comm );

    //! Constructor
    /*!
     * @param fullMesh. Original mesh, before partitioning
     * @param localMesh. Local mesh of the current proc
     * @param map. Original map without overlapping
     * @param comm. Communicator
     */
    GhostHandler ( meshPtr_Type fullMesh,
                   meshPtr_Type localMesh,
                   mapPtr_Type map,
                   commPtr_Type const& comm );

    //! Destructor
    ~GhostHandler() {}

    //@}

    //! @name Get Methods
    //@{

    //! Full mesh getter
    mesh_Type const& fullMesh()
    {
        return *M_fullMesh;
    }

    //! Local mesh getter
    mesh_Type const& localMesh()
    {
        return *M_localMesh;
    }

    //! Standard map getter
    map_Type const& map()
    {
        return *M_map;
    }

    //! List of point neighbors to a point (identified by the global ID)
    neighborList_Type const& pointPointNeighborsList()
    {
        ASSERT ( !M_pointPointNeighborsList.empty(), "M_pointPointNeighborsList is empty" );
        return M_pointPointNeighborsList;
    }

    //! List of edge neighbors to a point (identified by the global ID)
    neighborList_Type const& pointEdgeNeighborsList()
    {
        ASSERT ( !M_pointEdgeNeighborsList.empty(), "M_pointEdgeNeighborsList is empty" );
        return M_pointEdgeNeighborsList;
    }

    //! List of element neighbors to a point (identified by the global ID)
    neighborList_Type const& pointElementNeighborsList()
    {
        ASSERT ( !M_pointElementNeighborsList.empty(), "M_pointElementNeighborsList is empty" );
        return M_pointElementNeighborsList;
    }

    //@}

    //! @name Set Methods
    //@{

    //! Set verbosity
    /*!
     * @param verbose
     */
    void setVerbose ( const bool& verbose )
    {
        M_verbose = verbose && ( M_me == 0 );
    }

    //@}

    //! @name General Methods
    //@{

    //! Initialize neighbors list
    void setUpNeighbors ( NeighborType const neighborType = ALL_NEIGHBORS );

    //! Release pointers to full and local mesh
    void release();

    //! Clean up neighbor lists
    void clean ( NeighborType const neighborType = ALL_NEIGHBORS );

#ifdef HAVE_HDF5
    //! Export neighbor lists to an hdf5 file
    /*!
     * @param fileName. Name of the file to write
     * @param truncate. Must be true when the file already exists on disk
     */
    void exportToHDF5 ( std::string const& fileName = "ghostmap", bool const& truncate = true );

    //! Import neighbor lists to an hdf5 file
    /*!
     * @param fileName. Name of the file to write
     */
    void importFromHDF5 ( std::string const& fileName = "ghostmap" );
#endif // HAVE_HDF5

    //! Create point neighbors to points and store them in the NeighborMarker
    void createPointNeighbors();

    //! Create the list of point neighbors to points
    void createPointPointNeighborsList();

    //! Create the list of point neighbors to points that are in the given list of MarkerIDs
    /*!
     * @param flags. The list of MarkerIDs to restrict to.
     */
    void createPointPointNeighborsList (markerIDListSigned_Type const& flags);

    //! Create neighbors to a given point, with a specified number of generations
    /*!
     * @param globalID. ID of the point to be examined.
     * @param nCircles. Number of circles (generations) to consider.
     * @return the set of neighbors global IDs
     */
    neighbors_Type circleNeighbors ( UInt globalID, UInt nCircles = 1 );

    //! Create neighbors to a given point within a specified radius
    /*!
     * @param globalID. ID of the point to be examined.
     * @param radius. The value of the circle radius within which neighbors are included
     * @return the set of neighbors global IDs
     */
    neighbors_Type neighborsWithinRadius ( UInt globalID, Real radius );

    //! Create the list of edge neighbors to the points
    void createPointEdgeNeighborsList();

    //! Create the list of element neighbors to the points
    void createPointElementNeighborsList();

    //! Create an overlapped map on points
    /*! Create a map based on points, expanding it across suddomain interfaces
     * with overlap 1, using NeighborMarker.
     */
    map_Type& ghostMapOnPoints();

    //! Create an overlapped map on points
    /*! Create a map based on points, expanding it across suddomain interfaces
     * with generic overlap.
     *  @param overlap. Level of overlap between subdomains
     *  @return the overlapped map
     */
    map_Type& ghostMapOnPoints ( UInt overlap );

    //! Create an overlapped map on edges
    /*! Create a map based on edges, expanding it across suddomain interfaces
     * with generic overlap.
     *  @param overlap. Level of overlap between subdomains
     *  @return the overlapped map
     */
    map_Type& ghostMapOnEdges ( UInt overlap );

    //! Create an overlapped map on elements for Finite Volumes
    /*! Create a map based on elements, expanding it across suddomain interfaces.
     *  The elements added are only those that share a facet with the current subdomain.
     *  This type of map is typically used for Finite Volumes.
     *  @return the overlapped map
     */
    // ghostMapOnElementsCommonFacet
    map_Type& ghostMapOnElementsFV();

    //! Create an overlapped map on elements for Finite Elements
    /*! Create a map based on elements, expanding it across suddomain interfaces.
     *  The elements added are all those that share a point with the current subdomain.
     *  This type of map is typically used for Finite Elements.
     *  @param overlap. Level of overlap between subdomains
     */
    // ghostMapOnElementsCommonPoints
    map_Type& ghostMapOnElementsFE ( UInt overlap );

    //! Extend the subdomains graph of the given overlap.
    /*!
     * This method enriches each subdomain with the closest elements such that
     * the partitions have the required overlap.
     * \param elemGraph. The list of subdomain elements
     * \param entityPID. Info about proc ownership of each mesh entity.
     * \param overlap. Level of overlap between partitions.
     */
    void extendGraphFE ( graphPtr_Type elemGraph, idList_Type const& pointPID, UInt overlap );

    //! Extend the subdomains graph of the given overlap.
    /*!
     * This method enriches each subdomain with the closest elements such that
     * the partitions have the required overlap.
     * \param elemGraph. The list of subdomain elements
     * \param entityPID. Info about proc ownership of each mesh entity.
     * \param overlap. Level of overlap between partitions.
     */
    void extendGraphFE ( const vertexPartitionPtr_Type& elemGraph,
                         idList_Type const& pointPID,
                         UInt overlap,
                         UInt partIndex);

    //! showMe method
    void showMe ( bool const verbose = false, std::ostream& out = std::cout );

    //@}

protected:

    //! @name Ghost Maps
    //@{

    mapPtr_Type M_ghostMapOnPoints;
    mapPtr_Type M_ghostMapOnEdges;
    mapPtr_Type M_ghostMapOnElementsFV;
    mapPtr_Type M_ghostMapOnElementsFE;

    //@}

    //! @name Protected Members
    //@{

    meshPtr_Type M_fullMesh;
    meshPtr_Type M_localMesh;
    mapPtr_Type const M_map;
    commPtr_Type const M_comm;
    UInt const M_me;

    neighborList_Type M_pointPointNeighborsList;
    neighborList_Type M_pointEdgeNeighborsList;
    neighborList_Type M_pointElementNeighborsList;

    bool M_verbose;
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    std::ofstream M_debugOut;
#endif
    //@}
};

template <typename MeshType>
GhostHandler<MeshType>::GhostHandler ( commPtr_Type const& comm ) :
    M_fullMesh(),
    M_localMesh(),
    M_map(),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_pointPointNeighborsList(),
    M_pointEdgeNeighborsList(),
    M_pointElementNeighborsList(),
    M_verbose ( 0 )
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    , M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#endif
{
}

template <typename MeshType>
GhostHandler<MeshType>::GhostHandler ( meshPtr_Type fullMesh,
                                       meshPtr_Type localMesh,
                                       mapPtr_Type map,
                                       commPtr_Type const& comm ) :
    M_fullMesh ( fullMesh ),
    M_localMesh ( localMesh ),
    M_map ( map ),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_pointPointNeighborsList(),
    M_pointEdgeNeighborsList(),
    M_pointElementNeighborsList(),
    M_verbose ( 0 )
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    , M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#endif
{
}

template <typename MeshType>
GhostHandler<MeshType>::GhostHandler ( meshPtr_Type fullMesh,
                                       commPtr_Type const& comm ) :
    M_fullMesh ( fullMesh ),
    M_comm ( comm ),
    M_me ( comm->MyPID() ),
    M_pointPointNeighborsList(),
    M_pointEdgeNeighborsList(),
    M_pointElementNeighborsList(),
    M_verbose ( 0 )
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    , M_debugOut ( ( "gh." + ( comm->NumProc() > 1 ? boost::lexical_cast<std::string> ( M_me ) : "s" ) + ".out" ).c_str() )
#endif
{
}

template <typename MeshType>
void GhostHandler<MeshType>::setUpNeighbors ( NeighborType const neighborType )
{
    if ( (neighborType & POINT_NEIGHBORS) != 0 )
    {
        this->createPointPointNeighborsList();
    }
    if ( (neighborType & RIDGE_NEIGHBORS) != 0 )
    {
        this->createPointEdgeNeighborsList();
    }
    if ( (neighborType & ELEMENT_NEIGHBORS) != 0 )
    {
        this->createPointElementNeighborsList();
    }
}

template <typename MeshType>
void GhostHandler<MeshType>::release()
{
    M_fullMesh.reset();
    M_localMesh.reset();
}

template <typename MeshType>
void GhostHandler<MeshType>::clean ( NeighborType const neighborType )
{
    if ( (neighborType & POINT_NEIGHBORS) != 0 )
    {
        clearVector ( M_pointPointNeighborsList );
    }
    if ( (neighborType & RIDGE_NEIGHBORS) != 0 )
    {
        clearVector ( M_pointEdgeNeighborsList );
    }
    if ( (neighborType & ELEMENT_NEIGHBORS) != 0 )
    {
        clearVector ( M_pointElementNeighborsList );
    }
}

#ifdef HAVE_HDF5

namespace
{

void writeNeighborMap ( EpetraExt::HDF5& file, neighborList_Type& list, std::string const& name )
{
    // copy the map into vectors
    ASSERT ( list.size() > 0, "the map " + name + " is empty!" )
    std::vector<Int> offsets ( list.size() + 1 );
    std::vector<Int> values;
    Int sum ( 0 );
    for ( UInt i ( 0 ); i < list.size(); i++ )
    {
        sum += list[ i ].size();
        offsets[ i + 1 ] = sum;
        for ( neighbors_Type::const_iterator j ( list[ i ].begin() ); j != list[ i ].end(); ++j )
        {
            values.push_back ( *j );
        }
    }

    // Save the vectors into the file
    file.Write ( name, "offsetSize", static_cast<Int> ( offsets.size() ) );
    file.Write ( name, "offsets", H5T_NATIVE_INT, offsets.size(), &offsets[ 0 ] );
    file.Write ( name, "valueSize", static_cast<Int> ( values.size() ) );
    file.Write ( name, "values", H5T_NATIVE_INT, values.size(), &values[ 0 ] );

    //#ifdef LIFEV_GHOSTHANDLER_DEBUG
    //    std::cerr << name << std::endl;
    //    for ( UInt i ( 0 ); i < map.size(); i++ )
    //    {
    //        std::cerr << i << "> ";
    //        for ( neighborList_Type::const_iterator j ( map[ i ].begin() ); j != map[ i ].end(); ++j )
    //        {
    //            std::cerr << *j << " ";
    //        }
    //        std::cerr <<std::endl;
    //    }
    //#endif
}

void readNeighborMap ( EpetraExt::HDF5& file, neighborList_Type& list, std::string const& name )
{
    // Read the vectors from the file
    Int offsetSize;
    file.Read ( name, "offsetSize", offsetSize );
    std::vector<Int> offsets ( offsetSize );
    file.Read ( name, "offsets", H5T_NATIVE_INT, offsetSize, &offsets[ 0 ] );
    Int valueSize;
    file.Read ( name, "valueSize", valueSize );
    std::vector<Int> values ( valueSize );
    file.Read ( name, "values", H5T_NATIVE_INT, valueSize, &values[ 0 ] );

    // setup the map
    list.resize ( offsetSize - 1 );
    for ( Int i ( 0 ); i < offsetSize - 1; i++ )
    {
        for ( Int j ( offsets[ i ] ); j < offsets[ i + 1 ]; j++ )
        {
            list[ i ].insert ( values[ j ] );
        }
    }

    //#ifdef LIFEV_GHOSTHANDLER_DEBUG
    //    std::cerr << name << std::endl;
    //    for ( UInt i ( 0 ); i < map.size(); i++ )
    //    {
    //        std::cerr << i << "> ";
    //        for ( neighborList_Type::const_iterator j ( map[ i ].begin() ); j != map[ i ].end(); ++j )
    //        {
    //            std::cerr << *j << " ";
    //        }
    //        std::cerr <<std::endl;
    //    }
    //#endif
}

}

template <typename MeshType>
void GhostHandler<MeshType>::exportToHDF5 ( std::string const& fileName, bool const& truncate )
{
    EpetraExt::HDF5 HDF5 ( *M_comm );

    if ( truncate )
    {
        // Create and open the file / Truncate and open the file
        HDF5.Create ( ( fileName + ".h5" ).data() );
    }
    else
    {
        // Open an existing file without truncating it
        HDF5.Open ( ( fileName + ".h5" ).data() );
    }

    // Check if the file is created
    if ( !HDF5.IsOpen () )
    {
        std::cerr << "Unable to create " + fileName + ".h5";
        abort();
    }

    writeNeighborMap ( HDF5, M_pointPointNeighborsList, "pointPointNeighborsMap" );
    writeNeighborMap ( HDF5, M_pointEdgeNeighborsList, "pointEdgeNeighborsMap" );
    writeNeighborMap ( HDF5, M_pointElementNeighborsList, "pointElementNeighborsMap" );

    // Close the file
    HDF5.Close();
}

template <typename MeshType>
void GhostHandler<MeshType>::importFromHDF5 ( std::string const& fileName )
{
    EpetraExt::HDF5 HDF5 ( *M_comm );

    // Open an existing file
    HDF5.Open ( ( fileName + ".h5" ).data() );

    // Check if the file is created
    if ( !HDF5.IsOpen () )
    {
        std::cerr << "Unable to open " + fileName + ".h5";
        abort();
    }

    readNeighborMap ( HDF5, M_pointPointNeighborsList, "pointPointNeighborsMap" );
    readNeighborMap ( HDF5, M_pointEdgeNeighborsList, "pointEdgeNeighborsMap" );
    readNeighborMap ( HDF5, M_pointElementNeighborsList, "pointElementNeighborsMap" );

    // Close the file
    HDF5.Close();
}

#endif // HAVE_HDF5

//! this routine generates point neighbors for the given mesh
/*! the routine assumes that the mesh is not yet partitioned or reordered
 *  (i.e. the local id and the global id are the same).
 *  if this is not true the method should be changed to use a more
 *  expensive STL find on the mesh points to get the correct point that has
 *  the given global id or construct a globalToLocal map beforehand.
 */
template <typename MeshType>
void GhostHandler<MeshType>::createPointNeighbors()
{
    // @TODO: ASSERT_COMPILE_TIME that MeshType::pointMarker == NeighborMarker
    // this guarantees that the pointNeighbors structure is available.

    // generate point neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        // fill fullMesh points neighboring
        M_fullMesh->point ( id0 ).pointNeighbors().insert ( id1 );
        M_fullMesh->point ( id1 ).pointNeighbors().insert ( id0 );
    }

    // update localMesh points
    for ( UInt ip = 0; ip < M_localMesh->numPoints(); ip++ )
    {
        M_localMesh->point ( ip ).pointNeighbors() = M_fullMesh->point ( M_localMesh->point ( ip ).id() ).pointNeighbors();
    }
}

template <typename MeshType>
void GhostHandler<MeshType>::createPointPointNeighborsList()
{
    M_pointPointNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate point neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        M_pointPointNeighborsList[ id0 ].insert ( id1 );
        M_pointPointNeighborsList[ id1 ].insert ( id0 );
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "M_pointPointNeighborsList on proc " << M_me << std::endl;
    for ( UInt i = 0; i < M_pointPointNeighborsList.size(); i++ )
    {
        M_debugOut << i << ": ";
        for ( neighbors_Type::const_iterator it = M_pointPointNeighborsList[ i ].begin();
                it != M_pointPointNeighborsList[ i ].end(); ++it )
        {
            M_debugOut << *it << " ";
        }
        M_debugOut << std::endl;
    }
#endif
}

namespace
{

inline bool isInside ( markerID_Type const& pointMarker, std::vector<int> const& markerIDList )
{
    for ( UInt i = 0; i < markerIDList.size(); ++i)
        if ( pointMarker == markerIDList[i] )
        {
            return true;
        }
    return false;
}

}

template <typename MeshType>
void GhostHandler<MeshType>::createPointPointNeighborsList (markerIDListSigned_Type const& flags)
{
    M_pointPointNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate point neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        if ( isInside (M_fullMesh->edge ( ie ).point ( 1 ).markerID(), flags) )
        {
            M_pointPointNeighborsList[ id0 ].insert ( id1 );
        }

        if ( isInside (M_fullMesh->edge ( ie ).point ( 0 ).markerID(), flags) )
        {
            M_pointPointNeighborsList[ id1 ].insert ( id0 );
        }

        //if(M_fullMesh->edge( ie ).point( 1 ).markerID()==20 || M_fullMesh->edge( ie ).point( 1 ).markerID()==1)
        //    M_pointPointNeighborsList[ id0 ].insert( id1 );

        //if(M_fullMesh->edge( ie ).point( 0 ).markerID()==20 || M_fullMesh->edge( ie ).point( 0 ).markerID()==1)
        //    M_pointPointNeighborsList[ id1 ].insert( id0 );
    }
}

template <typename MeshType>
neighbors_Type GhostHandler<MeshType>::circleNeighbors ( UInt globalID, UInt nCircles )
{
    neighbors_Type neighbors;
    neighbors_Type newEntries;
    neighbors_Type newNeighbors;

    neighbors = this->pointPointNeighborsList() [ globalID ];

    for (UInt i = 0; i < nCircles - 1; ++i)
    {
        for (neighbors_Type::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            newEntries = this->pointPointNeighborsList() [*it];
            for (neighbors_Type::iterator ii = newEntries.begin(); ii != newEntries.end(); ++ii)
            {
                newNeighbors.insert (*ii);
            }
        }
        neighbors = newNeighbors;
    }

    return neighbors;
}


template <typename MeshType>
neighbors_Type GhostHandler<MeshType>::neighborsWithinRadius ( UInt globalID, Real radius )
{
    neighbors_Type neighbors;
    neighbors_Type newEntries;
    neighbors_Type newNeighbors;
    bool isInside = true;
    Real d = 0;
    UInt mysize = 0;

    neighbors = this->pointPointNeighborsList() [globalID];

    typename mesh_Type::point_Type const& p = M_fullMesh->point (globalID);

    while (isInside)
    {
        for (neighbors_Type::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            newEntries = this->pointPointNeighborsList() [*it];
            for (neighbors_Type::iterator ii = newEntries.begin(); ii != newEntries.end(); ++ii)
            {
                typename mesh_Type::point_Type const& n = M_fullMesh->point (*ii);
                d = std::sqrt ( ( n.x() - p.x() ) * ( n.x() - p.x() ) +
                                ( n.y() - p.y() ) * ( n.y() - p.y() ) +
                                ( n.z() - p.z() ) * ( n.z() - p.z() ) );
                if (d < radius)
                {
                    newNeighbors.insert (*ii);
                }
            }
        }

        neighbors = newNeighbors;

        if (neighbors.size() == mysize)
        {
            isInside = false;
        }

        mysize = neighbors.size();

    }

    return neighbors;

}

template <typename MeshType>
void GhostHandler<MeshType>::createPointEdgeNeighborsList()
{
    M_pointEdgeNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate point neighbors by watching edges
    // note: this can be based also on faces or volumes
    for ( UInt ie = 0; ie < M_fullMesh->numEdges(); ie++ )
    {
        ID id0 = M_fullMesh->edge ( ie ).point ( 0 ).id();
        ID id1 = M_fullMesh->edge ( ie ).point ( 1 ).id();

        ASSERT ( M_fullMesh->point ( id0 ).id() == id0 && M_fullMesh->point ( id1 ).id() == id1,
                 "the mesh has been reordered, the point must be found" );

        M_pointEdgeNeighborsList[ id0 ].insert ( ie );
        M_pointEdgeNeighborsList[ id1 ].insert ( ie );
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "M_pointEdgeNeighborsList on proc " << M_me << std::endl;
    for ( UInt i = 0; i < M_pointEdgeNeighborsList.size(); i++ )
    {
        M_debugOut << i << ": ";
        for ( neighbors_Type::const_iterator it = M_pointEdgeNeighborsList[ i ].begin();
                it != M_pointEdgeNeighborsList[ i ].end(); ++it )
        {
            M_debugOut << *it << " ";
        }
        M_debugOut << std::endl;
    }
#endif
}

template <typename MeshType>
void GhostHandler<MeshType>::createPointElementNeighborsList()
{
    M_pointElementNeighborsList.resize ( M_fullMesh->numGlobalPoints() );
    // generate element neighbors by cycling on elements
    for ( UInt ie = 0; ie < M_fullMesh->numElements(); ie++ )
    {
        ASSERT ( M_fullMesh->element ( ie ).id() == ie,
                 "the mesh has been reordered, the point must be found" );

        for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
        {
            ID id ( M_fullMesh->element ( ie ).point ( k ).id() );
            M_pointElementNeighborsList[ id ].insert ( ie );
        }
    }
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnPoints()
{
    // if the map has already been created, return it
    if ( M_ghostMapOnPoints )
    {
        return *M_ghostMapOnPoints;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnPoints()" << std::endl;
    }

    // check that the pointNeighbors have been created
    if ( M_localMesh->point ( 0 ).pointNeighbors().empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the pointNeighbors are empty, will be generated now" << std::endl;
        }
        this->createPointNeighbors();
    }

    // create map
    M_ghostMapOnPoints.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnPoints );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    // use a set to avoid duplicates
    std::set<Int> myGlobalElementsSet;

    // iterate on local mesh points
    // @todo: this can start from the repeated map and add only neighbors for SUBDOMAIN_INTERFACE marked points
    for ( UInt k = 0; k < M_localMesh->numPoints(); k++ )
    {
        // iterate on each point neighborhood
        for ( typename mesh_Type::PointMarker::neighborConstIterator_Type neighborIt = M_localMesh->point ( k ).pointNeighbors().begin();
                neighborIt != M_localMesh->point ( k ).pointNeighbors().end(); ++neighborIt )
        {
            myGlobalElementsSet.insert ( *neighborIt );
        }
    }

    std::vector<Int> myGlobalElements ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::mapPtr_Type repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    return *M_ghostMapOnPoints;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnPoints ( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnPoints )
    {
        return *M_ghostMapOnPoints;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnPoints( UInt )" << std::endl;
    }

    // check that the pointPointNeighborsMap has been created
    if ( M_pointPointNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the pointPointNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createPointPointNeighborsList();
    }

    // create map
    M_ghostMapOnPoints.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnPoints );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    Int*          pointer;
    std::set<Int> myGlobalElementsSet, myOriginalElementsSet;;

    // get all elements from the repeated map
    pointer = M_map->map ( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map->map ( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        myOriginalElementsSet.insert ( *pointer );
    }

    // todo: optimize this!!
    // 1: work only on the boundary
    // 2: copy back only if necessary
    // repeat on actual points to expand overlap
    for ( UInt i = 0; i < overlap; i++ )
    {
        // iterate on points adding all neighbors
        for ( std::set<Int>::const_iterator pointIt = myOriginalElementsSet.begin();
                pointIt != myOriginalElementsSet.end(); ++pointIt )
        {
            // iterate on each point neighborhood
            for ( neighbors_Type::const_iterator neighborIt = M_pointPointNeighborsList[ *pointIt ].begin();
                    neighborIt != M_pointPointNeighborsList[ *pointIt ].end(); ++neighborIt )
            {
                myGlobalElementsSet.insert ( *neighborIt );
            }
        }
        myOriginalElementsSet = myGlobalElementsSet;
    }

    std::vector<Int> myGlobalElements ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::mapPtr_Type repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    return *M_ghostMapOnPoints;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnEdges ( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnEdges )
    {
        return *M_ghostMapOnEdges;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnEdges()" << std::endl;
    }

    // check that the pointEdgeNeighborsMap has been created
    if ( M_pointEdgeNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the pointEdgeNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createPointEdgeNeighborsList();
    }

    // set up Unique (first) and Repeated edges based on the GHOST flag
    MapEpetraData mapData;
    mapData.unique.reserve   ( M_localMesh->numEdges() );
    mapData.repeated.reserve ( M_localMesh->numEdges() );
    // loop on local mesh edges
    for ( ID ie = 0; ie < M_localMesh->numEdges(); ie++ )
    {
        mapData.repeated.push_back ( M_localMesh->edge ( ie ).id() );
        if ( M_localMesh->edge ( ie ).isOwned() )
        {
            mapData.unique.push_back ( M_localMesh->edge ( ie ).id() );
        }
    }

    M_ghostMapOnEdges.reset ( new map_Type ( mapData, M_comm ) );

    if ( overlap > 0 )
    {
        map_Type& ghostMap ( *M_ghostMapOnEdges );

        std::set<Int> myGlobalElementsSet;
        std::set<Int> addedElementsSet;

        for (  UInt k ( 0 ); k < mapData.repeated.size(); k++ )
        {
            typename mesh_Type::edge_Type const& edge = M_fullMesh->edge ( mapData.repeated[ k ] );
            for ( UInt edgePoint = 0; edgePoint < mesh_Type::edge_Type::S_numPoints; edgePoint++ )
            {
                addedElementsSet.insert ( edge.point ( edgePoint ).id() );
            }
        }
        ( mapData.repeated.begin(), mapData.repeated.end() );


        // @todo: optimize this!!
        // 1: work only on the boundary
        // 2: copy back only if necessary
        // repeat on actual points to expand overlap
        for ( UInt i = 0; i < overlap; i++ )
        {
            // iterate on points adding all neighbors
            for ( std::set<Int>::const_iterator pointIt = addedElementsSet.begin();
                    pointIt != addedElementsSet.end(); ++pointIt )
            {
                // iterate on each point neighborhood
                for ( neighbors_Type::const_iterator neighborIt = M_pointEdgeNeighborsList[ *pointIt ].begin();
                        neighborIt != M_pointEdgeNeighborsList[ *pointIt ].end(); ++neighborIt )
                {
                    std::pair<std::set<Int>::iterator, bool> isInserted = myGlobalElementsSet.insert ( *neighborIt );
                    if ( isInserted.second )
                    {
                        typename mesh_Type::EdgeType const& edge = M_fullMesh->edge ( *neighborIt );
                        for ( UInt edgePoint = 0; edgePoint < mesh_Type::EdgeType::S_numPoints; edgePoint++ )
                        {
                            // TODO exclude already included points
                            //                            if ( edge.point( edgePoint ).id() != *pointIt )
                            addedElementsSet.insert ( edge.point ( edgePoint ).id() );
                        }
                    }
                }
            }
        }

        mapData.repeated.assign ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

        map_Type::mapPtr_Type repeatedMap ( new Epetra_Map ( -1,
                                                             mapData.repeated.size(),
                                                             &mapData.repeated[0],
                                                             0,
                                                             *M_comm ) );
        ghostMap.setMap ( repeatedMap, Repeated );
    }

    return *M_ghostMapOnEdges;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnElementsFV()
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsFV )
    {
        return *M_ghostMapOnElementsFV;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsFV()" << std::endl;
    }

    // create the map
    M_ghostMapOnElementsFV.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnElementsFV );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    Int*          pointer;
    std::set<Int> map;

    // get all elements from the repeated map
    pointer = M_map->map ( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map->map ( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        map.insert ( *pointer );
    }

    // add all facing elements
    typedef typename mesh_Type::facet_Type const* facetPtr_Type;
    std::vector<facetPtr_Type>
    facetsOnSubdInt = M_localMesh->facetList().extractElementsWithFlag (
                          EntityFlags::SUBDOMAIN_INTERFACE, &Flag::testOneSet );
    for ( typename std::vector<facetPtr_Type>::const_iterator facetIt = facetsOnSubdInt.begin();
            facetIt != facetsOnSubdInt.end(); ++facetIt )
    {
        map.insert ( (*facetIt)->secondAdjacentElementIdentity() );
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<Int> myGlobalElements ( map.begin(), map.end() );

    // generate map
    map_Type::mapPtr_Type repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    return *M_ghostMapOnElementsFV;
}

template <typename MeshType>
typename GhostHandler<MeshType>::map_Type& GhostHandler<MeshType>::ghostMapOnElementsFE ( UInt overlap )
{
    // if the map has already been created, return it
    if ( M_ghostMapOnElementsFE )
    {
        return *M_ghostMapOnElementsFE;
    }

    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsFE()" << std::endl;
    }

    // check that the pointElementNeighborsMap has been created
    if ( M_pointElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the pointElementNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createPointElementNeighborsList();
    }

    // create the map
    M_ghostMapOnElementsFE.reset ( new map_Type() );
    map_Type& ghostMap ( *M_ghostMapOnElementsFE );

    // use the same Unique map and comm of the original map
    ghostMap.setMap ( M_map->map ( Unique ), Unique );
    ghostMap.setComm ( M_comm );

    Int*          pointer;
    std::set<Int> myGlobalElementsSet;

    // get all elements from the repeated map
    pointer = M_map->map ( Repeated )->MyGlobalElements();
    for ( Int ii = 0; ii < M_map->map ( Repeated )->NumMyElements(); ++ii, ++pointer )
    {
        myGlobalElementsSet.insert ( *pointer );
    }

    // add all elements with a point on SUBDOMAIN_INTERFACE
    typedef typename mesh_Type::point_Type const* pointPtr_Type;
    std::vector<pointPtr_Type>
    pointsOnSubdInt = M_localMesh->pointList.extractElementsWithFlag (
                          EntityFlags::SUBDOMAIN_INTERFACE, &Flag::testOneSet );
    // must work on global IDs since added elements are not on localMesh
    std::vector<ID> pointIDOnSubdInt ( pointsOnSubdInt.size() );
    for ( UInt i = 0; i < pointsOnSubdInt.size(); i++)
    {
        pointIDOnSubdInt[i] = pointsOnSubdInt[i]->id();
    }

    std::vector<ID> addedPoints;

    for ( UInt n = 0; n < overlap; n++ )
    {
        for ( std::vector<ID>::const_iterator globalId = pointIDOnSubdInt.begin();
                globalId != pointIDOnSubdInt.end(); ++globalId )
        {
            // iterate on each point neighborhood
            for ( neighbors_Type::const_iterator neighborIt = M_pointElementNeighborsList[ *globalId ].begin();
                    neighborIt != M_pointElementNeighborsList[ *globalId ].end(); ++neighborIt )
            {
                std::pair<std::set<Int>::iterator, bool> isInserted = myGlobalElementsSet.insert ( *neighborIt );
                if ( isInserted.second )
                {
                    typename mesh_Type::element_Type const& elem = M_fullMesh->element ( *neighborIt );
                    for ( UInt elemPoint = 0; elemPoint < mesh_Type::element_Type::S_numPoints; elemPoint++ )
                    {
                        // TODO exclude already included points
                        addedPoints.push_back ( elem.point ( elemPoint ).id() );
                    }
                }
            }
        }
        if ( overlap > 1 )
        {
            pointIDOnSubdInt = addedPoints;
        }
    }

    // convert unique list to vector to assure continuity in memorization
    std::vector<Int> myGlobalElements ( myGlobalElementsSet.begin(), myGlobalElementsSet.end() );

    // generate map
    map_Type::mapPtr_Type repeatedMap ( new Epetra_Map ( -1, myGlobalElements.size(), &myGlobalElements[0], 0, *M_comm ) );
    ghostMap.setMap ( repeatedMap, Repeated );

    return *M_ghostMapOnElementsFE;
}

template <typename MeshType>
void GhostHandler<MeshType>::extendGraphFE ( graphPtr_Type elemGraph,
                                             idList_Type const& pointPID,
                                             UInt overlap )
{
    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP1( graph )" << std::endl;
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChronoManager<> timeMgr ( M_comm );
    LifeChrono timeNL;
    timeMgr.add ( "point-element ngbr list", &timeNL );
    timeNL.start();
#endif
    // check that the pointElementNeighborsMap has been created
    if ( M_pointElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the pointElementNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createPointElementNeighborsList();
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeNL.stop();
#endif

    std::vector<int>& myElems = (*elemGraph) [M_me];

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    // show own elements
    M_debugOut << "own elements on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myElems.size(); i++ )
    {
        M_debugOut << myElems[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timePG;
    timeMgr.add ( "point graph", &timePG );
    timePG.start();
#endif
    // generate graph of points
    // check: parallel algorithm seems to be faster for this
    graph_Type pointGraph ( M_comm->NumProc() );

    std::set<int> localPointsSet;
    for ( UInt e = 0; e < (*elemGraph) [ M_me ].size(); e++ )
    {
        // point block
        for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
        {
            const ID& pointID = M_fullMesh->element ( (*elemGraph) [ M_me ][ e ] ).point ( k ).id();
            localPointsSet.insert ( pointID );
        }
    }
    pointGraph[ M_me ].assign ( localPointsSet.begin(), localPointsSet.end() );

    std::vector<Int> pointGraphSize ( M_comm->NumProc(), -1 );
    pointGraphSize[ M_me ] = pointGraph[ M_me ].size();
    for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
    {
        M_comm->Broadcast ( &pointGraphSize[ p ], 1, p );
    }

    for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
    {
        pointGraph[ p ].resize ( pointGraphSize[ p ] );
    }

    // communicate other proc point graphs
    for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
    {
        M_comm->Broadcast ( &pointGraph[p][0], pointGraph[p].size(), p );
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timePG.stop();
#endif

    std::vector<int> const& myPoints = pointGraph[ M_me ];

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own points on proc " << M_me << std::endl;
    for ( UInt i = 0; i < myPoints.size(); i++ )
    {
        M_debugOut << myPoints[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeSI;
    timeMgr.add ( "find SUBD_INT", &timeSI );
    timeSI.start();
#endif
    // initialize a bool vector that tells if an element is in the current partition
    std::vector<bool> isInPartition ( M_fullMesh->numElements(), false );
    for ( UInt e = 0; e < myElems.size(); e++ )
    {
        isInPartition[ myElems[ e ] ] = true;
    }

    // find subdomain interface points
    std::set<Int> mySubdIntPoints;
    for ( UInt k = 0; k < myPoints.size(); k++ )
    {
        int const& currentPoint = myPoints[ k ];
        // mark as SUBD_INT point only if the point is owned by current process
        if ( pointPID[ currentPoint ] == M_me )
        {
            // check if all element neighbors are on this proc
            for ( neighbors_Type::const_iterator neighborIt = M_pointElementNeighborsList[ currentPoint ].begin();
                    neighborIt != M_pointElementNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                // add the point if a neighbor is missing
                if ( !isInPartition[ *neighborIt ] )
                {
                    mySubdIntPoints.insert ( currentPoint );
                    break;
                }
            }
        }
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeSI.stop();
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own SUBDOMAIN_INTERFACE points on proc " << M_me << std::endl;
    for ( std::set<Int>::const_iterator i = mySubdIntPoints.begin(); i != mySubdIntPoints.end(); ++i )
    {
        M_debugOut << *i << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeOP;
    timeMgr.add ( "overlapping points", &timeOP );
    timeOP.start();
#endif

    std::vector<int> workingPoints ( mySubdIntPoints.begin(), mySubdIntPoints.end() );
    std::set<int> newPoints;
    std::set<int> augmentedElemsSet ( myElems.begin(), myElems.end() );

    for ( UInt o = 0; o < overlap; o++ )
    {

#ifdef LIFEV_GHOSTHANDLER_DEBUG
        M_debugOut << "workingPoints" << std::endl;
        for ( UInt i = 0; i < workingPoints.size(); i++ )
        {
            M_debugOut << workingPoints[ i ] << std::endl;
        }
#endif
        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            const int& currentPoint = workingPoints[ k ];
            // iterate on point neighborhood
            for ( neighbors_Type::const_iterator neighborIt = M_pointElementNeighborsList[ currentPoint ].begin();
                    neighborIt != M_pointElementNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                std::pair<std::set<Int>::iterator, bool> isInserted = augmentedElemsSet.insert ( *neighborIt );
                if ( isInserted.second )
                {
                    // if the element is inserted in the list, we add its points to the ones
                    // to be checked for next overlap value
                    for ( UInt j = 0; j < mesh_Type::element_Type::S_numPoints; j++ )
                    {
                        newPoints.insert ( M_fullMesh->element ( *neighborIt ).point ( j ).id() );
                    }
                }
            }
        }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
        M_debugOut << "augmentedElemsSet" << std::endl;
        for ( std::set<Int>::const_iterator i = augmentedElemsSet.begin(); i != augmentedElemsSet.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif
        // clean up newPoints from already analized points
        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            newPoints.erase ( newPoints.find ( workingPoints[ k ] ) );
        }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
        M_debugOut << "newPoints" << std::endl;
        for ( std::set<int>::const_iterator i = newPoints.begin(); i != newPoints.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif
        // set up workingPoints if we are not exiting
        if ( o + 1 < overlap  )
        {
            workingPoints.assign ( newPoints.begin(), newPoints.end() );
        }
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeOP.stop();
#endif

    // assign the augmentedElems to the element graph
    myElems.assign ( augmentedElemsSet.begin(), augmentedElemsSet.end() );

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeMgr.print();
#endif
}

template <typename MeshType>
void GhostHandler<MeshType>::extendGraphFE ( const vertexPartitionPtr_Type& elemGraph,
                                             idList_Type const& pointPID,
                                             UInt overlap,
                                             UInt partIndex)
{
    if ( M_verbose )
    {
        std::cout << " GH- ghostMapOnElementsP1( graph )" << std::endl;
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChronoManager<> timeMgr ( M_comm );
    LifeChrono timeNL;
    timeMgr.add ( "node-element ngbr list", &timeNL );
    timeNL.start();
#endif
    // check that the pointElementNeighborsMap has been created
    if ( M_pointElementNeighborsList.empty()  )
    {
        if ( M_verbose )
        {
            std::cerr << "the pointElementNeighborsList is empty, will be generated now" << std::endl;
        }
        this->createPointElementNeighborsList();
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeNL.stop();
#endif

    std::vector<int>& myElems = * (elemGraph->at (partIndex) );

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    // show own elements
    M_debugOut << "own elements on proc " << partIndex << std::endl;
    for ( UInt i = 0; i < myElems.size(); i++ )
    {
        M_debugOut << myElems[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timePG;
    timeMgr.add ( "point graph", &timePG );
    timePG.start();
#endif
    // generate graph of points
    // check: parallel algorithm seems to be faster for this
    Int numParts = elemGraph->size();
    graph_Type pointGraph ( numParts );

    std::set<int> localPointsSet;
    for ( UInt e = 0; e < elemGraph->at (partIndex)->size(); e++ )
    {
        // point block
        for ( UInt k = 0; k < mesh_Type::element_Type::S_numPoints; k++ )
        {
            const ID& pointID = M_fullMesh->element ( elemGraph->at (partIndex)->at (e) ).point ( k ).id();
            localPointsSet.insert ( pointID );
        }
    }
    pointGraph[ partIndex ].assign ( localPointsSet.begin(), localPointsSet.end() );

    std::vector<Int> pointGraphSize ( numParts, -1 );
    pointGraphSize[ partIndex ] = pointGraph[ partIndex ].size();
    if (M_comm->NumProc() > 1)
    {
        for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
        {
            M_comm->Broadcast ( &pointGraphSize[ p ], 1, p );
        }

        for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
        {
            pointGraph[ p ].resize ( pointGraphSize[ p ] );
        }

        // communicate other proc point graphs
        for ( UInt p = 0; p < static_cast<UInt> ( M_comm->NumProc() ); p++ )
        {
            M_comm->Broadcast ( &pointGraph[p][0], pointGraph[p].size(), p );
        }
    }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timePG.stop();
#endif

    std::vector<int> const& myPoints = pointGraph[ partIndex ];

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own points on proc " << partIndex << std::endl;
    for ( UInt i = 0; i < myPoints.size(); i++ )
    {
        M_debugOut << myPoints[ i ] << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeSI;
    timeMgr.add ( "find SUBD_INT", &timeSI );
    timeSI.start();
#endif
    // initialize a bool vector that tells if an element is in the current partition
    std::vector<bool> isInPartition ( M_fullMesh->numElements(), false );
    for ( UInt e = 0; e < myElems.size(); e++ )
    {
        isInPartition[ myElems[ e ] ] = true;
    }

    // find subdomain interface points
    std::set<Int> mySubdIntPoints;
    for ( UInt k = 0; k < myPoints.size(); k++ )
    {
        int const& currentPoint = myPoints[ k ];
        // mark as SUBD_INT point only if the point is owned by current process
        if ( pointPID[ currentPoint ] == partIndex )
        {
            // check if all element neighbors are on this proc
            for ( neighbors_Type::const_iterator neighborIt = M_pointElementNeighborsList[ currentPoint ].begin();
                    neighborIt != M_pointElementNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                // add the point if a neighbor is missing
                if ( !isInPartition[ *neighborIt ] )
                {
                    mySubdIntPoints.insert ( currentPoint );
                    break;
                }
            }
        }
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeSI.stop();
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    M_debugOut << "own SUBDOMAIN_INTERFACE points on proc " << partIndex << std::endl;
    for ( std::set<Int>::const_iterator i = mySubdIntPoints.begin(); i != mySubdIntPoints.end(); ++i )
    {
        M_debugOut << *i << std::endl;
    }
#endif

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    LifeChrono timeOP;
    timeMgr.add ( "overlapping points", &timeOP );
    timeOP.start();
#endif

    std::vector<int> workingPoints ( mySubdIntPoints.begin(), mySubdIntPoints.end() );
    std::set<int> newPoints;
    std::set<int> augmentedElemsSet ( myElems.begin(), myElems.end() );

    for ( UInt o = 0; o < overlap; o++ )
    {

#ifdef LIFEV_GHOSTHANDLER_DEBUG
        M_debugOut << "workingPoints" << std::endl;
        for ( UInt i = 0; i < workingPoints.size(); i++ )
        {
            M_debugOut << workingPoints[ i ] << std::endl;
        }
#endif
        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            const int& currentPoint = workingPoints[ k ];
            // iterate on point neighborhood
            for ( neighbors_Type::const_iterator neighborIt = M_pointElementNeighborsList[ currentPoint ].begin();
                    neighborIt != M_pointElementNeighborsList[ currentPoint ].end(); ++neighborIt )
            {
                std::pair<std::set<Int>::iterator, bool> isInserted = augmentedElemsSet.insert ( *neighborIt );
                if ( isInserted.second )
                {
                    // if the element is inserted in the list, we add its points to the ones
                    // to be checked for next overlap value
                    for ( UInt j = 0; j < mesh_Type::element_Type::S_numPoints; j++ )
                    {
                        newPoints.insert ( M_fullMesh->element ( *neighborIt ).point ( j ).id() );
                    }
                }
            }
        }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
        M_debugOut << "augmentedElemsSet" << std::endl;
        for ( std::set<Int>::const_iterator i = augmentedElemsSet.begin(); i != augmentedElemsSet.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif
        // clean up newPoints from already analized points
        for ( UInt k = 0; k < workingPoints.size(); k++ )
        {
            newPoints.erase ( newPoints.find ( workingPoints[ k ] ) );
        }

#ifdef LIFEV_GHOSTHANDLER_DEBUG
        M_debugOut << "newPoints" << std::endl;
        for ( std::set<int>::const_iterator i = newPoints.begin(); i != newPoints.end(); ++i )
        {
            M_debugOut << *i << std::endl;
        }
#endif
        // set up workingPoints if we are not exiting
        if ( o + 1 < overlap  )
        {
            workingPoints.assign ( newPoints.begin(), newPoints.end() );
        }
    }
#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeOP.stop();
#endif

    // assign the augmentedElems to the element graph
    myElems.assign ( augmentedElemsSet.begin(), augmentedElemsSet.end() );

#ifdef LIFEV_GHOSTHANDLER_DEBUG
    timeMgr.print();
#endif
}

template <typename MeshType>
void GhostHandler<MeshType>::showMe ( bool const /*verbose*/, std::ostream& out )
{
    out << "GhostHandler::showMe()" << std::endl;
    out << "M_pointPointNeighborsMap" << std::endl;
    for ( UInt i = 0; i < M_pointPointNeighborsList.size(); i++ )
    {
        out << i << " > ";
        for ( neighbors_Type::const_iterator nIt = M_pointPointNeighborsList[ i ].begin();
                nIt != M_pointPointNeighborsList[ i ].end(); ++nIt )
        {
            out << *nIt << " ";
        }
        out << std::endl;
    }
    out << "M_pointElementNeighborsMap" << std::endl;
    for ( UInt i = 0; i < M_pointPointNeighborsList.size(); i++ )
    {
        out << i << " > ";
        for ( neighbors_Type::const_iterator nIt = M_pointPointNeighborsList[ i ].begin();
                nIt != M_pointPointNeighborsList[ i ].end(); ++nIt )
        {
            out << *nIt << " ";
        }
        out << std::endl;
    }
}


}

#endif /* _GHOSTHANDLER_HPP_ */
