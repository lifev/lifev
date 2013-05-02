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
 *  @file
 *  @brief File containing Mesh Classes
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Miguel Fernandez
 *
 *  @contributor Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 *  @contributor Antonio Cervone <ant.cervone@gmail.com>
 *  @contributor Mauro Perego <mperego@fsu.edu>
 *  @mantainer Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 */

#ifndef _REGIONMESH_HH_
#define _REGIONMESH_HH_

#include <cstdlib>
#include <iomanip>
#include <fstream>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeDebug.hpp>
#include <lifev/core/mesh/MeshElementMarked.hpp>
#include <lifev/core/util/Switch.hpp>
#include <lifev/core/mesh/MeshElementBare.hpp>

#include <lifev/core/mesh/MeshEntityContainer.hpp>
#include <lifev/core/array/ArraySimple.hpp>
#include <lifev/core/mesh/ElementShapes.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace LifeV
{

/**
 *  @class RegionMesh
 *  @brief Class for 3D, 2D and 1D Mesh
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Miguel Fernandez
 *  @contributor Mauro Perego <mperego@fsu.edu>
 *
 *  This is the class that stores the mesh entities for a single region.
 *
 *  In a region elements are all of the same type.
 *
 *  Note: to provide data useful in a parallel setting some methods
 *  return either the number of entities in the current mesh region or
 *  the ones in the global mesh, before partitioning. The latter are identified
 *  by the keyword Global in the name, e.g. numGlobalFaces() versus numFaces()
 */
template <typename GeoShapeType, typename MCType = defaultMarkerCommon_Type >
class RegionMesh
    :
public MeshEntity
{
public:

    typedef GeoShapeType geoShape_Type;

    //! Common Markers
    typedef MCType markerCommon_Type;

    //! static const S_geoDimensions: the dimensions (1,2,3) of the geometry
    static const UInt S_geoDimensions = geoShape_Type::S_nDimensions;


    /** @name GeoDim Types
     *  @ingroup public_types
     *  dummy types that allows to select the correct method based on the dimension of the geometry.
     *
     *  @{
     */
    typedef GeoDim<1> oneD_Type;
    typedef GeoDim<2> twoD_Type;
    typedef GeoDim<3> threeD_Type;
    typedef GeoDim<S_geoDimensions> geoDim_Type;
    /** @} */ // End of group GeoDim Types

    /** @name Marker Types
     *  @ingroup public_types
     *  Markers for Point, Edge, Face, Volume and Region.
     *
     *  @{
     */

    //! Point Marker
    typedef typename markerCommon_Type::pointMarker_Type pointMarker_Type;
    //! Edge Marker
    typedef typename markerCommon_Type::edgeMarker_Type edgeMarker_Type;
    //! Face Marker
    typedef typename markerCommon_Type::faceMarker_Type faceMarker_Type;
    //! Volume Marker
    typedef typename markerCommon_Type::volumeMarker_Type volumeMarker_Type;
    //! Region Marker
    typedef typename markerCommon_Type::regionMarker_Type regionMarker_Type;
    /** @} */ // End of group Marker Types


    /** @name Geometric Element Types
     *  @ingroup public_types
     *  volume_Type, face_Type, edge_Type, point_Type, element_Type, facet_Type, ridge_Type, peak_Type.
     *  Correspondences table:
     *  @table
     *  |       | element_Type | facet_Type | ridge_Type | peak_Type  |
     *  --------------------------------------------------------------|
     *  |  3D   | volume_Type  | face_Type  | edge_Type  | point_Type |
     *  |  2D   | face_Type    | edge_Type  | point_Type |     -      |
     *  |  1D   | edge_Type    | point_Type |     -      |     -      |
     *  @endtable
     *  @{
     */


    typedef MeshElementMarked<3, S_geoDimensions, geoShape_Type, markerCommon_Type>  volume_Type;
    typedef MeshElementMarked<2, S_geoDimensions, geoShape_Type, markerCommon_Type> face_Type;
    typedef MeshElementMarked<1, S_geoDimensions, geoShape_Type, markerCommon_Type> edge_Type;
    typedef MeshElementMarked<0, S_geoDimensions, geoShape_Type, markerCommon_Type> point_Type;

    typedef MeshElementMarked<S_geoDimensions, S_geoDimensions, geoShape_Type, markerCommon_Type>  element_Type;
    typedef MeshElementMarked < S_geoDimensions - 1, S_geoDimensions, geoShape_Type, markerCommon_Type > facet_Type;
    typedef MeshElementMarked < S_geoDimensions - 2, S_geoDimensions, geoShape_Type, markerCommon_Type > ridge_Type;
    typedef MeshElementMarked < S_geoDimensions - 3, S_geoDimensions, geoShape_Type, markerCommon_Type >  peak_Type;

    /** @} */ // End of group Geometric Element Types


    /** @name Basic Element Shape Types
     *  @ingroup public_types
     *  volumeShape_Type, faceShape_Type, edgeShape_Type, facetShape_Type, ridgeShape_Type, ridgeShape_Type.
     *  @{
     */

    //! Element Shape.
    typedef typename volume_Type::geoShape_Type volumeShape_Type;
    typedef typename face_Type::geoShape_Type faceShape_Type;
    typedef typename edge_Type::geoShape_Type edgeShape_Type;

    typedef geoShape_Type elementShape_Type;
    typedef typename geoShape_Type::GeoBShape facetShape_Type;
    typedef typename facetShape_Type::GeoBShape ridgeShape_Type;
    /** @} */ // End of group Basic Element Shape Types

    /** @name Geometric Element Container Types
     *  @ingroup public_types
     *  Typedefs for STL complaint containers of mesh geometric entities.
     *  @{
     */
    typedef MeshEntityContainer<volume_Type > volumes_Type;
    typedef MeshEntityContainer<face_Type>    faces_Type;
    typedef MeshEntityContainer<edge_Type>    edges_Type;
    typedef MeshEntityContainer<point_Type>   points_Type;

    typedef MeshEntityContainer<element_Type>  elements_Type;
    typedef MeshEntityContainer<facet_Type>    facets_Type;
    typedef MeshEntityContainer<ridge_Type>    ridges_Type;
    typedef MeshEntityContainer<peak_Type>   peaks_Type;

    /** @} */ // End of group Geometric Element Container Types

    typedef boost::shared_ptr<Epetra_Comm> commPtr_Type;

    /** @name Constructors & Destructor
     *  Default and Copy Constructor for the class.
     *  @{
     */

    //! Default constructor
    RegionMesh();

    //! Constructor
    /**
     *  @param comm communicator to manage output
     */
    explicit RegionMesh ( commPtr_Type const& comm );

    //! Constructor
    /**
     *  @param comm communicator to manage output
     *  @param id markerId of the RegionMesh
     */
    RegionMesh ( UInt id, commPtr_Type const& comm );

    //! Destructor
    virtual ~RegionMesh();

    /** @} */ // End of group Constructors & Destructor


    /** @defgroup public_methods Public Methods
     *
     */

    /** @name Utilities
     *  @ingroup utilities
     *  Utilities for mesh checking and debugging.
     *
     *  @{
     */

    //! Display general information about the content of the class.
    /**
     *  Pretty output of information about the mesh.
     *
     *  @param verbose If true output is verbose, false otherwise (default);
     *  @param out Output stream (std::cout default);
     *  @return Output stream (for concatenation).
     */
    std::ostream& showMe ( bool verbose = false, std::ostream& out = std::cout ) const;

    //! Basic tests for mesh consistency.
    /**
     *  Check consistency of the mesh and fix errors.
     *
     *  @param level Indicates level.
     *  - 0 => minimal level, just sets the switches to the appropriate value
     *  - 1 => launch an external function from mesh_util.h which performs a series of
     *   sophisticated tests.
     *  @param fix If true and level=1 it fixes the mesh.
     *  @param verbose If true output is verbose, false otherwise (default);
     *  @param out Output stream (std::cerr default);
     *  @return Severity level. If different from 0 the mesh has problem. If positive the problems are such that the mesh may not
     *  work in any case. If less than zero it may work in some cases.
     */
    Int check ( Int level = 0, bool const fix = false, bool verbose = true, std::ostream& out = std::cerr );

    //! Display local to global mapping.
    /**
     *  @param os Output stream.
     */
    void printLtGMap (std::ostream& os);

    //! Return the handle to perform transormations on the mesh
    inline MeshUtility::MeshTransformer<RegionMesh<geoShape_Type, markerCommon_Type>, markerCommon_Type >& meshTransformer();

    //! Return the communicator
    commPtr_Type comm() const;

    //! Setter for the communicator
    void setComm ( commPtr_Type const& comm );

    /** @} */ // End of group Utilities

    /** @name Switches Methods
     *  @ingroup public_methods
     *  Switches are used to store the status of the RegionMesh The switches
     *  are used internally to control whether some data structures have been
     *  set up.
     *
     *  The possible Switches are:
     *  - \c HAS_ALL_FACETS
     *  - \c HAS_ALL_RIDGES
     *  - \c HAS_BOUNDARY_FACETS
     *  - \c HAS_BOUNDARY_RIDGES
     *  - \c HAS_ELEMENT_TO_FACETS
     *  - \c HAS_ELEMENT_TO_RIDGES
     *  - \c HAS_BEEN_CHECKED
     *  - \c FACETS_HAVE_ADIACENCY
     *
     *  @{
     */

    //! Get the number of switch which have been set.
    /**
     *  @return The number of switches which have been set.
     */
    const UInt& numSwitches() const
    {
        return switches.size();
    }

    //! Interrogate Switch.
    /**
     *  @param _s Name of the switch.
     *  @return true if the name is in the switch, false otherwise.
     */
    bool getLinkSwitch ( std::string const& _s ) const;

    //! Set a switch using a given name
    /**
     *  @param _s Name of the switch.
     */
    void setLinkSwitch ( std::string const& _s );

    //! Unset a switch.
    /**
     *  @param _s Name of the switch.
     */
    void unsetLinkSwitch ( std::string const& _s );
    //! Output switches contents
    /**
     *  @param verbose Verbosity of the output.
     *  @param out Output stream.
     *  @return Output stream for concatenation.
     */
    std::ostream& showLinkSwitch ( bool verbose = false, std::ostream& out = std::cout )
    {
        return switches.showMe ( verbose, out );
    }

    /** @} */ // End of group Switch Methods

    /** @name Generic Methods
     *  @ingroup public_methods
     *
     *  These are the generic methods to get information about the number of
     *  elements.
     *
     *  It is a generic interface common for all RegionMeshes (3D -- 1D).
     *
     *  @{
     */

    //! get mesh marker.
    /**
     *  @return regionmesh marker.
     */
    typename markerCommon_Type::regionMarker_Type& markerClass() const
    {
        return M_marker;
    }

    //! get only mesh marker id.
    /**
     *  @return regionmesh marker id.
     */
    markerID_Type markerID() const
    {
        return M_marker.markerID();
    }

    //! get M_isPartitioned bool
    /**
     * @brief isPartitioned
     * @return M_isPartitioned member
     */
    bool isPartitioned() const
    {
        return M_isPartitioned;
    }

    /** @} */ // End of group Generic Methods


    /** @name Volume Methods
     *  @anchor volume_methods
     *  @ingroup public_methods
     *
     *  Methods which operates on 3D elements.
     *
     *  There are different way of counting gemetry entities:
     *  - <strong>Counter</strong>: A counter stored the number of entities (element, faces
     *    etc) in the mesh. It does NOT necessarily correspond to the numer of entities
     *    actually stored. Indeed, I may not stor edges, yet having a counter with
     *    the number of edges in the mesh;
     *
     *  - <strong>Stored</strong>: The number of stored entities corresponds to the size()
     *    of the container.  A container may store all mesh entities (this is the
     *    default for the Points and the Volumes), only the boundary entities
     *    (this is the default with the Faces) or even none (this is the default
     *    with the edges);
     *
     *  - <strong>Capacity</strong>: The maximum number of entities that may stored
     *    before the container is resized.
     *
     *  @{
     */

    //! Returns Number of Volumes.
    /**
     *  Returns number of Volume elements in the mesh as given by the internal counter.
     *  @return Number of Volumes.
     */
    UInt numVolumes() const
    {
        return M_numVolumes;
    }

    //! Returns Global Number of Volumes
    /**
     *  Returns number of Global Volume elements in the mesh as given by the internal counter.
     *  @return Number of Global Volumes.
     */
    UInt numGlobalVolumes() const
    {
        return M_numGlobalVolumes;
    }

    //! Volumes actually stored in list.
    /**
     *  @return Number of stored Volumes.
     */
    UInt storedVolumes() const
    {
        return volumeList.size();
    }

    //! Current capacity of Volumes Container.
    /**
     *  @return how many volumes may be stored.
     */
    UInt maxNumVolumes() const
    {
        return maxNumItems ( volumeList );
    }


    //! set mesh marker.
    /**
     *  @param marker to be set.
     */
    void setMarkerClass ( typename MCType::regionMarker_Type const& marker )
    {
        M_marker = marker;
    }


    //! set only the mesh marker id.
    /**
     *  @param marker id to be set.
     */
    void setMarkerID ( markerID_Type const& markerId )
    {
        M_marker.setMarkerID ( markerId );
    }

    void setIsPartitioned ( bool const isPartitioned )
    {
        M_isPartitioned = isPartitioned;
    }

    //! Changes Current capacity of Volumes.
    /**
     *  Changes Current capacity of Volumes (Optionally sets internal counter).
     *
     *  @param n Maximum number of volumes.
     *  @param setcounter true to set the counter, false otherwise (default).
     */
    void setMaxNumVolumes      ( UInt const n, bool const setcounter = false );

    //! Set the number of global volumes.
    /**
     *  Set the number of global volumes.
     *
     *  @param n maximum number of global volumes.
     */
    void setMaxNumGlobalVolumes ( UInt const n )
    {
        M_numGlobalVolumes = n;
    }

    //! Set number of volumes.
    /**
     *  Set number of volumes in the mesh by changing internal counter.
     *  @param n number of volumes.
     */
    void setNumVolumes      ( UInt const n )
    {
        M_numVolumes = n;
    }

    //! Adds volume.
    /**
     *  Adds volume. Local and global ID computed automatically.
     *  @return Reference to added volume.
     */
    element_Type& addVolume();

    //! Adds volume
    /**
     *  Adds volume. Local ID computed automatically.
     *  Global ID is left unchanged
     *  @param v Volume to be added.
     *  @return Reference to the newly added volume.
     */
    element_Type& addVolume ( element_Type const& vol );

    //! Adds volume in a certain position.
    /**
     *  Adds volume to a specified position.
     *  @param v Volume to be added.
     *  @param pos Position of the volume.
     *  @return Reference to the newly added volume.
     */
    element_Type& setVolume ( element_Type const& v, UInt const pos );

    //! set numVolumes counter.
    void setVolumeCounter()
    {
        M_numVolumes = volumeList.size();
    }

    //! Reference to last volume stored in list.
    /**
     *  Reference to last volume stored in list.
     *  Useful for mesh readers.
     *  @return reference of the last volume in the list.
     */
    element_Type& lastVolume()
    {
        return volumeList.back();
    }

    //! i-th mesh 3D Element.
    /**
     *  @param i index of the mesh 3D Element.
     *  @return the i-th volume.
     */
    element_Type const& volume ( UInt const i ) const;

    //! i-th mesh 3D Element.
    /**
     *  @param i index of the mesh volume.
     *  @return reference to the ith mesh volume.
     */
    element_Type& volume ( UInt const i );

    /** @} */ // End of group Volume Methods


    /** @name Element Adjacency Methods
     *  @ingroup public_methods
     *  Methods to obtain the ID of Face and Edge belonging to an element.
     *
     *  Accessing this information requires that the appropriate data
     *  structures have been set by using the updateElementEdges() or
     *  updateElementFaces() methods.
     *
     *  It is NOT required to have the full information about edges and faces:
     *  The ID of the Face and Edge entities may be calculated without
     *  contructing the corresponding Edge of Face Object. This saves
     *  memory. However the  methods has the capability of
     *  building the actual internal Faces (the boundary GeoFaces must
     *  already exist!) of the Edges. In this case the MarkerFlag is inherited using the rules stated in the
     *  marker_traits class.
     *
     *  @{
     */

    //! Is the array for local facets set up?
    /**
     *  It does not use switches, but interrogates the container directly
     *
     *  @return true if the array for local facets in not empty.
     */
    bool hasLocalFacets() const
    {
        return ! M_ElemToFacet.empty();
    }

    //! Is the array for local faces set up?
    /**
     *  It does not use switches, but interrogates the container directly
     *
     *  @return true if array for local faces in not empty.
     */
    bool hasLocalFaces() const
    {
        return hasLocalFaces ( twoD_Type() );
    }

    //! Build localFacetId table and optionally fills the list of Facets.
    /**
     *  @param createFacets is set true if we want also to create the actual list
     *  of internal facets. There is another utility (in mesh_util.hpp), which
     *  might be used for the same purpose if we want just to create the faces
     *  and not also the LocaFaceID table.
     *  @param verbose if true, output is verbose.
     *  @param estimateFacetNumber is a guess provided by the user of the total
     *  number of faces. It is relevant only when createFaces=true. Setting it
     *  to a proper value helps in reducing time and memory.
     *
     *  @note Facets are renumbered so that boundary facets are stored first
     *  @pre The routine assumes that the boundary facets are properly set, if not use the
     *  methods in #include MeshChecks.hpp
     *
     */
    void updateElementFacets ( bool createFaces = false, bool verbose = false, UInt estimateFacetNumber = 0 );



    //! Build localFaceId table and optionally fills the list of Faces.
    /**
     *  @param createFacet is set true if we want also to create the actual list
     *  of internal faces. There is another utility (in mesh_util.hpp), which
     *  might be used for the same purpose if we want just to create the faces
     *  and not also the LocaFaceID table.
     *  @param verbose if true, output is verbose.
     *  @param estimateFaceNumber is a guess provided by the user of the total
     *  number of faces. It is relevant only when createFaces=true. Setting it
     *  to a proper value helps in reducing time and memory.
     *
     *  @pre The routine assumes that the boundary facets are properly set, if not use the
     *  methods in #include MeshChecks.hpp
     *
     */
    void updateElementFaces ( bool createFaces = false, const bool verbose = false, UInt estimateFaceNumber = 0 );

    //! Destroys element-to-facet container. Useful to save memory!
    void cleanElementFacets();

    //! Destroys element-to-face container. Useful to save memory!
    void cleanElementFaces();

    //! Local Facet Id.
    /** @param elemId Id of element.
     *  @param locF local facet number 0 \< LocF \< numLocalFacets().
     *  @return ID of the facet.
     */
    UInt localFacetId ( UInt const elemId, UInt const locF ) const;

    //! Local Face Id.
    /** @param volId Id of volume.
     *  @param locF local face number 0 \< LocF \< numLocalFaces().
     *  @return ID of the face.
     */
    UInt localFaceId ( UInt const volId, UInt const locF ) const;

    //! Local Face Id.
    /** @param ivol Reference to a volume.
     *  @param locF local face number 0 \< LocF \< numLocalFaces().
     *  @return ID of the face.
     */
    UInt localFaceId ( const volume_Type& vol, UInt const locF ) const
    {
        return localFacetId ( vol, locF );
    }

    //! Is the array for ridges set up?
    /**
     *  It does not use switches, but interrogates the container directly.
     *
     *  @return true if ridges information are available, false otherwise.
     */
    bool hasLocalRidges() const
    {
        return ! M_ElemToRidge.empty();
    }

    //! Is the array for local Edges set up?
    /**
     *  It does not use switches, but interrogates the container directly.
     *
     *  @return true if edges information are available, false otherwise.
     */
    bool hasLocalEdges() const
    {
        return hasLocalEdges ( edge_Type() );
    }

    //! Build localRidgeId table and optionally fills the list of Ridges
    /**
     *
     * @param createRdges is set true if we want also to create the actual list
     *  of edges. There is another utility (MeshChecks.hpp), which
     *  might be used for the same purpose if we want just to create the faces
     *  and not also the LocalRidgeID table.
     *  @param verbose If true, output is verbose.
     *  @param estimateRidgeNumber is a guess provided by the user of the total
     *  number of ridges. It is relevant only when createFacets=true. Setting it
     *  to a proper value helps in reducing time and memory.
     *  @param renumber Relevant only if createFacets=true.It makes sure that boundary edges are first
     *  if set to false possibly existing edges are never moved
     *
     *  @note This method does not assume that boundary edges are stores, since
     *  this condition is NOT a a paradigm for a RegionMesh.
     */
    void updateElementRidges ( bool createRidges = false, const bool verbose = false,
                               UInt estimateRidgeNumber = 0, bool renumber = true)
    {
        updateElementRidges ( M_geoDim, createRidges, verbose, estimateRidgeNumber, renumber);
    }

    //! Builds localEdgeId table and optionally fills the list of Edges
    /**
     *
     * @param createEdges is set true if we want also to create the actual
     * edges. There is another utility (MeshChecks.hpp),
     * might be used for the same purpose if we want just to create the edges
     * and not also the LocalEdgeID table.
     *  @param verbose If true, output is verbose.
     *  @param estimateEdgeNumber is a guess provided by the user of the total
     *  number of edges. It is relevant only when createEdges=true. Setting it
     *  to a proper value helps in reducing time and memory.
     *  @param renumber Relevant only if createFaces=true.It makes sure that boundary edges are first
     *  if set to false possibly existing edges are never moved
     *
     *  @note This method does not assume that boundary edges are already stored, since
     *  this condition is NOT an invariant of a RegionMesh3D.
     */
    void updateElementEdges ( bool createEdges = false, const bool verbose = false,
                              UInt estimateEdgeNumber = 0, bool renumber = true)
    {
        updateElementEdges ( edge_Type(), createEdges, verbose, estimateEdgeNumber, renumber);
    }

    //! Destroys Ridge-To-Facet lookup table.
    void cleanElementRidges();

    //! Destroys edge To facet lookup table.
    void cleanElementEdges()
    {
        cleanElementEdges ( edge_Type() );
    }

    //! Local Ridge ID of a ridge in an element stored in the mesh
    /** @param elemId local ID of the element
     *  @param locR local ridge number (elemental) 0 \<= locR \< numLocalRidges().
     *  @return local ID of the ridge.
     */
    ID localRidgeId ( UInt const elemId, UInt const locR ) const
    {
        return localRidgeId ( M_geoDim, elemId, locR );
    }


    //! Local Ridge.
    /** @param elem Reference of the element.
     *  @param locR local ridge number.
     *  @return local ID of the ridge.
     */
    ID localRidgeId ( const element_Type& elem, UInt const locR ) const
    {
        return localRidgeId ( M_geoDim, elem.localId(), locR );
    }


    //! Local Edge.
    /** @param elemId Id of element.
     *  @param locE local edge number 0 \<= LocE \< numLocalEdges().
     *  @return local ID of the edge.
     */
    ID localEdgeId ( UInt const elemId, UInt const locE ) const
    {
        return localEdgeId ( M_geoDim, elemId, locE );
    }

    //! Local Edge (specialization for 3D geometries).
    /** @param elem Reference of the element.
     *  @param locE local edge number 0 \<= LocE \< numLocalEdges().
     *  @return local ID of the edge.
     */
    ID localEdgeId ( const volume_Type& elem, UInt const locE ) const
    {
        return localRidgeId ( elem, elem.localId(), locE );
    }


    //! Local Edge (specialization for 2D geometries).
    /** @param elem Reference of the element.
     *  @param locE local edge number 0 \<= LocE \< numLocalEdges().
     *  @return local ID of the edge.
     */
    ID localEdgeId ( const face_Type& elem, UInt const locE ) const
    {
        return localFacetId ( elem.localId(), locE );
    }
    /** @} */ // End of group Element Adjacency Methods


    /** @name Faces Methods
     *  @anchor face_methods
     *  @ingroup public_methods
     *
     *  Methods to access/create/modify faces data
     *
     *  There are different way of counting Faces:
     *  - <strong>Number of Faces</strong>: is the declared number of total faces in the mesh.
     *  It uses an internal counter.
     *  @warning A value different from zero does NOT imply that the
     *  faces are actually stored. This counter declares the number of Faces
     *  that the mash have not those that are actually stored in the list of
     *  faces!
     *  - <strong>Number of Boundary Faces</strong>: is the declared number of boundary faces
     *  in the mesh. It uses an internal counter.
     *  @warning A value different from zero does NOT imply that the
     *  boundary faces are actually stored.
     *  - <strong>Number of Stored Faces</strong>: It is the number of Faces actually stored
     *  on the face container.
     *  - <strong>Maximum number of stored faces</strong>: The number of faces that may stored
     *  before the container is resized.
     *  @note This parameter has to be set BEFORE inserting faces
     *  in the container if we want that pointer into the container maintains
     *  their validity. The container will also have a better performance.
     *
     *  See also \ref volume_methods.
     *
     *  @{
     */

    //! Returns Number of Faces
    /**
     *  Returns number of Face elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Number of Faces.
     */
    UInt numFaces() const;


    //! Returns Number of Facets
    /**
     *  Returns number of Facets in the mesh
     *  as given by the internal counter.
     *
     *  @return Number of Facets.
     */
    UInt numFacets() const
    {
        return numFacets ( M_geoDim );
    }


    //! Returns Global Number of Faces
    /**
     *  Returns global number of Face elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global Number of Faces.
     */
    UInt numGlobalFaces() const;

    //! Number of Faces actually stored.
    /**
     *  Faces actually stored in list
     *  @return Number of faces actually stored in list.
     */
    UInt storedFaces() const;

    //! Current capacity of Faces Container.
    /**
     * How many elements may be stored.
     * @return Number of maximum element that can be stored.
     */
    UInt maxNumFaces() const;

    //! Changes Current capacity of Faces.
    /**
     *  Changes Current capacity of Faces (Optionally sets internal counter).
     *
     *  @param n Maximum number of faces.
     *  @param setcounter true to set the counter, false otherwise.
     */
    void setMaxNumFaces ( UInt const n, bool const setcounter = false );

    //! Changes Current capacity of Global Faces.
    /**
     *  Changes Current capacity of Global Faces (Optionally sets internal counter).
     *
     *  @param n maximum number of global faces.
     */
    void setMaxNumGlobalFaces ( UInt const n );

    //! Changes Current capacity of Global Facets.
    /**
     *  Changes Current capacity of Global Facets (Optionally sets internal counter).
     *
     *  @param n maximum number of global facets.
     */
    void setMaxNumGlobalFacets ( UInt const n )
    {
        setMaxNumGlobalFacets ( M_geoDim, n );
    }

    //! Adds a face.
    /**
     *  Adds a face. Id computed automatically.
     *  @return Reference to added face.
     */
    face_Type& addFace()
    {
        return addFace ( face_Type() );
    }

    //! Adds a face.
    /**
     *  Adds a face (optionally a boundary face). Local and global ID
     *  are computed automatically.
     *  @param boundary true if it's a boundary face.
     *  @return Reference to added face.
     */
    face_Type& addFace ( bool const boundary );

    //! Adds a face.
    /**
     *  Adds a face (optionally a boundary face). Local ID computed automatically.
     *  It assumes that all attributes of face f have been properly set
     *  @param f Face to be added.
     *  @return Reference to the newly added face.
     */
    face_Type& addFace ( face_Type const& f );

    //! Adds a face in a certain position.
    /**
     *  Add face to a specified position (optionally a boundary face).
     *  It assumes that all attributes of face f have been properly set a part the id
     *  which is set to pos
     *  @param f Face to add.
     *  @param pos Position of the face.
     *  @return Reference to the newly added face.
     *  @note If you add a face on the boundary you may need to reorder the list of faces
     */
    face_Type& setFace ( face_Type const& f, UInt pos);

    //! Reference to last face stored in list.
    /**
     *  Reference to last face stored in list.
     *  Useful for mesh readers.
     *  @return reference of the last face in the list.
     */
    face_Type& lastFace();

    //! i-th mesh Face.
    /**
     *  @param i index of the mesh face.
     *  @return the i-th face.
     */
    face_Type const& face ( UInt const i ) const;

    //! i-th mesh face.
    /**
     *  @param i index of the mesh face.
     *  @return reference to the ith mesh face.
     */
    face_Type& face ( UInt const i );

    //! i-th mesh boundary face.
    /**
     *  @param i index of the mesh boundary face.
     *  @return the i-th face.
     */
    face_Type const& boundaryFace ( UInt const i ) const;

    //! i-th mesh boundary face.
    /**
     *  @param i index of the mesh boundary face.
     *  @return the i-th face.
     */
    face_Type& boundaryFace ( UInt const i );

    //! Set counter of faces.
    /**
     *  @param n Number of faces.
     */
    void setNumFaces ( UInt const n );

    //! Set counter of boundary faces.
    /**
     *  @param n Number of boundary faces.
     */
    void setNumBFaces ( UInt const n ) ;

    //! Do I store mesh faces?
    bool hasFaces() const;

    //! Do I store also internal faces?
    bool hasInternalFaces() const;

    //! Number of Boundary Faces.
    /**
     *  @return Number of boundary faces.
     */
    UInt numBFaces() const ;

    //! Is face whose id is given on boundary?
    /**
     *  @param id Face Id.
     *  @return true if f in on boundary.
     */
    bool isBoundaryFace ( UInt const& id  ) const;

    //! Number of local faces for each 3Delement.
    /**
     *  @return Number of local faces.
     */
    UInt numLocalFaces() const ;

    //! Does this id corresponds to a full face?
    /** A FULL FACE is a 2DElement  that is actually
     *  stored in the Face container.
     *  @param id Face id.
     *  @return false if id does not corresponfd to a boundary face and internal faces are not stored.
     */
    bool isFullFace ( UInt const& id ) const ;

    //! Id of the Volume Element adjacent to a Face.
    /** The first element is the one <em>ORIENTED coherently with the
     *  face</em> (AS STORED in Faces). It means that the face orientation is
     *  OUTWARD with respect to the element. The second element is either null
     *  (boundary face) or indicates that the normal of the face appears INWARD
     *  with respect to that element.
     *
     *  @param faceId Id of the face.
     *  @param Pos is equal to 0 or 1 and indicates first or second element.
     *  @return Id of adjacent volume element or NotAnId if none.
     */
    UInt faceElement ( UInt const faceId, UInt const Pos ) const ;

    //! 3DElement adjacent to a FACE. Face reference given.
    /**
     *  @param f Face reference.
     *  @param Pos is equal to 0 or 1 and indicates first or second element.
     *  @return Id of adjacent volume element or 0 if none.
     */
    UInt faceElement ( facet_Type const& f, UInt const Pos ) const;

    /** @} */ // End of group Faces Methods


    /** @name Edge Methods
     *  @ingroup public_methods
     *  All methods which operates on 1D elements.
     *
     *  There are different point counters which may bi interrogated or set:
     *  - <strong>Number of Edges</strong>: is the declared number of total edges in the mesh.\n
     *   <strong>Beware</strong>: A value different from zero does NOT imply that the
     *   edges are actually stored.
     *  - <strong>Number of Boundary Edges</strong>: is the declared number of boundary edges
     *   in the mesh.\n
     *   <strong>Beware</strong>: A value different from zero does NOT imply that
     *   the boundary edges are actually stored.
     *  - <strong>Number of Stored Edges</strong>: It is the number of Edges actually stored
     *   on the edge container.
     *  - <strong>Maximum number of stored edges</strong>: The number of edges that may stored
     *   before the container is resized.\n
     *   <strong>Important</strong>: This parameter has to be
     *   set BEFORE inserting edges in the container if we want that pointer
     *   into the container maintains their validity. The container will also
     *   have a better performance.
     *
     *  @note To have more information on the Edge methods
     *  look at the documentation of th eanalogous Edges methods @ref face_methods and in @ref volume_methods.
     *
     *  @{
     */

    //! Number of Edges.
    /**
     *  Returns number of Edge elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Number of Edges.
     */
    UInt numEdges() const;

    //! Global number of Edges.
    /**
     *  Returns global number of Edge elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global number of Edges.
     */
    UInt numGlobalEdges() const;

    //! Number of local edges for each (3D) element.
    /**
     * @return Number of local edges for each (3D) element
     */
    UInt numLocalEdges() const;

    //! Number of stored Edges.
    /**
     *  Edges actually stored in list.
     *
     *  @return number of edges actually stored in list.
     */
    UInt storedEdges() const;

    //! Capacity of Edge Container.
    /**
     *  Current capacity of Edges Container,
     *  i.e. how many elements may be stored.
     *
     *  @return Edges Container Capacity.
     */
    UInt maxNumEdges() const;

    //! Changes Current capacity of Edges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of edges.
     *  @param setcounter true to set the counter, false otherwise.
     */
    void setMaxNumEdges ( UInt const n, bool const setcounter = false );

    //! Changes Current capacity of Global Edges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of global edges.
     */
    void setMaxNumGlobalEdges ( UInt const n );

    //! Adds an Edge.
    /**
     *  Adds an edge. Id computed automatically.
     *
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    edge_Type& addEdge ( bool const boundary );

    //! Adds an Edge.
    /**
     *  Adds an edge to the end of the list
     *  and adjourn its Id. The edge attributes (a part the id) should be
     *  correctly set.
     *
     *  @param f Edge to add.
     *  @return Reference to added edge.
     */
    edge_Type& addEdge ( edge_Type const& r );

    //! Add an Edge to specified position.
    /**
     *  Adds an edge and sets its Id to position.
     *
     *  @param e Edge to add.
     *  @param position Position of the edge.
     *  @return Reference to added edge.
     */
    edge_Type& setEdge ( edge_Type const& e, UInt position );

    //! Reference to last edge stored in list.
    /**
     *  Useful for mesh readers
     *  @return reference of the last edge in the list.
     */
    edge_Type& lastEdge();

    //! i-th mesh edge.
    /**
     *  Returns the i-th edge.
     *
     *  @param i Index of the mesh edge.
     *  @return The i-th edge.
     */
    edge_Type const& edge ( UInt const i ) const;

    //! i-th mesh edge reference.
    /**
     *  Returns a reference to the i-th mesh Edge.
     *
     *  @param i Index of the mesh 1D Edge.
     *  @return Reference to the i-th Edge.
     */
    edge_Type& edge ( UInt const i );

    //! i-th mesh 1D Boundary Edge.
    /**
     *  Returns the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return i-th Boundary Edge.
     */
    edge_Type const& boundaryEdge ( UInt const i ) const;

    //! i-th mesh 1D Boundary Edge Reference.
    /**
     *  Returns a reference to the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return Reference to the i-th Boundary edge.
     */
    edge_Type& boundaryEdge ( UInt const i );

    //! Set boundary Edge counter.
    /**
     *  Set the Boundary Edges counter to a given number.
     *
     *  @param n Count of Boundary Edge.
     */
    void setNumBEdges ( UInt const n );

    //! Do I store mesh edges?
    /**
     *  Returns true if edges are stored.
     *
     *  @return true if edge are stored, false otherwise.
     */
    bool hasEdges() const;

    //! Edge on boundary check.
    /**
     *  Is this edge on boundary?
     *
     *  @param e The Edge.
     *  @return true if the edge is on the boundary, false otherwise.
     */
    bool isBoundaryEdge ( ridge_Type const& e ) const
    {
        return e.boundary();
    }

    //! Edge on boundary check by id.
    /**
     *  Is this edge, of given id, on boundary?
     *
     *  @param id Id of the edge.
     *  @return true if the edge is on the boundary, false otherwise.
     */
    bool isBoundaryEdge ( UInt const& id ) const
    {
        return (id < M_numBEdges );
    }

    //! Full Edge check by id.
    /**
     *  Does this ID corresponds to a full edge?
     *
     *  A FULL EDGE is a 1D Element that is actually stored in the Edge container.
     *
     *  @param id Id of the edge.
     *  @return true if the edge is actually stored, false otherwise.
     */
    bool isFullEdge ( UInt const& id ) const;

    //! Number of Boundary Edges.
    /**
     *  Returns number of boundary edge elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Number of Boundary Edges.
     */
    UInt numBEdges() const;

    //! Set internal counter of number of Edges.
    /**
     *  @param n Number of Edges.
     */
    void setNumEdges ( UInt const n );

    //! Do I store internal Edges?
    /**
     *  @return true if internal edges are stored.
     */
    bool hasInternalEdges() const;

    //! Number of edges on each face.
    /**
     *  @return Number of edges for each face.
     */
    UInt numLocalEdgesOfFace() const ;

    /** @} */ // End of group Edges Methods


    /** @name Points Methods
     *  @ingroup public_methods
     *  All methods which operates on Point Elements.
     *
     *  There are different Point counters which may be interrogated or set:
     *  - <strong>Number of Points</strong>: is the declared number of total points in
     *   the mesh.
     *  - <strong>Number of Boundary Points</strong>: is the declared number of boundary
     *   points in the mesh.\n
     *   Does NOT imply that the boundary points are actually stored.
     *  - <strong>Number of Stored Points</strong>: It is the number of Points actually stored
     *   on the points container.
     *  - <strong>Maximum number of stored points</strong>: The number of points that may
     *   stored before the container is resized.\n
     *   <strong>Very Important</strong>: This parameter has to be set <strong>BEFORE</strong>
     *   inserting points in the container. Since MeshElementMarkeds will contain POINTERS
     *   into the container!\n
     *   This is a debatable point!
     *
     *   @note To have more information on the Point methods look at the
     *   documentation of the analogous Edges methods (\ref face_methods).
     *
     *  @{
     */

    //! Returns number of points in the mesh.
    /**
     *  Interrogation to the counter of points in the mesh.
     *
     *  @return Number of points in the mesh.
     */
    UInt numPoints() const;

    //! Returns a reference to number of points in the mesh.
    /**
     *  Returns a reference to the internal counter of the number of points in the mesh.
     *
     *  @return Reference to the internal counter of number of points in the mesh.
     */
    UInt& numPoints();

    //! Returns number of points in the mesh actually stored.
    /**
     *  Returns number of points in the mesh actually stored interrogating the internal counter.
     *
     *  @return Number of stored points in the mesh.
     */
    UInt storedPoints() const;

    //! Returns number of boundary points in the mesh actually stored.
    /**
     *  Returns number of boundary points in the mesh actually stored interrogating the internal counter.
     *
     *  @return Number of stored boundary points in the mesh.
     */
    UInt storedBPoints() const;

    //! Returns the number of storable points in the mesh.
    /**
     *  Returns the number of storable points in the mesh interrogating the internal counter.
     *
     *  @return The number of storable points.
     */
    UInt maxNumPoints() const;

    //! Set the number of storable points in the mesh.
    /**
     *  Set the internal counter of storable points in the mesh.
     *
     *  @param n Maximum number of storable points in the mesh.
     *  @param setcounter If true, it sets the internal counter, otherwise not (default).
     */
    void setMaxNumPoints ( UInt const n, bool const setcounter = false );

    //! Set the number of storable global points in the mesh.
    /**
     *  Set the internal counter of storable global points in the mesh.
     *
     *  @param n Maximum number of storable global points in the mesh.
     */
    void setMaxNumGlobalPoints ( UInt const n);

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    point_Type& addPoint ( bool const boundary, bool const vertices);

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh.
     *  Point will be add at the end of the list and the Id is computed automatically
     *
     *  @param p Point to be added.
     *  @return Reference to the newly added Point.
     */
    point_Type& addPoint ( point_Type const& p);

    //! Returns the first mesh Point.
    /**
     *  Returns the first Point in the mesh.
     *
     *  @return const reference to the first mesh Point.
     */
    point_Type const& firstPoint() const
    {
        return pointList.front();
    }

    //! Returns the last mesh Point.
    /**
     *  Returns the last Point in the mesh.
     *
     *  @return const reference to the last mesh Point.
     */
    point_Type const& lastPoint() const
    {
        return pointList.back();
    }

    //! Returns the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *
     *  @param i Id of the Point.
     *  @return i-th mesh Point.
     */
    point_Type const& point ( UInt const i ) const;

    //! Returns a reference to the i-th mesh point.
    /**
     *  Returns the i-th point in the mesh.
     *UInt const i
     *  @param i Id of the point.
     *  @return Reference i-th mesh point.
     */
    point_Type& point ( UInt const i );

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    point_Type const& boundaryPoint ( UInt const i ) const;

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    point_Type& boundaryPoint ( UInt const i );

    //! Returns the number of Boundary Points.
    /**
     *  Returns the counter of the number of Boundary Points.
     *
     *  @return Number Boundary Points.
     */
    UInt numBPoints() const;

    //! Sets the number of Boundary Points.
    /**
     *  Sets the counter of the number of Boundary Points.
     *
     *  @param n Number Boundary Points.
     */
    void setNumBPoints ( UInt const n );

    //! Is this point on boundary?
    /**
     *  Is this point on boundary?
     *
     *  @param p The Point.
     *  @return true if the point is on the boundary, false otherwise.
     */
    bool isBoundaryPoint ( point_Type const& p ) const ;

    //! Is this point on boundary?
    /**
     *  Is this point, of given id, on boundary?
     *
     *  @param id Id of the point.
     *  @return true if the point is on the boundary, false otherwise.
     */
    bool isBoundaryPoint ( UInt const& id ) const;

    //! List of points
    /**
     *  @param fct Function of three Real arguments.
     *  @param list_pts List of Points.
     *  @todo Move away, this can be done using the utility of the list of pts
     */
    void getListOfPoints ( bool ( *fct ) ( Real, Real, Real ), std::vector<UInt>& list_pts );

    /** @} */ // End of group Points Methods


    /** @name Vertices Methods
     *  @ingroup public_methods
     *  All methods which operates on Vertex Elements.
     *
     *  @{
     */

    //! Number of Vertices in Region.
    /**
     *  Returns the number of Vertices in Region.
     *
     *  @return Number of Vertices in Region.
     */
    UInt  numVertices () const;

    //! Number of Boundary Vertices in Region.
    /**
     *  Returns the number of Boundary Vertices in Region.
     *
     *  @return Number of Vertices in Region.
     */
    UInt numBVertices() const;

    //! Reference to the number of Boundary Vertices in Region.
    /**
     *  Returns a reference to the number of Boundary Vertices in Region.
     *
     *  @return Reference to the number of Vertices in Region.
     */
    UInt& numBVertices();

    //! Number of local vertices for each (3D) element.
    /**
     *  @return Number of local vertices for each (3D) element.
     */
    UInt numLocalVertices() const;

    //! Returns the global number of vertices in the mesh.
    /**
     *  Interrogation to the counter of vertices in the mesh.
     *
     *  @return Number of vertices in the mesh.
     */
    UInt numGlobalVertices() const;

    //! Returns the global number of points in the mesh.
    /**
     *  Interrogation to the counter of points in the mesh.
     *
     *  @return Number of points in the mesh.
     */
    UInt numGlobalPoints() const
    {
        return M_numGlobalPoints;
    }

    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param id Point's id.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex ( UInt const& id ) const;

    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param p A Point.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex ( point_Type const& p ) const;

    //! Changes number of Vertices.
    /**
     *  Allows to change number of vertices in Region.
     *
     *  @param n Number of vertices in the mesh.
     */
    void setNumVertices (UInt const n);

    //! Set the number of vertices in the mesh.
    /**
     *  Set the internal counter of vertices points in the mesh.
     *
     *  @param n Number of vertices in the mesh.
     */
    void setNumGlobalVertices ( UInt const n );

    //! Changes number of Boundary Vertices.
    /**
     *  Allows to change number of boundary vertices in Region.
     *
     *  @param n Number of boundary vertices in the mesh.
     */
    void setNumBVertices (UInt const n);

    /** @} */ // End of group Vertices Methods



    /** @name Polytope Methods
     *  @ingroup public_methods
     *
     *  These are the polytope methods to get information about the number of
     *  polytopes (elements, facets, ridges, peaks).
     *
     *  Interfaces common for all RegionMeshes (3D -- 1D).
     *
     *  @{
     */

    //! Number of elements in mesh.
    /**
     *  @return Number of elements in mesh.
     *  @sa numFaces
     *  @note Alias to numFaces()
     */
    UInt numElements() const
    {
        return numElements (M_geoDim);
    }

    //! Number of Boundary facets.
    /**
     *  @return Number of Boundary facets.
     */
    UInt numBoundaryFacets() const
    {
        return numBoundaryFacets ( M_geoDim);
    }

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    element_Type& element ( const UInt& i )
    {
        return element (M_geoDim, i );
    }

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    const element_Type& element ( const UInt& i ) const
    {
        return element (M_geoDim, i );
    }

    //! Get boundary facet at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary facet at index i
     */
    facet_Type& boundaryFacet ( const UInt& i )
    {
        return boundaryFacet (M_geoDim, i);
    }

    //! Get boundary facet at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary facet at index i
     */
    const facet_Type& boundaryFacet ( const UInt& i ) const
    {
        return boundaryFacet (M_geoDim, i);
    }

    //! Number of global elements.
    /**
     * @return Number of global Elements (Volumes).
     */
    UInt numGlobalElements() const
    {
        return numGlobalElements ( M_geoDim );
    }

    //! Current capacity of the container of Elements.
    /**
     *  @return how many elements may be stored.
     */
    UInt maxNumElements() const
    {
        return maxNumVolumes ( M_geoDim );
    }

    //! Set counter of elements.
    /**
     *  @param n Number of elements.
     */
    void setNumElements ( UInt const n )
    {
        setNumElements ( M_geoDim, n );
    }

    //! Changes Current capacity of the container of elements.
    /**
     *  Changes Current capacity of the container of elements (Optionally sets internal counter).
     *
     *  @param n Maximum number of elements.
     *  @param setcounter true to set the counter, false otherwise (default).
     */

    void setMaxNumElements   ( UInt const n, bool const setcounter = false )
    {
        setMaxNumElements (M_geoDim, n, setcounter);
    }

    //! Set the number of global elements.
    /**
     *  Set the number of global elements.
     *
     *  @param n maximum number of global elements.
     */
    void setMaxNumGlobalElements ( UInt const n )
    {
        setMaxNumGlobalElements (M_geoDim, n) ;
    }

    //! Adds element.
    /**
     *  Adds element. Id computed automatically.
     *  @return Reference to added element.
     */
    element_Type& addElement()
    {
        return addElement ( element_Type() );
    }

    //! Adds element in a certain position (specialization for 3D geometry).
    /**
     *  Adds element to a specified position.
     *  @param element element to be added.
     *  @param pos Position of the element.
     *  @return Reference to the newly added element.
     */
    element_Type& setElement ( volume_Type const& elem, UInt const pos )
    {
        return setVolume ( elem, pos );
    }

    //! Adds element in a certain position (specialization for 2D geometry).
    /**
     *  Adds element to a specified position.
     *  @param element element to be added.
     *  @param pos Position of the element.
     *  @return Reference to the newly added element.
     */
    element_Type& setElement ( face_Type const& elem, UInt const pos )
    {
        return setFace ( elem, pos );
    }

    //! Adds element in a certain position (specialization for 1D geometry).
    /**
     *  Adds element to a specified position.
     *  @param element element to be added.
     *  @param pos Position of the element.
     *  @return Reference to the newly added element.
     */
    element_Type& setElement ( edge_Type const& elem, UInt const pos )
    {
        return setEdge ( elem, pos );
    }

    //! Returns Global Number of Facets
    /**
     *  Returns global number of Facets in the mesh
     *  as given by the internal counter.
     *
     *  @return Global Number of Facets.
     */
    UInt numGlobalFacets() const
    {
        return numGlobalFacets ( M_geoDim );
    }

    //! Changes Current capacity of Facets.
    /**
     *  Changes Current capacity of Facets (Optionally sets internal counter).
     *
     *  @param n Maximum number of facets.
     *  @param setcounter true to set the counter, false otherwise.
     */
    void setMaxNumFacets ( UInt const n, bool const setcounter = false )
    {
        setMaxNumFacets ( M_geoDim, n, setcounter);
    }

    //! Adds a facet.
    /**
     *  Adds a facet (optionally a boundary facet). Id computed automatically.
     *  @param boundary true if it's a boundary facet.
     *  @return Reference to added facet.
     */
    facet_Type& addFacet ( bool const boundary )
    {
        return addFacet ( M_geoDim, boundary);
    }


    //! Adds a facet. 3D specialization
    /**
     *  Adds a facet (optionally a boundary face). Id computed automatically.
     *  It assumes that all attributes of facet f have been properly set
     *  @param f Facet to be added.
     *  @return Reference to the newly added facet.
     */
    facet_Type& addFacet ( face_Type const& facet )
    {
        return addFace (facet);
    }

    //! Adds a facet. 2D specialization
    /**
     *  Adds a facet (optionally a boundary face). Id computed automatically.
     *  It assumes that all attributes of facet f have been properly set
     *  @param f Facet to be added.
     *  @return Reference to the newly added facet.
     */
    facet_Type& addFacet ( edge_Type const& facet )
    {
        return addEdge (facet);
    }

    //! Adds a facet. 1D specialization
    /**
     *  Adds a facet (optionally a boundary face). Id computed automatically.
     *  It assumes that all attributes of facet f have been properly set
     *  @param f Facet to be added.
     *  @return Reference to the newly added facet.
     */
    facet_Type& addFacet ( point_Type const& facet )
    {
        return addPoint (facet);
    }

    //! i-th mesh Facet.
    /**
     *  @param i index of the mesh facet.
     *  @return the i-th facet.
     */
    facet_Type const& facet ( UInt const i ) const
    {
        return facet ( M_geoDim, i);
    }


    //! i-th mesh facet.
    /**
     *  @param i index of the mesh facet.
     *  @return reference to the ith mesh facet.
     */
    facet_Type& facet ( UInt const i )
    {
        return facet ( M_geoDim, i);
    }

    //! Set counter of facets.
    /**
     *  @param n Number of facets.
     */
    void setNumFacets ( UInt const n )
    {
        setNumFacets ( M_geoDim, n );
    }

    //! Set counter of boundary facets.
    /**
     *  @param n Number of boundary facets.
     */
    void setNumBoundaryFacets ( UInt const n )
    {
        setNumBoundaryFacets ( M_geoDim, n ) ;
    }

    //! Is facet whose id is given on boundary?
    /**
     *  @param id Facet Id.
     *  @return true if f in on boundary.
     */
    bool isBoundaryFacet ( UInt const& id ) const
    {
        return isBoundaryFacet ( M_geoDim, id );
    }

    //! Number of Ridges.
    /**
     *  Returns number of Ridges in the mesh
     *  as given by the internal counter.
     *
     *  @return Number of Ridges.
     */
    UInt numRidges() const
    {
        return numRidges (M_geoDim);
    }

    //! Global number of Ridges.
    /**
     *  Returns global number of Ridges in the mesh
     *  as given by the internal counter.
     *
     *  @return Global number of Ridges.
     */
    UInt numGlobalRidges() const
    {
        return numGlobalRidges (M_geoDim );
    }

    //! Changes Current capacity of Ridges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of ridges.
     *  @param setcounter true to set the counter, false otherwise.
     */
    void setMaxNumRidges ( UInt const n, bool const setcounter = false )
    {
        setMaxNumRidges (M_geoDim, n, setcounter);
    }

    //! Changes Current capacity of Global Ridges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of global ridges.
     */
    void setMaxNumGlobalRidges ( UInt const n )
    {
        setMaxNumGlobalRidges (M_geoDim, n);
    }

    //! Set counter of facets.
    /**
     *  @param n Number of facets.
     */
    void setNumRidges ( UInt const n )
    {
        setNumRidges ( M_geoDim, n );
    }

    //! Adds a Ridge.
    /**
     *  Adds a ridge. Id computed automatically.
     *
     *  @param boundary true if is on boundary.
     *  @return Reference to added ridge.
     */
    ridge_Type& addRidge ( bool const boundary )
    {
        return addRidge ( M_geoDim, boundary);
    }

    //! Adds a Ridge. Specialization for 3D geometries
    /**
     *  Adds a ridge to the end of the list
     *  and adjourn its Id. The ridge attributes (a part the id) should be
     *  correctly set.
     *
     *  @param e ridge to add.
     *  @return Reference to added ridge.
     */
    ridge_Type& addRidge ( edge_Type const& r )
    {
        return addEdge ( r );
    }

    //! Adds a Ridge. Specialization for 2D geometries
    /**
     *  Adds a ridge to the end of the list
     *  and adjourn its Id. The ridge attributes (a part the id) should be
     *  correctly set.
     *
     *  @param e ridge to add.
     *  @return Reference to added ridge.
     */
    ridge_Type& addRidge ( point_Type const& r )
    {
        return addPoint ( r );
    }

    //! i-th mesh ridge.
    /**
     *  Returns the i-th ridge.
     *
     *  @param i Index of the mesh ridge.
     *  @return The i-th ridge.
     */
    ridge_Type const& ridge ( UInt const i ) const
    {
        return ridge (M_geoDim, i);
    }

    //! i-th mesh ridge reference.
    /**
     *  Returns a reference to the i-th mesh ridge.
     *
     *  @param i Index of the mesh ridge.
     *  @return Reference to the i-th ridge.
     */
    ridge_Type& ridge ( UInt const i )
    {
        return ridge ( M_geoDim, i);
    }

    //! Set boundary ridge counter.
    /**
     *  Set the Boundary ridges counter to a given number.
     *
     *  @param n Count of Boundary Ridge.
     */
    void setNumBoundaryRidges ( UInt const n )
    {
        setNumBoundaryRidges ( M_geoDim, n );
    }

    //! Ridge on boundary check. Specialization for 3D geometries.
    /**
     *  Is this ridge on boundary?
     *
     *  @param e The ridge.
     *  @return true if thee ridge is on the boundary, false otherwise.
     */
    bool isBoundaryRidge ( edge_Type const& r ) const
    {
        return isBoundaryEdge (r);
    }

    //! Ridge on boundary check. Specialization for 2D geometries.
    /**
     *  Is this ridge on boundary?
     *
     *  @param e The ridge.
     *  @return true if thee ridge is on the boundary, false otherwise.
     */
    bool isBoundaryRidge ( point_Type const& r ) const
    {
        return isBoundaryPoint (r);
    }

    //! Ridge on boundary check by id.
    /**
     *  Is this ridge, of given id, on boundary?
     *
     *  @param id Id of the ridge.
     *  @return true if the ridge is on the boundary, false otherwise.
     */
    bool isBoundaryRidge ( UInt const& id ) const
    {
        return isBoundaryRidge ( M_geoDim, id );
    }

    //! Returns the i-th mesh Peak.
    /**
     *  Returns the i-th peak in the mesh.
     *
     *  @param i Id of the peak.
     *  @return i-th mesh peak.
     */
    peak_Type const& peak ( UInt const i ) const
    {
        return peak ( M_geoDim, i );
    }

    //! Returns a reference to the i-th mesh peak.
    /**
     *  Returns the i-th peak in the mesh.
     *UInt const i
     *  @param i Id of the peak.
     *  @return Reference i-th mesh peak.
     */
    //  peak_Type & peak( UInt const i ) {return peak( M_geoDim, i );}

    //! Returns the global number of peaks in the mesh.
    /**
     *  Interrogation to the counter of peaks in the mesh.
     *
     *  @return Number of peaks in the mesh.
     */
    UInt numGlobalPeaks() const
    {
        return numGlobalPeaks (M_geoDim );
    }

    /** @} */ // End of group Polytope Methods

    /** @defgroup public_attributes Public Attributes
     *
     */

    /** @name Region Containers
     *  @ingroup public_attributes
     *
     *  STL compliant containers for basic structure.
     *
     *  I expose them since they are standard containers.
     *
     *  @{
     */

    //! Container of mesh Points/Vertices.
    points_Type pointList;
    //! Container of mesh Faces.
    faces_Type faceList;

    //! Container of mesh Edges.
    edges_Type edgeList;
    //! Container of mesh 3D Elements
    volumes_Type volumeList;
    //! Boundary points list.
    std::vector<point_Type* > _bPoints;
    /** @} */ // End of group Region Containers

    /** @name Container Polytope Getters
     *  @ingroup public_attributes
     *  getters for region Containers, by polytope names
     *  @{
     */

    //! returns a reference to the elements' container
    elements_Type& elementList()
    {
        return elementList (M_geoDim );
    }

    //! returns a reference to the facets' container
    facets_Type& facetList()
    {
        return facetList (M_geoDim);
    }

    //! returns a reference to the ridges' container
    ridges_Type& ridgeList()
    {
        return ridgeList (M_geoDim );
    }

    /** @} */ // End of group Container Polytope Getters


    /** @name Switches
     *  @ingroup public_attributes
     *
     *  @{
     */

    //! Switches
    Switch switches;

    /** @} */ // End of group Switches

    //! Is the array for local Edges set up? Specialization for 3D geometries
    bool hasLocalEdges ( ridge_Type ) const
    {
        return hasLocalRidges();
    }

    //! Is the array for local Edges set up? Specialization for 2D geometries
    bool hasLocalEdges ( facet_Type) const
    {
        return hasLocalFacets();
    }


private:

    /*! Arrays containing the ids of Edges and Faces of each element
      I use a Define to use localto global array or directly the
      bareedges */
    ArraySimple<UInt> M_ElemToFacet;
    ArraySimple<UInt> M_ElemToRidge;

    UInt M_numVolumes;

    UInt M_numVertices;
    UInt M_numBVertices;

    UInt M_numPoints;
    UInt M_numBPoints;

    UInt M_numFaces;
    UInt M_numBFaces;

    UInt M_numEdges;
    UInt M_numBEdges;

    UInt M_numGlobalPoints;
    UInt M_numGlobalVertices;
    UInt M_numGlobalEdges;
    UInt M_numGlobalFaces;
    UInt M_numGlobalVolumes;

    bool M_isPartitioned;
    typename  markerCommon_Type::regionMarker_Type M_marker;
    MeshUtility::MeshTransformer<RegionMesh<geoShape_Type, markerCommon_Type>, markerCommon_Type > M_meshTransformer;

    // communicator
    commPtr_Type M_comm;

    //used to select the correct method specialization
    geoDim_Type M_geoDim;

    //fake ridge, used in 1D geometries to avoid compilation errors when a ridge& is required in output
    ridge_Type M_aRidge;

    //fake peak, used in 1D geometries to avoid compilation errors when a peak& is required in output
    peak_Type M_aPeak;

    UInt numElements (threeD_Type) const
    {
        return numVolumes();
    }
    UInt numElements (twoD_Type) const
    {
        return numFaces();
    }
    UInt numElements (oneD_Type) const
    {
        return numEdges();
    }

    //! Number of Boundary facets
    UInt numBoundaryFacets (threeD_Type) const
    {
        return numBFaces();
    }
    UInt numBoundaryFacets (twoD_Type) const
    {
        return numBEdges();
    }
    UInt numBoundaryFacets (oneD_Type) const
    {
        return numBVertices();
    }

    //! Get element at the i-th index.
    element_Type& element ( threeD_Type, const UInt& i )
    {
        return volume (i);
    }
    element_Type& element ( twoD_Type, const UInt& i )
    {
        return face (i);
    }
    element_Type& element ( oneD_Type, const UInt& i )
    {
        return edge (i);
    }

    //! Get element at the i-th index.
    const element_Type& element ( threeD_Type, const UInt& i ) const
    {
        return volume (i);
    }
    const element_Type& element ( twoD_Type, const UInt& i ) const
    {
        return face (i);
    }
    const element_Type& element ( oneD_Type, const UInt& i ) const
    {
        return edge (i);
    }

    //! Get boundary facet at the i-th index.
    facet_Type& boundaryFacet ( threeD_Type, const UInt& i )
    {
        return boundaryFace (i);
    }
    facet_Type& boundaryFacet ( twoD_Type, const UInt& i )
    {
        return boundaryEdge (i);
    }
    inline facet_Type& boundaryFacet ( oneD_Type, const UInt& i );

    //! Get boundary facet at the i-th index.
    const facet_Type& boundaryFacet ( threeD_Type, const UInt& i ) const
    {
        return boundaryFace (i);
    }
    const facet_Type& boundaryFacet ( twoD_Type, const UInt& i ) const
    {
        return boundaryEdge (i);
    }
    const facet_Type& boundaryFacet ( oneD_Type, const UInt& i ) const
    {
        return boundaryPoint (i);
    }

    //! Number of global elements.
    UInt numGlobalElements (threeD_Type) const
    {
        return numGlobalVolumes();
    }
    UInt numGlobalElements (twoD_Type) const
    {
        return numGlobalFaces();
    }
    UInt numGlobalElements (oneD_Type) const
    {
        return numGlobalEdges();
    }

    //! Current capacity of the container of Elements.
    UInt maxNumElements (threeD_Type) const
    {
        return maxNumVolumes();
    }
    UInt maxNumElements (twoD_Type) const
    {
        return maxNumFaces();
    }
    UInt maxNumElements (oneD_Type) const
    {
        return maxNumEdges();
    }

    //! Set counter of elements.
    void setNumElements ( threeD_Type, UInt const n )
    {
        setNumVolumes ( n );
    }
    void setNumElements ( twoD_Type, UInt const n )
    {
        setNumFaces ( n );
    }
    void setNumElements ( oneD_Type, UInt const n )
    {
        setNumEdges ( n );
    }

    //! Changes Current capacity of the container of elements.
    void setMaxNumElements   ( threeD_Type, UInt const n, bool const setcounter = false )
    {
        setMaxNumVolumes ( n, setcounter);
    }
    void setMaxNumElements   ( twoD_Type, UInt const n, bool const setcounter = false )
    {
        setMaxNumFaces ( n, setcounter);
    }
    void setMaxNumElements   ( oneD_Type, UInt const n, bool const setcounter = false )
    {
        setMaxNumEdges ( n, setcounter);
    }

    //! Set the number of global elements.
    void setMaxNumGlobalElements ( threeD_Type, UInt const n )
    {
        setMaxNumGlobalVolumes (n) ;
    }
    void setMaxNumGlobalElements ( twoD_Type, UInt const n )
    {
        setMaxNumGlobalFaces (n) ;
    }
    void setMaxNumGlobalElements ( oneD_Type, UInt const n )
    {
        setMaxNumGlobalPoints (n) ;
    }

    //! Adds element.
    element_Type& addElement ( volume_Type const& elem )
    {
        return addVolume (elem);
    }
    element_Type& addElement ( face_Type const& elem )
    {
        return addFace (elem);
    }
    element_Type& addElement ( edge_Type const& elem )
    {
        return addEdge (elem);
    }

    //! Is the array for local faces set up?
    bool hasLocalFaces (facet_Type) const
    {
        return hasLocalFacets();
    }
    bool hasLocalFaces (ridge_Type) const
    {
        return hasLocalRidges();
    }

    //! Build localEdgeId table and optionally fills the list of Edges
    void updateElementEdges ( ridge_Type, bool createEdges = false, const bool verbose = false,
                              UInt estimateEdgeNumber = 0, bool renumber = true)
    {
        updateElementRidges ( M_geoDim, createEdges, verbose, estimateEdgeNumber, renumber);
    }

    //! Build localEdgeId table and optionally fills the list of Edges
    void updateElementEdges ( facet_Type, bool createEdges = false, const bool verbose = false,
                              UInt estimateEdgeNumber = 0, bool /*renumber*/ = true)
    {
        updateElementFacets ( createEdges, verbose, estimateEdgeNumber);
    }

    //! Build localRidgeId table and optionally fills the list of Ridges
    void updateElementRidges ( threeD_Type, bool createRidges = false, const bool verbose = false,
                               UInt estimateRidgeNumber = 0, bool renumber = true);
    void updateElementRidges ( twoD_Type, bool, const bool, UInt, bool)
    {
        ERROR_MSG ("RegionMesh::updateElementRidges, It is not possible to use this method with 2D geometries.");
    }
    void updateElementRidges ( oneD_Type, bool, const bool, UInt, bool)
    {
        ERROR_MSG ("RegionMesh::updateElementRidges, It is not possible to use this method with 1D geometries.");
    }

    //! specializations for cleanElementEdges
    void cleanElementEdges (ridge_Type)
    {
        cleanElementRidges();
    }
    void cleanElementEdges (facet_Type)
    {
        cleanElementFacets();
    }

    //! Local Ridge.
    ID localRidgeId ( threeD_Type,  UInt const elemId, UInt const locR ) const;
    ID localRidgeId ( twoD_Type, UInt const elemId, UInt const locR ) const
    {
        return element (elemId).point (locR).localId();
    }
    ID localRidgeId ( oneD_Type, UInt const, UInt const) const
    {
        return NotAnId;
    }


    //! Local Edge (specialization for 3D geometries).
    ID localEdgeId ( const threeD_Type, UInt const elemId, UInt const locE ) const
    {
        return localRidgeId ( threeD_Type(), elemId, locE );
    }

    //! Local Edge (specialization for 2D geometries).
    ID localEdgeId ( const twoD_Type, UInt const elemId, UInt const locE ) const
    {
        return localFacetId ( elemId, locE );
    }

    //! Local Edge (specialization for 1D geometries). It calls an error.
    ID localEdgeId ( const oneD_Type, UInt const /*elemId*/, UInt const /*locE*/ ) const
    {
        ERROR_MSG ( "localEdgeId not implemented for oneD_Type." );
        return static_cast<ID> (0);
    }

    //! specializations for numFacets
    UInt numFacets (threeD_Type) const
    {
        return numFaces();
    }
    UInt numFacets (twoD_Type) const
    {
        return numEdges();
    }
    UInt numFacets (oneD_Type) const
    {
        return numVertices();
    }

    //! specializations numGlobalFacets
    UInt numGlobalFacets (threeD_Type) const
    {
        return numGlobalFaces();
    }
    UInt numGlobalFacets (twoD_Type) const
    {
        return numGlobalEdges();
    }
    UInt numGlobalFacets (oneD_Type) const
    {
        return numGlobalVertices();
    }

    //! Changes Current capacity of Facets.
    void setMaxNumFacets ( threeD_Type, UInt const n, bool const setcounter = false )
    {
        setMaxNumFaces ( n, setcounter);
    }
    void setMaxNumFacets ( twoD_Type, UInt const n, bool const setcounter = false )
    {
        setMaxNumEdges ( n, setcounter);
    }
    void setMaxNumFacets ( oneD_Type, UInt const n, bool const setcounter = false )
    {
        setMaxNumPoints ( n, setcounter);
    }

    //! Changes Current capacity of Global Facets.
    void setMaxNumGlobalFacets ( threeD_Type, UInt const n )
    {
        setMaxNumGlobalFaces ( n );
    }
    void setMaxNumGlobalFacets ( twoD_Type, UInt const n )
    {
        setMaxNumGlobalEdges ( n );
    }
    void setMaxNumGlobalFacets ( oneD_Type, UInt const n )
    {
        setMaxNumGlobalPoints ( n );
    }

    //! Adds a facet.
    facet_Type& addFacet ( threeD_Type, bool const boundary )
    {
        return addFace (boundary);
    }
    facet_Type& addFacet ( twoD_Type, bool const boundary )
    {
        return addEdge (boundary);
    }
    facet_Type& addFacet ( oneD_Type, bool const boundary )
    {
        return addPoint (boundary);
    }

    //! i-th mesh Facet.
    facet_Type const& facet ( threeD_Type, UInt const i ) const
    {
        return face (i);
    }
    facet_Type const& facet ( twoD_Type, UInt const i ) const
    {
        return edge (i);
    }
    facet_Type const& facet ( oneD_Type, UInt const i ) const
    {
        return point (i);
    }

    //! i-th mesh facet.
    facet_Type& facet ( threeD_Type, UInt const i )
    {
        return face (i);
    }
    facet_Type& facet ( twoD_Type, UInt const i )
    {
        return edge (i);
    }
    facet_Type& facet ( oneD_Type, UInt const i )
    {
        return point (i);
    }

    //! Set counter of facets.
    void setNumFacets ( threeD_Type, UInt const n )
    {
        setNumFaces ( n );
    }
    void setNumFacets ( twoD_Type, UInt const n )
    {
        setNumEdges ( n );
    }
    void setNumFacets ( oneD_Type, UInt const n )
    {
        setNumVertices ( n );
    }

    //! Set counter of boundary facets.
    void setNumBoundaryFacets ( threeD_Type, UInt const n )
    {
        setNumBFaces ( n ) ;
    }
    void setNumBoundaryFacets ( twoD_Type, UInt const n )
    {
        setNumBEdges ( n ) ;
    }
    void setNumBoundaryFacets ( oneD_Type, UInt const n )
    {
        setNumBVertices ( n ) ;
    }

    //! Is facet whose id is given on boundary?
    bool isBoundaryFacet ( threeD_Type, UInt const& id ) const
    {
        return isBoundaryFace ( id );
    }
    bool isBoundaryFacet ( twoD_Type, UInt const& id ) const
    {
        return isBoundaryEdge ( id );
    }
    bool isBoundaryFacet ( oneD_Type, UInt const& id ) const
    {
        return isBoundaryPoint ( id );
    }

    //! Number of Ridges.
    UInt numRidges (threeD_Type) const
    {
        return numEdges();
    }
    UInt numRidges (twoD_Type) const
    {
        return numPoints();
    }
    UInt numRidges (oneD_Type) const
    {
        return 0;
    }

    //! Global number of Ridges.
    UInt numGlobalRidges (threeD_Type) const
    {
        return numGlobalEdges();
    }
    UInt numGlobalRidges (twoD_Type) const
    {
        return numGlobalVertices();
    }
    UInt numGlobalRidges (oneD_Type) const
    {
        return 0;
    }

    //! Changes Current capacity of Ridges.
    void setMaxNumRidges ( threeD_Type, UInt const n, bool const setcounter )
    {
        setMaxNumEdges ( n, setcounter);
    }
    void setMaxNumRidges ( twoD_Type, UInt const n, bool const setcounter )
    {
        setMaxNumPoints ( n, setcounter);
    }
    void setMaxNumRidges ( oneD_Type, UInt const, bool const )
    {
        ERROR_MSG ("RegionMesh::setMaxNumRidges, No ridges in 1D");
    }

    //! Changes Current capacity of Global Ridges.
    void setMaxNumGlobalRidges ( threeD_Type, UInt const n )
    {
        setMaxNumGlobalEdges ( n );
    }
    void setMaxNumGlobalRidges ( twoD_Type, UInt const n )
    {
        setMaxNumGlobalPoints ( n );
    }
    void setMaxNumGlobalRidges ( oneD_Type, UInt const )
    {
        ERROR_MSG ("RegionMesh::setMaxNumGlobalRidges, No ridges in 1D");
    }

    //! Adds a Ridge.
    ridge_Type& addRidge ( threeD_Type, bool const boundary )
    {
        return addEdge ( boundary);
    }
    ridge_Type& addRidge ( twoD_Type, bool const boundary )
    {
        return addPoint ( boundary, false );
    }
    ridge_Type& addRidge ( oneD_Type, bool const )
    {
        ERROR_MSG ("RegionMesh::addRidge, No ridges in 1D"); ;
        return M_aRidge;
    }

    //! i-th mesh ridge.
    ridge_Type const& ridge ( threeD_Type, UInt const i ) const
    {
        return edge (i);
    }
    ridge_Type const& ridge ( twoD_Type, UInt const i ) const
    {
        return point (i);
    }
    ridge_Type const& ridge ( oneD_Type, UInt const ) const
    {
        ERROR_MSG ("RegionMesh::ridge, No ridges in 1D");
        return M_aRidge;
    }

    //! i-th mesh ridge reference.
    ridge_Type& ridge ( threeD_Type, UInt const i )
    {
        return edge (i);
    }
    ridge_Type& ridge ( twoD_Type, UInt const i )
    {
        return point (i);
    }
    ridge_Type& ridge ( oneD_Type, UInt const )
    {
        ERROR_MSG ("RegionMesh::ridge, No ridges in 1D");
        return M_aRidge;
    }

    //! Set counter of ridges.
    void setNumRidges ( threeD_Type, UInt const n )
    {
        setNumEdges ( n );
    }

    //! Set boundary ridge counter.
    void setNumBoundaryRidges ( threeD_Type, UInt const n )
    {
        setNumBEdges ( n );
    }
    void setNumBoundaryRidges ( twoD_Type, UInt const n )
    {
        setNumBPoints ( n );
    }
    void setNumBoundaryRidges ( oneD_Type, UInt const )
    {
        ERROR_MSG ("RegionMesh::setNumBoundaryRidges, No ridges in 1D");
    }

    //! Ridge on boundary check by id.
    bool isBoundaryRidge ( threeD_Type, UInt const& id ) const
    {
        return isBoundaryEdge ( id );
    }
    bool isBoundaryRidge ( twoD_Type, UInt const& id ) const
    {
        return isBoundaryPoint ( id );
    }
    bool isBoundaryRidge ( oneD_Type, UInt const& ) const
    {
        ERROR_MSG ("RegionMesh::isBoundaryRidge, No ridges in 1D");
        return bool();
    }

    //! Returns the i-th mesh Peak.
    peak_Type const& peak ( threeD_Type, UInt const i ) const
    {
        return point (i);
    }
    peak_Type const& peak ( twoD_Type, UInt const ) const
    {
        ERROR_MSG ("RegionMesh::peak, No peak in 2D");
        return M_aPeak;
    }
    peak_Type const& peak ( oneD_Type, UInt const ) const
    {
        ERROR_MSG ("RegionMesh::peak, No peak in 1D");
        return M_aPeak;
    }

    //! Returns a reference to the i-th mesh peak.
    peak_Type& peak ( threeD_Type, UInt const i )
    {
        return point (i);
    }
    peak_Type& peak ( twoD_Type, UInt const )
    {
        ERROR_MSG ("RegionMesh::peak, No peak in 2D");
        return M_aPeak;
    }
    peak_Type& peak ( oneD_Type, UInt const )
    {
        ERROR_MSG ("RegionMesh::peak, No peak in 2D");
        return M_aPeak;
    }

    //! Returns the global number of peaks in the mesh.
    UInt numGlobalPeaks (threeD_Type) const
    {
        return numGlobalVertices();
    }
    UInt numGlobalPeaks (twoD_Type) const
    {
        return 0;
    }
    UInt numGlobalPeaks (oneD_Type) const
    {
        return 0;
    }

    //! returns a reference to the elements' container
    elements_Type& elementList (threeD_Type)
    {
        return volumeList;
    }
    elements_Type& elementList (twoD_Type)
    {
        return faceList;
    }
    elements_Type& elementList (oneD_Type)
    {
        return edgeList;
    }

    //! returns a reference to the facets' container
    facets_Type& facetList (threeD_Type)
    {
        return faceList;
    }
    facets_Type& facetList (twoD_Type)
    {
        return edgeList;
    }
    facets_Type& facetList (oneD_Type)
    {
        return pointList;
    }

    //! returns a reference to the ridges' container
    ridges_Type& ridgeList (threeD_Type)
    {
        return edgeList;
    }
    ridges_Type& ridgeList (twoD_Type)
    {
        return pointList;
    }
    ridges_Type& ridgeList (oneD_Type)
    {
        ERROR_MSG ("RegionMesh::ridgeList, no RidgeList in 1D");
        return ridges_Type();
    }
}; // End of class RegionMesh


// =================================================== //
// =================================================== //
//                    IMPLEMENTATION                   //
// =================================================== //
// =================================================== //


void set_switches_for_regionmesh ( Switch& sw );

template <typename GeoShapeType, typename MCType>
inline RegionMesh<GeoShapeType, MCType>::RegionMesh() :
    MeshEntity(),
    switches(),
    M_numVolumes ( 0 ),
    M_numVertices ( 0 ),
    M_numBVertices ( 0 ),
    M_numPoints ( 0 ),
    M_numBPoints ( 0 ),
    M_numFaces ( 0 ),
    M_numBFaces ( 0 ),
    M_numEdges ( 0 ),
    M_numBEdges ( 0 ),
    M_isPartitioned ( false ),
    M_meshTransformer ( *this ),
    M_comm()
{
    set_switches_for_regionmesh ( switches );
}

template <typename GeoShapeType, typename MCType>
inline RegionMesh<GeoShapeType, MCType>::RegionMesh ( commPtr_Type const& comm ) :
    MeshEntity(),
    switches(),
    M_numVolumes ( 0 ),
    M_numVertices ( 0 ),
    M_numBVertices ( 0 ),
    M_numPoints ( 0 ),
    M_numBPoints ( 0 ),
    M_numFaces ( 0 ),
    M_numBFaces ( 0 ),
    M_numEdges ( 0 ),
    M_numBEdges ( 0 ),
    M_isPartitioned ( false ),
    M_meshTransformer ( *this ),
    M_comm ( comm )
{
    //Modif Miguel:11/2002
    set_switches_for_regionmesh ( switches );
}


template <typename GeoShapeType, typename MCType>
inline RegionMesh<GeoShapeType, MCType>::RegionMesh ( UInt id, commPtr_Type const& comm ) :
    MeshEntity ( id ),
    switches(),
    M_numVolumes ( 0 ),
    M_numVertices ( 0 ),
    M_numBVertices ( 0 ),
    M_numPoints ( 0 ),
    M_numBPoints ( 0 ),
    M_numFaces ( 0 ),
    M_numBFaces ( 0 ),
    M_numEdges ( 0 ),
    M_numBEdges ( 0 ),
    M_isPartitioned ( false ),
    M_meshTransformer ( *this ),
    M_comm ( comm )
{
    //Modif Miguel:11/2002
    set_switches_for_regionmesh ( switches );
}

template <typename GeoShapeType, typename MCType>
inline RegionMesh<GeoShapeType, MCType>::~RegionMesh()
{
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numLocalVertices() const
{
    return element_Type::S_numLocalVertices;
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numLocalFaces() const
{
    return element_Type::S_numLocalFaces;
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numLocalEdges() const
{
    return element_Type::S_numLocalEdges;
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numLocalEdgesOfFace() const
{
    return face_Type::S_numLocalEdges;
}

// ***************************** VOLUMES

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setMaxNumVolumes ( UInt const n, bool const setcounter )
{
    volumeList.setMaxNumItems (n);
    if ( setcounter )
    {
        M_numVolumes = n;
    }
}

template <typename GeoShapeType, typename MCType>
typename RegionMesh<GeoShapeType, MCType>::element_Type&
RegionMesh<GeoShapeType, MCType>::addVolume()
{
    // I need to set the global ID
    element_Type aVolume;
    return addVolume ( aVolume );
}

template <typename GeoShapeType, typename MCType>
inline typename RegionMesh<GeoShapeType, MCType>::element_Type&
RegionMesh<GeoShapeType, MCType>::addVolume ( element_Type const& v )
{
    volumeList.push_back ( v );
    volume_Type& thisVolume (volumeList.back() );
    thisVolume.setLocalId ( volumeList.size() - 1 );
    return thisVolume;
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::element_Type&
RegionMesh<GeoShapeType, MCType>::setVolume ( element_Type const& v, UInt const pos )
{
    ASSERT_PRE ( pos < volumeList.capacity() , "position requested exceed capacity" <<
                 pos << " " << volumeList.capacity() ) ;
    volumeList ( pos ) = v;
    volumeList ( pos ).setLocalId ( pos );
    return volumeList ( pos );
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::element_Type const&
RegionMesh<GeoShapeType, MCType>::volume ( UInt const i ) const
{
    ASSERT_BD ( i < volumeList.size() ) ;
    return volumeList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::element_Type&
RegionMesh<GeoShapeType, MCType>::volume ( UInt const i )
{
    ASSERT_BD ( i < volumeList.size() ) ;
    return volumeList[ i ];
}

// ************************* FACES ******************************
template <typename GeoShapeType, typename MCType>
inline UInt /*const*/
RegionMesh<GeoShapeType, MCType>::numFaces() const
{
    return M_numFaces;
}

template <typename GeoShapeType, typename MCType>
inline UInt /*const*/
RegionMesh<GeoShapeType, MCType>::numGlobalFaces() const
{
    return M_numGlobalFaces;
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::storedFaces() const
{
    return faceList.size();
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::maxNumFaces() const
{
    return faceList.maxNumItems();
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setMaxNumFaces ( UInt const n, bool const setcounter )
{
    faceList.setMaxNumItems (n);
    if ( setcounter )
    {
        M_numFaces = n;
    }
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setMaxNumGlobalFaces ( UInt const n )
{
    M_numGlobalFaces = n;
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type&
RegionMesh<GeoShapeType, MCType>::addFace ( bool const boundary )
{
    face_Type aFace;
    aFace.setBoundary ( boundary );
    return this->addFace ( aFace );
}


template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type&
RegionMesh<GeoShapeType, MCType>::addFace ( face_Type const& f )
{
    faceList.push_back ( f );
    face_Type& thisFace = faceList.back();
    thisFace.setLocalId ( faceList.size() - 1 );
    return thisFace;
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type&
RegionMesh<GeoShapeType, MCType>::setFace ( face_Type const& f, UInt position)
{
    ASSERT_PRE ( position < faceList.capacity(), "Face list size exceeded" <<
                 position << " " << faceList.capacity() ) ;
    faceList ( position ) = f;
    faceList ( position ).setLocalId ( position );
    return faceList ( position );
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type&
RegionMesh<GeoShapeType, MCType>::lastFace()
{
    return faceList.back();
}


template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type const&
RegionMesh<GeoShapeType, MCType>::face ( UInt const i ) const
{
    ASSERT_BD ( i < faceList.size() ) ;
    return faceList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type&
RegionMesh<GeoShapeType, MCType>::face ( UInt const i )
{
    ASSERT_BD ( i < faceList.size() ) ;
    return faceList[ i ];
}


template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type const&
RegionMesh<GeoShapeType, MCType>::boundaryFace ( UInt const i ) const
{
    ASSERT_PRE ( faceList.size() != 0, "Boundary Faces not stored" ) ;
    //ASSERT_BD( i < this->numBFaces()) ; //TODO TO BE FIXED! LF
    ASSERT_BD ( i < faceList.size() ) ;
    return faceList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::face_Type&
RegionMesh<GeoShapeType, MCType>::boundaryFace ( UInt const i )
{
    ASSERT_PRE ( faceList.size() != 0, "Boundary Faces not stored" ) ;
    //    ASSERT_BD( i < this->numBFaces()) ;
    ASSERT_BD ( i < faceList.size() ) ;
    return faceList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numBFaces() const
{
    return M_numBFaces;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumFaces ( UInt const n )
{
    M_numFaces = n;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumBFaces ( UInt const n )
{
    M_numBFaces = n;
}

// ************************* EDGES ******************************

template <typename GeoShapeType, typename MCType>
inline UInt /*const*/
RegionMesh<GeoShapeType, MCType>::numEdges() const
{
    return M_numEdges;
}

template <typename GeoShapeType, typename MCType>
inline UInt /*const*/
RegionMesh<GeoShapeType, MCType>::numGlobalEdges() const
{
    return M_numGlobalEdges;
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::storedEdges() const
{
    return edgeList.size();
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::maxNumEdges() const
{
    return edgeList.maxNumItems();
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setMaxNumEdges ( UInt const n, bool const setcounter )
{
    edgeList.setMaxNumItems (n);
    if ( setcounter )
    {
        M_numEdges = n;
    }
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setMaxNumGlobalEdges ( UInt const n )
{
    M_numGlobalEdges = n;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumEdges ( UInt const n)
{
    M_numEdges = n;
}

template <typename GeoShapeType, typename MCType>
inline typename RegionMesh<GeoShapeType, MCType>::edge_Type&
RegionMesh<GeoShapeType, MCType>::addEdge ( bool const boundary )
{
    edge_Type anEdge;
    anEdge.setBoundary ( boundary );
    return addEdge ( anEdge );
}

template <typename GeoShapeType, typename MCType>
inline typename RegionMesh<GeoShapeType, MCType>::edge_Type&
RegionMesh<GeoShapeType, MCType>::addEdge ( edge_Type const& f)
{
    edgeList.push_back ( f );
    edge_Type& thisEdge = edgeList.back();
    thisEdge.setLocalId ( edgeList.size() - 1 );
    return thisEdge;
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::edge_Type&
RegionMesh<GeoShapeType, MCType>::setEdge ( edge_Type const& f, UInt position)
{
    ASSERT_PRE ( position < edgeList.capacity(), "Edge list size exceeded" <<
                 position << " " << edgeList.capacity() ) ;
    edgeList ( position ) = f;
    edge_Type& thisEdge (edgeList ( position ) );
    thisEdge.setLocalId ( position );
    return thisEdge;
}


template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::edge_Type&
RegionMesh<GeoShapeType, MCType>::lastEdge()
{
    return edgeList.back();
}



template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::edge_Type const&
RegionMesh<GeoShapeType, MCType>::edge ( UInt const i ) const
{
    ASSERT_BD ( i < edgeList.size() ) ;
    return edgeList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::edge_Type&
RegionMesh<GeoShapeType, MCType>::edge ( UInt const i )
{
    ASSERT_BD ( i < edgeList.size() ) ;
    return edgeList[ i ];
}


template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::edge_Type const&
RegionMesh<GeoShapeType, MCType>::boundaryEdge ( UInt const i ) const
{
    ASSERT_PRE ( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD ( i < edgeList.size() ) ;
    return edgeList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::edge_Type&
RegionMesh<GeoShapeType, MCType>::boundaryEdge ( UInt const i )
{
    ASSERT_PRE ( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD ( i < edgeList.size() ) ;
    return edgeList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numBEdges() const
{
    return M_numBEdges;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumBEdges ( UInt const n )
{
    M_numBEdges = n;
}

// ************************ Points/Vertices

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numPoints() const
{
    return M_numPoints;
}

template <typename GeoShapeType, typename MCType>
inline UInt&
RegionMesh<GeoShapeType, MCType>::numPoints()
{
    return M_numPoints;
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::storedPoints() const
{
    return pointList.size();
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::storedBPoints() const
{
    return _bPoints.size();
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::maxNumPoints() const
{
    return pointList.maxNumItems();
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setMaxNumPoints ( UInt const n, bool const setcounter )
{
    pointList.setMaxNumItems (n);
    if ( setcounter )
    {
        M_numPoints = n;
    }
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumVertices ( UInt const n )
{
    M_numVertices = n;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumGlobalVertices ( UInt const n )
{
    M_numGlobalVertices = n;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumBVertices ( UInt const n )
{
    M_numBVertices = n;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setMaxNumGlobalPoints ( UInt const n )
{
    M_numGlobalPoints = n;
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::point_Type&
RegionMesh<GeoShapeType, MCType>::addPoint ( bool const boundary, bool const vertex )
{
    point_Type aPoint;
    aPoint.setBoundary ( boundary );
    if ( vertex )
    {
        aPoint.setFlag ( EntityFlags::VERTEX );
    }
    return addPoint ( aPoint );
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::point_Type&
RegionMesh<GeoShapeType, MCType>::addPoint ( point_Type const& p)
{
    ASSERT_PRE ( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
                 pointList.size() + 1 << " " << pointList.capacity() ) ;

    pointList.push_back ( p );
    point_Type& thisPoint ( pointList.back() );
    thisPoint.setLocalId ( pointList.size() - 1 );
    //todo This is bug prone!
    if ( thisPoint.boundary() )
    {
        ASSERT_PRE ( _bPoints.size() < _bPoints.capacity(), "Boundary Point list size exceeded" <<
                     _bPoints.size() + 1 << " " << _bPoints.capacity() ) ;
        _bPoints.push_back ( & pointList.back() );
    }
    return thisPoint;
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::point_Type const&
RegionMesh<GeoShapeType, MCType>::point ( UInt const i ) const
{
    ASSERT_BD ( i < pointList.size() ) ;

    return pointList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::point_Type&
RegionMesh<GeoShapeType, MCType>::point ( UInt const i )
{
    ASSERT_BD ( i < pointList.size() ) ;

    return pointList[ i ];
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::point_Type const&
RegionMesh<GeoShapeType, MCType>::boundaryPoint ( UInt const i ) const
{
    ASSERT_PRE ( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD ( i < _bPoints.size() ) ;
    return * ( _bPoints[ i ] );
}

template <typename GeoShapeType, typename MCType>
inline
typename RegionMesh<GeoShapeType, MCType>::point_Type&
RegionMesh<GeoShapeType, MCType>::boundaryPoint ( UInt const i )
{
    ASSERT_PRE ( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD ( i < _bPoints.size() ) ;
    return * ( _bPoints[ i ] );
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numBPoints() const
{
    return M_numBPoints;
}

template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::setNumBPoints ( UInt const n )
{
    M_numBPoints = n;
    _bPoints.reserve ( n );
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numVertices() const
{
    return M_numVertices;
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numGlobalVertices() const
{
    return M_numGlobalVertices;
}


template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::numBVertices() const
{
    return M_numBVertices;
}

template <typename GeoShapeType, typename MCType>
inline UInt&
RegionMesh<GeoShapeType, MCType>::numBVertices()
{
    return M_numBVertices;
}

// *************** GENERAL *******************

template <typename GeoShapeType, typename MCType>
inline std::ostream&
RegionMesh<GeoShapeType, MCType>::showMe ( bool verbose, std::ostream& out ) const
{
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << "                      RegionMesh                " << std::endl;
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << " Global ID: " << this->id() << " Marker Flag: " << this->markerID() << std::endl;
    //  out <<"Faces local to  volumes stored: "<<hasLocalFaces()<<std::endl;
    //out <<"Edges local to  volumes stored: "<<hasLocalEdges()<<std::endl;
    //out <<"Edges local to  faces   stored:"<<switches.test("FACEtoEDGE")<<std::endl;
    //out <<"Volumes adjacent to Faces stored: "<<switches.test("FACEtoVOLUME")<<std::endl;
    //out <<"Faces adjacent to Edges stored: "<<switches.test("EDGEtoFACE")<<std::endl<<std::endl;
    //out <<"Internal Faces Stored: " << hasInternalFaces()<<std::endl;
    //out <<"Edges Stored: " << hasEdges()<<" Internal: "
    //    << hasInternalEdges()<<std::endl;
    out << "**************** COUNTERS ********************" << std::endl;
    out << "NumPoints=" << numPoints() << "  " << "numBPoints=" << numBPoints() << std::endl;
    out << "NumVertices=" << numVertices() << "  " << "numBVerices=" << numBVertices() << std::endl;
    out << "NumVolumes=" << numVolumes() << "  " << "numFaces=" << numFaces() << std::endl;
    out << "NumBFaces=" << numBFaces() << "  " << "numEdges=" << numEdges() << std::endl;
    out << "NumBEdges=" << numBEdges() << std::endl;
    out << "**************************************************" << std::endl;
    out << "************ACTUALLY STORED  ********************" << std::endl;
    out << "Points=" << pointList.size() << "  " << "Edges=  " << edgeList.size() << std::endl;
    out << "Faces= " << faceList.size() << "  " << "Volumes=" << volumeList.size() << std::endl;
    out << "**************************************************" << std::endl;
    out << "*******         STATUS OF SWITCHES    ************" << std::endl;
    switches.showMe ( verbose, out );
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    if ( verbose )
    {
        std::cout << "Verbose version not implemented yet";
    }
    return out;

}




template <typename GeoShapeType, typename MCType>
inline Int
RegionMesh<GeoShapeType, MCType>::check ( Int level, bool const fix, bool verb, std::ostream& out )
{
    verb = verb && ( M_comm->MyPID() == 0 );
    Int severity = 0;
    Switch testsw;
    if ( verb )
    {
        out << "**************************************************" << std::endl;
        out << "         Checkin  RegionMesh                " << std::endl;
        out << " ID: " << this->id() << std::endl;
        out << "**************************************************" << std::endl;
    }
    if ( level == 1 )
    {
        checkMesh3D ( *this, testsw, fix, verb, out, out, out );
        if ( verb )
        {
            out << "**********************************************" << std::endl <<
                "DETAILS OF EXTENDED  CHECK:" << std::endl;
            testsw.showMe ( true, out );
            out << "**********************************************" << std::endl << std::endl;
            if ( testsw.test ( "ABORT_CONDITION" ) )
            {
                return 1;
            }

        }
    }


    if ( pointList.size() != M_numPoints )
    {
        if ( verb ) out << " Point list size " << pointList.size() << " not equal to internal counter value "
                            << M_numPoints << std::endl;
        if ( fix )
        {
            M_numPoints = pointList.size();
            if ( verb )
            {
                out << "Fixed";
                out.flush();
            }
        }
    }

    if ( edgeList.size() == 0 )
    {
        if ( verb )
        {
            out << "WARNING: No Edges Stored" << std::endl;
            out << "MAY NOT WORK IF WE HAVE DOF ON EDGES AND ESSENTIAL BC!!!!" << std::endl;
        }
        severity = -1;
        unsetLinkSwitch ( "HAS_ALL_RIDGES" );
        unsetLinkSwitch ( "HAS_BOUNDARY_RIDGES" );
    }
    else if ( edgeList.size() == numBEdges() )
    {

        setLinkSwitch ( "HAS_BOUNDARY_RIDGES" );
        unsetLinkSwitch ( "HAS_ALL_RIDGES" );
        if ( verb )
        {
            out << "INFORMATION: Only Boundary Edges Stored" << std::endl;
        }
    }
    else if ( edgeList.size() == numEdges() )
    {
        setLinkSwitch ( "HAS_BOUNDARY_RIDGES" );
        setLinkSwitch ( "HAS_ALL_RIDGES" );
    }
    else
    {
        if ( verb )
        {
            out << "ERROR: Inconsistent numbering of edges" << std::endl;
            out << "Counters: Boundary=" << numBEdges() << " Total=" << numEdges() << std::endl;
            out << "Stored =" << edgeList.size() << std::endl;
            severity = 1;
        }

    }

    if ( faceList.size() == 0 )
    {
        if ( verb )
        {
            out << "ERROR: No Faces Stored: at least boundary faces are needed" << std::endl;
            out << "" << std::endl;
        }
        severity = 1;
        unsetLinkSwitch ( "HAS_ALL_FACETS" );
        unsetLinkSwitch ( "HAS_BOUNDARY_FACETS" );
    }
    else if ( faceList.size() == numBFaces() )
    {
        setLinkSwitch ( "HAS_BOUNDARY_FACETS" );
        unsetLinkSwitch ( "HAS_ALL_FACETS" );
        if ( verb )
        {
            out << "INFORMATION: Only Boundary Faces Stored" << std::endl;
        }
    }
    else if ( faceList.size() == numFaces() )
    {
        setLinkSwitch ( "HAS_BOUNDARY_FACETS" );
        setLinkSwitch ( "HAS_ALL_FACETS" );
    }
    else
    {
        if ( verb )
        {
            out << "ERROR: Inconsistent numbering of faces" << std::endl;
            out << "Counters: Boundary=" << numBFaces() << " Total=" << numFaces() << std::endl;
            out << "Stored =" << faceList.size() << std::endl;
            severity = 1;
        }
    }


    UInt count = 0;
    for ( typename points_Type::iterator i = pointList.begin(); i != pointList.end(); ++i )
        if ( i->boundary() )
        {
            ++count;
        }
    if ( count == 0 )
    {
        severity = 4;
    }
    if ( count != M_numBPoints )
    {
        if ( verb ) out << "Num Boundary points " << count << " not equal to internal counter value "
                            << M_numBPoints << std::endl;
        if ( ( count != 0 ) & fix )
        {
            M_numBPoints = count;
            if ( verb )
            {
                out << "Fixed Counter";
                out.flush();
            }
        }
    }

    if ( M_numVertices == 0 )
    {
        severity = 6;

        if ( verb )
        {
            out << " SEVERITY ERROR: internal Vertices Counter unset";
        }
    }

    if ( M_numPoints == 0 )
    {
        severity = 6;
        if ( verb )
        {
            out << " SEVERITY ERROR: internal Points Counter unset";
        }
    }
    if ( M_numPoints == 0 )
    {
        severity = 6;
        if ( verb )
        {
            out << " SEVERITY ERROR: internal Points Counter unset";
        }
    }
    if ( M_numBPoints == 0 )
    {
        severity = 6;
        if ( verb )
        {
            out << " SEVERITY ERROR: boundary Points Counter unset";
        }
    }
    if ( M_numBVertices == 0 )
    {
        severity = 6;
        if ( verb )
        {
            out << " SEVERITY ERROR: boundary Vertices Counter unset";
        }
    }

    if ( verb )
        out << "   Check Finished              " << std::endl <<
            "***********************************************" << std::endl;
    setLinkSwitch ( "HAS_BEEN_CHECKED" );

    return severity;

}


template <typename GeoShapeType, typename MCType>
inline
UInt
RegionMesh<GeoShapeType, MCType>::faceElement ( UInt const i, UInt const Pos ) const
{
    ASSERT_PRE ( i < faceList.size(), "Not enough faces stored" ) ;
    ASSERT_PRE ( Pos <= 1 , "Wrong position (0 or 1)" ) ;
    if ( Pos == 0 )
    {
        return ( faceList[ i ] ).firstAdjacentElementIdentity();
    }
    else
    {
        return ( faceList[ i ] ).secondAdjacentElementIdentity();
    }
}

template <typename GeoShapeType, typename MCType>
inline
UInt
RegionMesh<GeoShapeType, MCType>::faceElement ( facet_Type const& f, UInt const Pos ) const
{
    ASSERT_BD ( ! faceList.empty() ) ;
    ASSERT_PRE ( Pos <= 1 , "Wrong position (0 or 1)" );
    if ( Pos == 0 )
    {
        return f.firstAdjacentElementIdentity();
    }
    else
    {
        return f.secondAdjacentElementIdentity();
    }
}


template <typename GeoShapeType, typename MCType>
inline
void
RegionMesh<GeoShapeType, MCType>::setLinkSwitch ( std::string const& _s )
{
    bool check = switches.set ( _s );
    ASSERT0 ( check, std::stringstream ( "Switch named " + _s + " is not allowed" ).str().c_str() );
    LIFEV_UNUSED ( check );
}

template <typename GeoShapeType, typename MCType>
inline
void
RegionMesh<GeoShapeType, MCType>::unsetLinkSwitch ( std::string const& _s )
{
    bool check = switches.unset ( _s );
    ASSERT0 ( check, std::stringstream ( "Switch named " + _s + " is not allowed" ).str().c_str() );
    LIFEV_UNUSED ( check );
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::getLinkSwitch ( std::string const& _s ) const
{
    return switches.test ( _s );
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::hasInternalFaces() const
{
    return faceList.size() > M_numBFaces;
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::hasInternalEdges() const
{
    return edgeList.size() > M_numBEdges;
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::hasEdges() const
{
    return ! edgeList.empty();
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::hasFaces() const
{
    return ! faceList.empty();
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::isVertex ( point_Type const& p ) const
{
    return Flag::testOneSet (p.flag(), EntityFlags::VERTEX);
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::isVertex ( UInt const& id ) const
{
    return isVertex (pointList[id]);
}


template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::isBoundaryPoint ( UInt const& id ) const
{
    return point ( id ).boundary();
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::isBoundaryPoint ( point_Type const& p ) const
{
    return p.boundary();
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::isBoundaryFace ( UInt const& id ) const
{
    return this->faceList[id].boundary();
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::localFacetId ( UInt const elemId, UInt const locF ) const
{
    ASSERT_PRE ( !M_ElemToFacet.empty(), "Element to Facet array not  set" );
    ASSERT_BD ( elemId < numElements() );
    ASSERT_BD ( locF < element_Type::S_numLocalFacets );

    return M_ElemToFacet.operator() ( locF, elemId );
}

template <typename GeoShapeType, typename MCType>
inline UInt
RegionMesh<GeoShapeType, MCType>::localFaceId ( UInt const volId, UInt const locF ) const
{
    ASSERT_PRE (S_geoDimensions == 3, "RegionMesh::localFaceId, It is not possible to use this method with 2D and 1D geometries.");
    return localFacetId ( volId, locF );
}

template <typename GeoShapeType, typename MCType>
inline
UInt
RegionMesh<GeoShapeType, MCType>::localRidgeId (threeD_Type,  UInt const elemId, UInt const locR )
const
{
    ASSERT_PRE ( !M_ElemToRidge.empty(), "Volume to Edges array not  set" );
    ASSERT_BD ( elemId < numElements() );
    ASSERT_BD ( locR < element_Type::S_numLocalRidges );
    return M_ElemToRidge ( locR, elemId );
}

template <typename GeoShapeType, typename MCType>
inline
bool
RegionMesh<GeoShapeType, MCType>::isFullEdge ( UInt const& id ) const
{
    return edgeList.size() > id;
}

template <typename GeoShapeType, typename MCType>
bool
RegionMesh<GeoShapeType, MCType>::isFullFace ( UInt const& id ) const
{
    return faceList.size() > id;
}

template <typename GeoShapeType, typename MCType>
typename RegionMesh<GeoShapeType, MCType>::facet_Type&
RegionMesh<GeoShapeType, MCType>::boundaryFacet ( oneD_Type, const UInt& i )
{
    ASSERT_BD (i < 2);
    return point (i * (numPoints() - 1) );
}

/********************************************************************************
                     ELEMENT3D:GLOBAL FACES/EDGES
 *******************************************************************************/

// Forward Declarations
template <typename GeoShapeType, typename MCType>
void
RegionMesh<GeoShapeType, MCType>::updateElementRidges (threeD_Type, bool ce, bool verb, UInt ee, bool renumber )
{
    bool verbose = verb && ( M_comm->MyPID() == 0 );

    if (S_geoDimensions != 3)
    {
        ERROR_MSG ("RegionMesh::updateElementRidges, It is not possible to use this method with 2D and 1D geometries.");
    }

    // If the counter is set we trust it! Otherwise we use Euler formula
    // this is ok for domains with at most 1 hole!

    if (verbose)
    {
        std::cout << "     Updating element ridges ... " << std::flush;
    }

    renumber = renumber && ce && !  this->ridgeList().empty();
    if ( ce && ee == 0 )
    {
        ee = M_numEdges > M_numBEdges ? M_numEdges : ( GeoShapeType::S_numFaces / 2 - 1 ) * numVolumes() + M_numBFaces / 2 + numVertices();
    }


    if ( ce )
    {
        // We want to create the edges, we need to reserve space
        ridgeList().setMaxNumItems (ee);
    }
    MeshElementBareHandler<BareEdge> bareEdge;
    std::pair<UInt, bool> e;
    M_ElemToRidge.reshape ( numLocalEdges(), numVolumes() ); // DIMENSION ARRAY

    UInt elemLocalID, i1, i2;
    std::pair<BareEdge, bool> _edge;
    GeoShapeType ele;
    facetShape_Type bele;
    // First We check if we have already Edges stored
    if ( ! ridgeList().empty() )
    {
        // dump first the existing edges, to maintain the correct numbering
        // if everything is correct the numbering in the bareedge
        // structure will reflect the actual edge numbering
        std::pair<UInt, bool> _check;
        for ( UInt j = 0; j < ridgeList().size(); ++j )
        {
            i1 = ( ridge ( j ).point ( 0 ) ).localId();
            i2 = ( ridge ( j ).point ( 1 ) ).localId();

            _edge  = makeBareEdge ( i1, i2 );
            _check = bareEdge.addIfNotThere ( _edge.first );
        }
    }

    ridge_Type edg;

    for ( typename faces_Type::iterator ifa = faceList.begin();
            ifa != faceList.begin() + M_numBFaces; ++ifa )
    {
        for ( UInt j = 0; j < numLocalEdgesOfFace(); j++ )
        {
            i1 = bele.edgeToPoint ( j, 0 );
            i2 = bele.edgeToPoint ( j, 1 );
            // go to global
            i1 = ( ifa->point ( i1 ) ).localId();
            i2 = ( ifa->point ( i2 ) ).localId();

            _edge = makeBareEdge ( i1, i2 );

            e = bareEdge.addIfNotThere ( _edge.first );

            if ( ce && e.second )
            {
                //
                for ( UInt k = 0; k < 2 + facetShape_Type::S_numPointsPerEdge; k++ )
                {
                    UInt inode = bele.edgeToPoint (j, k);
                    edg.setPoint ( k, ifa->point ( inode ) );
                }
                MeshUtility::inheritPointsWeakerMarker ( edg );
                edg.setBoundary ( true );
                edg.setId ( ridgeList().size() );
                addRidge ( edg );
            }
        }

    }

    if ( ce )
    {
        M_numBEdges = ridgeList().size();
        setLinkSwitch ( "HAS_BOUNDARY_RIDGES" );
    }

    for ( typename elements_Type::iterator elemIt = elementList().begin();
            elemIt != elementList().end(); ++elemIt )
    {
        elemLocalID = elemIt->localId();

        for ( UInt j = 0; j < numLocalEdges(); j++ )
        {
            i1 = ele.edgeToPoint ( j, 0 );
            i2 = ele.edgeToPoint ( j, 1 );
            // go to global
            i1 = ( elemIt->point ( i1 ) ).localId();
            i2 = ( elemIt->point ( i2 ) ).localId();
            _edge = makeBareEdge ( i1, i2 );

            e = bareEdge.addIfNotThere ( _edge.first );
            M_ElemToRidge.operator() ( j, elemLocalID ) = e.first;
            if ( ce && e.second )
            {
                for ( UInt k = 0; k < 2 + geoShape_Type::S_numPointsPerEdge; k++ )
                {
                    UInt inode = ele.edgeToPoint (j, k);
                    edg.setPoint ( k, elemIt->point ( inode ) );
                }
                MeshUtility::inheritPointsWeakerMarker ( edg );
                edg.setBoundary ( true );
                edg.setId ( ridgeList().size() );
                addRidge ( edg );
            }
        }
    }

    if ( ce )
    {
        M_numEdges = ridgeList().size();
        this->M_numBEdges =
            ridgeList().countElementsWithFlag (EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);
        setLinkSwitch ( "HAS_ALL_RIDGES" );
        if (this->M_numGlobalEdges == 0)
        {
            this->M_numGlobalEdges = M_numEdges;
        }
    }

    if (renumber && !ridgeList().empty() )
    {
        ridgeList().reorderAccordingToFlag (EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);
        std::vector<ID>newToOld = ridgeList().resetId(); //reset the ids so that they are in accord with position in the container.
        //Unfortunately I need oldToNew!
        std::vector<ID> oldToNew ( newToOld.size() );
        for (UInt j = 0; j < newToOld.size(); ++j)
        {
            oldToNew[ newToOld[j] ] = j;
        }
        // Save some memory annihilating newToOld
        std::vector<ID>().swap (newToOld);
        // Fix element to ridge array to reflect new ridge numbering
        // M_ElemToRidge is in fact a vector!
        std::vector<UInt> tmp ( M_ElemToRidge.size() );
        std::vector<UInt>::iterator tmp_it = tmp.begin();
        for (std::vector<UInt>::iterator it = M_ElemToRidge.begin(); it < M_ElemToRidge.end(); ++it, ++tmp_it)
        {
            *tmp_it = oldToNew[*it];
        }
        std::copy (tmp.begin(), tmp.end(), M_ElemToRidge.begin() );
    }

    UInt n = bareEdge.maxId();

    if (!ce)
    {
        if ( M_numEdges == 0 || M_numEdges == M_numBEdges )
        {
            M_numEdges = n;
        }
    }

    if (verbose)
    {
        std::cout << n << " edges found";
    }
    ASSERT_POS ( n == M_numEdges , "#Edges found is not equal to that in RegionMesh" << n << " " << M_numEdges ) ;
    setLinkSwitch ( std::string ( "HAS_ELEMENT_TO_RIDGES" ) );

    if (verbose)
    {
        std::cout << " done." << std::endl;
    }
}


//
// Update element faces
//

template <typename GeoShapeType, typename MCType>
void
RegionMesh<GeoShapeType, MCType>::updateElementFacets ( bool cf, bool verbose, UInt ef )
{
    verbose = verbose && ( M_comm->MyPID() == 0 );

    typedef BareEntitySelector<typename facetShape_Type::BasRefSha> bareEntitySelector_Type;
    typedef typename bareEntitySelector_Type::bareEntity_Type bareFacet_type;

    if (verbose)
    {
        std::cout << "     Updating element facets ... " << std::flush;
    }

    ASSERT0 ( ! cf || numBoundaryFacets() > 0, std::stringstream ( std::string ("Boundary Facets Must have been set") +
                                                                   std::string ("in order to call updateElementFacets with createFacets=true") +
                                                                   std::string ("\nUse buildBoundaryFacets(..) from mesh_util.h") ).str().c_str() );
    // If the counter is set we trust it! Otherwise we use Euler formula

    if ( cf && ef == 0 )
    {
        ef = numFacets() > numBoundaryFacets() ? numFacets() : ( geoShape_Type::S_numFacets * numElements() + numBoundaryFacets() ) / 2;
    }

    ASSERT ( cf || numFacets() > 0 , "Mesh is not properly set!" );

    if ( cf )
    {
        facetList().setMaxNumItems ( ef );
    }



    facet_Type aFacet;

    MeshElementBareHandler<bareFacet_type> bareFacet;
    // Extra map for facets stored which are not boundary facets
    MeshElementBareHandler<bareFacet_type> extraBareFacet;
    std::pair<UInt, bool> e;
    M_ElemToFacet.reshape ( element_Type::S_numLocalFacets, numElements() ); // DIMENSION ARRAY

    UInt elemLocalID;
    std::pair<bareFacet_type, bool>_facet;

    GeoShapeType ele;
    // If we have all facets and the facets store all adjacency info
    // everything is easier
    if ( (facetList().size() == numFacets() ) && getLinkSwitch ( "FACETS_HAVE_ADIACENCY" ) && getLinkSwitch ( "HAS_ALL_FACETS" ) )
    {
        for ( typename facets_Type::iterator itf = facetList().begin(); itf != facetList().end(); ++itf )
        {
            if ( itf->firstAdjacentElementPosition() != NotAnId && itf->firstAdjacentElementIdentity() != NotAnId)
            {
                M_ElemToFacet ( itf->firstAdjacentElementPosition() , itf->firstAdjacentElementIdentity() ) = itf->localId();
            }
            if ( itf->secondAdjacentElementPosition() != NotAnId && itf->secondAdjacentElementIdentity() != NotAnId)
            {
                M_ElemToFacet ( itf->secondAdjacentElementPosition(), itf->secondAdjacentElementIdentity() ) = itf->localId();
            }
        }
        // we finish here
        setLinkSwitch ( "HAS_ELEMENT_TO_FACETS" );
        if (verbose)
        {
            std::cout << " done." << std::endl;
        }

        return ;
    }

    // If I have only boundary facets I need to process them first to keep the correct numbering

    // First We check if we have already Facets stored
    UInt _numOriginalStoredFacets = facetList().size();
    ID points[facetShape_Type::S_numVertices];
    if ( ! facetList().empty() )
    {
        // dump all facets in the container, to maintain the correct numbering
        // if everything is correct the numbering in the bareFacet structure
        // will reflect the actual facet numbering. However, if I want to create
        // the internal facets I need to make sure that I am processing only the
        // boundary ones in a special way.
        std::pair<UInt, bool> _check;
        for ( UInt j = 0; j < facetList().size(); ++j )
        {
            for (UInt k = 0; k < facetShape_Type::S_numVertices; k++)
            {
                points[k] = ( facet ( j ).point ( k ) ).localId();
            }
            _facet = bareEntitySelector_Type::makeBareEntity ( points );
            _check = bareFacet.addIfNotThere ( _facet.first );
            if ( ! ( this->facet ( j ).boundary() ) )
            {
                extraBareFacet.addIfNotThere ( _facet.first, j);
            }
        }
    }
    UInt numFoundBoundaryFacets = bareFacet.size();
    UInt facetCount = numFoundBoundaryFacets;
    for ( typename elements_Type::iterator elemIt = elementList().begin();
            elemIt != elementList().end(); ++elemIt )
    {
        elemLocalID = elemIt->localId();
        for ( UInt j = 0; j < element_Type::S_numLocalFacets; j++ )
        {
            for (UInt k = 0; k < facetShape_Type::S_numVertices; k++)
            {
                UInt id = ele.facetToPoint ( j, k );
                points[k] = elemIt->point ( id ).localId();
            }
            _facet = bareEntitySelector_Type::makeBareEntity ( points );

            e = bareFacet.addIfNotThere ( _facet.first );
            M_ElemToFacet ( j, elemLocalID ) = e.first;
            bool _isBound = e.first < numFoundBoundaryFacets;
            // Is the facet an extra facet (not on the boundary but originally included in the list)?
            bool _isExtra = (e.first >= numFoundBoundaryFacets && e.first < _numOriginalStoredFacets);
            if ( _isBound )
            {
                facet_Type& _thisFacet (facet (e.first) );
                _thisFacet.firstAdjacentElementIdentity()   = elemLocalID;
                _thisFacet.firstAdjacentElementPosition()   = j;
                _thisFacet.secondAdjacentElementIdentity()  = NotAnId;
                _thisFacet.secondAdjacentElementPosition()  = NotAnId;
            }
            else if (_isExtra)
            {
                // This is not a bfacets and I need to set up all info about adjacency properly
                facet_Type& _thisFacet (facet (e.first) );
                // I need to check if it is the first time I meet it. Then I delete it from the
                // map: if it as there it means that it is the first time I am treating this face
                if (extraBareFacet.deleteIfThere (_facet.first) )
                {
                    // I need to be sure about orientation, the easiest thing is to rewrite the facet points
                    for ( UInt k = 0; k < facet_Type::S_numPoints; ++k )
                    {
                        _thisFacet.setPoint ( k, elemIt->point ( ele.facetToPoint ( j, k ) ) );
                    }
                    _thisFacet.firstAdjacentElementIdentity()  = elemLocalID;
                    _thisFacet.firstAdjacentElementPosition()  = j;

                }
                else
                {
                    _thisFacet.secondAdjacentElementIdentity()  = elemLocalID;
                    _thisFacet.secondAdjacentElementPosition()  = j;
                }
            }
            else if ( cf ) // A facet not contained in the original list.
                // I process it only if requested!
            {
                if ( e.second )
                {
                    // a new facet It must be internal.
                    for ( UInt k = 0; k < facet_Type::S_numPoints; ++k )
                    {
                        aFacet.setPoint ( k, elemIt->point ( ele.facetToPoint ( j, k ) ) );
                    }

                    aFacet.firstAdjacentElementIdentity()  = elemLocalID;
                    aFacet.firstAdjacentElementPosition() = j;

                    // gets the marker from the RegionMesh
                    aFacet.setMarkerID ( NotAnId );
                    aFacet.setBoundary (false);
                    aFacet.setId ( facetCount++ );
                    addFacet ( aFacet); //The id should be correct
                }
                else
                {
                    facet ( e.first ).secondAdjacentElementIdentity() = elemLocalID;
                    facet ( e.first ).secondAdjacentElementPosition() = j;
                }
            }
        }
    }

    UInt n = bareFacet.maxId();
    // LF Fix _numfacets. This part has to be checked. One may want to use
    // this method on a partitioned mesh, in which case the Global facets are there
    setNumFacets (n); // We have found the right total number of facets in the mesh
    if (numGlobalFacets() == 0)
    {
        setMaxNumGlobalFacets (n);    // If not already set fix it.
    }

    if (verbose)
    {
        std::cout << n << " facets ";
    }
    //ASSERT_POS( n == M_numFacets , "#Facets found inconsistent with that stored in RegionMesh" ) ;
    setLinkSwitch ( "HAS_ELEMENT_TO_FACETS" );
    if ( cf )
    {
        setLinkSwitch ( "HAS_ALL_FACETS" );
    }
    //if ( cf ) Facets have adjacency in any case!
    setLinkSwitch ( "FACETS_HAVE_ADIACENCY" );
    if (verbose)
    {
        std::cout << " done." << std::endl;
    }
}

template <typename GeoShapeType, typename MCType>
void RegionMesh<GeoShapeType, MCType>::updateElementFaces ( bool createFaces, const bool verbose, UInt estimateFaceNumber)
{
    ASSERT_PRE (S_geoDimensions == 3, "RegionMesh::updateElementFaces, It is not possible to use this method with 2D and 1D geometries.");
    updateElementFacets ( createFaces , verbose, estimateFaceNumber );
}

template <typename GeoShapeType, typename MCType>
void
RegionMesh<GeoShapeType, MCType>::cleanElementFacets()
{
    M_ElemToFacet.clearArray();

    unsetLinkSwitch ( "HAS_ELEMENT_TO_FACETS" );
}

template <typename GeoShapeType, typename MCType>
void
RegionMesh<GeoShapeType, MCType>::cleanElementFaces()
{
    ASSERT_PRE (S_geoDimensions == 3, "RegionMesh::cleanElementFaces(), It is not possible to use this method with 2D and 1D geometries.");
    cleanElementFacets();
}

template <typename GeoShapeType, typename MCType>
void
RegionMesh<GeoShapeType, MCType>::cleanElementRidges()
{
    M_ElemToRidge.clearArray();
    unsetLinkSwitch ( "HAS_ELEMENT_TO_RIDGES" );
}


template <typename GeoShapeType, typename MCType>
inline void
RegionMesh<GeoShapeType, MCType>::
getListOfPoints ( bool ( *fct ) ( Real, Real, Real ), std::vector<UInt>& list_pts )
{
    for ( UInt i = 0; i < M_numPoints; i++ )
    {
        MeshVertex& pt = pointList ( i );
        if ( fct ( pt.x(), pt.y(), pt.z() ) )
        {
            list_pts.push_back ( i );
        }
    }
}


template <typename GeoShapeType, typename MCType>
inline MeshUtility::MeshTransformer<RegionMesh<GeoShapeType, MCType>, MCType >&
RegionMesh<GeoShapeType, MCType>::meshTransformer()
{
    return this->M_meshTransformer;
}

template <typename GeoShapeType, typename MCType>
inline typename RegionMesh<GeoShapeType, MCType>::commPtr_Type
RegionMesh<GeoShapeType, MCType>::comm() const
{
    return this->M_comm;
}

template <typename GeoShapeType, typename MCType>
void inline RegionMesh<GeoShapeType, MCType>::setComm ( commPtr_Type const& comm )
{
    M_comm = comm;
}

} // End of namespace LifeV

#endif //REGIONMESH_H

