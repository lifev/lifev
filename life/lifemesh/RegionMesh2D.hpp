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
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with LifeV. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************
 */
//@HEADER

/*!
 *  @file
 *  @brief File containing 2D Mesh Classes
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *
 *  @contributor Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 *  @mantainer Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 */

#ifndef _REGIONMESH2D_HH_
#define _REGIONMESH2D_HH_

#include <cstdlib>
#include <iomanip>
#include <fstream>

#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/LifeDebug.hpp>
#include <life/lifemesh/MeshElementMarked.hpp>
#include <life/lifecore/Switch.hpp>
#include <life/lifemesh/MeshElementBare.hpp>

#include <life/lifearray/MeshEntityContainer.hpp>
#include <life/lifearray/ArraySimple.hpp>
#include <life/lifemesh/ElementShapes.hpp>
#include <life/lifemesh/MeshUtility.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#ifdef HAVE_MPI
//headers useful only for reordering:
#include "mpi.h"
#include <parmetis.h>
#endif

namespace LifeV
{

/**
 *  @class RegionMesh2D
 *  @brief Class for 2D Mesh
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *
 *  This is the class that stores the mesh entities for a single 2D region.
 *
 *  In a region elements are all of the same type.
 *
 *  Note: to provide data useful in a parallel setting some methods
 *  return either the number of entities in the current mesh region or
 *  the ones in the global mesh, before partitioning. The latter are identified
 *  by the keyword Global in the name, e.g. numGlobalFaces() versus numFaces()
 */
template <typename GEOSHAPE, typename MC = defaultMarkerCommon_Type >
class RegionMesh2D
:
public MeshEntity,
public MC::regionMarker_Type
{
public:
    /** @name Marker Types
     *  @ingroup public_types
     *  Markers for Point, Edge, Face, Volume and Region.
     *
     *  @{
     */

	static const Int geoDimensions = GEOSHAPE::S_nDimensions;

    //! Common Markers
    typedef MC MarkerCommon;
    //! Point Marker
    typedef typename MC::pointMarker_Type pointMarker_Type;
    //! Edge Marker
    typedef typename MC::edgeMarker_Type edgeMarker_Type;
    //! Face Marker
    typedef typename MC::faceMarker_Type faceMarker_Type;
    //! Volume Marker
    typedef typename MC::volumeMarker_Type volumeMarker_Type;
    //! Region Marker
    typedef typename MC::regionMarker_Type regionMarker_Type;
    //! Region Marker (obsolete)
    typedef typename MC::regionMarker_Type  Marker;
    //! Region Marker (generic name)
    typedef typename MC::regionMarker_Type  marker_Type;

    /** @} */ // End of group Marker Types




    /** @name Geometric Element Types
     *  @ingroup public_types
     *  Volumes, Faces, Edges and Points.
     *  @{
     */

    //! Volume Element (3D)
    typedef MeshElementMarked<3, geoDimensions, GEOSHAPE, MC>  volume_Type;
    typedef MeshElementMarked<geoDimensions, geoDimensions, GEOSHAPE, MC>  element_Type;
    
    //! Face Element (2D)
    typedef MeshElementMarked<geoDimensions-1, geoDimensions, GEOSHAPE, MC> facet_Type;
    typedef MeshElementMarked<2, geoDimensions, GEOSHAPE, MC> face_Type;

    //! Edge Element (1D)
    typedef MeshElementMarked<geoDimensions-2, geoDimensions, GEOSHAPE, MC> ridge_Type;
    typedef MeshElementMarked<1, geoDimensions, GEOSHAPE, MC> edge_Type;
    //! Point Element (0D)
    typedef MeshElementMarked<0, geoDimensions, GEOSHAPE, MC>	point_Type;
    typedef MeshElementMarked<geoDimensions-3, geoDimensions, GEOSHAPE, MC>  peak_Type;

    /** @} */ // End of group Geometric Element Types


    /** @name Basic Element Shape Types
     *  @ingroup public_types
     *  Volume, Face and Edge geometric shapes.
     *  @{
     */

    //! Element Shape.
    typedef typename volume_Type::geoShape_Type volumeShape_Type;
    typedef GEOSHAPE elementShape_Type;

    //! Facet Shape (Boundary Facet).
    typedef typename face_Type::geoShape_Type faceShape_Type;
    typedef typename GEOSHAPE::GeoBShape facetShape_Type;

    //! Ridge Shape (Boundary of Boundary Facet)
    typedef typename edge_Type::geoShape_Type edgeShape_Type;
    typedef typename facetShape_Type::GeoBShape ridgeShape_Type;

    /** @} */ // End of group Basic Element Shape Types

    /** @name Geometric Element Container Types
     *  @ingroup public_types
     *  Typedefs for STL compliant containers of mesh geometric entities.
     *  @{
     */
    //! Points Container.
    typedef MeshEntityContainer<point_Type>   points_Type;
    typedef MeshEntityContainer<peak_Type>   peaks_Type;
    
    //! Elements Container.
    typedef MeshEntityContainer<volume_Type > volumes_Type;
    typedef MeshEntityContainer<element_Type>  elements_Type;
    
    //! Facets Container: it may contain only Boundary facets.
    typedef MeshEntityContainer<face_Type>    faces_Type;
    typedef MeshEntityContainer<facet_Type>    facets_Type;
    
    //! Ridges Container: it may be empty.
    typedef MeshEntityContainer<edge_Type>    edges_Type;
    typedef MeshEntityContainer<ridge_Type>    ridges_Type;


    /** @} */ // End of group Geometric Element Container Types



    /** @name Constructors & Destructor
     *  Default and Copy Constructor for the class.
     *  @{
     */

    //! Default constructor
    explicit RegionMesh2D();

    //! Default constructor
    /**
     *  @param id marker of the RegionMesh2D
     */
    explicit RegionMesh2D( UInt id );

    //! Destructor
    virtual ~RegionMesh2D<GEOSHAPE, MC>();

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
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;

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
    int check( int level = 0, bool const fix = false, bool const verbose = true, std::ostream & out = std::cerr );

    //! Display local to global mapping.
    /**
     *  @param os Output stream.
     */
    void printLtGMap(std::ostream & os);

    //! Return the handle to perform transformations on the mesh
    inline MeshUtility::MeshTransformer<RegionMesh2D<GEOSHAPE, MC> > & meshTransformer();

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
    bool getLinkSwitch( std::string const & _s ) const;

    //! Set a switch using a given name
    /**
     *  @param _s Name of the switch.
     */
    void setLinkSwitch( std::string const & _s );

    //! Unset a switch.
    /**
     *  @param _s Name of the switch.
     */
    void unsetLinkSwitch( std::string const & _s );
    //! Output switches contents
    /**
     *  @param verbose Verbosity of the output.
     *  @param out Output stream.
     *  @return Output stream for concatenation.
     */
    std::ostream & showLinkSwitch( bool verbose = false, std::ostream & out = std::cout )
    {
        return switches.showMe( verbose, out );
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

    //! Number of elements in mesh.
    /**
     *  @return Number of elements in mesh.
     *  @sa numFaces
     *  @note Alias to numFaces()
     */
    UInt numElements() const;

    //! Number of elements in mesh.
    /**
     *  @return Number of elements in mesh.
     *  @sa numFaces
     *  @note Alias to numElements()
     */
    UInt & numElements();

    //! Number of Boundary facets.
    /**
     *  @return Number of Boundary facets.
     */
    UInt numBFacets() const;

    //! Number of boundary faces.
    /**
     *  @return Number of boundary facets.
     */
    UInt& numBFacets();

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    element_Type& element( const UInt& i );

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    const element_Type& element( const UInt& i ) const;

    //! Get boundary facet at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary facet at index i
     */
    facet_Type& bFacet( const UInt& i );

    //! Get boundary facet at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary facet at index i
     */
    const facet_Type& bFacet( const UInt& i ) const;
    const facet_Type& bFace( const UInt& i ) const { return bFacet(i); }

    //! Number of global elements.
    /**
     * @return Number of global Elements (Volumes).
     */
    UInt numGlobalElements() const;

    //! Access to number of global elements.
    /**
     * @return Access to number of global Elements (Volumes).
     */
    UInt & numGlobalElements();

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
    UInt numVolumes()       const;

    //! Returns Global Number of Volumes
    /**
     *  Returns number of Global Volume elements in the mesh as given by the internal counter.
     *  @return Number of Global Volumes.
     */
    UInt numGlobalVolumes() const;

    //! Volumes actually stored in list.
    /**
     *  @return Number of stored Volumes.
     */
    UInt storedVolumes() const {return 0;};

    //! Current capacity of Volumes Container.
    /**
     *  @return how many elements may be stored.
     */
    UInt maxNumVolumes() const;
    inline UInt maxNumElements() const {return maxNumFaces();}

    //! Changes Current capacity of Volumes.
    /**
     *  Changes Current capacity of Volumes (Optionally sets internal counter).
     *
     *  @param n Maximum number of volumes.
     *  @param setcounter true to set the counter, false otherwise (default).
     */
    void setMaxNumVolumes( UInt const /*n*/, bool const /*setcounter*/ = false) {};
    inline void setMaxNumElements   ( UInt const n, bool const setcounter = false ) {setMaxNumFaces( n, setcounter);}

    //! Changes Current capacity of Global Volumes.
    /**
     *  Changes Current capacity of Global Volumes (Optionally sets internal counter).
     *
     *  @param n maximum number of global volumes.
     */
    void setMaxNumGlobalVolumes( UInt const n ) {};
    inline void setMaxNumGlobalElements( UInt const n ) {setMaxNumGlobalFaces(n) ;}

    //! Set Number of Volumes.
    /**
     *  Set number of Volume elements in the mesh by changing internal counter.
     *  @param n Number of volumes.
     */
    void setNumVolumes      ( UInt const n );
    inline void setNumElements      ( UInt const n ) {setNumFaces(n);}

    //! Adds volumes.
    /**
     *  Adds volume. Id computed automatically.
     *  @return Reference to added volume.
     */
    inline volume_Type & addVolume(){};
    inline element_Type & addElement() {return addFace();}

    //! Adds volumes.
    /**
     *  Adds volume. Id computed automatically.
     *  @param v Volume to be added.
     *  @return Reference to the newly added volume.
     */
    element_Type & addVolume( element_Type const & v );
    inline element_Type & addElement( element_Type const & v ) {return addFace(v);}

    //! Adds volume in a certain position.
    /**
     *  Adds volume to a specified position.
     *  @param v Volume to be added.
     *  @param pos Position of the volume.
     *  @return Reference to the newly added volume.
     */
    element_Type & setVolume( element_Type const & v, UInt const pos );
    inline element_Type & setElement( element_Type const & elem, UInt const pos ) {return setFace( elem, pos );}

    //! set numVolumes counter.
    void setVolumeCounter();
    inline void setElementCounter() {setVolumeCounter();}

    //! Reference to last volume stored in list.
    /**
     *  Reference to last volume stored in list.
     *  Useful for mesh readers.
     *  @return reference of the last volume in the list.
     */
    element_Type & lastVolume();

    //! i-th mesh 3D Element.
    /**
     *  @param i index of the mesh 3D Element.
     *  @return the i-th volume.
     */
    element_Type const & volume( UInt const i ) const;

    //! i-th mesh 3D Element.
    /**
     *  @param i index of the mesh volume.
     *  @return reference to the ith mesh volume.
     */
    element_Type & volume( UInt const i );

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

    //! Is the array for local Faces set up?
    /**
     *  It does not use switches, but interrogates the container directly
     *
     *  @return true if array for local Faces in set up.
     */
    bool hasLocalFaces() const;
    inline bool hasLocalFacets() const {return hasLocalEdges();}

    //! Build localFacetId table and optionally fills the list of Faces.
    /**
     *  @param createFaces is set true if we want also to create the actual list
     *  of internal faces. There is another utility (in mesh_util.h), which
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
    void updateElementFaces( bool /*createFaces*/ = false, const bool /*verbose*/ = false, UInt /*estimateFaceNumber*/ = 0 ){};
    void updateElementFacets( bool createFacets = false, const bool verbose = false, UInt estimateFacetNumber = 0 )
     		{updateElementEdges( createFacets, verbose, estimateFacetNumber );}

    //! Destroys element-to-face container. Useful to save memory!
    void cleanElementFacets();
    inline void cleanElementFaces(){cleanElementFacets();}

    //! Local Face Id.
    /** @param volId Id of volume (element).
     *  @param locF local face number 0 \< LocF \< numLocalFaces().
     *  @return ID of the face.
     */
    UInt localFacetId( const UInt elemId, const UInt locE ) const {return localEdgeId( elemId, locE );}
    inline UInt localFaceId( UInt const volId, UInt const locF ) const {return UInt();}

    //! Local Face Id.
    /** @param iv Reference to a volume (element).
     *  @param locF local face number 0 \< LocF \< numLocalFaces().
     *  @return ID of the face.
     */
    UInt localFacetId( const element_Type & iEl, UInt const locF ) const;
    inline UInt localFaceId( const element_Type & iv, UInt const locF ) const 
    	{return localFacetId( iv, locF ); }
    

    //! Is the array for local Edges set up?
    /**
     *  It does not use switches, but interrogates the container directly.
     *
     *  @return true if edges information are available, false otherwise.
     */
    bool hasLocalEdges() const;
    bool hasLocalRidges() const { return true;}

    //! Build localEdgeId table and optionally fills the list of Edges
    /**
     *
     * @param createEdges is set true if we want also to create the actual list
     *  of edges. There is another utility (MeshChecks.hpp), which
     *  might be used for the same purpose if we want just to create the faces
     *  and not also the LocalEdgeID table.
     *  @param verbose If true, output is verbose.
     *  @param estimateEdgeNumber is a guess provided by the user of the total
     *  number of edges. It is relevant only when createFaces=true. Setting it
     *  to a proper value helps in reducing time and memory.
     *  @param renumber Relevant only if createFaces=true.It makes sure that boundary edges are first
     *  if set to false possibly existing edges are never moved
     *
     *  @note This method does not assume that boundary edges are stores, since
     *  this condition is NOT a a paradigm for a RegionMesh2D.
     */
    void updateElementEdges( bool createEdges = false, const bool verbose = false,
                             UInt estimateEdgeNumber = 0, bool renumber=true);
    inline void updateElementRidges( bool createEdges = false, const bool verbose = false,
                             UInt estimateEdgeNumber = 0, bool renumber=true){};
       //    updateElementEdges( createEdges, verbose, estimateEdgeNumber, renumber); }

    //! Destroys Edge-To-Face lookup table.
    void cleanElementRidges(){};
    inline void cleanElementEdges() {cleanElementRidges();}

    //! Local Edge.
    /** @param volId Id of volume (element).
     *  @param locE local edge number 0 \< LocE \< numLocalEdges().
     *  @return ID of the edge.
     */
    UInt localEdgeId( const UInt elemId, const UInt locE ) const;
    inline UInt localRidgeId( UInt const elemId, UInt const locE )
		const {return element(elemId).point(locE).localId();}

    //! Local Edge.
    /** @param iv Reference of the volume.
     *  @param locE local edge number 0 \< LocE \< numLocalEdges().
     *  @return ID of the edge.
     */
    UInt localEdgeId( const element_Type& iElem, const UInt locE ) const;
    inline UInt localRidgeId( const element_Type & elem, UInt const locE )
		const {return element(elem.id()).point(locE).localId();}

    //! Returns Global-to-Local Node map.
    /**
     *  @return Global-to-Local map.
     */
    std::map<int,int> & globalToLocalNode() {return M_globalToLocalNode;}
    std::map<int,int> & globalToLocalPeak() {return M_globalToLocalNode;}
    //! Returns Local-to-Global Node map.
    /**
     *  @return Local-to-Global map.
     */
    std::map<int,int> & localToGlobalNode() {return M_localToGlobalNode;}
    std::map<int,int> & localToGlobalPeak() {return M_localToGlobalNode;}

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
    void setNumFaces( UInt n) {M_numFaces = n; }
    inline UInt numFacets() const {return numEdges();}

    //! Returns Global Number of Faces
    /**
     *  Returns global number of Face elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global Number of Faces.
     */
    UInt numGlobalFaces() const;
    inline UInt numGlobalFacets() const {return numGlobalEdges();}

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
    void setMaxNumFaces( UInt const n, bool const setcounter = false );
    void setMaxNumFacets( UInt const n, bool const setcounter = false ) {setMaxNumEdges( n, setcounter );}

    //! Changes Current capacity of Global Faces.
    /**
     *  Changes Current capacity of Global Faces (Optionally sets internal counter).
     *
     *  @param n maximum number of global faces.
     */
    void setMaxNumGlobalFaces( UInt const n );
    void setMaxNumGlobalFacets( UInt const n ) {setMaxNumGlobalEdges( n );}

    //! Adds faces.
    /**
     *  Adds a face (optionally a boundary face). Id computed automatically.
     *  @param boundary true if it's a boundary face.
     *  @return Reference to added face.
     */
    element_Type & addFace();

    //! Adds faces.
    /**
     *  Adds faces. Id computed automatically.
     *  @param v Face to add.
     *  @return Reference to the newly added face.
     */
    element_Type & addFace( element_Type const & v );

    //! Add face in a certain position.
    /**
     *  Add face to a specified position.
     *  @param v Face to add.
     *  @param pos Position of the face.
     *  @return Reference to the newly added face.
     */
    face_Type & setFace( face_Type const & f, UInt pos);

    //! Reference to last face stored in list.
    /**
     *  Reference to last face stored in list.
     *  Useful for mesh readers.
     *  @return reference of the last face in the list.
     */
    element_Type & lastFace();

    //! i-th mesh 2Delement.
    /**
     *  @param i index of the mesh 2Delement.
     *  @return the i-th face.
     */
    element_Type const & face( UInt const i ) const;

    //! i-th mesh 2Delement.
    /**
     *  @param i index of the mesh face.
     *  @return reference to the ith mesh face.
     */
    element_Type & face( UInt const i );

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
     *  look at the documentation of th eanalogous Edges methods.
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
    UInt numRidges() const {return numVertices(); }

    //! Global number of Edges.
    /**
     *  Returns global number of Edge elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global number of Edges.
     */
    UInt numGlobalEdges() const;
    UInt numGlobalRidges() const {return numGlobalVertices();}

    //! Number of local edges for each (2D) element.
    /**
     * @return Number of local edges for each (2D) element
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
    void setMaxNumEdges( UInt const n, bool const setcounter = false );
    void setMaxNumRidges( UInt const /*n*/, bool const /*setcounter*/ = false ) {};
    //! Changes Current capacity of Global Edges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of global edges.
     */
    void setMaxNumGlobalEdges( UInt const n );
    void setMaxNumGlobalRidges( UInt const /*n*/ ) {};

    //! Adds Edges.
    /**
     *  Adds edges. Id computed automatically.
     *
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    facet_Type & addEdge( bool const boundary = false );
    facet_Type & addFacet( bool const boundary = false ) {return addEdge(boundary);}

    //! Adds Edges.
    /**
     *  Adds a edge (optionally a boundary edge) to the end of the list
     *  and adjourn its ID.
     *
     *  @param f Edge to add.
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    facet_Type & addEdge( facet_Type const & f );
    facet_Type & addFacet( facet_Type const & f ){return addEdge(f);}

    //! Adds Edges to specified position.
    /**
     *  Adds a edge (optionally a boundary edge) and adjourn its Id.
     *
     *  @param f Edge to add.
     *  @param position Position of the edge.
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    facet_Type & setEdge( facet_Type const & f, UInt position, bool const boundary = false );

    //! Reference to last edge stored in list.
    /**
     *  Useful for mesh readers
     *  @return reference of the last edge in the list.
     */
    facet_Type & lastEdge();

    //! i-th mesh Edge.
    /**
     *  Returns the i-th Edge.
     *
     *  @param i Index of the mesh Edge.
     *  @return The i-th Edge.
     */
    facet_Type const & edge( UInt const i ) const;
    facet_Type const & facet( UInt const i ) const {return edge(i);}

    //! i-th mesh 1D Edge reference.
    /**
     *  Returns a reference to the i-th mesh Edge.
     *
     *  @param i Index of the mesh 1D Edge.
     *  @return Reference to the i-th Edge.
     */
    facet_Type & edge( UInt const i );
    facet_Type & facet( UInt const i ) {return edge(i);}

    //! i-th mesh 1D Boundary Edge.
    /**
     *  Returns the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return i-th Boundary Edge.
     */
    facet_Type const & boundaryEdge( UInt const i ) const;

    //! i-th mesh 1D Boundary Edge Reference.
    /**
     *  Returns a reference to the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return Reference to the i-th Boundary edge.
     */
    facet_Type & boundaryEdge( UInt const i );

    //! Set boundary Edge counter.
    /**
     *  Set the Boundary Edges counter to a given number.
     *
     *  @param n Count of Boundary Edge.
     */
    void setNumBEdges( UInt const n );
    void setNumBFacets( UInt const n ) { setNumBEdges(n);}
    void setNumBRidges( UInt const /*n*/ ) {}

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
     *  @param f The Edge.
     *  @return true if the edge is on the boundary, false otherwise.
     */
    bool isBoundaryFacet( facet_Type const & f ) const;

    //! Edge on boundary check by id.
    /**
     *  Is this edge, of given id, on boundary?
     *
     *  @param id Id of the edge.
     *  @return true if the edge is on the boundary, false otherwise.
     */
    bool isBoundaryFacet( UInt const & id ) const;
    bool isBoundaryRidge( UInt const & id ) const {return isBoundaryPoint(id);}

    //! Full Edge check by id.
    /**
     *  Does this ID corresponds to a full edge?
     *
     *  A FULL EDGE is a 1D Element that is actually stored in the Edge container.
     *
     *  @param id Id of the edge.
     *  @return true if the edge is actually stored, false otherwise.
     */
    bool isFullEdge( UInt const & id ) const;

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
    void setNumEdges ( UInt const n ) {M_numEdges = n;};
    void setNumRidges ( UInt const n ) {setNumVertices ( n );}

    //! Do I store internal edges?
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
    UInt & numPoints();

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
    void setMaxNumPoints( UInt const n, bool const setcounter = false );

    //! Set the number of storable global points in the mesh.
    /**
     *  Set the internal counter of storable global points in the mesh.
     *
     *  @param n Maximum number of storable global points in the mesh.
     */
    void setMaxNumGlobalPoints( UInt const n);

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    point_Type & addPoint( bool const boundary, bool const vertices = false);
    point_Type & addRidge( bool const boundary, bool const vertices = false){return addPoint(boundary, vertices);}

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param p Point to be added.
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    point_Type & addPoint( point_Type const & p );
    point_Type & addRidge( point_Type const & p ) {return addPoint(p);}

    //! Adds a Point in the mesh giving an id.
    /**
     *  Adds a Point inside the mesh giving an id, eventually specifing if it's a boundary point or a vertex.
     *
     *  @param p Point to be added.
     *  @param position Desired id.
     *  @param boundary If true, it's a boundary point, otherwise not (default).
     *  @param vertices If true, it's a vertex, otherwise not (default).
     *  @return Reference to the newly added Point.
     */
    point_Type & setPoint( point_Type const & p, UInt const position );

    //! Chang boundary flag of a given point
    /**
     * This method is required because in the present implementation of RegionMesh we keep
     * track of the boundary points, so just doing point(pos).setBoundary(boundary) will not
     * work since the internal list will  not be updated. We remind, however, that MeshEntityContainer
     * methods allow to extract all entities with a given flag set, so you may use that
     * method to extract boundary points irrespectively of the information stored in the
     * mesh boundary point list.
     *
     * @param position The position in the list of the point to be changed
     * @param boundary true or false if the point is or is not on the boundary
     * @return a reference to the point
     */
    point_Type & changePointBoundaryFlag(UInt const & position, bool const boundary);
    //! Returns the last mesh Point.
    /**
     *  Returns the last Point in the mesh.
     *
     *  @return Reference to the last mesh Point.
     */
    point_Type & lastPoint();

    //! Returns the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *
     *  @param i Id of the Point.
     *  @return i-th mesh Point.
     */
    point_Type const & point( UInt const i ) const;
    point_Type const & ridge( UInt const i ) const {return point(i);}
    inline point_Type const & peak( UInt const i ) const {return point(i);}

    //! Returns a reference to the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *
     *  @param i Id of the Point.
     *  @return Reference i-th mesh Point.
     */
    point_Type & point( UInt const i );
    point_Type & ridge( UInt const i ) {return point(i);}
    inline peak_Type & peak( UInt const i ) {return peak_Type();}

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    point_Type const & boundaryPoint( UInt const i ) const;

    //! Returns a reference to the i-th mesh Boundary Point.
    /**
     *  Returns the i-th Boundary Point in the mesh.
     *
     *  @param i Id of the Boundary Point.
     *  @return Reference i-th mesh Boundary Point.
     */
    point_Type & boundaryPoint( UInt const i );

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
    void setNumBPoints( UInt const n );

    //! Is this point on boundary?
    /**
     *  Is this point on boundary?
     *
     *  @param p The Point.
     *  @return true if the point is on the boundary, false otherwise.
     */
    bool isBoundaryPoint( point_Type const & p ) const ;

    //! Is this point on boundary?
    /**
     *  Is this point, of given id, on boundary?
     *
     *  @param id Id of the point.
     *  @return true if the point is on the boundary, false otherwise.
     */
    bool isBoundaryPoint( UInt const & id ) const;

    //! List of points
    /**
     *  @param fct Function of three double arguments.
     *  @param list_pts List of Points.
     */
    void getListOfPoints( bool ( *fct ) ( double, double, double ), std::vector<UInt>& list_pts );

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
    UInt & numBVertices();

    //! Number of local vertices for each (2D) element.
    /**
     *  @return Number of local vertices for each (2D) element.
     */
    UInt numLocalVertices() const;

    //! Returns the global number of vertices in the mesh.
    /**
     *  Interrogation to the counter of vertices in the mesh.
     *
     *  @return Number of vertices in the mesh.
     */
    UInt numGlobalVertices() const;
    inline UInt numGlobalPeaks() const { return numGlobalVertices();}


    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param id Point's id.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex ( UInt const & id ) const;

    //! Vertex check.
    /**
     *  Is this Point a Vertex?
     *
     *  @param p A Point.
     *  @return true if the Point is a Vertex, false otherwise.
     */
    bool isVertex ( point_Type const & p ) const;

    //! Changes number of Vertices.
    /**
     *  Allows to change number of vertices in Region.
     *
     *  @param n Number of vertices in the mesh.
     */
    void setNumVertices(UInt const n);

    //! Set the number of vertices in the mesh.
    /**
     *  Set the internal counter of vertices points in the mesh.
     *
     *  @param n Number of vertices in the mesh.
     */
    void setNumGlobalVertices( UInt const n ) { M_numGlobalVertices = n; }

    //! Changes number of Boundary Vertices.
    /**
     *  Allows to change number of boundary vertices in Region.
     *
     *  @param n Number of boundary vertices in the mesh.
     */
    void setNumBVertices(UInt const n);

    /** @} */ // End of group Vertices Methods


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
    //! Boundary points list.
    std::vector<point_Type * > _bPoints;
    /** @} */ // End of group Region Containers


    /** @name Switches
     *  @ingroup public_attributes
     *
     *  @{
     */

    //! Switches
    Switch switches;

    /** @} */ // End of group Switches

    elements_Type& elementList() {return faceList;}

    facets_Type& facetList() {return edgeList;}

    ridges_Type& ridgeList() {return pointList;}

    int inline dimension() const {return 2;}


private:

    /*! Arrays containing the ids of Edges and Faces of each element
      I use a Define to use localto global array or directly the
      bareedges */
    ArraySimple<UInt> M_FToE;

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

    std::map<int, int>      M_globalToLocalNode;
    std::map<int, int>      M_localToGlobalNode;
    std::map<int, int>      M_globalToLocalEdge;
    std::map<int, int>      M_globalToLocalFace;
    std::map<int, int>      M_globalToLocalVolume;

    MeshUtility::MeshTransformer<RegionMesh2D<GEOSHAPE, MC> > M_meshTransformer;

}; // End of class RegionMesh2D


// =================================================== //
// =================================================== //
//                    IMPLEMENTATION                   //
// =================================================== //
// =================================================== //


void set_switches_for_regionmesh( Switch & sw );


template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::RegionMesh2D() :
        MeshEntity(),
        MC::regionMarker_Type(),
        switches(),
        M_numVertices( 0 ),
        M_numBVertices( 0 ),
        M_numPoints( 0 ),
        M_numBPoints( 0 ),
        M_numEdges( 0 ),
        M_numBEdges( 0 ),
        M_numFaces( 0 ),
        M_globalToLocalNode(),
        M_localToGlobalNode(),
        M_globalToLocalEdge(),
        M_globalToLocalFace(),
        M_globalToLocalVolume(),
        M_meshTransformer(*this)
{
	set_switches_for_regionmesh( switches );
}


template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::RegionMesh2D( UInt id ) :
        MeshEntity( id ),
        MC::regionMarker_Type(),
        switches(),
        M_numVertices( 0 ),
        M_numBVertices( 0 ),
        M_numPoints( 0 ),
        M_numBPoints( 0 ),
        M_numFaces( 0 ),
        M_numEdges( 0 ),
        M_numBEdges( 0 ),
        M_meshTransformer(*this)
{
	set_switches_for_regionmesh( switches );
}

template <typename GEOSHAPE, typename MC>
RegionMesh2D<GEOSHAPE, MC>::~RegionMesh2D() {}

template <typename GEOSHAPE, typename MC>
inline
void
RegionMesh2D<GEOSHAPE, MC>::setLinkSwitch( std::string const & _s )
{
    ASSERT0( switches.set( _s ), std::stringstream( "Switch named " + _s + " is not allowed" ).str().c_str() );
};

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::getLinkSwitch( std::string const & _s ) const
{
    return switches.test( _s );
};

template <typename GEOSHAPE, typename MC>
inline
void
RegionMesh2D<GEOSHAPE, MC>::unsetLinkSwitch( std::string const & _s )
{
    std::ostringstream _err_msg;
    _err_msg << "Switch named " << _s << " is not allowed";
    ASSERT0( switches.unset( _s ), _err_msg.str().c_str() );
};


template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numLocalVertices() const
{
    return element_Type::S_numLocalVertices;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numLocalEdges() const
{
    return element_Type::S_numLocalEdges;
}

// ************** Generic Methods
template <typename GEOSHAPE, typename MC>
UInt RegionMesh2D<GEOSHAPE, MC>::numElements() const
{
    return M_numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt &RegionMesh2D<GEOSHAPE, MC>::numElements()
{
    return M_numFaces;
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh2D<GEOSHAPE, MC>::numGlobalElements() const
{
    return M_numGlobalFaces;
}

template <typename GEOSHAPE, typename MC>
inline UInt & RegionMesh2D<GEOSHAPE, MC>::numGlobalElements()
{
    return M_numGlobalFaces;
}

template <typename GEOSHAPE, typename MC>
UInt RegionMesh2D<GEOSHAPE, MC>::numBFacets() const
{
    return M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt &RegionMesh2D<GEOSHAPE, MC>::numBFacets()
{
    return M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::element_Type &
RegionMesh2D<GEOSHAPE, MC>::element( UInt const & i )
{
    return face( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::element_Type const &
RegionMesh2D<GEOSHAPE, MC>::element( UInt const & i ) const
{
    return face( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type &
RegionMesh2D<GEOSHAPE, MC>::bFacet( UInt const & i )
{
    return boundaryEdge( i );
}

template <typename GEOSHAPE, typename MC>
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type const &
RegionMesh2D<GEOSHAPE, MC>::bFacet( UInt const & i ) const
{
    return boundaryEdge( i );
}


//-------------------------------------------------------------------

// ***************************** FACES
template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numFaces() const
{
    return M_numFaces;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numGlobalFaces() const
{
    return M_numGlobalFaces;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::maxNumFaces() const
{
    return maxNumItems( faceList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedFaces() const
{
    return numItems( faceList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumFaces( UInt const n, bool const setcounter )
{
	faceList.setMaxNumItems(n);
    if ( setcounter )
        M_numFaces = n;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumGlobalFaces( UInt const n)
{
    M_numGlobalFaces = n;
}

// \todo use addItem
template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::element_Type &
RegionMesh2D<GEOSHAPE, MC>::addFace()
{
    return addFace( element_Type() );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::element_Type &
RegionMesh2D<GEOSHAPE, MC>::addFace( element_Type const & v )
{
    ASSERT_PRE( faceList.size() < faceList.capacity() , "Face list size exceeded" <<
                faceList.size() + 1 << " " << faceList.capacity() ) ;
    faceList.push_back( v );
    ( faceList.back() ).setId(faceList.size() - 1);
    return faceList.back();
}
// \todo Use setItem

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::face_Type &
RegionMesh2D<GEOSHAPE, MC>::setFace( face_Type const & v, UInt const pos )
{
    ASSERT_PRE( pos < faceList.capacity() , "position requested exceed capacity" <<
                pos << " " << faceList.capacity() ) ;
    faceList( pos ) = v;
    faceList( pos ).setId(pos);
    return faceList( pos );
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::element_Type &
RegionMesh2D<GEOSHAPE, MC>::lastFace()
{
    return faceList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::element_Type const &
RegionMesh2D<GEOSHAPE, MC>::face( UInt const i ) const
{
    ASSERT_BD( i < faceList.size() ) ;
    return faceList( i );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::element_Type &
RegionMesh2D<GEOSHAPE, MC>::face( UInt const i )
{
    ASSERT_BD( i < faceList.size() ) ;
    return faceList( i );
}

// ************************* EDGES ******************************
template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numEdges() const
{
    return M_numEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numGlobalEdges() const
{
    return M_numGlobalEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedEdges() const
{
    return numItems( edgeList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::maxNumEdges() const
{
    return maxNumItems( edgeList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumEdges( UInt const n, bool const setcounter )
{
	edgeList.setMaxNumItems(n);
	if ( setcounter )
		M_numEdges = n;
}


template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumGlobalEdges( UInt const n)
{
    M_numGlobalEdges = n;
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type &
RegionMesh2D<GEOSHAPE, MC>::addEdge( bool const boundary )
{
	facet_Type aEdge;
	aEdge.setBoundary(boundary);
    return this->addEdge( aEdge );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type &
RegionMesh2D<GEOSHAPE, MC>::addEdge( facet_Type const & f )
{
    ASSERT_PRE( edgeList.size() < edgeList.capacity(), "Edge list size exceeded" <<
                edgeList.size() + 1 << " " << edgeList.capacity() ) ;
    edgeList.push_back( f );
    edgeList.back().setId( edgeList.size()-1 );

    return edgeList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type &
RegionMesh2D<GEOSHAPE, MC>::setEdge( facet_Type const & f, UInt position, bool const boundary )
{
    ASSERT_PRE( position < edgeList.capacity(), "Edge list size exceeded" <<
                position << " " << edgeList.capacity() ) ;
    edgeList( position ) = f;
    edgeList( position ).setId( position );
#ifdef NOT_BDATA_FIRST

    if ( boundary )
    {
        ASSERT_PRE( position < _bEdges.capacity(), "Boundary Edge list size exceeded" <<
                    _bEdges.size() << " " << bEdges.capacity() ) ;
        _bEdges.push_back( &( edgeList( position ) ) );
    }
#endif
    return edgeList( position );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type &
RegionMesh2D<GEOSHAPE, MC>::lastEdge()
{
    return edgeList.back();
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type const &
RegionMesh2D<GEOSHAPE, MC>::edge( UInt const i ) const
{
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList( i );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type &
RegionMesh2D<GEOSHAPE, MC>::edge( UInt const i )
{
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList( i );
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type const &
RegionMesh2D<GEOSHAPE, MC>::boundaryEdge( UInt const i ) const
{
#ifdef NOT_BDATA_FIRST
    ASSERT_PRE( _bEdges.size() != 0 " Boundary Edges not Stored" ) ;
    ASSERT_BD( i < _bEdges.size() ) ;
    return *( _bEdges( i ) );
#else

    ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList( i );
#endif
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::facet_Type &
RegionMesh2D<GEOSHAPE, MC>::boundaryEdge( UInt const i )
{
#ifdef NOT_BDATA_FIRST
    ASSERT_PRE( _bEdges.size() != 0 " Boundary Edges not Stored" ) ;
    ASSERT_BD( i < _bEdges.size() ) ;
    return *( _bEdges( i ) );
#else

    ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList( i );
#endif
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setNumBEdges( UInt const n )
{
    M_numBEdges = n;
#ifdef NOT_BDATA_FIRST

    M_bEdges.reserve( n );
#endif
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::hasEdges() const
{
    return ! edgeList.empty();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::hasInternalEdges() const
{
    return edgeList.size() > M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBEdges() const
{
    return M_numBEdges;
}




template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryFacet( facet_Type const & e ) const
{
#ifdef NOT_BDATA_FIRST
    //ASSERT(false,"In this version Boundary edges must be stored first");
    bool isboundary = true;
    for ( UInt k = 0; k < facet_Type::S_numVertices; ++k )
    {
        isboundary = isboundary & e.point( k ).boundary();
    }
    return isboundary;
#else

    return e.id() < M_numBEdges;
#endif
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryFacet( UInt const & id ) const
{
    return isBoundaryFacet( edge( id ) );
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isFullEdge( UInt const & id ) const
{
    return edgeList.size() >= id;
}

/*
template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh2D<GEOSHAPE, MC>::edgeElement( UInt const i, UInt const Pos ) const
{
    ASSERT_PRE( i < edgeList.size(), "Not enough faces stored" ) ;
    ASSERT_BD( i > 0 ) ;
    return edgeElement( edge( i ), Pos );
};

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh2D<GEOSHAPE, MC>::edgeElement( facet_Type const & f, UInt const Pos ) const
{
    ASSERT_BD( ! edgeList.empty() ) ;
    ASSERT_PRE( Pos <= 1 , "Wrong position (0 or 1)" ) ;
    if ( Pos == 0 )
    {
        return f.firstAdjacentElementIdentity();
    }
    else
    {
        return f.secondAdjacentElementIdentity();
    }
};
*/

// ************************ Points/Vertices

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numPoints() const
{
    return M_numPoints;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numPoints()
{
    return M_numPoints;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedPoints() const
{
    return numItems( pointList );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::storedBPoints() const
{
    return _bPoints.size();
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::maxNumPoints() const
{
    return maxNumItems( pointList );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumPoints( UInt const n, bool const setcounter )
{
    pointList.setMaxNumItems(n);
    if ( setcounter )
        M_numPoints = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh2D<GEOSHAPE, MC>::setNumVertices( UInt const n )
{
    M_numVertices = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh2D<GEOSHAPE, MC>::setNumBVertices( UInt const n )
{
    M_numBVertices = n;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setMaxNumGlobalPoints( UInt const n )
{
    M_numGlobalPoints = n;
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type &
RegionMesh2D<GEOSHAPE, MC>::addPoint( bool const boundary, bool const vertex )
{
    point_Type aPoint;
    aPoint.setBoundary(boundary);

    if(vertex) aPoint.replaceFlag(aPoint.flag() | EntityFlags::VERTEX);
    return addPoint(aPoint);
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type &
RegionMesh2D<GEOSHAPE, MC>::addPoint( point_Type const & p )
{
    ASSERT_PRE( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
            pointList.size() + 1 << " " << pointList.capacity() ) ;

    pointList.push_back( p );

    point_Type * pp = & pointList.back();
    pp->setId( pointList.size() - 1 );

    if ( pp->boundary() )
    {
       ASSERT_PRE( _bPoints.size() < _bPoints.capacity(), "Boundary Point list size exceeded" <<
                  _bPoints.size() + 1 << " " << _bPoints.capacity() ) ;
        _bPoints.push_back( pp );
        pp->setBoundary(true);
    }
    return pointList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type &
RegionMesh2D<GEOSHAPE, MC>::setPoint
( point_Type const & p, UInt position )
{
    ASSERT_PRE( position < pointList.size(), "Position  exceed lpoint list size" <<
                position << " " << pointList.size() ) ;

    bool setToBoundary=p.boundary();
    bool originalBoundary=pointList[position].boundary();

    pointList [position]=p;
    point_Type * pp = & pointList[position];
    if (setToBoundary!=originalBoundary){
        if(setToBoundary){
            // add to list of boundary points
            _bPoints.push_back( pp );
        }
        else
        {
            // This is rather complex, since I do not know a priori
            // if and where point was already stored in the list!
            // No way to avoid it, sorry
            typename std::vector<point_Type *>::iterator bp;
            for (bp = _bPoints.begin(); bp != _bPoints.end(); ++bp )
            {
                if ( ( *bp ) ->id() == position )
                {
                    _bPoints.erase(bp);
                    break;
                }
            }
        }
        return pointList[position];
    }
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type &
RegionMesh2D<GEOSHAPE, MC>::lastPoint()
{
    return pointList.back();
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type const &
RegionMesh2D<GEOSHAPE, MC>::point( UInt const i ) const
{
    ASSERT_BD( i < pointList.size() ) ;
    return pointList( i );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type &
RegionMesh2D<GEOSHAPE, MC>::point( UInt const i )
{
    ASSERT_BD( i < pointList.size() ) ;
    return pointList( i );
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type const &
RegionMesh2D<GEOSHAPE, MC>::boundaryPoint( UInt const i ) const
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i < _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh2D<GEOSHAPE, MC>::point_Type &
RegionMesh2D<GEOSHAPE, MC>::boundaryPoint( UInt const i )
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i < _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numBPoints() const
{
    return M_numBPoints;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::setNumBPoints( UInt const n )
{
    M_numBPoints = n;
    _bPoints.reserve( n );
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numVertices() const
{
    return M_numVertices;
}

template <typename GEOSHAPE, typename MC>
UInt
RegionMesh2D<GEOSHAPE, MC>::numGlobalVertices() const
{
    return M_numGlobalVertices;
}

template <typename GEOSHAPE, typename MC>
UInt &
RegionMesh2D<GEOSHAPE, MC>::numBVertices()
{
    return M_numBVertices;
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isVertex( point_Type const & p ) const
{
    return p.id() < M_numVertices;
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isVertex( UInt const & id ) const
{
    return id < M_numVertices;
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryPoint( UInt const & id ) const
{
    return point( id ).boundary();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryPoint( point_Type const & p ) const
{
    return p.boundary();
}

/*
template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryRidge( UInt const & id ) const
{
    return point( id ).boundary();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::isBoundaryRidge( point_Type const & p ) const
{
    return p.boundary();
}
*/
/********************************************************************************
                        Element Adiacency Methods
*******************************************************************************/

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh2D<GEOSHAPE, MC>::hasLocalEdges() const
{
    return ! M_FToE.empty();
}

#ifdef SAVEMEMORY
class BareEdge;

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( const element_Type & ifac, UInt const locE ) const
{
    ASSERT_PRE( !M_FToE.empty(), "Face to Edges array not  set" );
    ASSERT_BD( locE < element_Type::S_numLocalEdges );
    std::pair<BareEdge, bool> it;
    UInt i1, i2;
    i1 = GEOSHAPE::edgeToPoint( locE, 0 );
    i2 = GEOSHAPE::edgeToPoint( locE, 1 );
    i1 = ( ifac.point( i1 ) ).id();
    i2 = ( ifac.point( i2 ) ).id();
    it = makeBareEdge( i1, i2 );
    return M_FToE.id( it.first );
}

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( UInt const facId, UInt const locE ) const
{
    ASSERT_BD( facId < M_numFaces );
    return localEdgeId( face( facId ), locE );
}

#else

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( const element_Type & ifac, UInt const locE )
const
{
    return M_FToE( ifac.id(), locE );
}


template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh2D<GEOSHAPE, MC>::localEdgeId( UInt const facId, UInt const locE )
const
{
    ASSERT_PRE( !M_FToE.empty(), "Face to Edges array not  set" );
    ASSERT_BD( facId < M_numFaces );
    ASSERT_BD( locE < element_Type::S_numLocalEdges );
    return M_FToE( locE, facId );
}


#endif

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::updateElementEdges( bool ce, const bool verbose, UInt ee, bool /*verbose*/ )
{
	if (verbose)
		std::cout << "     Updating element edges ... " << std::flush;

    ASSERT0( ! ce || M_numBEdges > 0, std::stringstream( std::string("Boundary Edges Must have been set") +
                                                         std::string("in order to call updateElementEdges with createEdges=true") +
                                                         std::string("\nUse buildBoundaryEdges(..) from mesh_util.h") ).str().c_str() );
    // If the counter is set we trust it! Otherwise we use Euler formula

    if ( ce && ee == 0 )
        //check!!!
        ee = M_numEdges > M_numBEdges ? M_numEdges : M_numFaces + M_numVertices -1;
    ASSERT( ce || numEdges() > 0 , "Mesh is not properly set!" );

    if ( ce )
        edgeList.reserve( ee );

    facet_Type edg;

    MeshElementBareHandler<BareEdge> _be;
    // Extra map for faces stored which are not boundary faces
    MeshElementBareHandler<BareEdge> _extraEdges;
    std::pair<UInt, bool> e;
    M_FToE.reshape( numLocalEdges(), numFaces() ); // DIMENSION ARRAY

    UInt fid, i1, i2;
    std::pair<BareEdge, bool>_edge;
    GEOSHAPE ele;
    // If we have all edges and the edges store all adjacency info
    // everything is easier
    if ( (edgeList.size() == numEdges()) & getLinkSwitch( "FACETS_HAVE_ADIACENCY" ) & getLinkSwitch( "HAS_ALL_FACETS" ) )
    {
        for ( typename edges_Type::iterator ite = edgeList.begin(); ite != edgeList.end(); ++ite )
        {
        	if ( ite->firstAdjacentElementPosition() != NotAnId && ite->firstAdjacentElementIdentity() != NotAnId)
                M_FToE( ite->firstAdjacentElementPosition() , ite->firstAdjacentElementIdentity() ) = ite->localId();

            if ( ite->secondAdjacentElementPosition() != NotAnId && ite->secondAdjacentElementIdentity() != NotAnId)
                M_FToE( ite->secondAdjacentElementPosition(), ite->secondAdjacentElementIdentity() ) = ite->localId();
        }
        // we finish here
        setLinkSwitch( "HAS_ELEMENT_TO_FACETS" );
        if (verbose)
        	std::cout << " done." << std::endl;
        return ;
    }


    // If I have only boundary faces I need to process them first to keep the correct numbering

    // First We check if we have already Faces stored
    UInt _numOriginalStoredEdges=edgeList.size();
    if ( ! edgeList.empty() )
    {
        // dump all faces in the container, to maintain the correct numbering
        // if everything is correct the numbering in the bareface structure
        // will reflect the actual face numbering However, if I want to create
        // the internal faces I need to make sure that I am processing only the
        // boundary ones. So I resize the container!
        std::pair<UInt, bool> _check;
        for ( UInt j = 0; j < edgeList.size(); ++j )
        {
            i1 = ( edgeList[ j ].point( 0 ) ).localId();
            i2 = ( edgeList[ j ].point( 1 ) ).localId();
            _edge = makeBareEdge( i1, i2);
            _check = _be.addIfNotThere( _edge.first );
            if (j>=this->M_numBFaces)_extraEdges.addIfNotThere( _edge.first, j);
        }
    }


    for ( typename faces_Type::iterator iface = faceList.begin();
            iface != faceList.end(); ++iface )
    {
        fid = iface->localId();
        for ( UInt j = 0; j < numLocalEdges(); j++ )
        {
            i1 = ele.edgeToPoint( j, 0 );
            i2 = ele.edgeToPoint( j, 1 );
            i1 = ( iface->point( i1 ) ).localId();
            i2 = ( iface->point( i2 ) ).localId();
            _edge = makeBareEdge( i1, i2);
            e = _be.addIfNotThere( _edge.first );
            M_FToE( j, fid ) = e.first;

            bool _isBound=e.first<this->M_numBEdges;
            // Is the edge an extra edge (not on the boundary but originally included in the list)?
            bool _isExtra = (e.first >=this->M_numBEdges  && e.first < _numOriginalStoredEdges);
            if (_isBound)
            {
                edge_Type & _thisEdge(edgeList[e.first]);
                _thisEdge.firstAdjacentElementIdentity()   = fid;
                _thisEdge.firstAdjacentElementPosition()   = j;
                _thisEdge.secondAdjacentElementIdentity()  = NotAnId;
                _thisEdge.secondAdjacentElementPosition()  = NotAnId;
            }
            else if (_isExtra)
            {
                // This is not a bfaces and I need to set up all info about adjacency properly
                edge_Type & _thisEdge(edgeList[e.first]);
                // I need to check if it is the first time I meet it. Then I delete it from the
                // map: if it as there it means that it is the first time I am treating this face
                if(_extraEdges.deleteIfThere(_edge.first)){
                    // I need to be sure about orientation, the easiest thing is to rewrite the face points
                    for ( UInt k = 0; k < edge_Type::S_numPoints; ++k )
                        _thisEdge.setPoint( k, iface->point( ele.edgeToPoint( j, k ) ) );
                    _thisEdge.firstAdjacentElementIdentity()  = fid;
                    _thisEdge.firstAdjacentElementPosition()  = j;

                }else{
                    _thisEdge.secondAdjacentElementIdentity()  = fid;
                    _thisEdge.secondAdjacentElementPosition()  = j;
                }
            }
            else if ( ce ) // An edge not contained in the original list.
            {
                if ( e.second )
                {
                    // a new face It must be internal.
                    for ( UInt k = 0; k < facet_Type::S_numPoints; ++k )
                        edg.setPoint( k, iface->point( ele.edgeToPoint( j, k ) ) );
                    edg.firstAdjacentElementIdentity()  = fid;
                    edg.firstAdjacentElementPosition() = j;

                    // gets the marker from the RegionMesh

                    edg.setMarker( this->marker() );
                    //        inheritPointsWeakerMarker( edg );

                    edg.setBoundary(false);
                    addEdge( edg ); //The id should be correct
                }
                else
                {
					edgeList( e.first ).secondAdjacentElementIdentity() = fid;
					edgeList( e.first ).secondAdjacentElementPosition() = j;

                }
            }
        }
    }

    UInt n = _be.maxId();
    // Fix _numfaces if it was not set or set to just the # of Bfaces

    M_numEdges = n;  // We have found the right total number of edges in the mesh
    if(M_numGlobalEdges==0) M_numGlobalEdges=n; // If not already set fix it.

    if (verbose) std::cout << n << " edges ";
    //ASSERT_POS( n == M_numEdges , "#Edges found inconsistent with that stored in RegionMesh" ) ;
    setLinkSwitch( "HAS_ELEMENT_TO_FACETS" );
    if ( ce )
    	setLinkSwitch( "HAS_ALL_FACETS" );
  //  if ( ce ) Edges have adjacency in any case!
        setLinkSwitch( "FACETS_HAVE_ADIACENCY" );
    if (verbose)
    	std::cout << " done." << std::endl;
}


/*


template <typename GEOSHAPE, typename MC>
template <typename function>
void RegionMesh2D<GEOSHAPE, MC>::transformMesh( const function& meshMapping)
{
	for ( unsigned int i = 0; i < pointList.size();++i )
	{
		point_Type& p = pointList[ i ];
		meshMapping(p.coordinate(0),p.coordinate(1),p.coordinate(2));
	}
}


*/

template <typename GEOSHAPE, typename MC>
void
RegionMesh2D<GEOSHAPE, MC>::cleanElementFacets()
{
    M_FToE.clear();
    unsetLinkSwitch( "HAS_ELEMENT_TO_FACETS" );
}


// *************** GENERAL INFO *******************

template <typename GEOSHAPE, typename MC>
std::ostream & RegionMesh2D<GEOSHAPE, MC>::showMe( bool verbose, std::ostream & out ) const
{
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << "                      RegionMesh2D                " << std::endl;
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << " ID: " << this->id() << std::endl;
    out << "Edges local to  faces stored: " << hasLocalEdges() << std::endl;
    //out <<"Edges local to  faces   stored:"<<switches.test("FACEtoEDGE")<<std::endl;
    //out <<"Faces adjacent to Edges stored: "<<switches.test("EDGEtoFACE")<<std::endl<<std::endl;
    out << "Edges Stored: " << hasEdges() << " Internal: "
    << hasInternalEdges() << std::endl;
    out << "**************************************************" << std::endl;
    out << "numPoints=" << numPoints() << "  " << "numBPoints=" << numBPoints() << std::endl;
    out << "numVertices=" << numVertices() << "  " << "numBVerices=" << numBVertices() << std::endl;
    out << "numFaces=" << numFaces() << std::endl;
    out << "numEdges=" << numEdges() << "  " << "numBEdges=" << numBEdges() << std::endl;
    out << "**************************************************" << std::endl;
    switches.showMe( verbose, out );
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    if ( verbose )
    {
        std::cout << "Verbose version not implemented yet" << std::endl;
    }
    return out;

}
template <typename GEOSHAPE, typename MC>
int
RegionMesh2D<GEOSHAPE, MC>::check( int /*level*/, bool const fix, bool const verb, std::ostream & out )
{
    int severity = 0;
    if ( verb )
    {
        out << "**************************************************" << std::endl;
        out << "         Checkin  RegionMesh2D                " << std::endl;
        out << " ID: " << this->id() << std::endl;
        out << "**************************************************" << std::endl;
    }
    if ( pointList.size() != M_numPoints )
    {
        out << " Point list size " << pointList.size() << " not equal to internal counter value "
        << M_numPoints << std::endl;
        if ( fix )
        {
            M_numPoints = pointList.size();
            out << "Fixed";
            out.flush();
        }
    }
    if ( edgeList.size() == 0 )
    {
        if ( verb )
            out << "Warning: No Edges Stored" << std::endl;
        severity = -1;
    }
    if ( faceList.size() == 0 )
    {
        if ( verb )
            out << "Warning: No Faces Stored" << std::endl;
        severity = 1;
    }
    UInt count = 0;
    for ( typename points_Type::iterator i = pointList.begin(); i != pointList.end(); ++i )

        if ( i->boundary() )
            ++count;
    if ( !count ) severity = 4;
    if ( count != M_numBPoints )
    {
        out << " Num Boundary points " << count << " not equal to internal counter value "
        << M_numBPoints << std::endl;
        if ( fix )
        {
            M_numBPoints = count;
            out << "Fixed";
            out.flush();
        }
    }
    UInt badid = 0;
    for ( UInt i = 0; i < storedPoints(); ++i )
        if ( point( i ).id() != i )
            ++badid;
    if ( badid != 0 )
    {
        out << " SEVERITY ERROR:" << badid << "Points ids are wrong";
        severity = 5;
    }

    badid = 0;
    for ( UInt i = 0; i < storedEdges(); ++i )
        if ( edge( i ).id() != i )
            ++badid;
    if ( badid != 0 )
    {
        out << " SEVERITY ERROR:" << badid << "Edges ids are wrong";
        severity = 5;
    }

    badid = 0;
    for ( UInt i = 0; i < storedFaces(); ++i )
        if ( face( i ).id() != i )
            ++badid;
    if ( badid != 0 )
    {
        out << " SEVERITY ERROR:" << badid << "Faces ids are wrong";
        severity = 5;
    }

    badid = 0;

    if ( M_numVertices == 0 )
    {
        out << " SEVERITY ERROR: internal Vertices Counter unset";
        severity = 6;
    }
    if ( M_numPoints == 0 )
    {
        out << " SEVERITY ERROR: internal Points Counter unset";
        severity = 6;
    }
    if ( M_numPoints == 0 )
    {
        out << " SEVERITY ERROR: internal Points Counter unset";
        severity = 6;
    }
    if ( M_numBPoints == 0 )
    {
        out << " SEVERITY ERROR: boundary Points Counter unset";
        severity = 6;
    }
    if ( M_numBVertices == 0 )
    {
        out << " SEVERITY ERROR: boundary Vertices Counter unset";
        severity = 6;
    }

    if ( verb )
        out << "   Check Finished              " << std::endl <<
        "***********************************************" << std::endl;

    return severity;
}

template <typename GEOSHAPE, typename MC>
inline MeshUtility::MeshTransformer<RegionMesh2D<GEOSHAPE, MC> > &
RegionMesh2D<GEOSHAPE, MC>::meshTransformer()
{
    return this->M_meshTransformer;
}

} // End of namespace LifeV

#endif //REGIONMESH2D_H
