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
 *  @brief File containing 3D Mesh Classes
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Miguel Fernandez
 *
 *  @contributor Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 *  @mantainer Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 */

#ifndef _REGIONMESH_HH_
#define _REGIONMESH_HH_

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
 *  @class RegionMesh
 *  @brief Class for 3D Mesh
 *
 *  @author Luca Formaggia <luca.formaggia@polimi.it>
 *  @author Miguel Fernandez
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
class RegionMesh
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
    explicit RegionMesh();

    //! Default constructor
    /**
     *  @param id marker of the RegionMesh
     */
    explicit RegionMesh( UInt id );

    //! Destructor
    virtual ~RegionMesh<GEOSHAPE, MC>();

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
    inline MeshUtility::MeshTransformer<RegionMesh<GEOSHAPE, MC> > & meshTransformer();

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
    UInt numElements() const {return numElements(element_Type() );}
    UInt numElements(volume_Type) const {return numVolumes();}
    UInt numElements(face_Type) const  {return numFaces(); }
    UInt numElements(edge_Type) const  {return numEdges(); }

    //! Number of Boundary facets.
    /**
     *  @return Number of Boundary facets.
     */
    UInt numBFacets() const { return numBFacets( facet_Type() );}
    UInt numBFacets(face_Type) const { return numBFaces();}
    UInt numBFacets(edge_Type) const { return numBEdges();}

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    element_Type& element( const UInt& i ) {return element(element_Type(), i );}
    element_Type& element( volume_Type, const UInt& i ) {return volume(i); }
    element_Type& element( face_Type, const UInt& i ) {return face(i); }
    element_Type& element( edge_Type, const UInt& i ) {return edge(i); }

    //! Get element at the i-th index.
    /**
     * @param i Index of the element
     * @return Element at index i
     */
    const element_Type& element( const UInt& i ) const {return element(element_Type(), i );}
    const element_Type& element( volume_Type, const UInt& i ) const {return volume(i); }
    const element_Type& element( face_Type, const UInt& i ) const {return face(i); }
    const element_Type& element( edge_Type, const UInt& i ) const {return edge(i); }



    //! Get boundary facet at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary facet at index i
     */
    facet_Type& bFacet( const UInt& i ){ return boundaryFacet(facet_Type(), i); }
    facet_Type& boundaryFacet( const UInt& i ){ return boundaryFacet(facet_Type(), i); }
    facet_Type& boundaryFacet( face_Type, const UInt& i ){ return boundaryFace(i); }
    facet_Type& boundaryFacet( edge_Type, const UInt& i ){ return boundaryEdge(i); }

    //! Get boundary facet at the i-th index.
    /**
     * @param i Index of the element
     * @return Boundary facet at index i
     */
    const facet_Type& bFacet( const UInt& i ) const { return boundaryFacet(facet_Type(), i); }
    const facet_Type& boundaryFacet( const UInt& i ) const{ return boundaryFacet(facet_Type(), i); }
    const facet_Type& boundaryFacet( face_Type, const UInt& i ) const { return boundaryFace(i); }
    const facet_Type& boundaryFacet( edge_Type, const UInt& i ) const { return boundaryEdge(i); }
    //! Number of global elements.
    /**
     * @return Number of global Elements (Volumes).
     */
    UInt numGlobalElements() const {return numGlobalElements( element_Type() );}
    UInt numGlobalElements(volume_Type) const {return numGlobalVolumes();}
    UInt numGlobalElements(face_Type) const {return numGlobalFaces();}

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
    UInt storedVolumes() const;

    //! Current capacity of Volumes Container.
    /**
     *  @return how many elements may be stored.
     */
    UInt maxNumVolumes() const;
    inline UInt maxNumElements() const {return maxNumVolumes( element_Type() );}
    inline UInt maxNumElements(volume_Type) const {return maxNumVolumes();}
    inline UInt maxNumElements(face_Type) const {return maxNumVolumes();}
    inline UInt maxNumElements(edge_Type) const {return maxNumVolumes();}

    //! Changes Current capacity of Volumes.
    /**
     *  Changes Current capacity of Volumes (Optionally sets internal counter).
     *
     *  @param n Maximum number of volumes.
     *  @param setcounter true to set the counter, false otherwise (default).
     */
    void setMaxNumVolumes      ( UInt const n, bool const setcounter = false );
    inline void setMaxNumElements   ( UInt const n, bool const setcounter = false ) {setMaxNumElements(element_Type(), n, setcounter);}
    inline void setMaxNumElements   ( volume_Type, UInt const n, bool const setcounter = false ) {setMaxNumVolumes( n, setcounter);}
    inline void setMaxNumElements   ( face_Type, UInt const n, bool const setcounter = false ) {setMaxNumFaces( n, setcounter);}
    inline void setMaxNumElements   ( edge_Type, UInt const n, bool const setcounter = false ) {setMaxNumEdges( n, setcounter);}

    //! Changes Current capacity of Global Volumes.
    /**
     *  Changes Current capacity of Global Volumes (Optionally sets internal counter).
     *
     *  @param n maximum number of global volumes.
     */
    void setMaxNumGlobalVolumes( UInt const n );
    inline void setMaxNumGlobalElements( UInt const n ) {setMaxNumGlobalElements(element_Type(), n) ;}
    inline void setMaxNumGlobalElements( volume_Type, UInt const n ) {setMaxNumGlobalVolumes(n) ;}
    inline void setMaxNumGlobalElements( face_Type, UInt const n ) {setMaxNumGlobalFaces(n) ;}
    inline void setMaxNumGlobalElements( edge_Type, UInt const n ) {setMaxNumGlobalPoints(n) ;}


    //! Set Number of Volumes.
    /**
     *  Set number of Volume elements in the mesh by changing internal counter.
     *  @param n Number of volumes.
     */
    void setNumVolumes      ( UInt const n );
    inline void setNumElements      ( UInt const n ) {setNumVolumes(n);}

    //! Adds volumes.
    /**
     *  Adds volume. Id computed automatically.
     *  @return Reference to added volume.
     */
    inline volume_Type & addVolume() { return addVolume( volume_Type() ); }
    inline element_Type & addElement() {return addElement( element_Type() );}

    //! Adds volumes.
    /**
     *  Adds volume. Id computed automatically.
     *  @param v Volume to be added.
     *  @return Reference to the newly added volume.
     */
    element_Type & addVolume( element_Type const & vol );
    inline element_Type & addElement( volume_Type const & elem ) {return addVolume(elem);}
    inline element_Type & addElement( face_Type const & elem ) {return addFace(elem);}
    inline element_Type & addElement( edge_Type const & elem ) {return addEdge(elem);}


    //! Adds volume in a certain position.
    /**
     *  Adds volume to a specified position.
     *  @param v Volume to be added.
     *  @param pos Position of the volume.
     *  @return Reference to the newly added volume.
     */
    element_Type & setVolume( element_Type const & v, UInt const pos );
    inline element_Type & setElement( volume_Type const & elem, UInt const pos ) {return setVolume( elem, pos );}
    inline element_Type & setElement( face_Type const & elem, UInt const pos ) {return setFace( elem, pos );}
    inline element_Type & setElement( edge_Type const & elem, UInt const pos ) {return setEdge( elem, pos );}

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
    bool hasLocalFacets() const;
    inline bool hasLocalFaces() const {return hasLocalFaces( face_Type() );}
    inline bool hasLocalFaces(facet_Type) const {return hasLocalFacets();}
    inline bool hasLocalFaces(ridge_Type) const {return hasLocalRidges();}

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
    void updateElementFacets( bool createFaces = false, const bool verbose = false, UInt estimateFaceNumber = 0 );
    void updateElementFaces( bool createFaces = false, const bool verbose = false, UInt estimateFaceNumber = 0 )
		{updateElementFaces( face_Type(), createFaces , verbose, estimateFaceNumber );}
    void updateElementFaces( facet_Type, bool createFaces = false, const bool verbose = false, UInt estimateFaceNumber = 0 )
		{updateElementFacets( createFaces , verbose, estimateFaceNumber );}
    void updateElementFaces( ridge_Type, bool createFaces = false, const bool verbose = false, UInt estimateFaceNumber = 0 )
		{updateElementRidges( createFaces , verbose, estimateFaceNumber );}

    //! Destroys element-to-face container. Useful to save memory!
    void cleanElementFacets();
    inline void cleanElementFaces(){cleanElementFaces( face_Type() );}
    inline void cleanElementFaces(facet_Type){cleanElementFacets();}
    inline void cleanElementFaces(ridge_Type){cleanElementRidges();}

    //! Local Face Id.
    /** @param volId Id of volume (element).
     *  @param locF local face number 0 \< LocF \< numLocalFaces().
     *  @return ID of the face.
     */
    UInt localFacetId( UInt const volId, UInt const locF ) const;
    inline UInt localFaceId( UInt const volId, UInt const locF ) const {return localFacetId( volId, locF );}

    //! Local Face Id.
    /** @param iv Reference to a volume (element).
     *  @param locF local face number 0 \< LocF \< numLocalFaces().
     *  @return ID of the face.
     */
    inline UInt localFaceId( const element_Type & iv, UInt const locF ) const 
    	{return localFacetId( iv, locF ); }
    UInt localFacetId( const volume_Type & vol, UInt const locF ) const;
    UInt localFacetId( const face_Type & face, UInt const locF ) const;
    UInt localFacetId( const point_Type & point, UInt const locF ) const;

    

    //! Is the array for local Edges set up?
    /**
     *  It does not use switches, but interrogates the container directly.
     *
     *  @return true if edges information are available, false otherwise.
     */
    bool hasLocalRidges() const;
    inline bool hasLocalEdges() const {return hasLocalEdges( edge_Type() );}
    inline bool hasLocalEdges( ridge_Type ) const {return hasLocalRidges();}
    inline bool hasLocalEdges( facet_Type) const {return hasLocalFacets();}

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
     *  this condition is NOT a a paradigm for a RegionMesh.
     */
    void updateElementRidges( bool createEdges = false, const bool verbose = false,
                                 UInt estimateEdgeNumber = 0, bool renumber=true)
		{updateElementRidges( ridge_Type(), createEdges, verbose, estimateEdgeNumber, renumber);}

    void updateElementRidges( edge_Type, bool createEdges = false, const bool verbose = false,
                                 UInt estimateEdgeNumber = 0, bool renumber=true);

    void updateElementRidges( point_Type, bool createEdges = false, const bool verbose = false,
            UInt estimateEdgeNumber = 0, bool renumber=true){;}

    inline void updateElementEdges( bool createEdges = false, const bool verbose = false,
                             UInt estimateEdgeNumber = 0, bool renumber=true){
           updateElementEdges( edge_Type(), createEdges, verbose, estimateEdgeNumber, renumber); }
    inline void updateElementEdges( ridge_Type, bool createEdges = false, const bool verbose = false,
                                 UInt estimateEdgeNumber = 0, bool renumber=true){
           updateElementRidges( createEdges, verbose, estimateEdgeNumber, renumber); }
    inline void updateElementEdges( facet_Type, bool createEdges = false, const bool verbose = false,
                                     UInt estimateEdgeNumber = 0, bool renumber=true){
           updateElementFacets( createEdges, verbose, estimateEdgeNumber, renumber); }

    //! Destroys Edge-To-Face lookup table.
    void cleanElementRidges();
    inline void cleanElementEdges() {cleanElementEdges(edge_Type());}
    inline void cleanElementEdges(ridge_Type) {cleanElementRidges();}
    inline void cleanElementEdges(facet_Type) {cleanElementFacets();}

    //! Local Edge.
    /** @param volId Id of volume (element).
     *  @param locE local edge number 0 \< LocE \< numLocalEdges().
     *  @return ID of the edge.
     */
    inline UInt localRidgeId( UInt const elemId, UInt const locE ) const {return localRidgeId( ridge_Type(), elemId, locE );}
    inline UInt localRidgeId( point_Type, UInt const elemId, UInt const locR ) const {return element(elemId).point(locR).localId();}
    UInt localRidgeId( edge_Type,  UInt const elemId, UInt const locR ) const;
    inline UInt localEdgeId( UInt const elemId, UInt const locE ) const {return localRidgeId( elemId, locE );}
    inline UInt localEdgeId( ridge_Type, UInt const elemId, UInt const locE ) const {return localRidgeId( ridge_Type(), elemId, locE );}
    inline UInt localEdgeId( facet_Type, UInt const elemId, UInt const locE ) const {return localFacetId( elemId, locE );}

    //! Local Edge.
    /** @param iv Reference of the volume.
     *  @param locE local edge number 0 \< LocE \< numLocalEdges().
     *  @return ID of the edge.
     */
    inline UInt localRidgeId( const element_Type & elem, UInt const locE ) const {return localRidgeId( elem.id(), locE );}
    UInt localEdgeId( const volume_Type & elem, UInt const locE ) const { return localRidgeId( elem.id(), locE );}
    UInt localEdgeId( const face_Type & elem, UInt const locE ) const { return localFacetId( elem.id(), locE );}



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
    inline UInt numFacets() const {return numFacets( facet_Type() );}
    inline UInt numFacets(face_Type) const {return numFaces();}
    inline UInt numFacets(edge_Type) const {return numEdges();}
    inline UInt numFacets(point_Type) const {return numPoints();}

    //! Returns Global Number of Faces
    /**
     *  Returns global number of Face elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global Number of Faces.
     */
    UInt numGlobalFaces() const;
    inline UInt numGlobalFacets() const {return numGlobalFacets( facet_Type() );}
    inline UInt numGlobalFacets(face_Type) const {return numGlobalFaces();}
    inline UInt numGlobalFacets(edge_Type) const {return numGlobalEdges();}

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
    void setMaxNumFacets( UInt const n, bool const setcounter = false ) {setMaxNumFacets( facet_Type(), n, setcounter);}
    void setMaxNumFacets( face_Type, UInt const n, bool const setcounter = false ) {setMaxNumFaces( n, setcounter);}
    void setMaxNumFacets( edge_Type, UInt const n, bool const setcounter = false ) {setMaxNumEdges( n, setcounter);}
    void setMaxNumFacets( point_Type, UInt const n, bool const setcounter = false ) {setMaxNumPoints( n, setcounter);}

    //! Changes Current capacity of Global Faces.
    /**
     *  Changes Current capacity of Global Faces (Optionally sets internal counter).
     *
     *  @param n maximum number of global faces.
     */
    void setMaxNumGlobalFaces( UInt const n );
    void setMaxNumGlobalFacets( UInt const n ) {setMaxNumGlobalFacets( facet_Type(), n );}
    void setMaxNumGlobalFacets( face_Type, UInt const n ) {setMaxNumGlobalFaces( n );}
    void setMaxNumGlobalFacets( edge_Type, UInt const n ) {setMaxNumGlobalEdges( n );}
    void setMaxNumGlobalFacets( point_Type, UInt const n ) {setMaxNumGlobalPoints( n );}

    //! Adds a face.
    /**
     *  Adds a face (optionally a boundary face). Id computed automatically.
     *  @param boundary true if it's a boundary face.
     *  @return Reference to added face.
     */

    face_Type & addFace() {return addFace( face_Type() ); }
    face_Type & addFace( bool const boundary );
    facet_Type & addFacet( bool const boundary ) {return addFacet( facet_Type(), boundary);}
    facet_Type & addFacet( face_Type, bool const boundary ) {return addFace(boundary);}
    facet_Type & addFacet( edge_Type, bool const boundary ) {return addEdge(boundary);}
    facet_Type & addFacet( point_Type, bool const boundary ) {return addPoint(boundary);}

    //! Adds a face.
    /**
     *  Adds a face (optionally a boundary face). Id computed automatically.
     *  It assumes that all attributes of face f have been properly set
     *  @param f Face to be added.
     *  @return Reference to the newly added face.
     */
    face_Type & addFace( face_Type const & f );
    facet_Type & addFacet( face_Type const & facet )  {return addFace(facet);}
    facet_Type & addFacet( edge_Type const & facet )  {return addEdge(facet);}
    facet_Type & addFacet( point_Type const & facet )  {return addPoint(facet);}

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
    face_Type & setFace( face_Type const & f, UInt pos);

    //! Reference to last face stored in list.
    /**
     *  Reference to last face stored in list.
     *  Useful for mesh readers.
     *  @return reference of the last face in the list.
     */
    face_Type & lastFace();

    //! i-th mesh Face.
    /**
     *  @param i index of the mesh face.
     *  @return the i-th face.
     */
   face_Type const & face( UInt const i ) const;
   facet_Type const & facet( UInt const i ) const {return facet( facet_Type(), i);};
   facet_Type const & facet( face_Type, UInt const i ) const {return face(i);};
   facet_Type const & facet( edge_Type, UInt const i ) const {return edge(i);};
   facet_Type const & facet( point_Type, UInt const i ) const {return point(i);};


    //! i-th mesh face.
    /**
     *  @param i index of the mesh face.
     *  @return reference to the ith mesh face.
     */
    face_Type & face( UInt const i );
    facet_Type & facet( UInt const i ) {return facet( facet_Type(), i);};
    facet_Type & facet( face_Type, UInt const i ) {return face(i);};
    facet_Type & facet( edge_Type, UInt const i ) {return edge(i);};
    facet_Type & facet( point_Type, UInt const i ) {return point(i);};


    //! i-th mesh boundary face.
    /**
     *  @param i index of the mesh boundary face.
     *  @return the i-th face.
     */
    face_Type const & boundaryFace( UInt const i ) const;

    //! i-th mesh boundary face.
    /**
     *  @param i index of the mesh boundary face.
     *  @return the i-th face.
     */
    face_Type & boundaryFace( UInt const i );

    //! Set counter of faces.
    /**
     *  @param n Number of faces.
     */
    void setNumFaces( UInt const n ) ;
    void setNumFacets( UInt const n ) {setNumFacets( facet_Type(), n );}
    void setNumFacets( face_Type, UInt const n ) {setNumFaces( n );}
    void setNumFacets( edge_Type, UInt const n ) {setNumEdges( n );}

    //! Set counter of boundary faces.
    /**
     *  @param n Number of boundary faces.
     */
    void setNumBFaces( UInt const n ) ;
    void setNumBFacets( UInt const n ) {setNumBFacets( facet_Type(), n ) ;}
    void setNumBFacets( face_Type, UInt const n ) {setNumBFaces( n ) ;}
    void setNumBFacets( edge_Type, UInt const n ) {setNumBEdges( n ) ;}
    void setNumBFacets( point_Type, UInt const n ) {setNumBPoints( n ) ;}


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
    bool isBoundaryFace( UInt const & id  ) const;
    bool isBoundaryFacet( UInt const & id ) const { return isBoundaryFacet( facet_Type(), id );};
    bool isBoundaryFacet( face_Type, UInt const & id ) const { return isBoundaryFace( id );};
    bool isBoundaryFacet( edge_Type, UInt const & id ) const { return isBoundaryEdge( id );};
    bool isBoundaryFacet( point_Type, UInt const & id ) const { return isBoundaryPoint( id );};


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
    bool isFullFace( UInt const & id ) const ;

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
    UInt faceElement( UInt const faceId, UInt const Pos ) const ;

    //! 3DElement adjacent to a FACE. Face reference given.
    /**
     *  @param f Face reference.
     *  @param Pos is equal to 0 or 1 and indicates first or second element.
     *  @return Id of adjacent volume element or 0 if none.
     */
    UInt faceElement( facet_Type const & f, UInt const Pos ) const;

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
    UInt numRidges() const {return numRidges(ridge_Type()); }
    UInt numRidges(edge_Type) const {return numEdges(); }
    UInt numRidges(point_Type) const {return numPoints(); }

    //! Global number of Edges.
    /**
     *  Returns global number of Edge elements in the mesh
     *  as given by the internal counter.
     *
     *  @return Global number of Edges.
     */
    UInt numGlobalEdges() const;
    UInt numGlobalRidges() const {return numGlobalRidges(element_Type() );}
    UInt numGlobalRidges(volume_Type) const {return numGlobalEdges();}
    UInt numGlobalRidges(face_Type) const {return numGlobalVertices();}
    UInt numGlobalRidges(edge_Type) const {return 0;}

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
    void setMaxNumEdges( UInt const n, bool const setcounter = false );
    void setMaxNumRidges( UInt const n, bool const setcounter = false ) {setMaxNumRidges(ridge_Type(), n, setcounter);}
    void setMaxNumRidges( edge_Type, UInt const n, bool const setcounter ) {setMaxNumEdges( n, setcounter);}
    void setMaxNumRidges( point_Type, UInt const n, bool const setcounter ) {setMaxNumPoints( n, setcounter);}

    //! Changes Current capacity of Global Edges.
    /**
     *  Optionally sets internal counter.
     *
     *  @param n Maximum number of global edges.
     */
    void setMaxNumGlobalEdges( UInt const n );
    void setMaxNumGlobalRidges( UInt const n ){ setMaxNumGlobalRidges(ridge_Type(), n);}
    void setMaxNumGlobalRidges( edge_Type, UInt const n ){ setMaxNumGlobalEdges( n );}
    void setMaxNumGlobalRidges( point_Type, UInt const n ){ setMaxNumGlobalPoints( n );}

    //! Adds an Edge.
    /**
     *  Adds an edge. Id computed automatically.
     *
     *  @param boundary true if is on boundary.
     *  @return Reference to added edge.
     */
    edge_Type & addEdge( bool const boundary );
    ridge_Type & addRidge( bool const boundary ) {return addRidge( ridge_Type(), boundary);}
    ridge_Type & addRidge( edge_Type, bool const boundary ) {return addEdge( boundary);}
    ridge_Type & addRidge( point_Type, bool const boundary ) {return addPoint( boundary);}

    //! Adds an Edge.
    /**
     *  Adds an edge to the end of the list
     *  and adjourn its Id. The edge attributes (a part the id) should be
     *  correctly set.
     *
     *  @param f Edge to add.
     *  @return Reference to added edge.
     */
    edge_Type & addEdge( edge_Type const & r );
    ridge_Type & addRidge( edge_Type const & e ){return addEdge( e );}
    ridge_Type & addRidge( point_Type const & p ){return addPoint( p );}

    //! Add an Edge to specified position.
    /**
     *  Adds an edge and sets its Id to position.
     *
     *  @param f Edge to add.
     *  @param position Position of the edge.
     *  @return Reference to added edge.
     */
    edge_Type & setEdge( edge_Type const & f, UInt position );

    //! Reference to last edge stored in list.
    /**
     *  Useful for mesh readers
     *  @return reference of the last edge in the list.
     */
    edge_Type & lastEdge();

    //! i-th mesh Edge.
    /**
     *  Returns the i-th Edge.
     *
     *  @param i Index of the mesh Edge.
     *  @return The i-th Edge.
     */
    edge_Type const & edge( UInt const i ) const;
    ridge_Type const & ridge( UInt const i ) const {return ridge(edge_Type(), i);}
    ridge_Type const & ridge( edge_Type, UInt const i ) const {return edge(i);}
    ridge_Type const & ridge( point_Type, UInt const i ) const {return point(i);}

    //! i-th mesh 1D Edge reference.
    /**
     *  Returns a reference to the i-th mesh Edge.
     *
     *  @param i Index of the mesh 1D Edge.
     *  @return Reference to the i-th Edge.
     */
    edge_Type & edge( UInt const i );
    ridge_Type & ridge( UInt const i ) {return ridge( ridge_Type(), i);}
    ridge_Type & ridge( edge_Type, UInt const i ) {return edge(i);}
    ridge_Type & ridge( point_Type, UInt const i ) {return point(i);}

    //! i-th mesh 1D Boundary Edge.
    /**
     *  Returns the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return i-th Boundary Edge.
     */
    edge_Type const & boundaryEdge( UInt const i ) const;

    //! i-th mesh 1D Boundary Edge Reference.
    /**
     *  Returns a reference to the i-th mesh Boundary Edge.
     *
     *  @param i Index of the mesh 1D Boundary Edge.
     *  @return Reference to the i-th Boundary edge.
     */
    edge_Type & boundaryEdge( UInt const i );

    //! Set boundary Edge counter.
    /**
     *  Set the Boundary Edges counter to a given number.
     *
     *  @param n Count of Boundary Edge.
     */
    void setNumBEdges( UInt const n );

    void setNumBRidges( UInt const n ) {setNumBRidges( ridge_Type(), n );}
    void setNumBRidges( edge_Type, UInt const n ) {setNumBEdges( n );}
    void setNumBRidges( point_Type, UInt const n ) {setNumBPoints( n );}

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
    bool isBoundaryEdge( ridge_Type const & e ) const {return (e.id() < M_numBEdges);}
    bool isBoundaryRidge( edge_Type const & e ) const {return isBoundaryEdge(e); }
    bool isBoundaryRidge( point_Type const & p ) const {return isBoundaryPoint(p); }


    //! Edge on boundary check by id.
    /**
     *  Is this edge, of given id, on boundary?
     *
     *  @param id Id of the edge.
     *  @return true if the edge is on the boundary, false otherwise.
     */
    bool isBoundaryEdge( UInt const & id ) const {return (id < M_numBEdges );}
    bool isBoundaryRidge( UInt const & id ) const {return isBoundaryRidge( ridge_Type(), id );}
    bool isBoundaryRidge( edge_Type, UInt const & id ) const {return isBoundaryEdge( id );}
    bool isBoundaryRidge( point_Type, UInt const & id ) const {return isBoundaryPoint( id );}



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
    void setNumEdges ( UInt const n );
    void setNumRidges ( UInt const n ) {setNumEdges ( n );}

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
    point_Type & addPoint( bool const boundary, bool const vertices);

    //! Adds a Point in the mesh.
    /**
     *  Adds a Point inside the mesh.
     *  Point will be add at the end of the list and the Id is computed automatically
     *
     *  @param p Point to be added.
     *  @return Reference to the newly added Point.
     */
    point_Type & addPoint( point_Type const & p);


    //! Adds a Point in the mesh giving an id.
    /**
     *  Stores a point at the specified position. The id is set equal to the position
     *  All other Point attributes must be set beforehand
     *
     *  @param p Point to be added.
     *  @param position Desired id.
     *  @return Reference to the newly added Point.
     */
    point_Type & setPoint( point_Type const & p, UInt const position);

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
    inline peak_Type const & peak( UInt const i ) const {return peak( peak_Type(), i );}
    inline peak_Type const & peak( point_Type, UInt const i ) const {return point(i);}

    //! Returns a reference to the i-th mesh Point.
    /**
     *  Returns the i-th Point in the mesh.
     *UInt const i
     *  @param i Id of the Point.
     *  @return Reference i-th mesh Point.
     */
    point_Type & point( UInt const i );
    inline peak_Type & peak( UInt const i ) {return peak( peak_Type(), i );}
    inline peak_Type & peak( point_Type, UInt const i ) {return point(i);}

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

    inline UInt numGlobalPeaks() const { return numGlobalPeaks(element_Type() );}
    inline UInt numGlobalPeaks(volume_Type) const { return numGlobalVertices();}
    inline UInt numGlobalPeaks(face_Type) const { return 0;}
    inline UInt numGlobalPeaks(edge_Type) const { return 0;}


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
    void setNumGlobalVertices( UInt const n );

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
    //! Container of mesh 3D Elements
    volumes_Type volumeList;
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

   // elements_Type& elementList() {return volumeList;}

    inline facets_Type& facetList() {return facetList(facet_Type());}
    inline facets_Type& facetList(face_Type) {return faceList;}
    inline facets_Type& facetList(edge_Type) {return edgeList;}


    inline ridges_Type& ridgeList() {return ridgeList(ridge_Type() );}
    inline ridges_Type& ridgeList(edge_Type) {return edgeList;}
    inline ridges_Type& ridgeList(point_Type) {return pointList;}


    inline elements_Type& elementList(){return elementList(element_Type() );};
    inline elements_Type& elementList(volume_Type) {return volumeList;};
    inline elements_Type& elementList(face_Type) {return faceList;};



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

    MeshUtility::MeshTransformer<RegionMesh<GEOSHAPE, MC> > M_meshTransformer;

}; // End of class RegionMesh


// =================================================== //
// =================================================== //
//                    IMPLEMENTATION                   //
// =================================================== //
// =================================================== //


void set_switches_for_regionmesh( Switch & sw );

template <typename GEOSHAPE, typename MC>
inline RegionMesh<GEOSHAPE, MC>::RegionMesh() :
MeshEntity(),
MC::regionMarker_Type(),
switches(),
M_numVolumes( 0 ),
M_numVertices( 0 ),
M_numBVertices( 0 ),
M_numPoints( 0 ),
M_numBPoints( 0 ),
M_numFaces( 0 ),
M_numBFaces( 0 ),
M_numEdges( 0 ),
M_numBEdges( 0 ),
M_meshTransformer(*this)
{ //Modif Miguel:11/2002
    set_switches_for_regionmesh( switches );
}


template <typename GEOSHAPE, typename MC>
inline RegionMesh<GEOSHAPE, MC>::RegionMesh( UInt id ) :
MeshEntity( id ),
MC::RegionMarker(),
switches(),
M_numVolumes( 0 ),
M_numVertices( 0 ),
M_numBVertices( 0 ),
M_numPoints( 0 ),
M_numBPoints( 0 ),
M_numFaces( 0 ),
M_numBFaces( 0 ),
M_numEdges( 0 ),
M_numBEdges( 0 ),
M_meshTransformer(*this)
{ //Modif Miguel:11/2002
    set_switches_for_regionmesh( switches );
}

template <typename GEOSHAPE, typename MC>
inline RegionMesh<GEOSHAPE, MC>::~RegionMesh()
{
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numLocalVertices() const
{
    return element_Type::S_numLocalVertices;
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numLocalFaces() const
{
    return element_Type::S_numLocalFaces;
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numLocalEdges() const
{
    return element_Type::S_numLocalEdges;
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numLocalEdgesOfFace() const
{
    return facet_Type::S_numLocalEdges;
}

// ***************************** VOLUMES

template <typename GEOSHAPE, typename MC>
inline UInt /*const*/
RegionMesh<GEOSHAPE, MC>::numVolumes() const
{
    return M_numVolumes;
}

template <typename GEOSHAPE, typename MC>
inline UInt /*const*/
RegionMesh<GEOSHAPE, MC>::numGlobalVolumes() const
{
    return M_numGlobalVolumes;
}

//     template <typename GEOSHAPE, typename MC>
//     UInt &
//     RegionMesh<GEOSHAPE, MC>::numVolumes()
//     {
//         return M_numVolumes;
//     }

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::maxNumVolumes() const
{
    return maxNumItems( volumeList );
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::storedVolumes() const
{
    return volumeList.numItems();
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumVolumes( UInt const n, bool const setcounter )
{
    volumeList.setMaxNumItems(n);
    if ( setcounter )
        M_numVolumes = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumVolumes( UInt const n )
{
    M_numVolumes = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumGlobalVolumes( UInt const n )
{
    M_numGlobalVolumes = n;
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::element_Type &
RegionMesh<GEOSHAPE, MC>::addVolume( element_Type const & v )
{
    ASSERT_PRE( volumeList.size() < volumeList.capacity() , "Volume list size exceeded" <<
                volumeList.size() + 1 << " " << volumeList.capacity() ) ;
    volumeList.push_back( v );

    ( volumeList.back() ).setId( volumeList.size() - 1 );
    return volumeList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::element_Type &
RegionMesh<GEOSHAPE, MC>::setVolume( element_Type const & v, UInt const pos )
{
    ASSERT_PRE( pos < volumeList.capacity() , "position requested exceed capacity" <<
                pos << " " << volumeList.capacity() ) ;
    volumeList( pos ) = v;
    volumeList( pos ).setId( pos );
    return volumeList( pos );
}

template <typename GEOSHAPE, typename MC>
inline
void
RegionMesh<GEOSHAPE, MC>::setVolumeCounter()
{
    M_numVolumes = volumeList.size();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::element_Type &
RegionMesh<GEOSHAPE, MC>::lastVolume()
{
    return volumeList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::element_Type const &
RegionMesh<GEOSHAPE, MC>::volume( UInt const i ) const
{
    ASSERT_BD( i < volumeList.size() ) ;
    return volumeList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::element_Type &
RegionMesh<GEOSHAPE, MC>::volume( UInt const i )
{
    ASSERT_BD( i < volumeList.size() ) ;
    return volumeList[ i ];
}

// ************************* FACES ******************************
template <typename GEOSHAPE, typename MC>
inline UInt /*const*/
RegionMesh<GEOSHAPE, MC>::numFaces() const
{
    return M_numFaces;
}

template <typename GEOSHAPE, typename MC>
inline UInt /*const*/
RegionMesh<GEOSHAPE, MC>::numGlobalFaces() const
{
    return M_numGlobalFaces;
}

//     template <typename GEOSHAPE, typename MC>
//     UInt &
//     RegionMesh<GEOSHAPE, MC>::numFaces()
//     {
//         return M_numFaces;
//     }

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::storedFaces() const
{
    return faceList.numItems();
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::maxNumFaces() const
{
    return faceList.maxNumItems();
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumFaces( UInt const n, bool const setcounter )
{
    faceList.setMaxNumItems(n);
    if ( setcounter )
        M_numFaces = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumGlobalFaces( UInt const n )
{
    M_numGlobalFaces = n;
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type &
RegionMesh<GEOSHAPE, MC>::addFace( bool const boundary )
{
    face_Type aFace;
    aFace.setBoundary(boundary);
    return this->addFace( aFace);
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type &
RegionMesh<GEOSHAPE, MC>::addFace( face_Type const & f )
{
    ASSERT_PRE( faceList.size() < faceList.capacity(), "Face list size exceeded" <<
                faceList.size() + 1 << " " << faceList.capacity() ) ;
    faceList.push_back( f );
    faceList.back().setId( faceList.size() -1 );

    return faceList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type &
RegionMesh<GEOSHAPE, MC>::setFace( face_Type const & f, UInt position)
{
    ASSERT_PRE( position < faceList.capacity(), "Face list size exceeded" <<
                position << " " << faceList.capacity() ) ;
    faceList( position ) = f;
    faceList( position ).setId( position );
    return faceList( position );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type &
RegionMesh<GEOSHAPE, MC>::lastFace()
{
    return faceList.back();
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type const &
RegionMesh<GEOSHAPE, MC>::face( UInt const i ) const
{
    ASSERT_BD( i < faceList.size() ) ;
    return faceList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type &
RegionMesh<GEOSHAPE, MC>::face( UInt const i )
{
    ASSERT_BD( i < faceList.size() ) ;
    return faceList[ i ];
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type const &
RegionMesh<GEOSHAPE, MC>::boundaryFace( UInt const i ) const
{
    ASSERT_PRE( faceList.size() != 0, "Boundary Faces not stored" ) ;
    //ASSERT_BD( i < this->numBFaces()) ; //TODO TO BE FIXED! LF
    ASSERT_BD( i < faceList.size() ) ;
    return faceList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::face_Type &
RegionMesh<GEOSHAPE, MC>::boundaryFace( UInt const i )
{
    ASSERT_PRE( faceList.size() != 0, "Boundary Faces not stored" ) ;
    //    ASSERT_BD( i < this->numBFaces()) ;
    ASSERT_BD( i < faceList.size() ) ;
    return faceList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numBFaces() const
{
    return M_numBFaces;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumFaces( UInt const n )
{
    M_numFaces = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumBFaces( UInt const n )
{
    M_numBFaces = n;
}

// ************************* EDGES ******************************

template <typename GEOSHAPE, typename MC>
inline UInt /*const*/
RegionMesh<GEOSHAPE, MC>::numEdges() const
{
    return M_numEdges;
}

template <typename GEOSHAPE, typename MC>
inline UInt /*const*/
RegionMesh<GEOSHAPE, MC>::numGlobalEdges() const
{
    return M_numGlobalEdges;
}

//     template <typename GEOSHAPE, typename MC>
//     UInt &
//     RegionMesh<GEOSHAPE, MC>::numEdges()
//     {
//         return M_numEdges;
//     }

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::storedEdges() const
{
    return edgeList.numItems();
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::maxNumEdges() const
{
    return edgeList.maxNumItems();
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumEdges( UInt const n, bool const setcounter )
{
    edgeList.setMaxNumItems(n);
    if ( setcounter )
        M_numEdges = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumGlobalEdges( UInt const n )
{
    M_numGlobalEdges = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumEdges( UInt const n)
{
    M_numEdges = n;
}

template <typename GEOSHAPE, typename MC>
inline typename RegionMesh<GEOSHAPE, MC>::edge_Type &
RegionMesh<GEOSHAPE, MC>::addEdge( bool const boundary )
{
    edge_Type anEdge;
    anEdge.setBoundary(boundary);
    return addEdge( anEdge);
}

template <typename GEOSHAPE, typename MC>
inline typename RegionMesh<GEOSHAPE, MC>::edge_Type &
RegionMesh<GEOSHAPE, MC>::addEdge( edge_Type const & f)
{
    ASSERT_PRE( edgeList.size() < edgeList.capacity(), "Edge list size exceeded" <<
                edgeList.size() + 1 << " " << edgeList.capacity() ) ;

    edgeList.push_back( f );
    ( edgeList.back() ).setId(edgeList.size() - 1);

    return edgeList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::edge_Type &
RegionMesh<GEOSHAPE, MC>::setEdge( edge_Type const & f, UInt position)
{
    ASSERT_PRE( position < edgeList.capacity(), "Edge list size exceeded" <<
                position << " " << edgeList.capacity() ) ;
    edgeList( position ) = f;
    edgeList( position ).setId( position );
    return edgeList( position );
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::edge_Type &
RegionMesh<GEOSHAPE, MC>::lastEdge()
{
    return edgeList.back();
}



template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::edge_Type const &
RegionMesh<GEOSHAPE, MC>::edge( UInt const i ) const
{
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::edge_Type &
RegionMesh<GEOSHAPE, MC>::edge( UInt const i )
{
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList[ i ];
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::edge_Type const &
RegionMesh<GEOSHAPE, MC>::boundaryEdge( UInt const i ) const
{
    ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::edge_Type &
RegionMesh<GEOSHAPE, MC>::boundaryEdge( UInt const i )
{
    ASSERT_PRE( edgeList.size() != 0, "Boundary Edges not stored" ) ;
    ASSERT_BD( i < edgeList.size() ) ;
    return edgeList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numBEdges() const
{
    return M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumBEdges( UInt const n )
{
    M_numBEdges = n;
}

// ************************ Points/Vertices

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numPoints() const
{
    return M_numPoints;
}

template <typename GEOSHAPE, typename MC>
inline UInt &
RegionMesh<GEOSHAPE, MC>::numPoints()
{
    return M_numPoints;
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::storedPoints() const
{
    return pointList.numItems();
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::storedBPoints() const
{
    return _bPoints.size();
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::maxNumPoints() const
{
    return pointList.maxNumItems();
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumPoints( UInt const n, bool const setcounter )
{
    pointList.setMaxNumItems(n);
    if ( setcounter )
        M_numPoints = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumVertices( UInt const n )
{
    M_numVertices = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumGlobalVertices( UInt const n )
{
    M_numGlobalVertices = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumBVertices( UInt const n )
{
    M_numBVertices = n;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setMaxNumGlobalPoints( UInt const n )
{
    M_numGlobalPoints = n;
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type &
RegionMesh<GEOSHAPE, MC>::addPoint( bool const boundary, bool const vertex )
{
    point_Type aPoint;
    aPoint.setBoundary(boundary);
    if(vertex) aPoint.replaceFlag(aPoint.flag() | EntityFlags::VERTEX);
    return addPoint(aPoint);
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type &
RegionMesh<GEOSHAPE, MC>::addPoint( point_Type const & p)
{
    ASSERT_PRE( pointList.size() < pointList.capacity(), "Point list size exceeded" <<
                pointList.size() + 1 << " " << pointList.capacity() ) ;

    pointList.push_back( p );

    point_Type * pp = & pointList.back();
    pp->setId(pointList.size()-1);

    if ( pp->boundary() )
    {
        ASSERT_PRE( _bPoints.size() < _bPoints.capacity(), "Boundary Point list size exceeded" <<
                    _bPoints.size() + 1 << " " << _bPoints.capacity() ) ;
        _bPoints.push_back( pp );
    }
    return pointList.back();
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type &
RegionMesh<GEOSHAPE, MC>::setPoint( point_Type const & p, UInt position)
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
inline typename RegionMesh<GEOSHAPE, MC>::point_Type &
RegionMesh<GEOSHAPE, MC>::
changePointBoundaryFlag(UInt const & position, bool const boundary)
{
    point_Type pp = pointList[position];
    if(pp.boundary()!=boundary)
    {
        pp.setBoundary(boundary);
        this->setPoint(pp,position);
    }
    return pointList[position];
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type &
RegionMesh<GEOSHAPE, MC>::lastPoint()
{
    return pointList.back();
}


template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type const &
RegionMesh<GEOSHAPE, MC>::point( UInt const i ) const
{
    ASSERT_BD( i < pointList.size() ) ;

    return pointList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type &
RegionMesh<GEOSHAPE, MC>::point( UInt const i )
{
    ASSERT_BD( i < pointList.size() ) ;

    return pointList[ i ];
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type const &
RegionMesh<GEOSHAPE, MC>::boundaryPoint( UInt const i ) const
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i < _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
inline
typename RegionMesh<GEOSHAPE, MC>::point_Type &
RegionMesh<GEOSHAPE, MC>::boundaryPoint( UInt const i )
{
    ASSERT_PRE( _bPoints.size() != 0, " Boundary Points not Stored" ) ;
    ASSERT_BD( i < _bPoints.size() ) ;
    return *( _bPoints( i ) );
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numBPoints() const
{
    return M_numBPoints;
}

template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::setNumBPoints( UInt const n )
{
    M_numBPoints = n;
    _bPoints.reserve( n );
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numVertices() const
{
    return M_numVertices;
}

template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numGlobalVertices() const
{
    return M_numGlobalVertices;
}


template <typename GEOSHAPE, typename MC>
inline UInt
RegionMesh<GEOSHAPE, MC>::numBVertices() const
{
    return M_numBVertices;
}

template <typename GEOSHAPE, typename MC>
inline UInt &
RegionMesh<GEOSHAPE, MC>::numBVertices()
{
    return M_numBVertices;
}

// *************** GENERAL *******************

template <typename GEOSHAPE, typename MC>
inline std::ostream &
RegionMesh<GEOSHAPE, MC>::showMe( bool verbose, std::ostream & out ) const
{
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << "                      RegionMesh                " << std::endl;
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    out << " ID: " << this->id() << " Marker Flag: " << this->marker() << std::endl;
    //  out <<"Faces local to  volumes stored: "<<hasLocalFacets()<<std::endl;
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
    switches.showMe( verbose, out );
    out << "**************************************************" << std::endl;
    out << "**************************************************" << std::endl;
    if ( verbose )
    {
        std::cout << "Verbose version not implemented yet";
    }
    return out;

}




template <typename GEOSHAPE, typename MC>
inline int
RegionMesh<GEOSHAPE, MC>::check( int level, bool const fix, bool const verb, std::ostream & out )
{
    int severity = 0;
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
        checkMesh3D( *this, testsw, fix, verb, out, out, out );
        if ( verb )
        {
            out << "**********************************************" << std::endl <<
                            "DETAILS OF EXTENDED  CHECK:" << std::endl;
            testsw.showMe( true, out );
            out << "**********************************************" << std::endl << std::endl;
            if ( testsw.test( "ABORT_CONDITION" ) )
                return 1;

        }
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
            out << "WARNING: No Edges Stored" << std::endl;
        if ( verb )
            out << "MAY NOT WORK IF WE HAVE DOF ON EDGES AND ESSENTIAL BC!!!!" << std::endl;
        severity = -1;
        unsetLinkSwitch( "HAS_ALL_RIDGES" );
        unsetLinkSwitch( "HAS_BOUNDARY_RIDGES" );
    }
    else if ( edgeList.size() == numBEdges() )
    {

        setLinkSwitch( "HAS_BOUNDARY_RIDGES" );
        unsetLinkSwitch( "HAS_ALL_RIDGES" );
        if ( verb )
            out << "INFORMATION: Only Boundary Edges Stored" << std::endl;
    }
    else if ( edgeList.size() == numEdges() )
    {
        setLinkSwitch( "HAS_BOUNDARY_RIDGES" );
        setLinkSwitch( "HAS_ALL_RIDGES" );
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
            out << "ERROR: No Faces Stored: at least boundary faces are needed" << std::endl;
        if ( verb )
            out << "" << std::endl;
        severity = 1;
        unsetLinkSwitch( "HAS_ALL_FACETS" );
        unsetLinkSwitch( "HAS_BOUNDARY_FACETS" );
    }
    else if ( faceList.size() == numBFaces() )
    {
        setLinkSwitch( "HAS_BOUNDARY_FACETS" );
        unsetLinkSwitch( "HAS_ALL_FACETS" );
        if ( verb )
            out << "INFORMATION: Only Boundary Faces Stored" << std::endl;
    }
    else if ( faceList.size() == numFaces() )
    {
        setLinkSwitch( "HAS_BOUNDARY_FACETS" );
        setLinkSwitch( "HAS_ALL_FACETS" );
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
            ++count;
    if ( count == 0 )
        severity = 4;
    if ( count != M_numBPoints )
    {
        out << "Num Boundary points " << count << " not equal to internal counter value "
                        << M_numBPoints << std::endl;
        if ( ( count != 0 ) & fix )
        {
            M_numBPoints = count;
            out << "Fixed Counter";
            out.flush();
        }
    }

    if ( M_numVertices == 0 )
    {
        severity = 6;

        out << " SEVERITY ERROR: internal Vertices Counter unset";
    }

    if ( M_numPoints == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: internal Points Counter unset";
    }
    if ( M_numPoints == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: internal Points Counter unset";
    }
    if ( M_numBPoints == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: boundary Points Counter unset";
    }
    if ( M_numBVertices == 0 )
    {
        severity = 6;
        out << " SEVERITY ERROR: boundary Vertices Counter unset";
    }

    if ( verb )
        out << "   Check Finished              " << std::endl <<
        "***********************************************" << std::endl;
    setLinkSwitch( "HAS_BEEN_CHECKED" );

    return severity;

}


template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh<GEOSHAPE, MC>::faceElement( UInt const i, UInt const Pos ) const
{
    ASSERT_PRE( i < faceList.size(), "Not enough faces stored" ) ;
    ASSERT_PRE( Pos <= 1 , "Wrong position (0 or 1)" ) ;
    if ( Pos == 0 )
    {
        return ( faceList[ i ] ).firstAdjacentElementIdentity();
    }
    else
    {
        return ( faceList[ i ] ).secondAdjacentElementIdentity();
    }
}

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh<GEOSHAPE, MC>::faceElement( facet_Type const & f, UInt const Pos ) const
{
    ASSERT_BD( ! faceList.empty() ) ;
    ASSERT_PRE( Pos <= 1 , "Wrong position (0 or 1)" );
    if ( Pos == 0 )
    {
        return f.firstAdjacentElementIdentity();
    }
    else
    {
        return f.secondAdjacentElementIdentity();
    }
}


template <typename GEOSHAPE, typename MC>
inline
void
RegionMesh<GEOSHAPE, MC>::setLinkSwitch( std::string const & _s )
{
    ASSERT0( switches.set( _s ), std::stringstream( "Switch named " + _s + " is not allowed" ).str().c_str() );
}

template <typename GEOSHAPE, typename MC>
inline
void
RegionMesh<GEOSHAPE, MC>::unsetLinkSwitch( std::string const & _s )
{
    ASSERT0( switches.unset( _s ), std::stringstream( "Switch named " + _s + " is not allowed" ).str().c_str() );
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::getLinkSwitch( std::string const & _s ) const
{
    return switches.test( _s );
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::hasInternalFaces() const
{
    return faceList.size() > M_numBFaces;
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::hasInternalEdges() const
{
    return edgeList.size() > M_numBEdges;
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::hasEdges() const
{
    return ! edgeList.empty();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::hasFaces() const
{
    return ! faceList.empty();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::isVertex( point_Type const & p ) const
{
    return Flag::testOneSet(p.flag(),EntityFlags::VERTEX);
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::isVertex( UInt const & id ) const
{
    return isVertex(pointList[id]);
}


template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::isBoundaryPoint( UInt const & id ) const
{
    return point( id ).boundary();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::isBoundaryPoint( point_Type const & p ) const
{
    return p.boundary();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::isBoundaryFace( UInt const & id ) const
{
    return this->faceList[id].boundary();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::hasLocalFacets() const
{
    return ! M_ElemToFacet.empty();
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::hasLocalRidges() const
{
    return ! M_ElemToRidge.empty();
}

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh<GEOSHAPE, MC>::localFacetId( UInt const elemId, UInt const locF ) const
{
    ASSERT_PRE( !M_ElemToFacet.empty(), "Volume to Face array not  set" );
    ASSERT_BD( elemId < numElements() );
    ASSERT_BD( locF < element_Type::S_numLocalFacets );
    return M_ElemToFacet.operator() ( locF, elemId );
}

template <typename GEOSHAPE, typename MC>
inline
UInt
RegionMesh<GEOSHAPE, MC>::localRidgeId(edge_Type,  UInt const elemId, UInt const locR )
const
{
    ASSERT_PRE( !M_ElemToRidge.empty(), "Volume to Edges array not  set" );
    ASSERT_BD( elemId < numElements() );
    ASSERT_BD( locR < element_Type::S_numLocalRidges );
    return M_ElemToRidge( locR, elemId );
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::isFullEdge( UInt const & id ) const
{
    return edgeList.size() > id;
}

template <typename GEOSHAPE, typename MC>
inline
bool
RegionMesh<GEOSHAPE, MC>::isFullFace( UInt const & id ) const
{
    return faceList.size() > id;
}


/********************************************************************************
                     ELEMENT3D:GLOBAL FACES/EDGES
 *******************************************************************************/
// Forward Declarations
template <typename GEOSHAPE, typename MC>
void
RegionMesh<GEOSHAPE, MC>::updateElementRidges(edge_Type, bool ce, bool verbose, UInt ee, bool renumber )
{
    // If the counter is set we trust it! Otherwise we use Euler formula
    // this is ok for domains with at most 1 hole!

    if (verbose)
        std::cout << "     Updating element edges ... " << std::flush;

    renumber=renumber && ce && !  this->edgeList.empty();
    if ( ce && ee == 0 )
        ee = M_numEdges > M_numBEdges ? M_numEdges : ( GEOSHAPE::S_numFaces / 2 - 1 ) * numVolumes() + M_numBFaces / 2 + numVertices();


    if ( ce )
    {
        // We want to create the edges, we need to reserve space
        edgeList.setMaxNumItems(ee);
    }
    MeshElementBareHandler<BareEdge> _be;
    std::pair<UInt, bool> e;
    M_ElemToRidge.reshape( numLocalEdges(), numVolumes() ); // DIMENSION ARRAY

    UInt elemLocalID, i1, i2;
    std::pair<BareEdge, bool> _edge;
    GEOSHAPE ele;
    facetShape_Type bele;
    // First We check if we have already Edges stored
    if ( ! edgeList.empty() )
    {
        // dump first the existing edges, to maintain the correct numbering
        // if everything is correct the numbering in the bareedge
        // structure will reflect the actual edge numbering
        std::pair<UInt, bool> _check;
        for ( UInt j = 0; j < edgeList.size(); ++j )
        {
            i1 = ( edgeList[ j ].point( 0 ) ).id();
            i2 = ( edgeList[ j ].point( 1 ) ).id();

            _edge  = makeBareEdge( i1, i2 );
            _check = _be.addIfNotThere( _edge.first );
        }
    }

    ridge_Type edg;

    for ( typename faces_Type::iterator ifa = faceList.begin();
                    ifa != faceList.begin() + M_numBFaces; ++ifa )
    {
        for ( UInt j = 0; j < numLocalEdgesOfFace(); j++ )
        {
            i1 = bele.edgeToPoint( j, 0 );
            i2 = bele.edgeToPoint( j, 1 );
            // go to global
            i1 = ( ifa->point( i1 ) ).id();
            i2 = ( ifa->point( i2 ) ).id();

            _edge = makeBareEdge( i1, i2 );

            e = _be.addIfNotThere( _edge.first );

            if ( ce && e.second )
            {
                //
                for ( UInt k = 0; k < 2 + facetShape_Type::S_numPointsPerEdge; k++ )
                {
                    UInt inode = bele.edgeToPoint(j, k);
                    edg.setPoint( k, ifa->point( inode ) );
                }
                MeshUtility::inheritPointsWeakerMarker( edg );
                edg.setBoundary(true);
                addEdge( edg);
            }
        }

    }

    if ( ce )
    {
        M_numBEdges = edgeList.size();
        setLinkSwitch( "HAS_BOUNDARY_RIDGES" );
    }

    for ( typename elements_Type::iterator elemIt = volumeList.begin();
            elemIt != volumeList.end(); ++elemIt )
    {
        elemLocalID = elemIt->localId();

        for ( UInt j = 0; j < numLocalEdges(); j++ )
        {
            i1 = ele.edgeToPoint( j, 0 );
            i2 = ele.edgeToPoint( j, 1 );
            // go to global
            i1 = ( elemIt->point( i1 ) ).id();
            i2 = ( elemIt->point( i2 ) ).id();
            _edge = makeBareEdge( i1, i2 );

            e = _be.addIfNotThere( _edge.first );
            M_ElemToRidge.operator() ( j, elemLocalID ) = e.first;
            if ( ce && e.second )
            {
                for ( UInt k = 0; k < 2 + elementShape_Type::S_numPointsPerEdge; k++ )
                {
                    UInt inode = ele.edgeToPoint(j, k);
                    edg.setPoint( k, elemIt->point( inode ) );
                }
                MeshUtility::inheritPointsWeakerMarker( edg );
                edg.setBoundary(true);;
                addEdge( edg);
            }
        }
    }

    if ( ce )
    {
        M_numEdges = edgeList.size();
        this->M_numBEdges=
        edgeList.countElementsWithFlag(EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);
        setLinkSwitch( "HAS_ALL_RIDGES" );
    }

    if(renumber && !edgeList.empty())
    {
        edgeList.reorderAccordingToFlag(EntityFlags::PHYSICAL_BOUNDARY,&Flag::testOneSet);
        std::vector<ID>newToOld=edgeList.resetId(); //reset the ids so that they are in accord with position in the container.
        //Unfortunately I need oldToNew!
        std::vector<ID> oldToNew( newToOld.size() );
        for (UInt j=0;j<newToOld.size();++j) oldToNew[ newToOld[j] ]=j;
        // Save some memory annihilating newToOld
        std::vector<ID>().swap(newToOld);
        // Fix volume to edge array to reflect new edge numbering
        // M_ElemToRidge is in fact a vector!
        std::vector<UInt> tmp( M_ElemToRidge.size() );
        std::vector<UInt>::iterator tmp_it=tmp.begin();
        for (std::vector<UInt>::iterator it=M_ElemToRidge.begin();it<M_ElemToRidge.end();++it,++tmp_it)
            *tmp_it = oldToNew[*it];
        std::copy(tmp.begin(),tmp.end(),M_ElemToRidge.begin());
    }

    UInt n = _be.maxId();

    if (!ce)
    {
        if ( M_numEdges == 0 || M_numEdges == M_numBEdges )
            M_numEdges = n;
    }

    if (verbose)
        std::cout << n << " edges found";
    ASSERT_POS( n == M_numEdges , "#Edges found is not equal to that in RegionMesh" << n << " " << M_numEdges ) ;
    setLinkSwitch( std::string( "HAS_ELEMENT_TO_RIDGES" ) );

    if (verbose)
        std::cout << " done." << std::endl;
}


//
// Update element faces
//

template <typename GEOSHAPE, typename MC>
void
RegionMesh<GEOSHAPE, MC>::updateElementFacets( bool cf, const bool verbose, UInt ef )
{

	typedef BareEntitySelector<facetShape_Type::S_numVertices> bareEntitySelector_Type;
	typedef typename bareEntitySelector_Type::bareEntity_Type bareFacet_type;

    if (verbose)
        std::cout << "     Updating element facets ... " << std::flush;

    ASSERT0( ! cf || numBFacets() > 0, std::stringstream( std::string("Boundary Facets Must have been set") +
                                                         std::string("in order to call updateElementFacets with createFacets=true") +
                                                         std::string("\nUse buildBoundaryFacets(..) from mesh_util.h") ).str().c_str() );
    // If the counter is set we trust it! Otherwise we use Euler formula

    if ( cf && ef == 0 )
        ef = numFacets() > numBFacets() ? numFacets() : ( elementShape_Type::S_numFacets * numElements() + numBFacets() ) / 2;

    ASSERT( cf || numFacets() > 0 , "Mesh is not properly set!" );

    if ( cf )
        facetList().setMaxNumItems( ef );



    facet_Type aFacet;

    MeshElementBareHandler<bareFacet_type> _be;
    // Extra map for facets stored which are not boundary facets
    MeshElementBareHandler<bareFacet_type> _extraFacets;
    std::pair<UInt, bool> e;
    M_ElemToFacet.reshape( element_Type::S_numLocalFacets, numElements() ); // DIMENSION ARRAY

    UInt elemLocalID;
    std::pair<bareFacet_type, bool>_facet;

    GEOSHAPE ele;
    // If we have all facets and the facets store all adjacency info
    // everything is easier
    if ( (facetList().size() == numFacets()) & getLinkSwitch( "FACETS_HAVE_ADIACENCY" ) & getLinkSwitch( "HAS_ALL_FACETS" ) )
    {
        for ( typename facets_Type::iterator itf = facetList().begin(); itf != facetList().end(); ++itf )
        {
            if ( itf->firstAdjacentElementPosition() != NotAnId && itf->firstAdjacentElementIdentity() != NotAnId)
                M_ElemToFacet( itf->firstAdjacentElementPosition() , itf->firstAdjacentElementIdentity() ) = itf->localId();
            if ( itf->secondAdjacentElementPosition() != NotAnId && itf->secondAdjacentElementIdentity() != NotAnId)
                M_ElemToFacet( itf->secondAdjacentElementPosition(), itf->secondAdjacentElementIdentity() ) = itf->localId();
        }
        // we finish here
        setLinkSwitch( "HAS_ELEMENT_TO_FACETS" );
        if (verbose)
            std::cout << " done." << std::endl;

        return ;
    }

    // If I have only boundary facets I need to process them first to keep the correct numbering

    // First We check if we have already Facets stored
    UInt _numOriginalStoredFacets=facetList().size();
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
        		points[k] = ( facet( j ).point( k ) ).localId();
           	_facet = bareEntitySelector_Type::makeBareEntity( points );
            _check = _be.addIfNotThere( _facet.first );
            if (j>=this->numBFacets() )_extraFacets.addIfNotThere( _facet.first, j);
        }
    }

    for ( typename elements_Type::iterator elemIt = elementList().begin();
    		elemIt != elementList().end(); ++elemIt )
    {
        elemLocalID = elemIt->localId();
        for ( UInt j = 0; j < element_Type::S_numLocalFacets; j++ )
        {
        	for (UInt k = 0; k < facetShape_Type::S_numVertices; k++)
        	{
        		UInt id = ele.facetToPoint( j, k );
        		points[k] = elemIt->point( id ).localId();
        	}
            _facet = bareEntitySelector_Type::makeBareEntity( points );

            e = _be.addIfNotThere( _facet.first );
            M_ElemToFacet( j, elemLocalID ) = e.first;
            bool _isBound=e.first<this->numBFacets();
            // Is the facet an extra facet (not on the boundary but originally included in the list)?
            bool _isExtra = (e.first >=this->numBFacets()  && e.first < _numOriginalStoredFacets);
            if (_isBound)
            {
                facet_Type & _thisFacet(facet(e.first) );
                _thisFacet.firstAdjacentElementIdentity()   = elemLocalID;
                _thisFacet.firstAdjacentElementPosition()   = j;
                _thisFacet.secondAdjacentElementIdentity()  = NotAnId;
                _thisFacet.secondAdjacentElementPosition()  = NotAnId;
            }
            else if (_isExtra)
            {
                // This is not a bfacets and I need to set up all info about adjacency properly
                facet_Type & _thisFacet(facet(e.first) );
                // I need to check if it is the first time I meet it. Then I delete it from the
                // map: if it as there it means that it is the first time I am treating this face
                if(_extraFacets.deleteIfThere(_facet.first)){
                    // I need to be sure about orientation, the easiest thing is to rewrite the facet points
                    for ( UInt k = 0; k < facet_Type::S_numPoints; ++k )
                        _thisFacet.setPoint( k, elemIt->point( ele.facetToPoint( j, k ) ) );
                    _thisFacet.firstAdjacentElementIdentity()  = elemLocalID;
                    _thisFacet.firstAdjacentElementPosition()  = j;

                }else{
                    _thisFacet.secondAdjacentElementIdentity()  = elemLocalID;
                    _thisFacet.secondAdjacentElementPosition()  = j;
                }
            }
            else if ( cf ) // A face not contained in the original list.
                // I process it only if requested!
            {
                if ( e.second )
                {
                    // a new face It must be internal.
                    for ( UInt k = 0; k < facet_Type::S_numPoints; ++k )
                        aFacet.setPoint( k, elemIt->point( ele.facetToPoint( j, k ) ) );

                    aFacet.firstAdjacentElementIdentity()  = elemLocalID;
                    aFacet.firstAdjacentElementPosition() = j;

                    // gets the marker from the RegionMesh

                    aFacet.setMarker( this->marker() );
                    aFacet.setBoundary(false);
                    addFacet( aFacet); //The id should be correct
                }
                else
                {
                    facet( e.first ).secondAdjacentElementIdentity() = elemLocalID;
                    facet( e.first ).secondAdjacentElementPosition() = j;
                }
            }
        }
    }

    UInt n = _be.maxId();
    // LF Fix _numfacets. This part has to be checked. One may want to use
    // this method on a partitioned mesh, in which case the Global facets are there
    setNumFacets(n); // We have found the right total number of facets in the mesh
    if(numGlobalFacets()==0) setMaxNumGlobalFacets(n); // If not already set fix it.

    if (verbose) std::cout << n << " facets ";
    //ASSERT_POS( n == M_numFacets , "#Facets found inconsistent with that stored in RegionMesh" ) ;
    setLinkSwitch( "HAS_ELEMENT_TO_FACETS" );
    if ( cf )
        setLinkSwitch( "HAS_ALL_FACETS" );
    //if ( cf ) Facets have adjacency in any case!
        setLinkSwitch( "FACETS_HAVE_ADIACENCY" );
    if (verbose)
        std::cout << " done." << std::endl;
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh<GEOSHAPE, MC>::cleanElementFacets()
{
    M_ElemToFacet.clearArray();

    unsetLinkSwitch( "HAS_ELEMENT_TO_FACETS" );
}

template <typename GEOSHAPE, typename MC>
void
RegionMesh<GEOSHAPE, MC>::cleanElementRidges()
{
    M_ElemToRidge.clearArray();
    unsetLinkSwitch( "HAS_ELEMENT_TO_RIDGES" );
}


template <typename GEOSHAPE, typename MC>
inline void
RegionMesh<GEOSHAPE, MC>::
getListOfPoints( bool ( *fct ) ( double, double, double ), std::vector<UInt>& list_pts )
{
    for ( UInt i = 0; i < M_numPoints; i++ )
    {
        MeshVertex& pt = pointList( i );
        if ( fct( pt.x(), pt.y(), pt.z() ) )
        {
            list_pts.push_back( i );
        }
    }
}


template <typename GEOSHAPE, typename MC>
inline MeshUtility::MeshTransformer<RegionMesh<GEOSHAPE, MC> > &
RegionMesh<GEOSHAPE, MC>::meshTransformer()
{
    return this->M_meshTransformer;
}

} // End of namespace LifeV

#endif //REGIONMESH_H

