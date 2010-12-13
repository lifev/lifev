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
    @brief EpetraMap

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 26-10-2006

    This class manages the distribution of elements of matrices or vectors on a parallel machine
 */

#ifndef _EPETRAMAP_
#define _EPETRAMAP_

#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Map.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Comm.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifefem/refFE.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/partitionMesh.hpp>


namespace LifeV
{

enum EpetraMapType {Unique = 0, Repeated};


//! EpetraMap - Wrapper for Epetra_Map
/*!
  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

  The EpetraMap class provides a general interface for the Epetra_Map class of Trilinos.

  Visit http://trilinos.sandia.gov for more informations about Epetra_Map.
 */
class EpetraMap
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_Map                                            map_type;
    typedef boost::shared_ptr<map_type>                           map_ptrtype;
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Export> > exporter_ptrtype;
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Import> > importer_ptrtype;
    typedef Epetra_Comm                                           comm_type;
    typedef boost::shared_ptr<comm_type>                          comm_ptrtype;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    EpetraMap();

    //! Constructor
    /*!
      To define a linear map, set MyGlobalElements = 0
      @param numGlobalElements Number of global elements
      @param numMyElements Number of local elements
      @param myGlobalElements Array of Id of the local element
      @param indexBase Starting index base (typically 0 or 1)
      @param commPtr Pointer to the communicator
    */
    EpetraMap( Int  numGlobalElements,
               Int  numMyElements,
               Int* myGlobalElements,
               Int  indexBase,
               const comm_ptrtype& commPtr );

    //! Constructor
    /*!
      Build a nearly equally distributed map.

      The map is equally distributed if NumGlobalElements % Rank = 0
      @param NumGlobalElements Total number of elements inside the map
      @param indexBase Starting index base (typically 0 or 1)
      @param CommPtr A pointer to the Epetra communicator
     */
    EpetraMap( const Int numGlobalElements,
               const Int indexBase,
               const comm_ptrtype& commPtr );

    //! Constructor
    /*!
      @param size Size of the map
      @param commPtr Pointer to the communicator
     */
    EpetraMap( const Int           size,
               const comm_ptrtype& commPtr );

    //! Constructor
    /*!
      Calls createImportExport from setUp()
      @param refFE Reference finite element
      @param meshPart Partition of the mesh
      @param commPtr Pointer to the communicator
     */
    template<typename Mesh>
    EpetraMap( const RefFE&               refFE,
               const partitionMesh<Mesh>& meshPart,
               const comm_ptrtype&        commPtr );

    //! Constructor
    /*!
      Calls createImportExport from setUp()
      @param refFE Reference finite element
      @param mesh Mesh
      @param commPtr Pointer to the communicator
    */
    template<typename Mesh>
    EpetraMap( const RefFE&         refFE,
               const Mesh&          mesh,
               const comm_ptrtype&  commPtr );

    //! Copy constructor
    /*!
      @param epetraMap An EpetraMap object
     */
    EpetraMap( const EpetraMap&  epetraMap );

    //! Constructor
    /*!
      Builds a submap of map _epetraMap with a given positive offset and
      the maximum id to consider

      e.g:
      <ol>
      <li> offset = 2,
      <li> maxid = 6
      <li> epetraMap = [ 0 2 5 7 8 10 1 ]
      <li> this = [ 0 3 5 6 ]
      </ol>

      if needed, indexBase may be changed (default values < 0 means "same as original map")
      @param blockMap Epetra_BlockMap
      @param offset Offset to be used to build the map
      @param maxId Maximum Id
      @param indexBase Starting index base (typically 0 or 1)
    */
    EpetraMap( const Epetra_BlockMap& blockMap,
               const Int offset,
               const Int maxId,
               Int indexBase = -1 );
private:
    //! Constructor from raw Epetra_Map
    /*!
      This constructor should be used only inside this class,
      therefore it is private
      @param map: underlying Epetra_Map
     */
    EpetraMap( const map_type map );

public:
    //! Destructor
    ~EpetraMap() {}

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
      The assignment operator will copy the pointers of the maps, exporter and importer
      @param epetraMap EpetraMap to be assigned to the current matrix
     */
    EpetraMap& operator  = ( const EpetraMap& epetraMap );

    //! Addition operator
    /*!
      The addition operator combines two map together
      @param epetraMap EpetraMap to be combined with the current map
     */
    EpetraMap& operator += ( const EpetraMap& epetraMap );

    //! Addition operator
    /*!
      The addition operator combines two map together to create a new map
      @param epetraMap EpetraMap to be combined with the current map
     */
    EpetraMap operator +  ( const EpetraMap& epetraMap );

    //! Addition operator
    /*!
      The addition operator create a map of size "size" and add it to the current map.
      @param size Size of the map to be added to the current map
     */
    EpetraMap& operator += ( Int const size );

    //! Addition operator
    /*!
      The addition operator create a map of size "size" and add it to the current map
      to create a new map
      @param size Size of the map to be added to the current map
     */
    EpetraMap operator +  ( Int const size );

    //@}

    //! @name Methods
    //@{

    //! This method creates a pointer to a EpetraMap that has points only on processor root
    /*!
      @param root processor on which to export all the points
     */
    boost::shared_ptr<EpetraMap> createRootMap( Int const root ) const;

    //! This method return true if both the unique map and the repeated map are identical
    bool MapsAreSimilar( EpetraMap const& epetraMap ) const;

    //! Show informations about the map
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Get Methods
    //@{

    //! Return the communicator
    comm_type const& Comm() const { return *M_commPtr; }

    //! Return a shared pointer on the communicator
    comm_ptrtype& CommPtr() { return M_commPtr; }

    //! Return a shared pointer on the internal Epetra_Map
    map_ptrtype const & getMap  ( EpetraMapType mapType ) const;

    //! Getter for the Epetra_Export
    Epetra_Export const& getExporter();

    //! Getter for the Epetra_Import
    Epetra_Import const& getImporter();

    //@}

private:

    //! @name Private Methods
    //@{

    //! Create a map
    /*!
      Note: createMap does not call createImportExport
      @param numGlobalElements Number of global elements of the map
      @param numMyElements number of local element
      @param myGlobalElements Array of Id of the global elements of the map
      @param indexBase Starting index base (typically 0 or 1)
      @param comm Communicator
    */
    void createMap( Int   numGlobalElements,
                    Int   numMyElements,
                    Int*  myGlobalElements,
                    Int   indexBase,
                    const comm_type& comm );

    //! Getter for the repeated map
    map_ptrtype const & getRepeatedMap() const { return M_repeatedEpetraMap; }

    //! Getter for the unique map
    map_ptrtype const & getUniqueMap()   const { return M_uniqueEpetraMap; }

    //! Reset the internal unique map and recompute it using the repeated map
    void  uniqueMap();

    //! Reset and rebuild the importer and exporter for the map
    void  createImportExport();

    //! Sort the element given using a bubble sort algorithm
    /*!
      @param elements Epetra_IntSerialDenseVector vector to be sorted
     */
    void  bubbleSort( Epetra_IntSerialDenseVector& elements );

    //! Setup a map using the finite element and using nodes, edges, faces and volumes numbering
    /*!
      Calls createImportExport
      @param refFE Reference finite element
      @param commPtr Pointer on the communicator
      @param repeatedNodeVector Vector containing the node ids
      @param repeatedEdgeVector Vector containing the edge ids
      @param repeatedFaceVector Vector containing the face ids
      @param repeatedVolumeVector Vector containing the volume ids
     */
    void setUp( const RefFE&        refFE,
                const comm_ptrtype& commPtr,
                std::vector<Int>& repeatedNodeVector,
                std::vector<Int>& repeatedEdgeVector,
                std::vector<Int>& repeatedFaceVector,
                std::vector<Int>& repeatedVolumeVector );

    //@}

    map_ptrtype        M_repeatedEpetraMap;
    map_ptrtype        M_uniqueEpetraMap;
    exporter_ptrtype   M_exporter;
    importer_ptrtype   M_importer;
    comm_ptrtype       M_commPtr;

};


// ===================================================
// Constructors & Destructor
// ===================================================
template<typename Mesh>
EpetraMap::
EpetraMap( const RefFE&               refFE,
           const partitionMesh<Mesh>& meshPart,
           const comm_ptrtype&        commPtr ):
        M_repeatedEpetraMap(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr( commPtr )
{
    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

    setUp( refFE,
           commPtr,
           const_cast<std::vector<Int>&>( meshPart.repeatedNodeVector() ),
           const_cast<std::vector<Int>&>( meshPart.repeatedEdgeVector() ),
           const_cast<std::vector<Int>&>( meshPart.repeatedFaceVector() ),
           const_cast<std::vector<Int>&>( meshPart.repeatedVolumeVector() ) );

}


template<typename Mesh>
EpetraMap::
EpetraMap( const RefFE&        refFE,
           const Mesh&         mesh,
           const comm_ptrtype& commPtr ):
        M_repeatedEpetraMap(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr( commPtr )
{
    std::vector<Int> repeatedNodeVector;
    std::vector<Int> repeatedEdgeVector;
    std::vector<Int> repeatedFaceVector;
    std::vector<Int> repeatedVolumeVector;

    if ( refFE.nbDofPerVertex() )
    {
        repeatedNodeVector.reserve(mesh.numPoints());
        for ( UInt ii = 1; ii <= mesh.numPoints(); ii++ )
            repeatedNodeVector.push_back( mesh.pointList(ii).id() );
    }

    if ( refFE.nbDofPerEdge() )
    {
        repeatedEdgeVector.reserve( mesh.numEdges() );

        for ( UInt ii = 1; ii <= mesh.numEdges(); ii++ )
            repeatedEdgeVector.push_back( mesh.edgeList(ii).id() );
    }

    if ( refFE.nbDofPerFace() )
    {
        repeatedFaceVector.reserve( mesh.numFaces() );

        for ( UInt ii = 1; ii <= mesh.numFaces(); ii++ )
            repeatedFaceVector.push_back( mesh.faceList(ii).id() );
    }

    if ( refFE.nbDofPerVolume() )
    {
        repeatedVolumeVector.reserve( mesh.numVolumes() );

        for ( UInt ii = 1; ii <= mesh.numVolumes(); ii++ )
            repeatedVolumeVector.push_back( mesh.volumeList(ii).id() );
    }


    setUp( refFE,
           commPtr,
           repeatedNodeVector,
           repeatedEdgeVector,
           repeatedFaceVector,
           repeatedVolumeVector );

    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

}


}


#endif

