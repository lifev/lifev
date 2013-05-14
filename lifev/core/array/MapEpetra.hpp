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
    @brief MapEpetra

    @author Gilles Fourestey <gilles.fourestey@epfl.ch>
    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 26-10-2006

    This class manages the distribution of elements of matrices or vectors on a parallel machine
 */

#ifndef _EPETRAMAP_
#define _EPETRAMAP_

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Map.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Comm.h>

#ifdef HAVE_HDF5
#include <EpetraExt_HDF5.h>
#endif

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/EnumMapEpetra.hpp>
#include <lifev/core/array/MapEpetraData.hpp>
#include <lifev/core/array/MapVector.hpp>

namespace LifeV
{


//! MapEpetra - Wrapper for Epetra_Map
/*!
  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
  @author Simone Deparis <simone.deparis@epfl.ch>
  @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

  The MapEpetra class provides a general interface for the Epetra_Map class of Trilinos.

  Visit http://trilinos.sandia.gov for more informations about Epetra_Map.
 */
class MapEpetra
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_Map                                            map_type;
    typedef boost::shared_ptr<map_type>                           map_ptrtype;

    typedef MapEpetraData                                         mapData_Type;

    /* Double shared_ptr are used here to ensure that all the similar MapEpetra
       point to the same exporter/importer. If double shared_ptr were not used, a
       map initialized without importer/exporter (i.e. with a shared_ptr
       pointing to 0) would not "gain" the importer/exporter created by another
       map.*/
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Export> > exporter_ptrtype;
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Import> > importer_ptrtype;


    typedef Epetra_Comm                                           comm_type;
    typedef boost::shared_ptr<comm_type>                          comm_ptrtype;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    MapEpetra();

    //! Constructor
    /*!
      To define a linear map, set MyGlobalElements = 0
      @param numGlobalElements Number of global elements
      @param numMyElements Number of local elements
      @param myGlobalElements Array of Id of the local element
      @param commPtr Pointer to the communicator
    */
    MapEpetra ( Int  numGlobalElements,
                Int  numMyElements,
                Int* myGlobalElements,
                const comm_ptrtype& commPtr );

    MapEpetra ( std::pair<std::vector<Int>, std::vector<Int> > myGlobalElements,
                const comm_ptrtype& commPtr );

    MapEpetra ( MapEpetraData const & mapData, comm_ptrtype const& commPtr );

    //! Constructor
    /*
      Build a nearly equally distributed map.

      The map is equally distributed if NumGlobalElements % Rank = 0
      @param NumGlobalElements Total number of elements inside the map
      @param CommPtr A pointer to the Epetra communicator
     */
    MapEpetra ( const Int numGlobalElements,
                const Int notUsed,
                const comm_ptrtype& commPtr );

    //! Constructor
    /*!
      @param size Size of the map
      @param commPtr Pointer to the communicator
     */
    MapEpetra ( const Int           size,
                const comm_ptrtype& commPtr );

    //! Copy constructor
    /*!
      @param epetraMap An MapEpetra object
     */
    MapEpetra ( const MapEpetra&  epetraMap );

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

      @param blockMap Epetra_BlockMap
      @param offset Offset to be used to build the map
      @param maxId Maximum Id
    */
    MapEpetra ( const Epetra_BlockMap& blockMap,
                const Int offset,
                const Int maxId);

private:
    //! Constructor from raw Epetra_Map
    /*!
      This constructor should be used only inside this class,
      therefore it is private
      @param map: underlying Epetra_Map
     */
    MapEpetra ( const map_type map );

public:
    //! Destructor
    ~MapEpetra() {}

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
      The assignment operator will copy the pointers of the maps, exporter and importer
      @param epetraMap MapEpetra to be assigned to the current matrix
     */
    MapEpetra& operator  = ( const MapEpetra& epetraMap );

    //! Addition operator
    /*!
      The addition operator combines two map together
      @param epetraMap MapEpetra to be combined with the current map
     */
    MapEpetra& operator += ( const MapEpetra& epetraMap );

    //! Addition operator
    /*!
      The addition operator combines two map together to create a new map
      @param epetraMap MapEpetra to be combined with the current map
     */
    MapEpetra operator +  ( const MapEpetra& epetraMap );

    //! Addition operator
    /*!
      The addition operator create a map of size "size" and add it to the current map.
      @param size Size of the map to be added to the current map
     */
    MapEpetra& operator += ( Int const size );

    //! Addition operator
    /*!
      The addition operator create a map of size "size" and add it to the current map
      to create a new map
      @param size Size of the map to be added to the current map
     */
    MapEpetra operator +  ( Int const size );

    //! Juxtaposition operator
    /*!
      This operator is used when block structures are used. Indeed, it creates
      from two different maps a MapVector that can be used to initialize
      block structures, such as matrices and vectors (see \ref BlockAlgebraPage "this page" for examples).
     */
    MapVector<MapEpetra> operator| (const MapEpetra& map) const
    {
        return MapVector<MapEpetra> (*this, map);
    }

    //@}

    //! @name Methods
    //@{

    //! This method creates a pointer to a MapEpetra that has points only on processor root
    /*!
      @param root processor on which to export all the points
     */
    boost::shared_ptr<MapEpetra> createRootMap ( Int const root ) const;

    //! This method return true if both the unique map and the repeated map are identical
    bool mapsAreSimilar ( MapEpetra const& epetraMap ) const;

#ifdef HAVE_HDF5
    //! Save the matrix into a HDF5 (.h5) file
    /*!
      @param fileName Name of the file where the map will be saved, without extension (.h5)
      @param mapName Name of the map in the HDF5 file
      @param truncate True if the file has to be truncated; False if the file already exist and should not be truncated
     */
    void exportToHDF5 ( std::string const& fileName, std::string const& mapName = "map", bool const& truncate = true );

    //! Read a matrix from a HDF5 (.h5) file
    /*!
      @param fileName Name of the file where the map will be saved, without extension (.h5)
      @param matrixName Name of the map in the HDF5 file
     */
    void importFromHDF5 ( std::string const& fileName, std::string const& mapName = "map" );
#endif

    //! Show informations about the map
    void showMe ( std::ostream& output = std::cout ) const;

    //! Getter for the global number of entries
    UInt mapSize() const
    {
        return map (Unique)->NumGlobalElements();
    }

    //! check if a global id is owned by the current partition
    bool isOwned ( const UInt globalId ) const
    {
        return ( M_uniqueMapEpetra->LID ( static_cast<int> (globalId) ) > -1 );
    }

    //@}

    //! @name Get Methods
    //@{

    //! Return the communicator
    comm_type const& comm() const
    {
        return *M_commPtr;
    }

    //! Return a shared pointer on the communicator
    comm_ptrtype const& commPtr() const
    {
        return M_commPtr;
    }
    comm_ptrtype& commPtr()
    {
        return M_commPtr;
    }

    //! Return a shared pointer on the internal Epetra_Map
    map_ptrtype const& map ( MapEpetraType mapType ) const;

    //! Getter for the Epetra_Export
    Epetra_Export const& exporter();

    //! Getter for the Epetra_Import
    Epetra_Import const& importer();
    //@}

    //! @name Set Methods
    //@{

    //! Set the communicator
    void setComm ( comm_ptrtype const& commPtr )
    {
        M_commPtr = commPtr;
    }

    //! set the internal Epetra_Maps
    void setMap ( map_ptrtype map, MapEpetraType mapType );

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
      @param comm Communicator
    */
    void createMap ( Int   numGlobalElements,
                     Int   numMyElements,
                     Int*  myGlobalElements,
                     const comm_type& comm );

    //! Getter for the repeated map
    map_ptrtype const& getRepeatedMap() const
    {
        return M_repeatedMapEpetra;
    }

    //! Getter for the unique map
    map_ptrtype const& getUniqueMap()   const
    {
        return M_uniqueMapEpetra;
    }

    //! Reset the internal unique map and recompute it using the repeated map
    void  uniqueMap();

    //! Reset and rebuild the importer and exporter for the map
    void  createImportExport();

    //! Sort the element given using a bubble sort algorithm
    /*!
      @param elements Epetra_IntSerialDenseVector vector to be sorted
     */
    void  bubbleSort ( Epetra_IntSerialDenseVector& elements );

    //@}

    map_ptrtype        M_repeatedMapEpetra;
    map_ptrtype        M_uniqueMapEpetra;
    exporter_ptrtype   M_exporter;
    importer_ptrtype   M_importer;
    comm_ptrtype       M_commPtr;
};

typedef MapVector<MapEpetra> MapEpetraVector;

} // end namespace LifeV

#endif

