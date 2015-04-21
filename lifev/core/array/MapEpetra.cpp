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

#include <Epetra_Util.h>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MapEpetra.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

MapEpetra::MapEpetra () :
    M_exporter (new boost::shared_ptr<Epetra_Export>()),
    M_importer (new boost::shared_ptr<Epetra_Import>())
{
    // Nothing to be done here
}

MapEpetra::MapEpetra ( Int  numGlobalElements,
                       Int  numMyElements,
                       Int* myGlobalElements,
                       const commPtr_Type& commPtr ) :
    M_commPtr ( commPtr )
{
    ASSERT (M_commPtr.get()!=0, "Error! The communicator pointer is not valid.\n");

    createMap ( numGlobalElements,
                numMyElements,
                myGlobalElements,
                *commPtr );
}

MapEpetra::MapEpetra ( mapData_Type const& mapData, commPtr_Type const& commPtr ) :
    M_exporter (new boost::shared_ptr<Epetra_Export>()),
    M_importer (new boost::shared_ptr<Epetra_Import>()),
    M_commPtr  ( commPtr )
{
    ASSERT (M_commPtr.get()!=0, "Error! The communicator pointer is not valid.\n");

    M_uniqueMapEpetra.reset ( new Epetra_Map ( -1,
                                               mapData.unique.size(),
                                               mapData.unique.data(),
                                               0,
                                               *M_commPtr ) );
    M_repeatedMapEpetra.reset ( new Epetra_Map ( -1,
                                                 mapData.repeated.size(),
                                                 mapData.repeated.data(),
                                                 0,
                                                 *M_commPtr ) );
}

MapEpetra::MapEpetra ( const Int numGlobalElements,
                       const Int /*notUsed*/,
                       const commPtr_Type& commPtr ) :
    M_exporter (new boost::shared_ptr<Epetra_Export>()),
    M_importer (new boost::shared_ptr<Epetra_Import>()),
    M_commPtr  ( commPtr )
{
    ASSERT (M_commPtr.get()!=0, "Error! The communicator pointer is not valid.\n");

    std::vector<Int> myGlobalElements ( numGlobalElements );

    for ( Int i = 0; i < numGlobalElements; ++i )
    {
        myGlobalElements[i] = i;
    }

    M_repeatedMapEpetra.reset ( new Epetra_Map ( -1, numGlobalElements, &myGlobalElements[0], 0, *commPtr ) );
    M_uniqueMapEpetra.reset ( new Epetra_Map ( numGlobalElements, 0, *commPtr ) );
}

MapEpetra::MapEpetra ( const Int           size,
                       const commPtr_Type& commPtr ) :
    M_exporter (new boost::shared_ptr<Epetra_Export>()),
    M_importer (new boost::shared_ptr<Epetra_Import>()),
    M_commPtr  ( commPtr )
{
    ASSERT (M_commPtr.get()!=0, "Error! The communicator pointer is not valid.\n");

    Int numGlobalElements ( size );
    Int numMyElements    ( numGlobalElements );
    std::vector<Int>  myGlobalElements ( size );

    for ( Int i (0); i < numGlobalElements; ++i )
    {
        myGlobalElements[i] = i;
    }
    M_repeatedMapEpetra.reset ( new Epetra_Map ( numGlobalElements,
                                                 numMyElements,
                                                 &myGlobalElements[0],
                                                 0,
                                                 *commPtr ) );

    if ( commPtr->MyPID() != 0 )
    {
        numMyElements = 0;
    }

    M_uniqueMapEpetra.reset ( new Epetra_Map ( numGlobalElements,
                                               numMyElements,
                                               &myGlobalElements[0],
                                               0,
                                               *commPtr ) );
}

MapEpetra::MapEpetra ( const map_Type map ) :
    M_repeatedMapEpetra ( new map_Type ( map ) ),
    M_commPtr(map.Comm().Clone())
{
    uniqueMap ();
}

MapEpetra::MapEpetra ( const Epetra_BlockMap& blockMap, const Int offset, const Int maxId) :
    M_commPtr(blockMap.Comm().Clone())
{
    std::vector<Int> myGlobalElements;
    Int* sourceGlobalElements ( blockMap.MyGlobalElements() );
    Int const startIdOrig ( offset );
    Int const endIdOrig  ( startIdOrig + maxId );
    const Int maxMyElements = std::min ( maxId, blockMap.NumMyElements() );
    myGlobalElements.reserve ( maxMyElements );

    //Sort MyGlobalElements to avoid a bug in Trilinos (9?) when multiplying two matrices (A * B^T)
    std::sort ( myGlobalElements.begin(), myGlobalElements.end() );

    // We consider that the source Map may not be ordered
    for ( Int i (0); i < blockMap.NumMyElements(); ++i )
        if ( sourceGlobalElements[i] < endIdOrig && sourceGlobalElements[i] >= startIdOrig )
        {
            myGlobalElements.push_back ( sourceGlobalElements[i] - offset );
        }

    createMap ( -1,
                myGlobalElements.size(),
                &myGlobalElements.front(),
                *M_commPtr );
}

// ===================================================
// Operators
// ===================================================

MapEpetra& MapEpetra::operator= (const MapEpetra& epetraMap)
{
    if ( this != &epetraMap )
    {
        M_repeatedMapEpetra = epetraMap.M_repeatedMapEpetra;
        M_uniqueMapEpetra   = epetraMap.M_uniqueMapEpetra;
        M_exporter          = epetraMap.M_exporter;
        M_importer          = epetraMap.M_importer;
        M_commPtr           = epetraMap.M_commPtr;
    }

    return *this;
}

MapEpetra& MapEpetra::operator+= ( const MapEpetra& epetraMap )
{
    if ( ! epetraMap.getUniqueMap() )
    {
        return *this;
    }

    if ( ! this->getUniqueMap() )
    {
        this->operator = ( epetraMap );
        return *this;
    }

    if (M_commPtr.get()==0)
    {
        // In case this map was created with default constructor
        M_commPtr = epetraMap.M_commPtr;
    }

    Int*             pointer;
    std::vector<Int> map;

    pointer = getRepeatedMap()->MyGlobalElements();
    for ( Int ii = 0; ii < getRepeatedMap()->NumMyElements(); ++ii, ++pointer )
    {
        map.push_back ( *pointer );
    }

    Int numGlobalElements = getUniqueMap()->NumGlobalElements();

    pointer = epetraMap.getRepeatedMap()->MyGlobalElements();
    for (Int ii = 0; ii < epetraMap.getRepeatedMap()->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back ( *pointer + numGlobalElements );
    }

    M_repeatedMapEpetra.reset ( new Epetra_Map (-1, map.size(), map.data(), 0, *M_commPtr ) );

    map.resize (0);
    pointer = getUniqueMap()->MyGlobalElements();

    for ( Int ii = 0; ii < getUniqueMap()->NumMyElements(); ++ii, ++pointer )
    {
        map.push_back ( *pointer );
    }

    pointer = epetraMap.getUniqueMap()->MyGlobalElements();
    for ( Int ii = 0; ii < epetraMap.getUniqueMap()->NumMyElements(); ++ii, ++pointer )
    {
        map.push_back ( *pointer + numGlobalElements );
    }

    M_uniqueMapEpetra.reset ( new Epetra_Map ( -1, map.size(), map.data(), 0, *M_commPtr ) );

    M_exporter.reset (new boost::shared_ptr<Epetra_Export>());
    M_importer.reset (new boost::shared_ptr<Epetra_Import>());

    return *this;
}

MapEpetra& MapEpetra::operator += ( Int const size )
{
    MapEpetra  lagrMap ( size, commPtr() );

    ASSERT ( this->getUniqueMap(), "operator+=(const Int) works only for an existing MapEpetra" );

    this->operator+= ( lagrMap );
    return *this;
}

// ===================================================
// Methods
// ===================================================

boost::shared_ptr<MapEpetra> MapEpetra::createRootMap (Int const root) const
{
    boost::shared_ptr<MapEpetra> rootMap ( new MapEpetra ( Epetra_Util::Create_Root_Map ( *getUniqueMap(), root ) ) );
    return rootMap;
}

bool MapEpetra::mapsAreSimilar ( MapEpetra const& epetraMap ) const
{
    if ( this == &epetraMap )
    {
        return true;
    }

    return ( getUniqueMap()->SameAs ( *epetraMap.getUniqueMap() ) &&
             getRepeatedMap()->SameAs ( *epetraMap.getRepeatedMap() ) );
}

#ifdef HAVE_HDF5

void MapEpetra::exportToHDF5 ( std::string const& fileName, std::string const& mapName, bool const truncate )
{
    ASSERT (M_commPtr.get()!=0, "Error! The stored communicator pointer is not valid.\n");
    ASSERT (M_uniqueMapEpetra.get()!=0 && M_repeatedMapEpetra.get()!=0, "Error! One (or both) the map pointers are not valid.\n");

    EpetraExt::HDF5 HDF5 ( *M_commPtr );

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

    // Save the maps into the file
    HDF5.Write ( ( mapName + "Unique" ).c_str(), *M_uniqueMapEpetra );
    HDF5.Write ( ( mapName + "Repeated" ).c_str(), *M_repeatedMapEpetra );

    // Close the file
    HDF5.Close();

} // exportToHDF5

void MapEpetra::importFromHDF5 ( std::string const& fileName, std::string const& mapName )
{
    ASSERT (M_commPtr.get()!=0, "Error! The stored communicator pointer is not valid.\n");

    EpetraExt::HDF5 HDF5 ( *M_commPtr );

    // Open an existing file
    HDF5.Open ( ( fileName + ".h5" ).data() );

    // Check if the file is created
    if ( !HDF5.IsOpen () )
    {
        std::cerr << "Unable to open " + fileName + ".h5";
        abort();
    }

    // Read the unique map from the file
    Epetra_Map* importedMap ( 0 );
    HDF5.Read ( ( mapName + "Unique" ).c_str(), importedMap );

    // Copy the loaded map to the member object
    M_uniqueMapEpetra.reset ( new map_Type ( *importedMap ) );

    // Read the repeated map from the file
    HDF5.Read ( ( mapName + "Repeated" ).c_str(), importedMap );

    // Copy the loaded matrix to the member object
    M_repeatedMapEpetra.reset ( new map_Type ( *importedMap ) );

    // Close the file
    HDF5.Close();

} // importFromHDF5

#endif // HAVE_HDF5

void MapEpetra::showMe ( std::ostream& output ) const
{
    output << "unique map:" << std::endl;
    output << *getUniqueMap();
    output << "repeated map:" << std::endl;
    output << *getRepeatedMap();
}

// ===================================================
// Get Methods
// ===================================================

const MapEpetra::mapPtr_Type& MapEpetra::map (MapEpetraType mapType) const
{
    switch ( mapType )
    {
        case Unique:
            return getUniqueMap();
        case Repeated:
            return getRepeatedMap();
    }
    return getUniqueMap();
}

const Epetra_Export& MapEpetra::exporter()
{
    ASSERT (M_uniqueMapEpetra.get()!=0 && M_repeatedMapEpetra.get()!=0, "Error! One (or both) the map pointers are not valid.\n");

    createImportExport();

    return **M_exporter;
}

const Epetra_Import& MapEpetra::importer()
{
    ASSERT (M_uniqueMapEpetra.get()!=0 && M_repeatedMapEpetra.get()!=0, "Error! One (or both) the map pointers are not valid.\n");

    createImportExport();

    return **M_importer;
}

// ===================================================
// Set Methods
// ===================================================

void MapEpetra::setComm (commPtr_Type const& commPtr)
{
    ASSERT (commPtr.get()!=0, "Error! The communicator pointer is not valid.\n");

    M_commPtr = commPtr;
}

void MapEpetra::setMap ( mapPtr_Type map, MapEpetraType mapType )
{
    switch ( mapType )
    {
        case Unique:
            M_uniqueMapEpetra = map;
        case Repeated:
            M_repeatedMapEpetra = map;
    }
}

// ===================================================
// Private Methods
// ===================================================

MapEpetra::MapEpetra (const MapEpetra& epetraMap)
{
    this->operator= ( epetraMap );
}

void MapEpetra::createMap (Int  numGlobalElements,
                           Int  numMyElements,
                           Int* myGlobalElements,
                           const comm_Type& comm)
{
    if ( numMyElements != 0 && myGlobalElements == 0 ) // linearMap
        M_repeatedMapEpetra.reset ( new Epetra_Map ( numGlobalElements,
                                                     numMyElements,
                                                     0,
                                                     comm ) );
    else // classic LifeV map
        M_repeatedMapEpetra.reset ( new Epetra_Map ( numGlobalElements,
                                                     numMyElements,
                                                     myGlobalElements,
                                                     0,
                                                     comm ) );

    uniqueMap();
}

void MapEpetra::uniqueMap()
{
    M_uniqueMapEpetra.reset ( new Epetra_Map ( Epetra_Util::Create_OneToOne_Map ( *getRepeatedMap(), false ) ) );

    M_exporter.reset (new boost::shared_ptr<Epetra_Export>());
    M_importer.reset (new boost::shared_ptr<Epetra_Import>());
}

void MapEpetra::createImportExport()
{
    if ( !getRepeatedMap() || !getUniqueMap() )
    {
        return;
    }

    if ( M_exporter->get() == 0 )
    {
        M_exporter->reset ( new Epetra_Export ( *getRepeatedMap(), *getUniqueMap() ) );
    }

    if ( M_importer->get() == 0 )
    {
        M_importer->reset ( new Epetra_Import ( *getRepeatedMap(), *getUniqueMap() ) );
    }
}

void
MapEpetra::bubbleSort (Epetra_IntSerialDenseVector& elements)
{
    Int hold;

    for ( Int pass (0); pass < elements.Length() - 1; pass++ )
        for ( Int j (0); j < elements.Length() - 1; j++ )
            if ( elements[j] > elements[j + 1] )
            {
                hold          = elements[j];
                elements[j]   = elements[j + 1];
                elements[j + 1] = hold;
            }
}

// ===================================================
// External operators
// ===================================================

MapEpetra operator+ (const MapEpetra& map1, const MapEpetra& map2)
{
    MapEpetra mapOut (map1);
    mapOut += map2;

    return mapOut;
}

MapEpetra operator+ (const MapEpetra& map, Int size)
{
    MapEpetra mapOut (map);
    mapOut += size;

    return mapOut;
}

MapVector<MapEpetra> operator| (const MapEpetra& map1, const MapEpetra& map2)
{
    return MapVector<MapEpetra> (map1, map2);
}

} // end namespace LifeV

