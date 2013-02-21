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

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Util.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MapEpetra.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MapEpetra::MapEpetra() :
    M_repeatedMapEpetra(),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer(),
    M_commPtr()
{}

MapEpetra::MapEpetra ( Int  numGlobalElements,
                       Int  numMyElements,
                       Int* myGlobalElements,
                       const comm_ptrtype& commPtr ) :
    M_repeatedMapEpetra(),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer(),
    M_commPtr ( commPtr )
{
    createMap ( numGlobalElements,
                numMyElements,
                myGlobalElements,
                *commPtr );
}

MapEpetra::MapEpetra ( std::pair< std::vector<Int>, std::vector<Int> > myGlobalElements,
                       const comm_ptrtype& commPtr ) :
    M_repeatedMapEpetra(),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer(),
    M_commPtr ( commPtr )
{
    std::vector<Int> const& myGlobalElementsUnique = myGlobalElements.first;
    std::vector<Int> const& myGlobalElementsRepeated = myGlobalElements.second;

    M_uniqueMapEpetra.reset ( new Epetra_Map ( -1,
                                               myGlobalElementsUnique.size(),
                                               &myGlobalElementsUnique[ 0 ],
                                               0,
                                               *M_commPtr ) );
    M_repeatedMapEpetra.reset ( new Epetra_Map ( -1,
                                                 myGlobalElementsRepeated.size(),
                                                 &myGlobalElementsRepeated[ 0 ],
                                                 0,
                                                 *M_commPtr ) );
}


MapEpetra::MapEpetra ( const Int numGlobalElements,
                       const Int /*notUsed*/,
                       const comm_ptrtype& commPtr ) :
    M_repeatedMapEpetra(),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer(),
    M_commPtr ( commPtr )
{
    std::vector<Int> myGlobalElements ( numGlobalElements );

    for ( Int i = 0; i < numGlobalElements; ++i )
    {
        myGlobalElements[i] = i;
    }

    M_repeatedMapEpetra.reset ( new Epetra_Map ( -1, numGlobalElements, &myGlobalElements[0], 0, *commPtr ) );
    M_uniqueMapEpetra.reset ( new Epetra_Map ( numGlobalElements, 0, *commPtr ) );
}

MapEpetra::MapEpetra ( const Int           size,
                       const comm_ptrtype& commPtr ) :
    M_repeatedMapEpetra(),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer(),
    M_commPtr ( commPtr )
{
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

MapEpetra::MapEpetra ( const map_type map ) :
    M_repeatedMapEpetra ( new map_type ( map ) ),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer(),
    M_commPtr()
{
    uniqueMap();
}

MapEpetra::MapEpetra ( const Epetra_BlockMap& blockMap, const Int offset, const Int maxId) :
    M_repeatedMapEpetra(),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer(),
    M_commPtr()
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
                blockMap.Comm() );
}

// ===================================================
// Operators
// ===================================================
MapEpetra&
MapEpetra::operator = ( const MapEpetra& epetraMap )
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


MapEpetra&
MapEpetra::operator += ( const MapEpetra& epetraMap )
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

    M_repeatedMapEpetra.reset ( new Epetra_Map (-1, map.size(), &map[0], 0, epetraMap.getRepeatedMap()->Comm() ) );

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

    M_uniqueMapEpetra.reset ( new Epetra_Map ( -1, map.size(), &map[0], 0, epetraMap.getRepeatedMap()->Comm() ) );

    M_exporter.reset();
    M_importer.reset();

    return *this;
}

MapEpetra
MapEpetra::operator + ( const MapEpetra& epetraMap )
{
    MapEpetra map ( *this );
    map += epetraMap;
    createImportExport();
    return map;
}

MapEpetra&
MapEpetra::operator += ( Int const size )
{
    MapEpetra  lagrMap ( size, commPtr() );

    ASSERT ( this->getUniqueMap(), "operator+=(const Int) works only for an existing MapEpetra" );

    this->operator+= ( lagrMap );
    return *this;
}

MapEpetra
MapEpetra::operator +  ( Int const size )
{
    MapEpetra map ( *this );
    map += size;
    createImportExport();
    return map;
}

// ===================================================
// Methods
// ===================================================
boost::shared_ptr<MapEpetra>
MapEpetra::createRootMap ( Int const root )   const
{
    boost::shared_ptr<MapEpetra> rootMap ( new MapEpetra ( Epetra_Util::Create_Root_Map ( *getUniqueMap(), root ) ) );
    return rootMap;
}

bool
MapEpetra::mapsAreSimilar ( MapEpetra const& epetraMap ) const
{
    if ( this == &epetraMap )
    {
        return true;
    }

    return ( getUniqueMap()->SameAs ( *epetraMap.getUniqueMap() ) &&
             getRepeatedMap()->SameAs ( *epetraMap.getRepeatedMap() ) );
}

#ifdef HAVE_HDF5

void MapEpetra::exportToHDF5 ( std::string const& fileName, std::string const& mapName, bool const& truncate )
{
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
    M_uniqueMapEpetra.reset ( new map_type ( *importedMap ) );

    // Read the repeated map from the file
    HDF5.Read ( ( mapName + "Repeated" ).c_str(), importedMap );

    // Copy the loaded matrix to the member object
    M_repeatedMapEpetra.reset ( new map_type ( *importedMap ) );

    // Close the file
    HDF5.Close();

} // importFromHDF5

#endif // HAVE_HDF5

void
MapEpetra::showMe ( std::ostream& output ) const
{
    output << "unique map:" << std::endl;
    output << *getUniqueMap();
    output << "repeated map:" << std::endl;
    output << *getRepeatedMap();
}

// ===================================================
// Get Methods
// ===================================================
MapEpetra::map_ptrtype const&
MapEpetra::map ( MapEpetraType mapType )   const
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

Epetra_Export const&
MapEpetra::exporter()
{
    createImportExport();
    return **M_exporter;
}

Epetra_Import const&
MapEpetra::importer()
{
    createImportExport();
    return **M_importer;
}

// ===================================================
// Set Methods
// ===================================================
void MapEpetra::setMap ( map_ptrtype map, MapEpetraType mapType )
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
MapEpetra::MapEpetra ( const MapEpetra& epetraMap ) :
    M_repeatedMapEpetra(),
    M_uniqueMapEpetra(),
    M_exporter(),
    M_importer()
{
    this->operator= ( epetraMap );
}


void
MapEpetra::createMap ( Int  numGlobalElements,
                       Int  numMyElements,
                       Int* myGlobalElements,
                       const comm_type& comm )
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


void
MapEpetra::uniqueMap()
{
    M_uniqueMapEpetra.reset ( new Epetra_Map ( Epetra_Util::Create_OneToOne_Map ( *getRepeatedMap(), false ) ) );
    M_exporter.reset();
    M_importer.reset();
    return;
}

void
MapEpetra::createImportExport()
{

    if ( !getRepeatedMap() || !getUniqueMap() )
    {
        return;
    }

    // The exporter is needed to import to a repeated vector
    if ( M_exporter.get() == 0 )
    {
        M_exporter.reset ( new boost::shared_ptr<Epetra_Export> );
    }

    if ( M_exporter->get() == 0 )
    {
        M_exporter->reset ( new Epetra_Export ( *getRepeatedMap(), *getUniqueMap() ) );
    }

    if ( M_importer.get() == 0 )
    {
        M_importer.reset ( new boost::shared_ptr<Epetra_Import> );
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



} // end namespace LifeV

