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

#include <lifeconfig.h>
#include <life/lifealg/EpetraMap.hpp>
#include <Epetra_Util.h>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
EpetraMap::EpetraMap():
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr()
{}

EpetraMap::EpetraMap(int                NumGlobalElements,
                     int                NumMyElements,
                     int*               MyGlobalElements,
                     int                IndexBase,
                     const comm_ptrtype&  CommPtr):
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr(CommPtr)
{

    //Sort MyGlobalElements to avoid a bug in Trilinos (9?) when multiplying two matrices (A * B^T)
    std::sort (MyGlobalElements, MyGlobalElements + NumMyElements);

    createMap( NumGlobalElements,
               NumMyElements,
               MyGlobalElements,
               IndexBase,
               *CommPtr );
}

/*
//! construct a map with entries [1:lagrangeMultipliers] distributed on all the processors
EpetraMap::EpetraMap(std::vector<int> const& lagrangeMultipliers,
                     const Epetra_Comm& Comm):
    M_repeatedEpetra_Map(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer()
{
    int NumGlobalElements(lagrangeMultipliers.size());
    int NumMyElements    (NumGlobalElements);
    std::vector<int>  MyGlobalElements(lagrangeMultipliers);
    int IndexBase = 1;

    for (int i(0); i < NumGlobalElements; ++i)
        MyGlobalElements[i] = i + 1;


    createMap( NumGlobalElements,
               NumMyElements,
               &MyGlobalElements[0],
               IndexBase,
               Comm);

}
*/

EpetraMap::EpetraMap( const int          NumGlobalElements,
                      const int          IndexBase,
                      const comm_ptrtype& CommPtr ) :
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr(CommPtr)
{
    std::vector<int> MyGlobalElements( NumGlobalElements );

    for ( int i = 0; i < NumGlobalElements; ++i )
        MyGlobalElements[i] = i + IndexBase;

    M_repeatedEpetra_Map.reset( new Epetra_Map( -1, NumGlobalElements, &MyGlobalElements[0], IndexBase, *CommPtr) );
    M_uniqueEpetraMap.reset( new Epetra_Map( NumGlobalElements, IndexBase, *CommPtr) );
}

EpetraMap::EpetraMap(const int          size,
                     const comm_ptrtype& CommPtr):
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr(CommPtr)
{
    int NumGlobalElements(size);
    int NumMyElements    (NumGlobalElements);
    std::vector<int>  MyGlobalElements(size);
    int IndexBase = 1;

    for (int i(0); i < NumGlobalElements; ++i)
        MyGlobalElements[i] = i + 1;

    /*
    createMap( NumGlobalElements,
               NumMyElements,
               &MyGlobalElements[0],
               IndexBase,
               Comm);
    */

    M_repeatedEpetra_Map.reset( new Epetra_Map(NumGlobalElements,
                                               NumMyElements,
                                               &MyGlobalElements[0],
                                               IndexBase,
                                               *CommPtr) );

    if (CommPtr->MyPID() != 0) NumMyElements = 0;

    M_uniqueEpetraMap.reset( new Epetra_Map(NumGlobalElements,
                                            NumMyElements,
                                            &MyGlobalElements[0],
                                            IndexBase,
                                            *CommPtr) );



}

EpetraMap::EpetraMap( const map_type map ):
        M_repeatedEpetra_Map(new map_type(map)),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr()
{
    uniqueMap();
}

/*! Builds a submap of map _epetraMap with a given positive offset and
  the maximum id to consider (in the new count)
  eg: offset = 2, maxid = 6;
  _epetraMap = [ 0 2 5 7 8 10 1]
  this  =      [   0 3 5 7 ]

  if needed, indexBase may be changed (default values < 0 means "same as original map")
*/
EpetraMap::EpetraMap(const Epetra_BlockMap& _blockMap, const int offset, const int maxid,
                     int indexbase)
        :
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr()
{

    if (indexbase < 0) indexbase = _blockMap.IndexBase();

    std::vector<int> MyGlobalElements;
    int* sourceGlobalElements(_blockMap.MyGlobalElements());
    int const startIdOrig(offset + indexbase );
    int const endIdOrig  (startIdOrig + maxid );
    const int maxMyElements = std::min(maxid, _blockMap.NumMyElements());
    MyGlobalElements.reserve(maxMyElements);

    //Sort MyGlobalElements to avoid a bug in Trilinos (9?) when multiplying two matrices (A * B^T)
    std::sort (MyGlobalElements.begin(), MyGlobalElements.end());

    // We consider that the source Map may not be ordered
    for (int i(0); i < _blockMap.NumMyElements(); ++i)
        if (sourceGlobalElements[i] < endIdOrig && sourceGlobalElements[i] >= startIdOrig)
            MyGlobalElements.push_back(sourceGlobalElements[i] - offset);

    createMap( -1,
               MyGlobalElements.size(),
               &MyGlobalElements.front(),
               indexbase,
               _blockMap.Comm() );


}

// ===================================================
// Operators
// ===================================================
EpetraMap &
EpetraMap::operator = (const EpetraMap& _epetraMap)
{

    if (this != &_epetraMap)
    {
        M_repeatedEpetra_Map = _epetraMap.M_repeatedEpetra_Map;
        M_uniqueEpetraMap    = _epetraMap.M_uniqueEpetraMap;
        M_exporter           = _epetraMap.M_exporter;
        M_importer           = _epetraMap.M_importer;
        M_commPtr            = _epetraMap.M_commPtr;
    }

    return *this;
}


EpetraMap &
EpetraMap::operator += (const EpetraMap& _epetraMap)
{
    if (! _epetraMap.getUniqueMap())
        return *this;

    if (! this->getUniqueMap())
    {
        this->operator = (_epetraMap);
        return *this;
    }

    int*             pointer;
    std::vector<int> map;

    pointer = getRepeatedMap()->MyGlobalElements();
    for (int ii = 0; ii < getRepeatedMap()->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back(*pointer);
    }

    int numGlobalElements = getUniqueMap()->NumGlobalElements()
                            + getUniqueMap()->IndexBase() - _epetraMap.getUniqueMap()->IndexBase();

//    std::cout << "NumGlobalElements = " << numGlobalElements << std::endl;

    pointer = _epetraMap.getRepeatedMap()->MyGlobalElements();
    for (int ii = 0; ii < _epetraMap.getRepeatedMap()->NumMyElements(); ++ii, ++pointer)
    {
//        std::cout << "pointer = " << *pointer << std::endl;
        map.push_back(*pointer + numGlobalElements);
    }

    int IndexBase = getRepeatedMap()->IndexBase();

    M_repeatedEpetra_Map.reset( new Epetra_Map(-1, map.size(), &map[0], IndexBase, _epetraMap.getRepeatedMap()->Comm()) );

    map.resize(0);
    pointer = getUniqueMap()->MyGlobalElements();

    for (int ii = 0; ii < getUniqueMap()->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back(*pointer);
    }

    pointer = _epetraMap.getUniqueMap()->MyGlobalElements();
    for (int ii = 0; ii < _epetraMap.getUniqueMap()->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back(*pointer + numGlobalElements);
    }

    M_uniqueEpetraMap.reset( new Epetra_Map(-1, map.size(), &map[0], IndexBase, _epetraMap.getRepeatedMap()->Comm()) );

    M_exporter.reset();
    M_importer.reset();

    return *this;
}

/*
EpetraMap &
EpetraMap::operator += (std::vector<int> const& lagrangeMultipliers)
{
    if ( lagrangeMultipliers.size() <= 0)
        return *this;

    ASSERT(this->getUniqueMap(), "operator+=(const int) works only for an existing EpetraMap");

    EpetraMap  lagrMap(lagrangeMultipliers, Comm());

    this->operator+=(lagrMap);

    return *this;
}
*/

EpetraMap &
EpetraMap::operator += (int const size)
{

    EpetraMap  lagrMap(size, CommPtr());

    ASSERT(this->getUniqueMap(), "operator+=(const int) works only for an existing EpetraMap");

//     if ( size <= 0)
//         return *this;

    this->operator+=(lagrMap);
    return *this;
}


// ===================================================
// Methods
// ===================================================
boost::shared_ptr<EpetraMap>
EpetraMap::createRootMap(int const root)   const
{

    boost::shared_ptr<EpetraMap> rootMap(new EpetraMap(Epetra_Util::Create_Root_Map(*getUniqueMap(), root) ));
    return rootMap;
}

bool
EpetraMap::MapsAreSimilar( EpetraMap const& _epetraMap) const
{

    if ( this == &_epetraMap )
        return true;

    return( getUniqueMap()->SameAs( *_epetraMap.getUniqueMap()) &&
            getRepeatedMap()->SameAs( *_epetraMap.getRepeatedMap()) );


}

// ===================================================
// Get Methods
// ===================================================
EpetraMap::map_ptrtype const &
EpetraMap::getMap( EpetraMapType maptype)   const
{
    switch (maptype)
    {
    case Unique:
        return getUniqueMap();
    case Repeated:
        return getRepeatedMap();
    }
    return getUniqueMap();
}

Epetra_Export const&
EpetraMap::getExporter()
{
    createImportExport();
    return **M_exporter;
}

Epetra_Import const&
EpetraMap::getImporter()
{
    createImportExport();
    return **M_importer;
}


// ===================================================
// Private Methods
// ===================================================
EpetraMap::EpetraMap(const EpetraMap& _epetraMap):
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer()
{
    this->operator=(_epetraMap);
}


void
EpetraMap::createMap(int  NumGlobalElements,
                     int  NumMyElements,
                     int* MyGlobalElements,
                     int  IndexBase,
                     const comm_type& Comm)
{

    if (NumMyElements !=0 && MyGlobalElements == 0) // linearMap
        M_repeatedEpetra_Map.reset( new Epetra_Map(NumGlobalElements,
                                                   NumMyElements,
                                                   IndexBase,
                                                   Comm) );
    else // classic LifeV map
        M_repeatedEpetra_Map.reset( new Epetra_Map(NumGlobalElements,
                                                   NumMyElements,
                                                   MyGlobalElements,
                                                   IndexBase,
                                                   Comm) );

    uniqueMap();
}


void
EpetraMap::uniqueMap()
{
    M_uniqueEpetraMap.reset( new Epetra_Map( Epetra_Util::Create_OneToOne_Map (*getRepeatedMap(), false) ) );
    M_exporter.reset();
    M_importer.reset();
    return;
}

void
EpetraMap::createImportExport()
{

    if ( !getRepeatedMap() || !getUniqueMap() ) return;

    // The exporter is needed to import to a repeated vector
    if ( M_exporter.get() == 0 )
        M_exporter.reset( new boost::shared_ptr<Epetra_Export> );

    if ( M_exporter->get() == 0 )
        M_exporter->reset( new Epetra_Export(*getRepeatedMap(), *getUniqueMap()) );

    if ( M_importer.get() == 0 )
        M_importer.reset( new boost::shared_ptr<Epetra_Import> );

    if ( M_importer->get() == 0 )
        M_importer->reset( new Epetra_Import(*getRepeatedMap(), *getUniqueMap()) );

}

void
EpetraMap::bubbleSort(Epetra_IntSerialDenseVector& Elements)
{
    int hold;

    for (int pass(0); pass < Elements.Length()-1; pass++)
        for (int j(0); j < Elements.Length()-1; j++)
            if (Elements[j] > Elements[j+1])
            {
                hold      = Elements[j];
                Elements[j]  = Elements[j+1];
                Elements[j+1]= hold;
            }
}

void
EpetraMap::setUp(const RefFE&               refFE,
                 const comm_ptrtype&       _commPtr,
                 std::vector<int>& repeatedNodeVector,
                 std::vector<int>& repeatedEdgeVector,
                 std::vector<int>& repeatedFaceVector,
                 std::vector<int>& repeatedVolumeVector)
{
    int indexBase = 1;

    if (refFE.nbDofPerVertex())
    {
        int numNode = repeatedNodeVector.size();
        EpetraMap repeatedNodeMap(-1, numNode, &repeatedNodeVector[0], indexBase, _commPtr);
        operator+=(repeatedNodeMap);
    }

    if (refFE.nbDofPerEdge())
    {
        int numEdge = repeatedEdgeVector.size();
        EpetraMap repeatedEdgeMap(-1, numEdge, &repeatedEdgeVector[0], indexBase, _commPtr);
        operator+=(repeatedEdgeMap);
    }

    if (refFE.nbDofPerFace())
    {
        int numFace = repeatedFaceVector.size();
        EpetraMap repeatedFaceMap(-1, numFace, &repeatedFaceVector[0], indexBase, _commPtr);
        operator+=(repeatedFaceMap);
    }

    if (refFE.nbDofPerVolume())
    {
        int numElem = repeatedVolumeVector.size();
        EpetraMap repeatedElemMap(-1, numElem, &repeatedVolumeVector[0], indexBase, _commPtr);
        operator+=(repeatedElemMap);
    }

    createImportExport();
}

} // end namespace LifeV

