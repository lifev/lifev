/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file EpetraMap.cpp
*/

#include <life/lifealg/EpetraMap.hpp>
#include <Epetra_Util.h>

namespace LifeV
{
EpetraMap::EpetraMap():
    M_repeatedEpetra_Map(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer()
{}



EpetraMap::EpetraMap(int                NumGlobalElements,
                     int                NumMyElements,
                     int*               MyGlobalElements,
                     int                IndexBase,
                     const Epetra_Comm& Comm):
    M_repeatedEpetra_Map(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer()
{
    createMap( NumGlobalElements,
               NumMyElements,
               MyGlobalElements,
               IndexBase,
               Comm);
}

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


EpetraMap::EpetraMap(const int          size,
                     const Epetra_Comm& Comm):
    M_repeatedEpetra_Map(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer()
{
    int NumGlobalElements(size);
    int NumMyElements    (NumGlobalElements);
    std::vector<int>  MyGlobalElements(size);
    int IndexBase = 1;

    for (int i(0); i < NumGlobalElements; ++i)
        MyGlobalElements[i] = i + 1;


    createMap( NumGlobalElements,
               NumMyElements,
               &MyGlobalElements[0],
               IndexBase,
               Comm);

}



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


EpetraMap::map_type
EpetraMap::getRootMap( int root)   const
{
    return Epetra_Util::Create_Root_Map(*getUniqueMap(), root);
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
    M_importer()
{

    if (indexbase < 0) indexbase = _blockMap.IndexBase();

    std::vector<int> MyGlobalElements;
    int* sourceGlobalElements(_blockMap.MyGlobalElements());
    int const startIdOrig(offset + indexbase );
    int const endIdOrig  (startIdOrig + maxid );
    const int maxMyElements = std::min(maxid, _blockMap.NumMyElements());
    MyGlobalElements.reserve(maxMyElements);



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


EpetraMap &
EpetraMap::operator = (const EpetraMap& _epetraMap)
{

   if (this != &_epetraMap)
   {
       M_repeatedEpetra_Map = _epetraMap.M_repeatedEpetra_Map;
       M_uniqueEpetraMap    = _epetraMap.M_uniqueEpetraMap;
       M_exporter           = _epetraMap.M_exporter;
       M_importer           = _epetraMap.M_importer;
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


//


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


//

EpetraMap &
EpetraMap::operator += (int const size)
{

    EpetraMap  lagrMap(size, Comm());

    ASSERT(this->getUniqueMap(), "operator+=(const int) works only for an existing EpetraMap");

//     if ( size <= 0)
//         return *this;

    this->operator+=(lagrMap);
    return *this;
}

//


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
                     const Epetra_Comm &Comm)
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


  int MyPID    ( getRepeatedMap()->Comm().MyPID() );
  int NumIDs   ( getRepeatedMap()->NumMyElements() );
  int indexBase( getRepeatedMap()->MinAllGID() );

  Epetra_IntSerialDenseVector GIDList  (NumIDs);
  getRepeatedMap()->MyGlobalElements(GIDList.Values());

  Epetra_IntSerialDenseVector PIDList  (NumIDs);
  Epetra_IntSerialDenseVector LIDList  (NumIDs);

  getRepeatedMap()->RemoteIDList(NumIDs, GIDList.Values(), PIDList.Values(), LIDList.Values());

// now use LIDList has a helping pointer

  int MyUniqueElements(0);

  for (int i(0); i< NumIDs; i++)
    if (PIDList[i] == MyPID)
      LIDList[MyUniqueElements++] = GIDList[i];

  LIDList.Resize(MyUniqueElements);

  bubbleSort(LIDList);

  M_uniqueEpetraMap.reset( new Epetra_Map(-1,
                                          LIDList.Length(),
                                          LIDList.Values(),
                                          indexBase,
                                          getRepeatedMap()->Comm()) );
  M_exporter.reset();
  M_importer.reset();


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
                 Epetra_Comm&               _comm,
                 std::vector<int>& repeatedNodeVector,
                 std::vector<int>& repeatedEdgeVector,
                 std::vector<int>& repeatedFaceVector,
                 std::vector<int>& repeatedVolumeVector)
{
    int indexBase = 1;

    if (refFE.nbDofPerVertex)
    {
        int numNode = repeatedNodeVector.size();
        EpetraMap repeatedNodeMap(-1, numNode, &repeatedNodeVector[0], indexBase,  _comm);
        operator+=(repeatedNodeMap);
    }

    if (refFE.nbDofPerEdge)
    {
        int numEdge = repeatedEdgeVector.size();
        EpetraMap repeatedEdgeMap(-1, numEdge, &repeatedEdgeVector[0], indexBase,  _comm);
        operator+=(repeatedEdgeMap);
    }

    if (refFE.nbDofPerFace)
    {
    	int numFace = repeatedFaceVector.size();
        EpetraMap repeatedFaceMap(-1, numFace, &repeatedFaceVector[0], indexBase,  _comm);
        operator+=(repeatedFaceMap);
    }

    if (refFE.nbDofPerVolume)
    {
    	int numElem = repeatedVolumeVector.size();
        EpetraMap repeatedElemMap(-1, numElem, &repeatedVolumeVector[0], indexBase,  _comm);
        operator+=(repeatedElemMap);
    }

    createImportExport();
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

bool
EpetraMap::MapsAreSimilar( EpetraMap const& _epetraMap) const
{

    if ( this == &_epetraMap )
        return true;

    return( getUniqueMap()->SameAs( *_epetraMap.getUniqueMap()) &&
            getRepeatedMap()->SameAs( *_epetraMap.getRepeatedMap()) );


}


} // end namespace LifeV

