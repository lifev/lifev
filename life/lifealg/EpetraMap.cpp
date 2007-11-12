/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file EpetraMap.cpp
*/

#include "EpetraMap.hpp"

namespace LifeV
{
EpetraMap::EpetraMap():
    M_epetraMap(0),
    M_uniqueEpetraMap(0)
{}



EpetraMap::EpetraMap(int                NumGlobalElements,
                     int                NumMyElements,
                     int*               MyGlobalElements,
                     int                IndexBase,
                     const Epetra_Comm& Comm):
    M_epetraMap(0),
    M_uniqueEpetraMap(0)
{
    createMap( NumGlobalElements,
               NumMyElements,
               MyGlobalElements,
               IndexBase,
               Comm);
}

EpetraMap::~EpetraMap()
{
    delete M_epetraMap;
    delete M_uniqueEpetraMap;
}


EpetraMap &
EpetraMap::operator = (const EpetraMap& _epetraMap)
{

   if (this != &_epetraMap)
   {
       if (!_epetraMap.getEpetra_Map() == 0)
       {
           if (M_epetraMap != 0)
           {
               *M_epetraMap = *_epetraMap.getEpetra_Map();
               *M_uniqueEpetraMap = *_epetraMap.getUniqueEpetra_Map();
           }
           else
            {
                M_epetraMap = new Epetra_Map(*_epetraMap.getEpetra_Map());
                M_uniqueEpetraMap = new Epetra_Map(*_epetraMap.getUniqueEpetra_Map());
            }
       }
   }

   return *this;
}


EpetraMap &
EpetraMap::operator += (const EpetraMap& _epetraMap)
{
    if (this->getEpetra_Map() == 0)
    {
        this->operator = (_epetraMap);
        return *this;
    }

    if (_epetraMap.getEpetra_Map() == 0)
        return *this;

    int*             pointer;
    std::vector<int> map;

    pointer = M_epetraMap->MyGlobalElements();
    for (int ii = 0; ii < M_epetraMap->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back(*pointer);
    }

    int numGlobalElements = getUniqueEpetra_Map()->NumGlobalElements();

//    std::cout << "NumGlobalElements = " << numGlobalElements << std::endl;

    pointer = _epetraMap.getEpetra_Map()->MyGlobalElements();
    for (int ii = 0; ii < _epetraMap.getEpetra_Map()->NumMyElements(); ++ii, ++pointer)
    {
//        std::cout << "pointer = " << *pointer << std::endl;
        map.push_back(*pointer + numGlobalElements);
    }

    int IndexBase = M_epetraMap->IndexBase();
    delete M_epetraMap;

    M_epetraMap       = new Epetra_Map(-1, map.size(), &map[0], IndexBase, _epetraMap.getEpetra_Map()->Comm());

    map.resize(0);
    pointer = M_uniqueEpetraMap->MyGlobalElements();

    for (int ii = 0; ii < M_uniqueEpetraMap->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back(*pointer);
    }

    pointer = _epetraMap.getUniqueEpetra_Map()->MyGlobalElements();
    for (int ii = 0; ii < _epetraMap.getUniqueEpetra_Map()->NumMyElements(); ++ii, ++pointer)
    {
        map.push_back(*pointer + numGlobalElements);
    }

    delete M_uniqueEpetraMap;

    M_uniqueEpetraMap       = new Epetra_Map(-1, map.size(), &map[0], IndexBase, _epetraMap.getEpetra_Map()->Comm());

//    uniqueMap();

    return *this;
}


EpetraMap::EpetraMap(const EpetraMap& _epetraMap):
    M_epetraMap      (0),
    M_uniqueEpetraMap(0)
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
    delete M_epetraMap;
    delete M_uniqueEpetraMap;

    M_epetraMap = new Epetra_Map(NumGlobalElements,
                                 NumMyElements,
                                 MyGlobalElements,
                                 IndexBase,
                                 Comm);
    uniqueMap();
}


void
EpetraMap::uniqueMap()
{
  int MyPID    ( M_epetraMap->Comm().MyPID() );
  int NumIDs   ( M_epetraMap->NumMyElements() );
  int indexBase( M_epetraMap->MinAllGID() );

  Epetra_IntSerialDenseVector GIDList  (NumIDs);
  M_epetraMap->MyGlobalElements(GIDList.Values());

  Epetra_IntSerialDenseVector PIDList  (NumIDs);
  Epetra_IntSerialDenseVector LIDList  (NumIDs);

  M_epetraMap->RemoteIDList(NumIDs, GIDList.Values(), PIDList.Values(), LIDList.Values());

// now use LIDList has a helping pointer

  int MyUniqueElements(0);

  for (int i(0); i< NumIDs; i++)
    if (PIDList[i] == MyPID)
      LIDList[MyUniqueElements++] = GIDList[i];

  LIDList.Resize(MyUniqueElements);

  bubbleSort(LIDList);

  delete M_uniqueEpetraMap;
  M_uniqueEpetraMap = new Epetra_Map(-1,
                                     LIDList.Length(),
                                     LIDList.Values(),
                                     indexBase,
                                     M_epetraMap->Comm());

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
                 std::vector<int>& repeatedElemVector)
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
        int numElem = repeatedElemVector.size();
        EpetraMap repeatedElemMap(-1, numElem, &repeatedElemVector[0], indexBase,  _comm);
        operator+=(repeatedElemMap);
    }
}

}
