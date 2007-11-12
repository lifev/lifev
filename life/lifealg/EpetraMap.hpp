/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  This file is part of the LifeV library

  Author(s): Gilles Fourestey gilles.fourestey@epfl.ch
       Date: 2004-10-26

  Copyright (C) 2004 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


*/
/** \file EpetraMap.hpp
*/


#ifndef _EPETRAMAP_
#define _EPETRAMAP_

#include <life/lifefem/refFE.hpp>

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "life/lifecore/life.hpp"
#include "life/lifemesh/partitionMesh.hpp"



namespace LifeV
{
////////////////////////////////////////////////////////////////
//
//  Epetra Matrix format Wrapper
//
///////////////////////////////////////////////////////////////


class EpetraMap
{
public:

    EpetraMap();
    EpetraMap(int                NumGlobalElements,
              int                NumMyElements,
              int*               MyGlobalElements,
              int                IndexBase,
              const Epetra_Comm& Comm);

    template<typename Mesh>
    EpetraMap(const RefFE&         refFE,
              const partitionMesh<Mesh>& meshPart,
              Epetra_Comm&         _comm);

    template<typename Mesh>
    EpetraMap(const RefFE&         refFE,
              const Mesh&          mesh,
              Epetra_Comm&         _comm);



    EpetraMap(const EpetraMap& _epetraMap);

    ~EpetraMap();

    EpetraMap&         operator  = (const EpetraMap& _epetraMap);
    EpetraMap&         operator += (const EpetraMap& _epetraMap);
    EpetraMap          operator +  (const EpetraMap& _epetraMap)
        {
            EpetraMap map( *this );
            return map += _epetraMap;
        }

    Epetra_Map const * getEpetra_Map()       const
        {return M_epetraMap;}
    Epetra_Map const * getRepeatedEpetra_Map() const
        {return M_epetraMap;}
    Epetra_Map const * getUniqueEpetra_Map()   const
        {return M_uniqueEpetraMap;}

    Epetra_Comm const& Comm() const { return M_uniqueEpetraMap->Comm(); }

//    Epetra_Map*        getEpetra_Map(){return M_epetraMap;}

    void               createMap(int   NumGlobalElements,
                                 int   NumMyElements,
                                 int*  MyGlobalElements,
                                 int   IndexBase,
                                 const Epetra_Comm &Comm)  ;

//    EpetraMap&          uniqueMap();

private:

    void               uniqueMap();
    void               bubbleSort(Epetra_IntSerialDenseVector& Elements);

    void setUp(const RefFE&               refFE,
               Epetra_Comm&               _comm,
               std::vector<int>& repeatedNodeVector,
               std::vector<int>& repeatedEdgeVector,
               std::vector<int>& repeatedFaceVector,
               std::vector<int>& repeatedElemVector);

    Epetra_Map*        M_epetraMap;
    Epetra_Map*        M_uniqueEpetraMap;

};

template<typename Mesh>
EpetraMap::
EpetraMap(const RefFE&               refFE,
          const partitionMesh<Mesh>& meshPart,
          Epetra_Comm&               _comm):
    M_epetraMap(0),
    M_uniqueEpetraMap(0)
{

    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

    setUp( refFE,
           _comm,
           const_cast<std::vector<int>&>(meshPart.repeatedNodeVector()),
           const_cast<std::vector<int>&>(meshPart.repeatedEdgeVector()),
           const_cast<std::vector<int>&>(meshPart.repeatedFaceVector()),
           const_cast<std::vector<int>&>(meshPart.repeatedElemVector()) );

}


template<typename Mesh>
EpetraMap::
EpetraMap(const RefFE&               refFE,
          const Mesh&                mesh,
          Epetra_Comm&               _comm):
    M_epetraMap(0),
    M_uniqueEpetraMap(0)
{

    std::vector<int> repeatedNodeVector;
    std::vector<int> repeatedEdgeVector;
    std::vector<int> repeatedFaceVector;
    std::vector<int> repeatedElemVector;

    if (refFE.nbDofPerVertex)
    {
        repeatedNodeVector.reserve(mesh.numPoints());
        for ( UInt ii = 1; ii <= mesh.numPoints(); ii++ )
            repeatedNodeVector.push_back(mesh.pointList( ii ).id());
    }

    if (refFE.nbDofPerEdge)
    {
        repeatedEdgeVector.reserve(mesh.numEdges());

        for ( UInt ii = 1; ii <= mesh.numEdges(); ii++ )
            repeatedEdgeVector.push_back(mesh.edgeList( ii ).id());
    }

    if (refFE.nbDofPerFace)
    {
        repeatedFaceVector.reserve(mesh.numFaces());

        for ( UInt ii = 1; ii <= mesh.numFaces(); ii++ )
            repeatedFaceVector.push_back(mesh.faceList( ii ).id());
    }

    if (refFE.nbDofPerVolume)
    {
        repeatedElemVector.reserve(mesh.numVolumes());

        for ( UInt ii = 1; ii <= mesh.numVolumes(); ii++ )
            repeatedElemVector.push_back(mesh.volumeList( ii ).id());
    }


    setUp( refFE,
           _comm,
           repeatedNodeVector,
           repeatedEdgeVector,
           repeatedFaceVector,
           repeatedElemVector );

    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

}


}


#endif

