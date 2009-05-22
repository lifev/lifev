/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  This file is part of the LifeV library

  Author(s): Gilles Fourestey gilles.fourestey@epfl.ch
             Simone Deparis   simone.deparis@epfl.ch
       Date: 2006-10-26

  Copyright (C) 2006 EPFL

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


#include <boost/shared_ptr.hpp>

#include <life/lifefem/refFE.hpp>

#include <Epetra_Map.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Comm.h>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/partitionMesh.hpp>



namespace LifeV
{
////////////////////////////////////////////////////////////////
//
//  Epetra Matrix format Wrapper
//
///////////////////////////////////////////////////////////////

enum EpetraMapType {Unique = 0, Repeated};

class EpetraMap
{
public:


    /** @name Typedefs
     */
    //@{
    typedef Epetra_Map map_type;
    typedef boost::shared_ptr<map_type> map_ptrtype;
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Export> > exporter_ptrtype;
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Import> > importer_ptrtype;
    //@}

    EpetraMap();
    // epetra map constructor. To define a linear map, set MyGlobalElements = 0
    EpetraMap(int                NumGlobalElements,
              int                NumMyElements,
              int*               MyGlobalElements,
              int                IndexBase,
              const Epetra_Comm& Comm);

    //! construct a map with entries lagrangeMultipliers.
    //! Repeated lagrange multipliers will be repeated on the repeatedmap
    //! Again: it is not necessary that the lagrangeMltiplier vector is the same on all
    //!       processors nor that it is different
    /*
    EpetraMap(std::vector<int> const& lagrangeMultipliers,
              const Epetra_Comm&      Comm);
    */

    //!
    EpetraMap(const int               size,
              const Epetra_Comm&      Comm);

    // Calls createImportExport from setUp()
    template<typename Mesh>
    EpetraMap(const RefFE&         refFE,
              const partitionMesh<Mesh>& meshPart,
              Epetra_Comm&         _comm);

    // Calls createImportExport from setUp()
    template<typename Mesh>
    EpetraMap(const RefFE&         refFE,
              const Mesh&          mesh,
              Epetra_Comm&         _comm);

    EpetraMap(const EpetraMap& _epetraMap);

    /*! Builds a submap of map _epetraMap with a given positive offset and
      the maximum id to consider
      eg: offset = 2, maxid = 6;
      _epetraMap = [ 0 2 5 7 8 10 1]
      this  =      [   0 3 5 7 ]

      if needed, indexBase may be changed (default values < 0 means "same as original map")
    */
    EpetraMap(const Epetra_BlockMap& _blockMap, const int offset, const int maxid,
              int indexbase = -1);

    ~EpetraMap() {}

    // The copy operator will copy the pointers of the maps, exporter and importer
    EpetraMap&         operator  = (const EpetraMap& _epetraMap);

    EpetraMap&         operator += (const EpetraMap& _epetraMap);
    EpetraMap          operator +  (const EpetraMap& _epetraMap)
        {
            EpetraMap map( *this );
            map += _epetraMap;
            createImportExport();
            return map;
        }

    /*
    EpetraMap&         operator += (std::vector<int> const&   lagrangeMultipliers);
    EpetraMap          operator +  (std::vector<int> const&   lagrangeMultipliers)
        {
            EpetraMap map( *this );
            map += lagrangeMultipliers;
            createImportExport();
            return map;
        }
    */

    EpetraMap&         operator += (int const size);
    EpetraMap          operator +  (int const size)
        {
            //int me =  M_uniqueEpetraMap->Comm().MyPID();
            EpetraMap map( *this );
            map += size;
            createImportExport();
            return map;
        }





    Epetra_Comm const& Comm() const { return M_uniqueEpetraMap->Comm(); }

//    Epetra_Map*        getRepeatedEpetra_Map(){return M_repeatedEpetra_Map;}

    // createMap does not call createImportExport
    void               createMap(int   NumGlobalElements,
                                 int   NumMyElements,
                                 int*  MyGlobalElements,
                                 int   IndexBase,
                                 const Epetra_Comm &Comm)  ;

    map_ptrtype const & getMap( EpetraMapType maptype)   const;
    map_type            getRootMap( int root)   const;

    Epetra_Export const& getExporter();
    Epetra_Import const& getImporter();


    bool MapsAreSimilar( EpetraMap const& _epetraMap) const;


//    EpetraMap&          uniqueMap();

private:

    /*
    Epetra_Map const * getRepeatedEpetra_Map()       const
        {return M_repeatedEpetra_Map;}
    Epetra_Map const * getUniqueEpetra_Map()   const
        {return M_uniqueEpetraMap;}
    */
    map_ptrtype const & getRepeatedMap() const { return M_repeatedEpetra_Map; }
    map_ptrtype const & getUniqueMap()   const { return M_uniqueEpetraMap; }


    void  uniqueMap();
    void  createImportExport();
    void  bubbleSort(Epetra_IntSerialDenseVector& Elements);

    // Calls createImportExport
    void setUp(const RefFE&               refFE,
               Epetra_Comm&               _comm,
               std::vector<int>& repeatedNodeVector,
               std::vector<int>& repeatedEdgeVector,
               std::vector<int>& repeatedFaceVector,
               std::vector<int>& repeatedVolumeVector);

    map_ptrtype        M_repeatedEpetra_Map;
    map_ptrtype        M_uniqueEpetraMap;
    exporter_ptrtype   M_exporter;
    importer_ptrtype   M_importer;


};

template<typename Mesh>
EpetraMap::
EpetraMap(const RefFE&               refFE,
          const partitionMesh<Mesh>& meshPart,
          Epetra_Comm&               _comm):
    M_repeatedEpetra_Map(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer()
{

    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

    setUp( refFE,
           _comm,
           const_cast<std::vector<int>&>(meshPart.repeatedNodeVector()),
           const_cast<std::vector<int>&>(meshPart.repeatedEdgeVector()),
           const_cast<std::vector<int>&>(meshPart.repeatedFaceVector()),
           const_cast<std::vector<int>&>(meshPart.repeatedVolumeVector()) );

}


template<typename Mesh>
EpetraMap::
EpetraMap(const RefFE&               refFE,
          const Mesh&                mesh,
          Epetra_Comm&               _comm):
    M_repeatedEpetra_Map(),
    M_uniqueEpetraMap(),
    M_exporter(),
    M_importer()
{

    std::vector<int> repeatedNodeVector;
    std::vector<int> repeatedEdgeVector;
    std::vector<int> repeatedFaceVector;
    std::vector<int> repeatedVolumeVector;

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
        repeatedVolumeVector.reserve(mesh.numVolumes());

        for ( UInt ii = 1; ii <= mesh.numVolumes(); ii++ )
            repeatedVolumeVector.push_back(mesh.volumeList( ii ).id());
    }


    setUp( refFE,
           _comm,
           repeatedNodeVector,
           repeatedEdgeVector,
           repeatedFaceVector,
           repeatedVolumeVector );

    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

}


}


#endif

