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

    //! @name Public Types
    //@{

    typedef Epetra_Map map_type;
    typedef boost::shared_ptr<map_type> map_ptrtype;
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Export> > exporter_ptrtype;
    typedef boost::shared_ptr< boost::shared_ptr<Epetra_Import> > importer_ptrtype;
    typedef Epetra_Comm comm_type;
    typedef boost::shared_ptr<comm_type> comm_ptrtype;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Default Constructor
    EpetraMap();
    // epetra map constructor. To define a linear map, set MyGlobalElements = 0
    EpetraMap(Int                NumGlobalElements,
              Int                NumMyElements,
              Int*               MyGlobalElements,
              Int                IndexBase,
              const comm_ptrtype&  CommPtr);

    //! construct a map with entries lagrangeMultipliers.
    //! Repeated lagrange multipliers will be repeated on the repeatedmap
    //! Again: it is not necessary that the lagrangeMltiplier vector is the same on all
    //!       processors nor that it is different
    /*
    EpetraMap(std::vector<Int> const& lagrangeMultipliers,
              comm_ptrtype&      CommPtr);
    */

    //! Build a nearly equally distributed map.
    /*!
     *  The map is equally distributed if NumGlobalElements % Rank = 0
     *  @param NumGlobalElements - Total number of elements inside the map
     *  @param indexBase - Starting index base (typically 0 or 1)
     *  @param CommPtr - a pointer to the Epetra communicator
     */
    EpetraMap( const Int NumGlobalElements, const Int IndexBase, const comm_ptrtype& CommPtr );

    EpetraMap(const Int               size,
              const comm_ptrtype&     CommPtr);

    // Calls createImportExport from setUp()
    template<typename Mesh>
    EpetraMap(const RefFE&               refFE,
              const partitionMesh<Mesh>& meshPart,
              const comm_ptrtype&       _commPtr);

    // Calls createImportExport from setUp()
    template<typename Mesh>
    EpetraMap(const RefFE&         refFE,
              const Mesh&          mesh,
              const comm_ptrtype& _commPtr);

    EpetraMap(const EpetraMap& _epetraMap);

    /*! Builds a submap of map _epetraMap with a given positive offset and
      the maximum id to consider
      eg: offset = 2, maxid = 6;
      _epetraMap = [ 0 2 5 7 8 10 1]
      this  =      [   0 3 5 6 ]

      if needed, indexBase may be changed (default values < 0 means "same as original map")
    */
    EpetraMap(const Epetra_BlockMap& _blockMap, const Int offset, const Int maxid,
              Int indexbase = -1);

    //! Constructor from raw Epetra_Map. This constructor should be used only inside this class,
    //! therefore it is private
    /*!
     * \param map: underlying Epetra_Map
     */
private:
    EpetraMap(const map_type map);

public:
    ~EpetraMap() {}

    //@}


    //! @name Operators
    //@{

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
    EpetraMap&         operator += (std::vector<Int> const&   lagrangeMultipliers);
    EpetraMap          operator +  (std::vector<Int> const&   lagrangeMultipliers)
        {
            EpetraMap map( *this );
            map += lagrangeMultipliers;
            createImportExport();
            return map;
        }
    */

    EpetraMap&         operator += (Int const size);
    EpetraMap          operator +  (Int const size)
    {
        //Int me =  M_uniqueEpetraMap->Comm().MyPID();
        EpetraMap map( *this );
        map += size;
        createImportExport();
        return map;
    }

    //@}

    //! @name Methods
    //@{

    //! This methods create a pointer to a EpetraMap that has points only on processor root
    /*!
     * \param this: EpetraMap that selects the relevant points
     * \param root: processor on which to export all the points
     */
    boost::shared_ptr<EpetraMap>         createRootMap( Int const     root)    const;

    bool MapsAreSimilar( EpetraMap const& _epetraMap) const;

    //EpetraMap&          uniqueMap();

    //@}

    //! @name Get Methods
    //@{

    //comm_type const& Comm() const { return M_uniqueEpetraMap->Comm(); }
    comm_type const& Comm() const { return *M_commPtr; }
    comm_ptrtype& CommPtr() { return M_commPtr; }

    //Epetra_Map*        getRepeatedEpetra_Map(){return M_repeatedEpetra_Map;}

    map_ptrtype const & getMap  ( EpetraMapType maptype) const;

    Epetra_Export const& getExporter();

    Epetra_Import const& getImporter();

    //@}

private:

    //! @name Private Methods
    //@{

    // createMap does not call createImportExport
    void               createMap(Int   NumGlobalElements,
                                 Int   NumMyElements,
                                 Int*  MyGlobalElements,
                                 Int   IndexBase,
                                 const comm_type& Comm)  ;

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
               const comm_ptrtype&       _commPtr,
               std::vector<Int>& repeatedNodeVector,
               std::vector<Int>& repeatedEdgeVector,
               std::vector<Int>& repeatedFaceVector,
               std::vector<Int>& repeatedVolumeVector);

    //@}

    map_ptrtype        M_repeatedEpetra_Map;
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
EpetraMap(const RefFE&               refFE,
          const partitionMesh<Mesh>& meshPart,
          const comm_ptrtype&        _commPtr):
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr(_commPtr)
{

    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

    setUp( refFE,
           _commPtr,
           const_cast<std::vector<Int>&>(meshPart.repeatedNodeVector()),
           const_cast<std::vector<Int>&>(meshPart.repeatedEdgeVector()),
           const_cast<std::vector<Int>&>(meshPart.repeatedFaceVector()),
           const_cast<std::vector<Int>&>(meshPart.repeatedVolumeVector()) );

}


template<typename Mesh>
EpetraMap::
EpetraMap(const RefFE&               refFE,
          const Mesh&                mesh,
          const comm_ptrtype&        _commPtr):
        M_repeatedEpetra_Map(),
        M_uniqueEpetraMap(),
        M_exporter(),
        M_importer(),
        M_commPtr(_commPtr)
{

    std::vector<Int> repeatedNodeVector;
    std::vector<Int> repeatedEdgeVector;
    std::vector<Int> repeatedFaceVector;
    std::vector<Int> repeatedVolumeVector;

    if (refFE.nbDofPerVertex())
    {
        repeatedNodeVector.reserve(mesh.numPoints());
        for ( UInt ii = 1; ii <= mesh.numPoints(); ii++ )
            repeatedNodeVector.push_back(mesh.pointList( ii ).id());
    }

    if (refFE.nbDofPerEdge())
    {
        repeatedEdgeVector.reserve(mesh.numEdges());

        for ( UInt ii = 1; ii <= mesh.numEdges(); ii++ )
            repeatedEdgeVector.push_back(mesh.edgeList( ii ).id());
    }

    if (refFE.nbDofPerFace())
    {
        repeatedFaceVector.reserve(mesh.numFaces());

        for ( UInt ii = 1; ii <= mesh.numFaces(); ii++ )
            repeatedFaceVector.push_back(mesh.faceList( ii ).id());
    }

    if (refFE.nbDofPerVolume())
    {
        repeatedVolumeVector.reserve(mesh.numVolumes());

        for ( UInt ii = 1; ii <= mesh.numVolumes(); ii++ )
            repeatedVolumeVector.push_back(mesh.volumeList( ii ).id());
    }


    setUp( refFE,
           _commPtr,
           repeatedNodeVector,
           repeatedEdgeVector,
           repeatedFaceVector,
           repeatedVolumeVector );

    // Epetra_Map is "badly" coded, in fact its constructor needs a non-constant pointer to indices, but it
    // never modify them

}


}


#endif

