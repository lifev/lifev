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
    @file MeshUtility.hpp
    @brief Functions to load meshes

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 16-01-2013
 */

#ifndef NSMESHUTILITY_HPP
#define NSMESHUTILITY_HPP

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>


namespace LifeV
{

namespace MeshUtility
{

//! setup and get the mesh data
/*!
  @param meshName name of the mesh file
  @param resourcesPath path to the mesh folder
  @param meshOrder order of the mesh elements
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
MeshData
getMeshData ( const std::string& meshName,
              const std::string& resourcesPath = "./",
              const std::string& meshOrder = "P1");

//! Print informations about the mesh
template< typename RegionMeshType>
void printMeshInfos ( boost::shared_ptr< RegionMeshType > mesh )
{
#ifdef HAVE_MPI
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif
    Displayer displayer ( Comm );
    MeshUtility::MeshStatistics::meshSize meshSize = MeshUtility::MeshStatistics::computeSize ( *mesh );
    displayer.leaderPrint ( "Mesh size (max): ", meshSize.maxH, "\n" );
    displayer.leaderPrint ( "Mesh size (min): ", meshSize.minH, "\n" );
    displayer.leaderPrint ( "Mesh size (av.): ", meshSize.meanH, "\n" );
}

//! Read and partitioned a *.mesh file
/*!
  @param meshLocal The partitioned mesh that we want to generate
  @param meshFull  The non partitioned mesh that we want to keep
  @param isPartitioned boolean to say if the mesh should be partitioned or just loaded
  @param meshName name of the mesh file
  @param resourcesPath path to the mesh folder
  @param meshOrder order of the mesh
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
template< typename RegionMeshType>
void fillWithMesh ( boost::shared_ptr< RegionMeshType >& meshLocal,
                    boost::shared_ptr< RegionMeshType >& meshFull ,
                    bool isPartitioned,
                    const std::string& meshName,
                    const std::string& resourcesPath = "./",
                    const std::string& meshOrder = "P1" )
{
    if (isPartitioned)
    {
        fillWithPartitionedMesh ( meshLocal, meshName, resourcesPath );
    }
    else
    {
        fillWithFullMesh ( meshLocal,  meshFull, meshName, resourcesPath, meshOrder );
    }
}

//! Read and partitioned a *.mesh file
/*!
  @param meshLocal The partitioned mesh that we want to generate
  @param isPartitioned boolean to say if the mesh should be partitioned or just loaded
  @param meshName name of the mesh file
  @param resourcesPath path to the mesh folder
  @param meshOrder order of the mesh
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
template< typename RegionMeshType>
void fillWithMesh ( boost::shared_ptr< RegionMeshType >& meshLocal,
                    bool isPartitioned,
                    const std::string& meshName,
                    const std::string& resourcesPath = "./",
                    const std::string& meshOrder = "P1" )
{
    boost::shared_ptr< RegionMeshType > tmpMeshFull ( new RegionMeshType );
    fillWithMesh ( meshLocal, tmpMeshFull, isPartitioned, meshName, resourcesPath, meshOrder );
}


//! Read and partitioned a *.mesh file
/*!
  @param meshLocal The partitioned mesh that we want to generate
  @param meshFull  The non partitioned mesh that we want to keep
  @param meshName name of the mesh file
  @param resourcesPath path to the mesh folder
  @param meshOrder order of the mesh
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
template< typename RegionMeshType>
void fillWithFullMesh (  boost::shared_ptr< RegionMeshType >& meshLocal,
                         boost::shared_ptr< RegionMeshType >& meshFull,
                         const std::string& meshName,
                         const std::string& resourcesPath = "./",
                         const std::string& meshOrder = "P1" )
{
#ifdef HAVE_MPI
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif
    Displayer displayer ( Comm );

    LifeChrono meshReadChrono;
    meshReadChrono.start();
    boost::shared_ptr<RegionMeshType > fullMesh ( new RegionMeshType );
    readMesh (*fullMesh, getMeshData (meshName, resourcesPath, meshOrder ) );
    printMeshInfos ( fullMesh );
    meshReadChrono.stop();
    displayer.leaderPrint ("Loading time: ", meshReadChrono.diff(), " s.\n");

    LifeChrono meshPartChrono;
    meshPartChrono.start();
    MeshPartitioner< RegionMeshType > meshPartitioner ( fullMesh, Comm );
    meshLocal = meshPartitioner.meshPartition();
    meshPartChrono.stop();
    displayer.leaderPrint ("Partitioning time: ", meshPartChrono.diff(), " s.\n");
    if ( meshFull )
    {
        meshFull = fullMesh;
    }
    else
    {
        fullMesh.reset();    //Freeing the global mesh to save memory
    }
}

//! Read and partitioned a *.mesh file
/*!
  @param meshLocal The partitioned mesh that we want to generate
  @param meshName name of the mesh file
  @param resourcesPath path to the mesh folder
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
template< typename RegionMeshType>
void fillWithPartitionedMesh ( boost::shared_ptr< RegionMeshType >& meshLocal,
                               const std::string& meshName,
                               const std::string& resourcesPath = "./" )
{
#ifdef HAVE_MPI
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif
    Displayer displayer ( Comm );

    LifeChrono meshReadChrono;
    meshReadChrono.start();
#ifdef LIFEV_HAS_HDF5
    PartitionIO< RegionMeshType > partitionIO ( ( resourcesPath + meshName ).data(), Comm);
    partitionIO.read (meshLocal);
#else
    ASSERT (false, "You must compile LifeV with HDF5 to load partitioned meshes");
#endif
    meshReadChrono.stop();
    displayer.leaderPrint ("Loading time: ", meshReadChrono.diff(), " s.\n");
}

//! Build a mesh from a partitioned mesh
/*!
  @param mesh The mesh that we want to generate
  @param regionFlag Flag of the region
  @param m Number of elements along the ( length, width, height )
  @param l length of the mesh ( length, width, height )
  @param t translation of the mesh along the (x,y,z)-axis
  @param verbose Verbose mode enabled/disabled
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
template< typename RegionMeshType>
void fillWithStructuredMesh ( boost::shared_ptr< RegionMeshType >& mesh,
                              boost::shared_ptr< RegionMeshType >& meshFull,
                              markerID_Type regionFlag,
                              const std::vector<UInt>& m,
                              bool verbose = false,
                              const std::vector<Real>& l = std::vector<Real> (3, 1),
                              const std::vector<Real>& t = std::vector<Real> (3, 0) )
{
#ifdef HAVE_MPI
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif
    Displayer displayer ( Comm );

    LifeChrono meshBuildChrono;
    meshBuildChrono.start();
    boost::shared_ptr< RegionMeshType > fullMesh ( new RegionMeshType ( Comm ) );
    if ( m.size() == 1)
    {
        //TODO structured mesh in 1D
    }
    else if ( m.size() == 2)
    {
        //TODO structured mesh in 2D
    }
    else
    {
        regularMesh3D ( *fullMesh,
                        regionFlag,
                        m[0], m[1], m[2],
                        verbose,
                        l[0], l[1], l[2],
                        t[0], t[1], t[2] );
    }
    printMeshInfos ( fullMesh );
    meshBuildChrono.stop();
    displayer.leaderPrint ("Building time: ", meshBuildChrono.diff(), " s.\n");

    LifeChrono meshPartChrono;
    meshPartChrono.start();
    MeshPartitioner<  RegionMeshType  > meshPartitioner ( fullMesh, Comm );
    mesh = meshPartitioner.meshPartition();
    meshPartChrono.stop();
    displayer.leaderPrint ("Partitioning time: ", meshPartChrono.diff(), " s.\n");
    if ( meshFull )
    {
        meshFull = fullMesh;
    }
    else
    {
        fullMesh.reset();    //Freeing the global mesh to save memory
    }
}

//! Build a mesh from a partitioned mesh
/*!
  @param mesh The mesh that we want to generate
  @param regionFlag Flag of the region
  @param m Number of elements along the ( length, width, height )
  @param l length of the mesh ( length, width, height )
  @param t translation of the mesh along the (x,y,z)-axis
  @param verbose Verbose mode enabled/disabled
*/
/*!
    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
 */
template< typename RegionMeshType>
void fillWithStructuredMesh ( boost::shared_ptr< RegionMeshType >& mesh,
                              markerID_Type regionFlag,
                              const std::vector<UInt>& m,
                              bool verbose = false,
                              const std::vector<Real>& l = std::vector<Real> (3, 1),
                              const std::vector<Real>& t = std::vector<Real> (3, 0) )
{
    boost::shared_ptr< RegionMeshType > tmpMeshFull ( new RegionMeshType );
    fillWithStructuredMesh ( mesh, tmpMeshFull, regionFlag, m, verbose, l, t );
}

} // namespace MeshUtility

} // namespace LifeV

#endif /* NSMESHUTILITY_HPP */
