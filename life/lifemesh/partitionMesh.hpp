/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file partitionMesh.hpp
*/
#ifndef PARTMESH
#define PARTMESH 1

//#include "Epetra_Map.h"
#include <life/lifecore/life.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefem/dof.hpp>
#include "Epetra_MpiComm.h"
#include <parmetis.h>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <set>

namespace LifeV
{
/*!
  \class partitionMesh

  This class implements the partitioning of a mesh using (par)Metis

*/
template<typename Mesh>
class partitionMesh
{

public:

    typedef boost::shared_ptr<Mesh> mesh_type;

    //constructors
    partitionMesh(Mesh           &_mesh,
                  Epetra_Comm    &_comm);
    //destructor
    ~partitionMesh();

    //! return a copy of vertexDist as a container
    const std::vector<int>& vertexDist() const { return M_vertexDist; };

    mesh_type          mesh(){return M_mesh;}


    const std::vector<int>&       repeatedNodeVector() const {return M_repeatedNodeVector;}
    const std::vector<int>&       repeatedEdgeVector() const {return M_repeatedEdgeVector;}
    const std::vector<int>&       repeatedFaceVector() const {return M_repeatedFaceVector;}
    const std::vector<int>&       repeatedElemVector() const {return M_repeatedElemVector;}

//    void               createFullMap(const RefFE& refFE);

private:

    mesh_type          M_mesh;

    std::vector<int>        M_vertexDist;
    std::vector<int>        M_iadj;
    std::vector<int>        M_jadj;

    std::vector<int>        M_localNodes;
    std::set<int>           M_localEdges;
    std::set<int>           M_localFaces;
    std::vector<int>        M_localElements;

    std::vector<int>        M_repeatedNodeVector;
    std::vector<int>        M_repeatedEdgeVector;
    std::vector<int>        M_repeatedFaceVector;
    std::vector<int>        M_repeatedElemVector;

    std::map<int, int>      M_globalToLocalNode;
    std::map<int, int>      M_globalToLocalEdge;
    std::map<int, int>      M_globalToLocalFace;
    std::map<int, int>      M_globalToLocalElem;

    Epetra_Comm*       M_comm;
    UInt               M_me;

    int*               M_localIndex;
}; // class partitionMesh



//
// IMPLEMENTATION
//



template<typename Mesh>
partitionMesh<Mesh>::partitionMesh( Mesh &_mesh, Epetra_Comm &_comm):
    M_mesh              (new Mesh),
    M_localNodes        (),
    M_localEdges        (),
    M_localFaces        (),
    M_localElements     (),
    M_repeatedNodeVector(),
    M_repeatedEdgeVector(),
    M_repeatedFaceVector(),
    M_repeatedElemVector(),
    M_globalToLocalNode (),
    M_globalToLocalEdge (),
    M_globalToLocalFace (),
    M_globalToLocalElem (),
    M_comm              (&_comm)
{

    UInt elementNodes;

    typedef typename Mesh::VolumeShape ElementShape;

    switch ( ElementShape::Shape )
    {
        case HEXA:
            elementNodes = 8;
            break;
        case TETRA:
            elementNodes = 4;
            break;
        default:
            ERROR_MSG( "Element shape not implement in partitionMesh" );
    }

    UInt elementFaces;

    typedef typename Mesh::FaceShape FaceShape;

    switch ( FaceShape::Shape )
    {
        case QUAD:
            elementFaces = 6;
            break;
        case TRIANGLE:
            elementFaces = 4;
            break;
        default:
            ERROR_MSG( "Face Shape not implemented in partitionMesh" );
    }


    UInt elementEdges = 6;

    int npes;

    npes = M_comm->NumProc();
    M_me = M_comm->MyPID();

    M_vertexDist.resize(npes + 1);

    M_vertexDist[0] = 0;

    UInt k = _mesh.numVolumes();

    for (int i = 0; i < npes; ++i)
    {
      UInt l = k/(npes - i);
      M_vertexDist[ i + 1 ] = M_vertexDist[ i ] + l;
      k -= l;
    }


    // building up the neighbour arrays


    UInt localStart = M_vertexDist[M_me] + 1;
    UInt localEnd   = M_vertexDist[M_me + 1] + 1;


    M_iadj.resize(0);
    M_iadj.push_back(0);

    UInt sum = 0;

    for ( UInt ie = localStart; ie < localEnd; ++ie )
    {
        for (UInt iface = 1; iface <= elementFaces; ++iface)
        {
            UInt face = _mesh.localFaceId( ie, iface );
            UInt vol  = _mesh.face(face).ad_first();
            if (vol == ie) vol = _mesh.face(face).ad_second();
            if (vol != 0)
                {
                    M_jadj.push_back(vol - 1);
                    ++sum;
                }
        }
        M_iadj.push_back(sum);
    }

    //MPI_Barrier(MPI_COMM_WORLD);

    // parMetis part

    std::vector<int> xadj;
    std::vector<int> adjncy;

    int*   vwgt    = 0;
    int*   adjwgt  = 0;

    int    wgtflag = 0;
    int    ncon    = 1;
    int    numflag = 0;

    int    edgecut;

    std::vector<int>  part(vertexDist()[M_me + 1] -
                      vertexDist()[M_me]);

    std::vector<int>  options(3,0);

    options[0] = 1;
    options[1] = 3;
    options[2] = 1;

    int nparts = M_comm->NumProc();

    std::vector<float> tpwgts(ncon*nparts, 1./nparts);
    std::vector<float> ubvec (ncon, 1.05);

    Epetra_MpiComm* mpiComm = dynamic_cast<Epetra_MpiComm*>(M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    int rank, nprocs;
    MPI_Comm_rank(MPIcomm, &rank);
    MPI_Comm_size(MPIcomm, &nprocs);

//     std::cout << "rank   = " << rank << std::endl;
//     std::cout << "nprocs = " << nprocs << std::endl;

    ParMETIS_V3_PartKway((int*) &M_vertexDist[0],
                         (int*) &M_iadj[0],
                         (int*) &M_jadj[0],
                         vwgt, adjwgt,
                         &wgtflag,
                         &numflag, &ncon, &nparts,
                         &tpwgts[0], &ubvec[0], &options[0],
                         &edgecut, &part[0], &MPIcomm);

    M_comm->Barrier();

    int nProc;
    nProc = M_comm->NumProc();

    std::vector< std::vector<int> > locProc(nProc, 0);

    std::vector<int>           locVolume(0);

    for (UInt ii = 0; ii < part.size(); ++ii)
    {
        locProc[part[ii]].push_back(ii + vertexDist()[M_me]);
    }

    for (int iproc = 0; iproc < nProc; ++iproc)
    {
        if (int(M_me) != iproc)
            {
                int size = locProc[iproc].size();
                MPI_Send(&size, 1, MPI_INT, iproc, 666, MPIcomm);
                MPI_Send(&locProc[iproc][0], size, MPI_INT, iproc, 666, MPIcomm);
            }
    }


    for (int iproc = 0; iproc < nProc; ++iproc)
    {
        if (int(M_me) != iproc)
        {
            int        size;
            MPI_Status status;
            MPI_Recv(&size, 1, MPI_INT, iproc, 666, MPIcomm, &status);
            std::vector<int> stack(size, 0);
            MPI_Recv(&stack[0], size, MPI_INT, iproc, 666, MPIcomm, &status);

            for (int jj = 0; jj < size; ++jj)
            {
                locProc[M_me].push_back(stack[jj]);
            }
        }
    }


    std::cout << M_me << " has " << locProc[M_me].size() << " elements." << std::endl;

    // local mesh construction

    if (!M_me) std::cout << "Building local mesh ... " << std::flush;

    std::map<int, int>::iterator  im;
    std::set<int>::iterator       is;

    int count = 1;
    UInt ielem;
    UInt inode;


    for (UInt jj = 0; jj < locProc[M_me].size(); ++jj)
    {
        ielem = locProc[M_me][jj];
        M_localElements.push_back(ielem);

        for (UInt ii = 1; ii <= elementNodes; ++ii)
        {
            inode = _mesh.volume(ielem + 1).point(ii).id();
            im    = M_globalToLocalNode.find(inode);

//            std::cout << inode << " ";

            if ( im == M_globalToLocalNode.end() )
            {
                M_globalToLocalNode.insert(std::make_pair(inode, count));
                ++count;
                M_localNodes.push_back(_mesh.volume(ielem + 1).point(ii).id());
            }
        }

        for (UInt ii = 1; ii <= elementEdges; ++ii)
        {
            M_localEdges.insert(_mesh.localEdgeId(ielem + 1, ii));
        }

        for (UInt ii = 1; ii <= elementFaces; ++ii)
        {
            M_localFaces.insert(_mesh.localFaceId(ielem + 1, ii));
        }
    }

    std::vector<int>::iterator it;

    // nodes construction

    UInt nBoundaryPoint = 0;
    M_mesh->pointList.reserve(M_localNodes.size());
    M_mesh->_bPoints.reserve(_mesh.numBPoints()*M_localNodes.size()/_mesh.numBPoints()); // guessing how many boundary points on this processor.

    inode = 1;

    typename Mesh::PointType * pp = 0;

    for (it = M_localNodes.begin(); it != M_localNodes.end(); ++it, ++inode)
    {

        typename Mesh::PointType point = 0;

        bool boundary = _mesh.isBoundaryPoint(*it);
        if (boundary)
        {
            ++nBoundaryPoint;
        }

        pp = &M_mesh->addPoint( boundary );
        pp->setMarker( _mesh.point(*it).marker() );

        pp->x() = _mesh.point(*it).x();
        pp->y() = _mesh.point(*it).y();
        pp->z() = _mesh.point(*it).z();

        UInt id = _mesh.point(*it).id();

        pp->setId(id);
        pp->setLocalId( inode );

        M_mesh->localToGlobalNode().insert(std::make_pair(id, inode));
        M_mesh->globalToLocalNode().insert(std::make_pair(inode, id));

//         point.setMarker( _mesh.point(*it).marker() );

//         point.x() = _mesh.point(*it).x();
//         point.y() = _mesh.point(*it).y();
//         point.z() = _mesh.point(*it).z();

//         std::cout << point.x() << " " << point.y() << " " << point.z() << " " << std::flush;
//         UInt id = _mesh.point(*it).id();

//         point.setId(id);
//         point.setLocalId( inode );

//         M_mesh->pointList[id - 1] = point;
//         std::cout << id << " " << inode << std::endl;


    }
    // volumes construction

    count = 1;

    typename Mesh::VolumeType * pv = 0;

    M_mesh->volumeList.reserve(M_localElements.size());

    for (it = M_localElements.begin(); it != M_localElements.end(); ++it, ++count)
    {
        pv = &M_mesh->addVolume();
//        std::cout << "volume " << _mesh.volume(*it + 1).id() << std::flush;
        pv->setId     ( _mesh.volume(*it + 1).id() );
        pv->setLocalId( count );

        M_globalToLocalElem.insert(make_pair( _mesh.volume(*it + 1).id(),  count ));

        for (ID id = 1; id <= elementNodes; ++id)
        {
            inode = _mesh.volume(*it + 1).point(id).id();
            im    = M_globalToLocalNode.find(inode);
            pv->setPoint( id, M_mesh->pointList( (*im).second ) );
        }

        int ibc = _mesh.volume(*it + 1).marker();

        pv->setMarker( EntityFlag( ibc ) );

    }

#if 0
    for (UInt ii = 0; ii < M_localElements.size(); ++ii)
    {
        for (UInt jj = 0; jj < 4; ++jj)
            std::cout << M_mesh->volume(ii + 1).point(jj + 1).id() << " " ;
        std::cout << std::endl;
    }

    for (UInt ii = 0; ii < M_localElements.size(); ++ii)
    {
        for (UInt jj = 0; jj < 4; ++jj)
            std::cout << M_mesh->volume(ii + 1).point(jj + 1).x() << " "
                      << M_mesh->volume(ii + 1).point(jj + 1).y() << " "
                      << M_mesh->volume(ii + 1).point(jj + 1).z() << " - ";
        std::cout << std::endl;
    }
#endif

    // edges construction

    typename Mesh::EdgeType * pe;

    count = 1;

    UInt nBoundaryEdges = 0;
    M_mesh->edgeList.reserve(M_localEdges.size());


    for (is = M_localEdges.begin(); is != M_localEdges.end(); ++is, ++count)
    {
        bool boundary = (_mesh.isBoundaryEdge(*is));
        if (boundary)
        {
            ++nBoundaryEdges;
        }

//        std::cout << "is = " << *is << " ";
//        pe = M_mesh->addEdge(_mesh.edge(*is), boundary);
        pe = &M_mesh->addEdge( boundary );

        pe->setId     (_mesh.edge(*is).id());
        pe->setLocalId( count );

        for (ID id = 1; id <= 2; ++id)
        {
            inode = _mesh.edge(*is).point(id).id();
//            std::cout << inode << std::endl;
            im = M_globalToLocalNode.find(inode);
            pe->setPoint(id, M_mesh->pointList((*im).second));
        }

        pe->setMarker( _mesh.edge(*is).marker() );
    }

    // faces construction

    typename Mesh::FaceType * pf = 0;

    count = 1;

    UInt nBoundaryFace = 0;
    M_mesh->faceList.reserve(M_localFaces.size());

    for (is = M_localFaces.begin(); is != M_localFaces.end(); ++is, ++count)
    {
        bool boundary = (_mesh.isBoundaryFace(*is));
        if (boundary)
        {
            ++nBoundaryFace;
        }

        pf =  &M_mesh->addFace( boundary );

//         std::cout << *is << " ";
//         std::cout << (_mesh.face(*is).id()) << std::endl;

        pf->setId     (_mesh.face(*is).id());
//        pf->setId     (*is);
        pf->setLocalId( count );

        int elem1 =  _mesh.face(*is).ad_first();
        int elem2 =  _mesh.face(*is).ad_second();

        im =  M_globalToLocalElem.find(elem1);

        int localElem1;

        if (im == M_globalToLocalElem.end()) localElem1 = 0;
        else
            localElem1 = (*im).second;

        im =  M_globalToLocalElem.find(elem2);

        int localElem2;

        if (im == M_globalToLocalElem.end()) localElem2 = 0;
        else
            localElem2 = (*im).second;


//        UInt localElem2 = (*im).second;


        if ((localElem1 == 0) && ! boundary) localElem1 = localElem2;
        if ((localElem2 == 0) && ! boundary) localElem2 = localElem1;

//         std::cout << me << " : ";
//         std::cout << elem1 << " -> ";
//         std::cout << localElem1 << " ";

        pf->ad_first() = localElem1;
        pf->pos_first() = _mesh.face(*is).pos_first();

//         std::cout << elem2 << " -> ";
//         std::cout << localElem2 << std::endl;

        pf->ad_second()  = localElem2;
        pf->pos_second() = _mesh.face(*is).pos_second();

        for (ID id = 1; id <= 3; ++id)
            {
                inode =_mesh.face(*is).point(id).id();
                im = M_globalToLocalNode.find(inode);
                pf->setPoint(id, M_mesh->pointList( (*im).second ));
            }
        pf->setMarker( _mesh.face(*is).marker());

        M_mesh->setLinkSwitch( "HAS_ALL_FACES" );
        M_mesh->setLinkSwitch( "FACES_HAVE_ADIACENCY" );

    }

// final setup


    UInt nVolumes = M_localElements.size();
    UInt nNodes   = M_localNodes.size();
    UInt nEdges   = M_localEdges.size();
    UInt nFaces   = M_localFaces.size();

    M_mesh->setMaxNumPoints ( nNodes, true );
    M_mesh->setMaxNumEdges  ( nEdges, true );
    M_mesh->setMaxNumFaces  ( nFaces, true );
    M_mesh->setMaxNumVolumes( nVolumes, true);

//    std::cout << "num points " << _mesh.numPoints() << std::endl;

    M_mesh->setMaxNumGlobalPoints ( _mesh.numPoints  () );
    M_mesh->setNumGlobalVertices  ( _mesh.numPoints  () );
    M_mesh->setMaxNumGlobalEdges  ( _mesh.numEdges   () );
    M_mesh->setMaxNumGlobalFaces  ( _mesh.numFaces   () );
    M_mesh->setMaxNumGlobalVolumes( _mesh.numVolumes () );

    M_mesh->setNumBPoints   ( nBoundaryPoint );
    M_mesh->setNumBFaces    ( nBoundaryFace  );
    M_mesh->setNumBEdges    ( nBoundaryEdges );

    M_mesh->setNumVertices ( nNodes );
    M_mesh->setNumBVertices( nBoundaryPoint );



//! barrier
//    MPI_Barrier(MPI_COMM_WORLD);
//     std::cout << "____________________________________" << std::endl;
//     for ( UInt j = 0; j < M_mesh->edgeList.size();++j )
//     {
//         std::cout << "edge " << j << " := ";
//         UInt i1 = ( M_mesh->edgeList[ j ].point( 1 ) ).id();
//         std::cout << i1 << " ";
//         UInt i2 = ( M_mesh->edgeList[ j ].point( 2 ) ).id();
//         std::cout << i2 << std::endl;
//     }
//     std::cout << "____________________________________" << std::endl;


    M_mesh->updateElementEdges();
    M_mesh->updateElementFaces();

//     Switch sw;
//     if ( !checkMesh3D( *M_mesh, sw, false, true, std::cout, std::cerr, std::cout ) )
//         abort();

//! barrier
//    MPI_Barrier(MPI_COMM_WORLD);

//    if (!M_me) std::cout << "ok. " << std::flush;


//     M_mesh->showMe(true);
//    _mesh.showMe(true);
//! barrier
//    MPI_Barrier(MPI_COMM_WORLD);


//      ostringstream myStream;
//      myStream << M_me << std::flush;

//      string filename = "mesh." + myStream.str() + ".mesh";

//      std::cout << M_me << " is wrinting the local mesh ... " << std::flush;
//      writeMesh(filename, *M_mesh);

//    wr_medit_ascii("mesh.mesh", *M_mesh);
//     std::cout << "done." << std::endl;
//    M_mesh->check();
//    M_mesh->showMe();
//    MPI_Barrier(MPI_COMM_WORLD);

    if (!M_me) std::cout << "ok." << std::endl;
    if (!M_me) std::cout << "Creating the map ... " << std::flush;

    // repeated map creation

    std::vector<int> elementList = locProc[M_me];

    // repeated element map creation

    for (UInt ii = 0; ii < elementList.size(); ++ii)
    {
        ielem = elementList[ii];
        M_repeatedElemVector.push_back(ielem + 1);
    }


    // repeated node map creation

    std::set<int>    repeatedNodeList;

    for (UInt ii = 0; ii < elementList.size(); ++ii)
    {
        ielem = elementList[ii];
        for (UInt jj = 1; jj <= elementNodes; ++jj)
        {
            inode = _mesh.volume(ielem + 1).point(jj).id();
            repeatedNodeList.insert(inode);
        }
    }

    M_repeatedNodeVector.reserve(repeatedNodeList.size());

    for (is = repeatedNodeList.begin(); is != repeatedNodeList.end(); ++is)
    {
        M_repeatedNodeVector.push_back(*is);
    }

    // repeated edge list creation

    std::set<int>    repeatedEdgeList;

    for (UInt ii = 0; ii < elementList.size(); ++ii)
    {
        ielem = elementList[ii];
        for (UInt jj = 0; jj < elementEdges; ++jj)
        {
            UInt iedge = _mesh.localEdgeId(ielem + 1, jj + 1);
            repeatedEdgeList.insert((int) iedge);
        }
    }

    M_repeatedEdgeVector.reserve(repeatedEdgeList.size());

    for (is = repeatedEdgeList.begin(); is != repeatedEdgeList.end(); ++is)
    {
//        std::cout << *is << std::endl;
        M_repeatedEdgeVector.push_back(*is);
    }

    // repeated face list creation

    std::set<int>    repeatedFaceList;

    for (UInt ii = 0; ii < elementList.size(); ++ii)
    {
        ielem = elementList[ii];
        for (UInt jj = 1; jj <= elementFaces; ++jj)
        {
            UInt iface = _mesh.localFaceId(ielem + 1, jj);
            repeatedFaceList.insert(iface);
        }
    }


    M_repeatedFaceVector.reserve(repeatedFaceList.size());

    for (is = repeatedFaceList.begin(); is != repeatedFaceList.end(); ++is)
    {
        M_repeatedFaceVector.push_back(*is);
    }

    if (!M_me) std::cout << "ok." << std::endl;

}


template<typename Mesh>
partitionMesh<Mesh>::~partitionMesh()
{
}


}
#endif
