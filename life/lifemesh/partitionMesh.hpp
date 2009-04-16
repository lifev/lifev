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
    // First of all, we want to know which kind of elements the mesh is built of:
    // How many nodes does each element hold?
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

    // How many faces does each element hold?
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

    // How many edges does each element hold?
    UInt elementEdges = 6;

    // ParMETIS is able to work in parallel: how many processors does it have at hand?
    int npes;

    npes = M_comm->NumProc();
    M_me = M_comm->MyPID();

    // CAREFUL: ParMetis works on a graph abstraction.
    // A graph is built over the data structure to be split, each vertex being a mesh element
    // so hereby a "vertex" is actually a _graph_ vertex, i. e. a mesh element
    M_vertexDist.resize(npes + 1);

    M_vertexDist[0] = 0;

    UInt k = _mesh.numVolumes();

    // Evenly distributed graph vertices
    for (int i = 0; i < npes; ++i)
    {
      UInt l = k/(npes - i);
      M_vertexDist[ i + 1 ] = M_vertexDist[ i ] + l;
      k -= l;
    }
    ASSERT( k == 0, "At this point we should have 0 volumes left" ) ;


    // Now each processor will take care of its own graph vertices (i. e. mesh elements).
    // Nothing guarantees about the neighbor elements distribution across the processors,
    // since as of now we just split the set of volumes based on IDs.
    // Here we building up the neighbor arrays.

    UInt localStart = M_vertexDist[M_me] + 1;
    UInt localEnd   = M_vertexDist[M_me + 1] + 1;


    M_iadj.resize(0);
    M_iadj.push_back(0);

    UInt sum = 0;

    for ( UInt ie = localStart; ie < localEnd; ++ie )
    {
        for (UInt iface = 1; iface <= elementFaces; ++iface)
        {
            // global ID of the iface-th face in element ie
            UInt face = _mesh.localFaceId( ie, iface );
            // first adjacent element to face "face"
            UInt vol  = _mesh.face(face).ad_first();
            if (vol == ie) vol = _mesh.face(face).ad_second();
            if (vol != 0)
                {
                    // this is the list of adjacency
                    // for each graph vertex, simply push back the ID of its neighbors
                    M_jadj.push_back(vol - 1);
                    ++sum;
                }
        }
        // this is the list of "keys" to access M_jadj
        // graph element i has neighbors M_jadj[ k ],
        // with M_iadj[i] <= k < M_iadj[i+1]
        M_iadj.push_back(sum);
    }

    //MPI_Barrier(MPI_COMM_WORLD);

    // **************
    // parMetis part
    // **************
    // Each processor has three lists: xadj, adjncy, vtxdist
    std::vector<int> xadj;
    std::vector<int> adjncy;

    // these two arrays are to be used for weighted graphs:
    // usually we will set them to NULL
    int*   vwgt    = 0;
    int*   adjwgt  = 0;

    int    wgtflag = 0; // 0 means the graphs is not weighted
    int    ncon    = 1; // number of weights attached to each graph vertex
    int    numflag = 0; // numbering scheme: 0 is C-style (arrays indexes starting from 0)

    int    edgecut; // here will be stored the number of edges cut in the partitioning process

    // This array's size is equal to the number of locally-stored vertices:
    // at the end of the partitioning process, "part" will contain the partitioning array:
    // part[m] = n; means that graph vertex m belongs to subdomain n
    std::vector<int>  part(vertexDist()[M_me + 1] -
                      vertexDist()[M_me]);

    std::vector<int>  options(3,0);

    // additional options
    options[0] = 1; // means that additional options are actually passed
    options[1] = 3; // level of information to be returned during execution (see ParMETIS's defs.h file)
    options[2] = 1; // random number seed for the ParMETIS routine

    // number of desired subdomains: can be different from the number of procs
    int nparts = M_comm->NumProc();

    // fraction of vertex weight to be distributed to each subdomain.
    // here we want the subdomains to be of the same size
    std::vector<float> tpwgts(ncon*nparts, 1./nparts);
    // imbalance tolerance for each vertex weight
    std::vector<float> ubvec (ncon, 1.05);

    Epetra_MpiComm* mpiComm = dynamic_cast<Epetra_MpiComm*>(M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    int rank, nprocs;
    MPI_Comm_rank(MPIcomm, &rank);
    MPI_Comm_size(MPIcomm, &nprocs);

//     std::cout << "rank   = " << rank << std::endl;
//     std::cout << "nprocs = " << nprocs << std::endl;
    /*
     (from ParMETIS v 3.1 manual)
     This routine is used to compute a k-way partitioning of a graph
     on p processors using the multilevel k-way multi-constraint
     partitioning algorithm.
     */
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

    // this is a vector of subdomains: each component is
    // the list of vertices belonging to the specific subdomain
    std::vector< std::vector<int> > locProc(nProc);

    std::vector<int>           locVolume(0);

    // cycling on locally stored vertices
    for (UInt ii = 0; ii < part.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        locProc[part[ii]].push_back(ii + vertexDist()[M_me]);
    }

    // sizeof(MPI_INT) is expressed in bytes
    // MPI_INT range is ( -2^(7*sizeof(MPI_INT), 2^(7*sizeof(MPI_INT) - 1 )
    // (MPI_INT is not unsigned int)
    //int _exp(7*sizeof(MPI_INT));
    //int _max_int( std::ldexp(1., _exp-1) ); // taking a little less: 2^(exp-1) instead of 2^exp - 1
    //int max_int (1000);

    // cycling on subdomains
    // TODO: Matteo please comment this part :)
    for (int iproc = 0; iproc < nProc; ++iproc)
    {
        // all processes other than me are sending vertices
        // belonging to my subdomain
        if (int(M_me) != iproc)
        {                                                  //start if
            int size = locProc[iproc].size();

            // tell me how many vertices belonging to me you have to send me
            MPI_Send(&size, 1, MPI_INT, iproc, 10, MPIcomm);

            // workaround for huge data to be passed
            if (size > 1000)
            {
                int incr = 1 ;
                int pos = 0;
                int size_part =size;

                // divide the whole data set into smaller packets
                while (size_part > 1000 )
                {
                    incr+=1;
                    size_part=size/incr;
                }

                MPI_Send(&incr, 1, MPI_INT, iproc, 20, MPIcomm);
                MPI_Send(&size_part, 1, MPI_INT, iproc, 30, MPIcomm);

                for ( int kk = 0; kk < incr; ++kk)
                {         //for
                    MPI_Send(&pos, 1, MPI_INT, iproc, 100+kk, MPIcomm);
                    MPI_Send(&locProc[iproc][pos], size_part, MPI_INT, iproc, 5000000+kk, MPIcomm);
                    pos = pos + size_part;
                }      //endfor

                int resto = size%incr;

                MPI_Send(&resto, 1, MPI_INT, iproc, 80, MPIcomm);

                if(resto!=0)
                {
                    MPI_Send(&pos, 1, MPI_INT, iproc, 40, MPIcomm);
                    MPI_Send(&locProc[iproc][pos],  resto, MPI_INT, iproc, 50, MPIcomm);
                }
            }
            else
                MPI_Send(&locProc[iproc][0], size, MPI_INT, iproc, 60, MPIcomm);
        }                    //end if
    }                  //end for


    for (int iproc = 0; iproc < nProc; ++iproc)
    {
        if (int(M_me) != iproc)
        {
            int        size;
            MPI_Status status;
            MPI_Recv(&size, 1, MPI_INT, iproc, 10, MPIcomm, &status);

            std::vector<int> stack(size, 0);

            if (size > 1000 )
            {
                int size_part, pos, incr;

                MPI_Recv(&incr, 1, MPI_INT, iproc, 20, MPIcomm, &status);
                MPI_Recv(&size_part, 1, MPI_INT, iproc, 30, MPIcomm, &status);

                for ( int kk = 0; kk < incr; ++kk)
                {
                    MPI_Recv(&pos, 1, MPI_INT, iproc, 100+kk, MPIcomm, &status);
                    MPI_Recv(&stack[pos], size_part , MPI_INT, iproc, 5000000+kk, MPIcomm, &status);
                }
                int resto=0;
                MPI_Recv(&resto, 1, MPI_INT, iproc, 80, MPIcomm, &status);

                if(resto!=0)
                {
                    MPI_Recv(&pos, 1, MPI_INT, iproc, 40, MPIcomm, &status);
                    MPI_Recv(&stack[pos],  resto, MPI_INT, iproc, 50, MPIcomm, &status);
                }
            }
            else
                MPI_Recv(&stack[0], size , MPI_INT, iproc, 60, MPIcomm, &status);

            for (int jj = 0; jj < size; ++jj)
            {
                locProc[M_me].push_back(stack[jj]);
            }
        }
    }


    std::cout << M_me << " has " << locProc[M_me].size() << " elements." << std::endl;

    // ***********************
    // local mesh construction
    // ***********************

    if (!M_me) std::cout << "Building local mesh ... \n" << std::flush;

    std::map<int, int>::iterator  im;
    std::set<int>::iterator       is;

    int count = 1;
    UInt ielem;
    UInt inode;

    // cycle on local element's ID
    for (UInt jj = 0; jj < locProc[M_me].size(); ++jj)
    {
        ielem = locProc[M_me][jj];
        M_localElements.push_back(ielem);

        // cycle on element's nodes
        for (UInt ii = 1; ii <= elementNodes; ++ii)
        {
            inode = _mesh.volume(ielem + 1).point(ii).id();
            im    = M_globalToLocalNode.find(inode);

//            std::cout << inode << " ";

            // if the node is not yet present in the list of local nodes, then add it
            // CAREFUL: also local numbering starts from 1 in RegionMesh
            if ( im == M_globalToLocalNode.end() )
            {
                M_globalToLocalNode.insert(std::make_pair(inode, count));
                ++count;
                // store here the global numbering of the node
                M_localNodes.push_back(_mesh.volume(ielem + 1).point(ii).id());
            }
        }

        // cycle on element's edges
        for (UInt ii = 1; ii <= elementEdges; ++ii)
        {
          // store here the global numbering of the edge
            M_localEdges.insert(_mesh.localEdgeId(ielem + 1, ii));
        }

        // cycle on element's faces
        for (UInt ii = 1; ii <= elementFaces; ++ii)
        {
          // store here the global numbering of the face
            M_localFaces.insert(_mesh.localFaceId(ielem + 1, ii));
        }
    }

    std::vector<int>::iterator it;

    // ******************
    // nodes construction
    // ******************

    UInt nBoundaryPoint = 0;
    M_mesh->pointList.reserve(M_localNodes.size());
    // guessing how many boundary points on this processor.
    M_mesh->_bPoints.reserve(_mesh.numBPoints()*M_localNodes.size()/_mesh.numBPoints());

    inode = 1;

    typename Mesh::PointType * pp = 0;

    // loop in the list of local nodes:
    // in this loop inode is the local numbering of the points
    for (it = M_localNodes.begin(); it != M_localNodes.end(); ++it, ++inode)
    {

        typename Mesh::PointType point = 0;

        // create a boundary point in the local mesh, if needed
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

        M_mesh->localToGlobalNode().insert(std::make_pair(inode, id));
        M_mesh->globalToLocalNode().insert(std::make_pair(id, inode));

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

    // ******************
    // volumes construction
    // ******************

    count = 1;

    typename Mesh::VolumeType * pv = 0;

    M_mesh->volumeList.reserve(M_localElements.size());

    // loop in the list of local elements
    // CAREFUL! in this loop inode is the global numbering of the points
    // We insert the local numbering of the nodes in the local volume list
    for (it = M_localElements.begin(); it != M_localElements.end(); ++it, ++count)
    {
        pv = &M_mesh->addVolume();
//        std::cout << "volume " << _mesh.volume(*it + 1).id() << std::flush;
        // CAREFUL! in ParMETIS data structures, numbering starts from 0
        pv->setId     ( _mesh.volume(*it + 1).id() );
        pv->setLocalId( count );

        M_globalToLocalElem.insert(make_pair( _mesh.volume(*it + 1).id(),  count ));

        for (ID id = 1; id <= elementNodes; ++id)
        {
            inode = _mesh.volume(*it + 1).point(id).id();
            // im is an iterator to a map element
            // im->first is the key (i. e. the global ID "inode")
            // im->second is the value (i. e. the local ID "count")
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

    // ******************
    // edges construction
    // ******************

    typename Mesh::EdgeType * pe;

    count = 1;

    UInt nBoundaryEdges = 0;
    M_mesh->edgeList.reserve(M_localEdges.size());

    // loop in the list of local edges
    for (is = M_localEdges.begin(); is != M_localEdges.end(); ++is, ++count)
    {
      // create a boundary edge in the local mesh, if needed
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
            // im is an iterator to a map element
            // im->first is the key (i. e. the global ID "inode")
            // im->second is the value (i. e. the local ID "count")
            im = M_globalToLocalNode.find(inode);
            pe->setPoint(id, M_mesh->pointList((*im).second));
        }

        pe->setMarker( _mesh.edge(*is).marker() );
    }

    // ******************
    // faces construction
    // ******************

    typename Mesh::FaceType * pf = 0;

    count = 1;

    UInt nBoundaryFace = 0;
    M_mesh->faceList.reserve(M_localFaces.size());

    // loop in the list of local faces
    for (is = M_localFaces.begin(); is != M_localFaces.end(); ++is, ++count)
    {
      // create a boundary face in the local mesh, if needed
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

        // find the mesh elements adjacent to the face
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

        // if this process does not own either of the adjacent elements
        // then the two adjacent elements coincide in the local mesh
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

    // ******************
    // final setup
    // ******************


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

    // *********************
    // repeated map creation
    // *********************

    std::vector<int> elementList = locProc[M_me];

    // repeated element map creation

    for (UInt ii = 0; ii < elementList.size(); ++ii)
    {
        ielem = elementList[ii];
        M_repeatedElemVector.push_back(ielem + 1);
    }


    // repeated node map creation

    // use a set to store each node only once
    std::set<int>    repeatedNodeList;
    // cycle on the elements nodes
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

    // use a set to store each node only once
    std::set<int>    repeatedEdgeList;

    // cycle on the elements edges
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

    // use a set to store each node only once
    std::set<int>    repeatedFaceList;

    // cycle on the elements faces
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
