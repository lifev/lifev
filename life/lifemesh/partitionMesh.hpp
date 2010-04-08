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
                  Epetra_Comm    &_comm,
                  Epetra_Map* interfaceMap = 0, Epetra_Map* interfaceMapRep = 0
                  );
    //destructor
    ~partitionMesh();

    //! return a copy of vertexDist as a container
    const std::vector<int>& vertexDist() const { return M_vertexDist; };

    const mesh_type     mesh() const {return M_mesh;}


    const std::vector<int>&       repeatedNodeVector() const {return M_repeatedNodeVector;}
    const std::vector<int>&       repeatedEdgeVector() const {return M_repeatedEdgeVector;}
    const std::vector<int>&       repeatedFaceVector() const {return M_repeatedFaceVector;}
    const std::vector<int>&       repeatedVolumeVector() const {return M_repeatedVolumeVector;}

//    void               createFullMap(const RefFE& refFE);

private:

    mesh_type          M_mesh;

    std::vector<int>        M_vertexDist;
    std::vector<int>        M_iadj;
    std::vector<int>        M_jadj;

    std::vector<int>        M_localNodes;
    std::set<int>           M_localEdges;
#ifdef TWODIM
    std::vector<int>        M_localFaces;
#elif defined THREEDIM
    std::set<int>           M_localFaces;
#endif
    std::vector<int>        M_localVolumes;

    std::vector<int>        M_repeatedNodeVector;
    std::vector<int>        M_repeatedEdgeVector;
    std::vector<int>        M_repeatedFaceVector;
    std::vector<int>        M_repeatedVolumeVector;

    std::map<int, int>      M_globalToLocalNode;
    std::map<int, int>      M_globalToLocalEdge;
    std::map<int, int>      M_globalToLocalFace;
    std::map<int, int>      M_globalToLocalVolume;

    Epetra_Comm*       M_comm;
    UInt               M_me;

}; // class partitionMesh



//
// IMPLEMENTATION
//



template<typename Mesh>
partitionMesh<Mesh>::partitionMesh( Mesh &_mesh, Epetra_Comm &_comm,
                                    Epetra_Map* interfaceMap,
                                    Epetra_Map* interfaceMapRep):
    M_mesh              (new Mesh),
    M_localNodes        (),
    M_localEdges        (),
    M_localFaces        (),
    M_localVolumes     (),
    M_repeatedNodeVector(),
    M_repeatedEdgeVector(),
    M_repeatedFaceVector(),
    M_repeatedVolumeVector(),
    M_globalToLocalNode (),
    M_globalToLocalEdge (),
    M_globalToLocalFace (),
    M_globalToLocalVolume (),
    M_comm              (&_comm)
{

    // First of all, we want to know which kind of elements the mesh is built of:
    // How many nodes, faces, edges does each element hold?
    UInt elementNodes, elementFaces, elementEdges;
    UInt faceNodes(0);

    typedef typename Mesh::ElementShape ElementShape;

    switch ( ElementShape::Shape )
    {
        case HEXA:
            elementNodes = 8;
            elementFaces = 6;
            elementEdges = 12;
            faceNodes    = 4;
            break;
        case TETRA:
            elementNodes = 4;
            elementFaces = 4;
            elementEdges = 6;
            faceNodes    = 3;
            break;
        case QUAD:
            elementNodes = 4;
            elementEdges = 4;
            elementFaces = 0;
            break;
        case TRIANGLE:
           	elementNodes = 3;
           	elementEdges = 3;
           	elementFaces = 0;
            break;
        default:
            ERROR_MSG( "Face Shape not implemented in partitionMesh" );
    }



    // ParMETIS is able to work in parallel: how many processors does it have at hand?
    int npes;

    npes = M_comm->NumProc();
    M_me = M_comm->MyPID();

    // CAREFUL: ParMetis works on a graph abstraction.
    // A graph is built over the data structure to be split, each vertex being a mesh element
    // so hereby a "vertex" is actually a _graph_ vertex, i. e. a mesh element
    M_vertexDist.resize(npes + 1);

    M_vertexDist[0] = 0;

    UInt k = _mesh.numElements();

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


    // this vector contains the weights for the edges of the graph,
    // it is set to null if it is not used.
    std::vector<int>        adjwgt;

    M_iadj.resize(0);
    M_iadj.push_back(0);

    UInt sum = 0;

#ifdef TWODIM

    for ( UInt ie = localStart; ie < localEnd; ++ie )
    {
        for (UInt iedge = 1; iedge <= elementEdges; ++iedge)
        {
        	// global ID of the iedge-th edge in element ie
            UInt edge = _mesh.localEdgeId( ie, iedge );
            // first adjacent element to edge "edge"
            UInt elem  = _mesh.edge(edge).ad_first();
            if (elem == ie) elem = _mesh.edge(edge).ad_second();
            if (elem != 0)
                {
            	    // this is the list of adjacency
            	    // for each graph vertex, simply push back the ID of its neighbors
                    M_jadj.push_back(elem - 1);
                    ++sum;
                }
        }
        // this is the list of "keys" to access M_jadj
        // graph element i has neighbors M_jadj[ k ],
        // with M_iadj[i] <= k < M_iadj[i+1]
        M_iadj.push_back(sum);
    }

#elif defined THREEDIM


    // In fluid-structure interaction:
    // *    If the solid mesh is not partitioned the following part won't be
    //      executed
    // *    If the solid mesh is partitioned:
    //      - The fluid is partitioned first
    //      - The solid mesh partition tries to follow the partition of the fluid
    //      This is achieved by specifying a weight to some edge of the graph.
    //      The interface between two processors is the set of the nodes that for
    //      at least one processor are on the repeated map and not on the unique map.
    //      That's why the constructor needs both the unique and repeated maps
    //      on the interface
    //////////////////// BEGIN OF SOLID PARTITION PART ////////////////////////
    boost::shared_ptr<std::vector<int> > repeatedFace;// used for the solid partitioning
    std::vector<int> myRepeatedFace;// used for the solid partitioning
    boost::shared_ptr<std::vector<int> >  myisOnProc;// used for the solid partitioning
    boost::shared_ptr<std::vector<int> >  isOnProc;// used for the solid partitioning
    if(interfaceMap)
        {
            myisOnProc.reset(new std::vector<int>(_mesh.numVolumes()));
            bool myFaceRep;
            bool myFace(false);
            short count;
            for(UInt h=0; h<_mesh.numVolumes(); ++h)
                {
                    //(*ordVolumeList)[h]=h;
                    (*myisOnProc)[h]=-1;
                }
            k=0;
            //order.reset(new std::vector<UInt>());
            //order->reserve(localEnd-localStart);

            // This loop is throughout the whole unpartitioned mesh,
            // it is expensive and not scalable.
            // Bad, this part should be done offline
        for ( UInt ie = 1; ie <= _mesh.numVolumes(); ++ie )
            {
                for (UInt iface = 1; iface <= elementFaces; ++iface)
                    {
                        UInt face = _mesh.localFaceId( ie, iface );
                        UInt vol  = _mesh.face(face).ad_first();
                        if (vol == ie) vol = _mesh.face(face).ad_second();
                        if (vol != 0)
                            {
                                myFace=false;
                                myFaceRep=false;
                                count=0;

                                for(int ipoint=1; ipoint<=(int)faceNodes; ++ipoint)//vertex-based dofs
                                    {
                                        myFaceRep = ((interfaceMap->LID(_mesh.face(face).point(ipoint).id()/*first is fluid*/) == -1)&&(interfaceMapRep->LID(_mesh.face(face).point(ipoint).id()/*first is fluid*/) != -1));//to avoid repeated
                                        myFace = myFace ||(interfaceMap->LID(_mesh.face(face).point(ipoint).id()) != -1);
                                        if(myFaceRep)
                                            {++count;}
                                    }
                                if(count > 1)
                                    {
                                        myRepeatedFace.push_back(1);
                                    }
                                else
                                    {
                                        myRepeatedFace.push_back(0);
                                    }
                            }
                        if(myFace)
                            {
                            (*myisOnProc)[ie-1]=M_me;
                            }
                    }
            }

        repeatedFace.reset(new std::vector<int>(myRepeatedFace.size()));
        isOnProc.reset(new std::vector<int>(*myisOnProc));

        // Lot of communication here!!
        MPI_Allreduce( &myRepeatedFace[0], &(*repeatedFace)[0], myRepeatedFace.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce( &(*myisOnProc)[0], &(*isOnProc)[0], myisOnProc->size(), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
         }
    //////////////////// END OF SOLID PARTITION PART ////////////////////////

    for ( UInt ie = localStart; ie < localEnd; ++ie )
    {
        for (UInt iface = 1; iface <= elementFaces; ++iface)
        {
            // global ID of the iface-th face in element ie
            UInt face = _mesh.localFaceId( ie, iface );
            // first adjacent element to face "face"
            UInt elem  = _mesh.face(face).ad_first();
            if (elem == ie) elem = _mesh.face(face).ad_second();
            if (elem != 0)
                {
                    // this is the list of adjacency
                    // for each graph vertex, simply push back the ID of its neighbors
                    M_jadj.push_back(elem - 1);
                    ++sum;
                    if(interfaceMap)// if I'm partitioning the solid in FSI
                        {
                            if((*repeatedFace)[sum])
                                {
                                    adjwgt.push_back(0);
                                }
                            else
                                adjwgt.push_back(10);
                        }
                }
        }
        // this is the list of "keys" to access M_jadj
        // graph element i has neighbors M_jadj[ k ],
        // with M_iadj[i] <= k < M_iadj[i+1]
        M_iadj.push_back(sum);
    }

#endif

    //MPI_Barrier(MPI_COMM_WORLD);

    // **************
    // parMetis part

    // this array is to be used for weighted nodes on the graph:
    // usually we will set it to NULL

    int*   vwgt    = 0;

    int    wgtflag;
    if(interfaceMap)
        wgtflag=1;
    else
        wgtflag=0;
    int    ncon    = 1;
    int    numflag = 0;

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
    int* adjwgtPtr(0);
    if (adjwgt.size() > 0) adjwgtPtr = (int*) &adjwgt[0];
    ParMETIS_V3_PartKway((int*) &M_vertexDist[0],
                         (int*) &M_iadj[0],
                         (int*) &M_jadj[0],
                         vwgt,  adjwgtPtr,
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

    // cycling on locally stored vertices
    for (UInt ii = 0; ii < part.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        locProc[part[ii]].push_back(ii + vertexDist()[M_me]);
    }

    //////////////// BEGIN OF SOLID PARTITION PART ////////////////

    if(interfaceMap)
        {
            int procOrder[nprocs];
            std::vector<std::vector<UInt> > mymatchesForProc(nprocs);
            std::vector<std::vector<UInt> > matchesForProc(nprocs);
            bool orderingError[nprocs];
            for(int i=0; i<nprocs ; ++i)
                {
                    orderingError[i]=false;
                    for(int j=0; j<nprocs ; ++j)
                        {
                            mymatchesForProc[i].push_back(0);
                            matchesForProc[i].push_back(0);
                        }
                }
            for(UInt kk=0; kk<part.size(); ++kk)
                {
                    if((*isOnProc)[kk+M_vertexDist[M_me]]!=-1)
                        ++mymatchesForProc[part[kk]][(*isOnProc)[kk+M_vertexDist[M_me]]];
                }
            for(UInt j=0; (int)j<nprocs; ++j)
                MPI_Allreduce( &mymatchesForProc[j][0], &matchesForProc[j][0], nprocs, MPI_INT, MPI_SUM, MPIcomm);
            M_comm->Barrier();
//             for(UInt i=0; i<nprocs ; ++i)
//                 for(UInt j=0; j<nprocs ; ++j)
//                 {
//                     std::cout<<M_me<<" matchesForProc ["<<j<<"] ["<< i << "] "<<matchesForProc[j][i]<<std::endl;
//                     std::cout<<M_me<<" mymatchesForProc ["<<j<<"] ["<< i << "] "<<mymatchesForProc[j][i]<<std::endl;
//                 }

            int suitableProcess=-1;
            UInt max=0;

            for(int ii=0; ii<nprocs; ++ii)
                {
                    if(matchesForProc[M_me][ii]>max)
                        {
                            suitableProcess=ii;
                            max=matchesForProc[M_me][ii];
                        }
                }


// // TO ASSIGN A BAD PARTITON ==> to keep commented

//             max=1000;
//             for(int ii=0; ii<nprocs; ++ii)
//                 {
//                     if(matchesForProc[M_me][ii]<max)
//                         {
//                             suitableProcess=ii;
//                             max=matchesForProc[M_me][ii];
//                         }
//                 }
//  // END OF BAD PARTITION COMMENTED PART

            //            int procOrder[nprocs];
            ASSERT(suitableProcess!=-1, "one partition is without interface nodes!");
            procOrder[M_me]= suitableProcess;

//             std::cout<<M_me<<" procorder(me) "<<procOrder[M_me]<<std::endl;
//             std::cout<<M_me<<" suitableProc "<<suitableProcess<<std::endl;
            M_comm->Barrier();

            std::vector<UInt> maxs(nprocs);
            maxs[M_me]=max;
            for(int j=0; j<nprocs ; ++j)//Allgather
                MPI_Bcast( &maxs[j], 1, MPI_INT, j, MPIcomm);//perhaps generates errors

            std::vector<pair<UInt, int> > procIndex(nprocs);
            for(int k=0; k<nprocs; ++k)
                {
                    procIndex[k]=std::make_pair( maxs[k], k);
                }
/*
struct booleanCondition
{
    bool reordering( std::pair<UInt, int>& a, std::pair<UInt, int>& b)
    {
        return a.first>b.first;
    }
};
*/

// std::vector<int> procIndex(nprocs);
            //std::vector<int>::const_iterator Iter1;
            //bool (*function)( int, int );
            //function = condition.condition;
            //boost::shared_ptr<bool(int, int)> fctprtr;
            //fctprtr.reset(&condition.condition);
            //bool (booleanCondition::*function)(int, int)=&booleanCondition::condition;
            //booleanCondition condition(maxs);

// transformation=std::sort(maxs.begin(), maxs.end());
 std::sort(procIndex.begin(), procIndex.end()/*, &booleanCondition::reordering*/);
 for(int l=0;l<nprocs;++l)
     //     std::cout<<M_me<< " :reordered maxs: " <<procIndex[l].first<<std::endl;
 for(int l=0;l<nprocs;++l)
     //     std::cout<<M_me<< " :reordered indices: " <<procIndex[l].second<<std::endl;
            //std::vector newMaxs(nprocs);
            //            max=0;
            // for(int k=0; k<nprocs; ++k)
//             for(int j=0; j<nprocs ; ++j)
//                 {
//                     if(maxs(j)>max)
//                         {
//                             max=maxs(j);
//                             //newMaxs[k]=max;
//                             loopIndex[k]=j;
//                         }
//                 }
            for(int j=0; j<nprocs ; ++j)//Allgather
                MPI_Bcast( &procOrder[j], 1, MPI_INT, j, MPIcomm);//perhaps generates errors
//             for(int j=0; j<nprocs ; ++j)
//                 {
//                     std::cout<<"me: "<<M_me<<" j: "<< j <<" procorder2 "<<procOrder[j]<<std::endl;
//                }
            std::vector< std::vector<int> > locProc2(locProc);
            for(int j=nprocs; j>0 ; --j)
                {
                    //if(procOrder[nprocs-1 - procIndex[j].second]>nprocs-1 - procIndex[j].second)
                    //if(orderingError[ procIndex[j-1].second]==false)
                    //locProc[procOrder[ procIndex[j-1].second]].swap(locProc[procIndex[j-1].second]);

                    if(orderingError[procOrder[procIndex[j-1].second]]==false)
                        locProc[procOrder[ procIndex[j-1].second]]=  locProc2[procIndex[j-1].second];
//                     if(j!=procOrder[j])
//                         {
//                             std::vector<int> tmp(locProc[j]);
//                             locProc[j]=locProc[procOrder[j]];
//                             locProc[procOrder[j]]=tmp;
//                         }
                    else
                        {
                            std::cout<<"Ordering error when assigning the processor"<<M_me<<" to the partition,"<<std::endl<<" parmetis did a bad job."<<std::endl;
                            // j not assigned
                            for(int i=nprocs; i>0 ; --i)
                                {
                                    if(orderingError[procIndex[i-1].second]==false)//means that i is the first proc not assigned
                                        {
                                            //locProc[procOrder[procIndex[j-1].second]].swap(locProc[/*nprocs-1 - procIndex[*/i/*].second*/]);
                                            procOrder[procIndex[j-1].second]=procIndex[i-1].second;
                                            locProc[procIndex[i-1].second]=locProc2[procIndex[j-1].second];
                                            break;
                                        }
                                }

                        }
                    orderingError[procOrder[procIndex[j-1].second]]=true;
                }
//             for(int j=0; j<nprocs ; ++j)
//                 {
//                     std::cout<<"me: "<<M_me<<" j: "<< j <<" procorder3 "<<procOrder[j]<<std::endl;
//                 }
        }
    ////////////////// END OF SOLID PARTITION PART /////////////////////

    // sizeof(MPI_INT) is expressed in bytes
    // MPI_INT range is ( -2^(7*sizeof(MPI_INT), 2^(7*sizeof(MPI_INT) - 1 )
    // (MPI_INT is not unsigned int)
    //int _exp(7*sizeof(MPI_INT));
    //int _max_int( std::ldexp(1., _exp-1) ); // taking a little less: 2^(exp-1) instead of 2^exp - 1
    int max_int (1000);

    int ssize[nProc];
    int rsize[nProc];
    // cycling on subdomains
    // TODO: Matteo please comment this part :)

    MPI_Status status;
    int size;

    MPI_Status  recv_status; //,  send_status;
    MPI_Request /*recv_request,*/ send_request;


    for (int iproc = 0; iproc < nProc; ++iproc)
    {
        // all processes other than me are sending vertices
        // belonging to my subdomain
        if (int(M_me) != iproc)
            {                                                  //start if
                size = locProc[iproc].size();

                // tell me how many vertices belonging to me you have to send me
                MPI_Isend(&size, 1, MPI_INT, iproc, 10, MPIcomm, &send_request);
                ssize[iproc] = size;

                //MPI_Wait(&send_request, &send_status);
                //            MPI_Wait(&recv_request,&recv_status);


            }
        else
            {
                for (int jproc = 0; jproc < nProc; ++jproc)
                    {
                        if ((int)M_me !=jproc)
                            {
                                MPI_Recv(&size, 1, MPI_INT, jproc, 10, MPIcomm, &recv_status);
                                rsize[jproc] = size;
                            }
                    }
            }
    }


    for (int iproc = 0; iproc < nProc; ++iproc)
        {

            if ((int)M_me != iproc)
                {
                    size = ssize[iproc];
                    // workaround for huge data to be passed
                    if (size > max_int)
                        {
                            int incr = 1 ;
                            int pos = 0;
                            int size_part =size;

                            // divide the whole data set into smaller packets
                            while (size_part > max_int )
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
                        if (size != 0)
                            MPI_Send(&locProc[iproc][0], size, MPI_INT, iproc, 60, MPIcomm);

//                 }
//         }

                }
            else
                {
                    for (int jproc = 0; jproc < nProc; ++jproc)
                        {
                            if (jproc != iproc)
                                {
                                    size = rsize[jproc];
                                    std::vector<int> stack(size, 0);

                                    if (size > max_int )
                                        {
                                            int size_part, pos, incr;

                                            MPI_Recv(&incr, 1, MPI_INT, jproc, 20, MPIcomm, &status);
                                            MPI_Recv(&size_part, 1, MPI_INT, jproc, 30, MPIcomm, &status);

                                            for ( int kk = 0; kk < incr; ++kk)
                                                {
                                                    MPI_Recv(&pos, 1, MPI_INT, jproc, 100+kk, MPIcomm, &status);
                                                    MPI_Recv(&stack[pos], size_part , MPI_INT, jproc, 5000000+kk, MPIcomm, &status);
                                                }
                                            int resto=0;
                                            MPI_Recv(&resto, 1, MPI_INT, jproc, 80, MPIcomm, &status);

                                            if(resto!=0)
                                                {
                                                    MPI_Recv(&pos, 1, MPI_INT, jproc, 40, MPIcomm, &status);
                                                    MPI_Recv(&stack[pos],  resto, MPI_INT, jproc, 50, MPIcomm, &status);
                                                }
                                        }
                                    else
                                        {
                                            if (size != 0)
                                                MPI_Recv(&stack[0], size , MPI_INT, jproc, 60, MPIcomm, &status);
                                        }
                                    for (int jj = 0; jj < size; ++jj)
                                        {
                                            locProc[M_me].push_back(stack[jj]);
                                        }

                                }
                        }
                }//end if
            //M_comm->Barrier();
        }                 //end for

#ifdef DEBUG
    Debug(4000) << M_me << " has " << locProc[M_me].size() << " elements.\n";
#endif
    // ***********************
    // local mesh construction
    // ***********************

    if (!M_me) std::cout << "Building local mesh ...                        " << std::flush;

    std::map<int, int>::iterator  im;
    std::set<int>::iterator       is;

    int count = 1;
    UInt ielem;
    UInt inode;

    // cycle on local element's ID
#ifdef TWODIM

    for (UInt jj = 0; jj < locProc[M_me].size(); ++jj)
    {
    	ielem = locProc[M_me][jj];
        M_localFaces.push_back(ielem);

        // cycle on element's nodes
        for (UInt ii = 1; ii <= elementNodes; ++ii)
        {
        	inode = _mesh.face(ielem + 1).point(ii).id();
        	im    = M_globalToLocalNode.find(inode);

        	// if the node is not yet present in the list of local nodes, then add it
        	// CAREFUL: also local numbering starts from 1 in RegionMesh
        	if ( im == M_globalToLocalNode.end() )
        	{
        		M_globalToLocalNode.insert(std::make_pair(inode, count));
        		++count;
        		// store here the global numbering of the node
        		M_localNodes.push_back(_mesh.face(ielem + 1).point(ii).id());
			}
		}

        for (UInt ii = 1; ii <= elementEdges; ++ii)
        	// store here the global numbering of the edge
        	M_localEdges.insert(_mesh.localEdgeId(ielem + 1, ii));

    }
    std::vector<int>::iterator it;

    // ******************
    // nodes construction
    // ******************

    UInt nBoundaryPoint = 0;
    M_mesh->pointList.reserve(M_localNodes.size());

    // guessing how many boundary points on this processor.
    M_mesh->_bPoints.reserve(M_localNodes.size()); // guessing how many boundary points on this processor.

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

    	UInt id = _mesh.point(*it).id();

    	pp->setId(id);
    	pp->setLocalId( inode );
    	M_mesh->localToGlobalNode().insert(std::make_pair(inode, id));
    	M_mesh->globalToLocalNode().insert(std::make_pair(id, inode));

    }

    // ******************
    // faces construction
    // ******************

    count = 1;

    typename Mesh::FaceType * pf = 0;
    M_mesh->faceList.reserve(M_localFaces.size());

    // loop in the list of local elements
    // CAREFUL! in this loop inode is the global numbering of the points
    // We insert the local numbering of the nodes in the local face list
    for (it = M_localFaces.begin(); it != M_localFaces.end(); ++it, ++count)
    {
    	pf = &M_mesh->addFace();
    	// CAREFUL! in ParMETIS data structures, numbering starts from 0
    	pf->setId     ( _mesh.face(*it + 1).id() );
    	pf->setLocalId( count );

    	M_globalToLocalFace.insert(make_pair( _mesh.face(*it + 1).id(),  count ));

    	for (ID id = 1; id <= elementNodes; ++id)
    	{
    		inode = _mesh.face(*it + 1).point(id).id();

    		// im is an iterator to a map element
    		// im->first is the key (i. e. the global ID "inode")
    		// im->second is the value (i. e. the local ID "count")
    		im    = M_globalToLocalNode.find(inode);
    		pf->setPoint( id, M_mesh->pointList( (*im).second ) );
		}

    	int ibc = _mesh.face(*it + 1).marker();

    	pf->setMarker( EntityFlag( ibc ) );
    }


    // ******************
    // edges construction
    // ******************

    typename Mesh::EdgeType * pe = 0;

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
        	// create a boundary edge in the local mesh, if needed
        	++nBoundaryEdges;
        }

        pe =  &M_mesh->addEdge( boundary );
        pe->setId     (_mesh.edge(*is).id());
        pe->setLocalId( count );

        int elem1 =  _mesh.edge(*is).ad_first();
        int elem2 =  _mesh.edge(*is).ad_second();

       // find the mesh elements adjacent to the edge
        im =  M_globalToLocalFace.find(elem1);

        int localElem1;

        if (im == M_globalToLocalFace.end())
        	localElem1 = 0;
        else
        	localElem1 = (*im).second;

		im =  M_globalToLocalFace.find(elem2);

		int localElem2;

		if (im == M_globalToLocalFace.end())
			localElem2 = 0;
		else
			localElem2 = (*im).second;

        // if this process does not own either of the adjacent elements
		// then the two adjacent elements and the respective edge positions coincide in the local mesh
		if ((localElem1 == 0) && ! boundary)
		{
			pe->ad_first() = localElem2;
			pe->pos_first()=_mesh.edge(*is).pos_second();
        }
		else
		{
			pe->ad_first() = localElem1;
			pe->pos_first() = _mesh.edge(*is).pos_first();
        }

		if ((localElem2 == 0) && ! boundary)
		{
			pe->ad_second() = localElem1;
			pe->pos_second()=_mesh.edge(*is).pos_first();
		}
		else
		{
			pe->ad_second()  = localElem2;
			pe->pos_second() = _mesh.edge(*is).pos_second();
		}


		for (ID id = 1; id <= _mesh.edge(*is).numLocalVertices; ++id)
		{
			inode =_mesh.edge(*is).point(id).id();
			im = M_globalToLocalNode.find(inode);
			pe->setPoint(id, M_mesh->pointList( (*im).second ));
		}
		pe->setMarker( _mesh.edge(*is).marker());

		M_mesh->setLinkSwitch( "HAS_ALL_EDGES" );
		M_mesh->setLinkSwitch( "EDGES_HAVE_ADIACENCY" );
	}

#elif defined THREEDIM

    for (UInt jj = 0; jj < locProc[M_me].size(); ++jj)
    {
        ielem = locProc[M_me][jj];
        M_localVolumes.push_back(ielem);

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
        	// store here the global numbering of the edge
          	M_localEdges.insert(_mesh.localEdgeId(ielem + 1, ii));


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

    M_mesh->volumeList.reserve(M_localVolumes.size());

    // loop in the list of local elements
    // CAREFUL! in this loop inode is the global numbering of the points
    // We insert the local numbering of the nodes in the local volume list
    for (it = M_localVolumes.begin(); it != M_localVolumes.end(); ++it, ++count)
    {
        pv = &M_mesh->addVolume();
//        std::cout << "volume " << _mesh.volume(*it + 1).id() << std::flush;
        // CAREFUL! in ParMETIS data structures, numbering starts from 0
        pv->setId     ( _mesh.volume(*it + 1).id() );
        pv->setLocalId( count );

        M_globalToLocalVolume.insert(make_pair( _mesh.volume(*it + 1).id(),  count ));

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
    for (UInt ii = 0; ii < M_localVolumes.size(); ++ii)
    {
        for (UInt jj = 0; jj < 4; ++jj)
            std::cout << M_mesh->volume(ii + 1).point(jj + 1).id() << " " ;
        std::cout << std::endl;
    }

    for (UInt ii = 0; ii < M_localVolumes.size(); ++ii)
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
      		// create a boundary edge in the local mesh, if needed
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
        im =  M_globalToLocalVolume.find(elem1);

        int localElem1;

        if (im == M_globalToLocalVolume.end()) localElem1 = 0;
        else
            localElem1 = (*im).second;

        im =  M_globalToLocalVolume.find(elem2);

        int localElem2;
        if (im == M_globalToLocalVolume.end()) localElem2 = 0;
        else
            localElem2 = (*im).second;


        // if ((localElem1 == 0) && ! boundary) localElem1 = localElem2;
        // if ((localElem2 == 0) && ! boundary) localElem2 = localElem1;

        // pf->ad_first() = localElem1;
        // pf->pos_first() = _mesh.face(*is).pos_first();

        // pf->ad_second()  = localElem2;
        // pf->pos_second() = _mesh.face(*is).pos_second();

        // if this process does not own either of the adjacent elements
        // then the two adjacent elements and the respective face positions coincide in the local mesh
        //possible bug fixed: not only the two adjacent elements face, but also the face positions should coincide.
        //otherwise it could happen that a pair(element, position) is associated to different faces.
        //This can lead to a wrong treatment of the dofPerFace (in 2D of the dofPerEdge, as occurred with P2)
        //
        if ((localElem1 == 0) && ! boundary)
        {
        	pf->ad_first() = localElem2;
        	pf->pos_first()=_mesh.face(*is).pos_second();
        }
        else
        {
        	pf->ad_first() = localElem1;
        	pf->pos_first() = _mesh.face(*is).pos_first();
        }

        if ((localElem2 == 0) && ! boundary)
        {
        	pf->ad_second() = localElem1;
        	pf->pos_second()=_mesh.face(*is).pos_first();
        }
        else
        {
        	pf->ad_second()  = localElem2;
        	pf->pos_second() = _mesh.face(*is).pos_second();
        }


        for (ID id = 1; id <= _mesh.face(*is).numLocalVertices; ++id)
        {
        	inode =_mesh.face(*is).point(id).id();
        	im = M_globalToLocalNode.find(inode);
        	pf->setPoint(id, M_mesh->pointList( (*im).second ));
		}

        pf->setMarker( _mesh.face(*is).marker());

        M_mesh->setLinkSwitch( "HAS_ALL_FACES" );
        M_mesh->setLinkSwitch( "FACES_HAVE_ADIACENCY" );

    }

#endif

    // ******************
    // final setup
    // ******************


    UInt nVolumes = M_localVolumes.size();
    UInt nNodes   = M_localNodes.size();
    UInt nEdges   = M_localEdges.size();
    UInt nFaces   = M_localFaces.size();

    M_mesh->setMaxNumPoints ( nNodes, true );
    M_mesh->setMaxNumEdges  ( nEdges, true );
    M_mesh->setMaxNumFaces  ( nFaces, true );
#ifndef TWODIM
    M_mesh->setMaxNumVolumes( nVolumes, true);
#endif

    M_mesh->setMaxNumGlobalPoints ( _mesh.numPoints  () );
    M_mesh->setNumGlobalVertices  ( _mesh.numPoints  () );
    M_mesh->setMaxNumGlobalEdges  ( _mesh.numEdges   () );
    M_mesh->setMaxNumGlobalFaces  ( _mesh.numFaces   () );
#ifndef TWODIM
    M_mesh->setMaxNumGlobalVolumes( _mesh.numVolumes () );
    M_mesh->setNumBFaces    ( nBoundaryFace  );
#endif
    M_mesh->setNumBPoints   ( nBoundaryPoint );
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

#ifndef TWODIM
    M_mesh->updateElementFaces();
#endif

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

    if (!M_me) std::cout << "done" << std::endl;
    if (!M_me) std::cout << "Creating the map ...                           " << std::flush;

    // *********************
    // repeated map creation
    // *********************

    std::vector<int> elementList = locProc[M_me];

    // repeated element map creation

#ifdef TWODIM

    // use sets to store each entity only once
    std::set<int>    repeatedNodeList;
    std::set<int>    repeatedEdgeList;
    for (UInt ii = 0; ii < elementList.size(); ++ii)
     {
         ielem = elementList[ii];
         M_repeatedFaceVector.push_back(ielem + 1);
         for (UInt jj = 1; jj <= elementNodes; ++jj)
         {
             inode = _mesh.element(ielem + 1).point(jj).id();
             repeatedNodeList.insert(inode);
         }
         for (UInt jj = 1; jj <= elementEdges; ++jj)
         {
        	 UInt iedge = _mesh.localEdgeId(ielem + 1, jj);
        	 repeatedEdgeList.insert(iedge);
         }
     }

     // repeated node map creation

     M_repeatedNodeVector.reserve(repeatedNodeList.size());

     for (is = repeatedNodeList.begin(); is != repeatedNodeList.end(); ++is)
     {
         M_repeatedNodeVector.push_back(*is);
     }

     // repeated edge list creation

     M_repeatedEdgeVector.reserve(repeatedEdgeList.size());

     for (is = repeatedEdgeList.begin(); is != repeatedEdgeList.end(); ++is)
     {
         M_repeatedEdgeVector.push_back(*is);
     }

#elif defined THREEDIM

    // use sets to store each entity only once
    std::set<int>    repeatedNodeList;
    std::set<int>    repeatedEdgeList;
    std::set<int>    repeatedFaceList;

    for (UInt ii = 0; ii < elementList.size(); ++ii)
    {
        ielem = elementList[ii];
        M_repeatedVolumeVector.push_back(ielem + 1);
        for (UInt jj = 1; jj <= elementNodes; ++jj)
        {
        	inode = _mesh.volume(ielem + 1).point(jj).id();
        	repeatedNodeList.insert(inode);
        }
        for (UInt jj = 1; jj <= elementEdges; ++jj)
        {
             UInt iedge = _mesh.localEdgeId(ielem + 1, jj);
             repeatedEdgeList.insert((int) iedge);
        }
        for (UInt jj = 1; jj <= elementFaces; ++jj)
        {
            UInt iface = _mesh.localFaceId(ielem + 1, jj);
            repeatedFaceList.insert(iface);
        }
    }


    // repeated node map creation

    M_repeatedNodeVector.reserve(repeatedNodeList.size());

    for (is = repeatedNodeList.begin(); is != repeatedNodeList.end(); ++is)
    {
        M_repeatedNodeVector.push_back(*is);
    }

    // repeated edge list creation
    M_repeatedEdgeVector.reserve(repeatedEdgeList.size());

    for (is = repeatedEdgeList.begin(); is != repeatedEdgeList.end(); ++is)
           M_repeatedEdgeVector.push_back(*is);

    // repeated face list creation
    M_repeatedFaceVector.reserve(repeatedFaceList.size());

    for (is = repeatedFaceList.begin(); is != repeatedFaceList.end(); ++is)
    {
        M_repeatedFaceVector.push_back(*is);
    }
#endif

    if (!M_me) std::cout << "done" << std::endl;
}


template<typename Mesh>
partitionMesh<Mesh>::~partitionMesh()
{
}


}
#endif
