/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*- */
//@HEADER

/*!
  \brief short description of file
  \file partitionMesh.hpp
  \author Gilles Fourestey gilles.fourestey@epfl.ch (Modified by Radu Popescu radu.popescu@epfl.ch)
*/

#ifndef PARTMESH
#define PARTMESH 1

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
  \brief short description of class
  \author Gilles Fourestey gilles.fourestey@epfl.ch (Modified by Radu Popescu radu.popescu@epfl.ch)

  This class implements the partitioning of a mesh using (par)Metis.
*/
template<typename Mesh>
class partitionMesh
{
public:
    //@{
    typedef boost::shared_ptr<Mesh> mesh_type;
    //@}
    //! \name Constructors & Destructors
    //@{
    //! Constructor
    /*!
      This is a "dummy" constructor. It takes as parameters the
      unpartitioned mesh (reference), the Epetra_Comm object in
      use (reference) and pointers to the Epetra interface and
      repeated interface maps. The constructor initializes the
      M_comm Epetra_Comm data member and calls a private method
      partitionMesh::execute which handles the mesh partitioning.
      \param _mesh - Mesh& - the unpartitioned mesh
      \param _comm - Epetra_Comm& - Epetra communicator object
      \param interfaceMap - Epetra_Map*
      \param interfaceMapRep - Epetra_Map*
     */
    partitionMesh(Mesh &_mesh, Epetra_Comm &_comm, Epetra_Map* interfaceMap = 0,
                  Epetra_Map* interfaceMapRep = 0);
    //! Empty destructor
    ~partitionMesh();
    //@}
    //! \name Get Methods
    //@{
    //! Return a reference to M_vertexDist
    const std::vector<int>& vertexDist()           const {return M_vertexDist;};
    //! Return a pointer to M_mesh
    const mesh_type         mesh()                 const {return M_mesh;}
    //! Return a reference to M_repeatedNodeVector
    const std::vector<int>& repeatedNodeVector()   const {return M_repeatedNodeVector;}
    //! Return a reference to M_repeatedEdgeVector
    const std::vector<int>& repeatedEdgeVector()   const {return M_repeatedEdgeVector;}
    //! Return a reference to M_repeatedFaceVector
    const std::vector<int>& repeatedFaceVector()   const {return M_repeatedFaceVector;}
    //! Return a reference to M_repeatedVolumeVector
    const std::vector<int>& repeatedVolumeVector() const {return M_repeatedVolumeVector;}
    //@}
private:
    //! Private Methods
    //@{
    //! Execute mesh partitioning using the configured MPI processes
    /*! 
      This function takes exactly the same arguments as the constructor. The body
      of the constructor has been moved inside this function. The is the first step
      towards having an empty constructor. Right now this doesn't change the
      interface of the class, but next step would be to make an empty constructor,
      make this function public and use it directly to partition the mesh.
      Sets current mesh element parameters: M_elementNodes, M_elementEdges,
      M_elementFaces, M_faceNodes
      Updates: M_locProc (indirectly)
      Other data members are changed indirectly by calling other private methods.
      \param _mesh - Mesh& - the unpartitioned mesh
      \param _comm - Epetra_Comm& - Epetra communicator object
      \param interfaceMap - Epetra_Map*
      \param interfaceMapRep - Epetra_Map*
    */
    void execute(Mesh &_mesh, Epetra_Comm &_comm, Epetra_Map* interfaceMap = 0,
                 Epetra_Map* interfaceMapRep = 0);
    //! Sets the element parameters according to the type of mesh element used.
    /*!
      Sets element parameters (nodes, faces, edges and number of nodes on each
      face according to the type of mesh element used (Mesh::ElementShape::Shape).
      Updates M_elementNodes, M_elementFaces, M_elementEdges, M_faceNodes.
     */
    void setElementParameters();
    //! Build the graph vertex distribution vector
    /*!
      Updates the member M_vertexDist according to the number of processors to be 
      used by ParMETIS (the number of processes started for MPI
      \param numElements - UInt - number of elements in the mesh
     */
    void distributeElements(UInt numElements);
    //! Find faces on the boundaries between domains (FSI)
    /*!
      Identifies the element faces that are common to both the fluid and the solid 
      meshes and creates a map between the faces that reside on the boundary between 
      the two meshes and the processors used. Updates the members M_repeatedFace and 
      M_isOnProc.
      \param _mesh - Mesh& - the unpartitioned mesh
      \param interfaceMap - Epetra_Map*
      \param interfaceMapRep - Epetra_Map*
     */
    void findRepeatedFacesFSI(Mesh &_mesh, Epetra_Map *interfaceMap,
                              Epetra_Map *interfaceMapRep);
    //! Partition the connectivity graph using ParMETIS
    /*!
      Partitions the connectivity graph using ParMETIS. The result is stored in the
      member M_part: this is a vector of integer values; its size is the number of
      elements in the unpartitioned mesh. Each value represents the number of the
      partition to which the element was assigned. Also creates M_locProc, the
      vector of elements in each subdomain.
      Updates: M_part, M_locProc
      \param _mesh - Mesh& - the unpartitioned mesh
      \param interfaceMap - Epetra_Map*
     */
    void partitionConnectivityGraph(Mesh &_mesh,
                                    Epetra_Map* interfaceMap);
    //! Updates the map between elements and processors in FSI
    /*!
      Updates M_locProc during FSI modeling.
     */
    void matchFluidPartitionsFSI();
    //! Redistribute elements among processes
    /*!
      Redistributes elements among processes, when needed, after the connectivity
      graph partitioning phase. Updates M_locProc
     */
    void redistributeElements();
    //! Construct local mesh
    /*!
      Constructs the data structures for the local mesh partition.
      Updates M_localNodes, M_localEdges, M_localFaces, M_localVolumes,
      M_globalToLocalNode.
      \param _mesh - Mesh& - the unpartitioned mesh
    */
    void constructLocalMesh(Mesh &_mesh);
    //! Construct nodes
    /*!
      Adds nodes to the partitioned mesh object. Updates M_nBoundaryPoints,
      M_mesh.
      \param _mesh - Mesh& - the unpartitioned mesh
     */
    void constructNodes(Mesh &_mesh);
    //! Construct volumes
    /*!
      Adds volumes to the partitioned mesh object. Updates M_globalToLocalVolume,
      M_mesh.
      \param _mesh - Mesh& - the unpartitioned mesh
     */
    void constructVolumes(Mesh &_mesh);
    //! Construct edges
    /*!
      Adds edges to the partitioned mesh object. Updates M_nBoundaryEdges,
      M_mesh.
      \param _mesh - Mesh& - the unpartitioned mesh
     */
    void constructEdges(Mesh &_mesh);
    //! Construct faces
    /*!
      Adds faces to the partitioned mesh object. Updates M_nBoundaryFaces,
      M_mesh.
      \param _mesh - Mesh& - the unpartitioned mesh
     */
    void constructFaces(Mesh &_mesh);
    //! Final setup of local mesh
    /*!
      Updates the partitioned mesh object data members after adding the mesh
      elements (nodes, edges, faces, volumes).
      Updates M_mesh.
      \param _mesh - Mesh& - the unpartitioned mesh
     */
    void finalSetup(Mesh &_mesh);
    //! Create repeated element map
    /*!
      Creates a map of the boundary elements (nodes, edges, faces, volumes).
      Updates M_repeatedNodeVector, M_repeatedEdgeVector, M_repeatedFaceVector,
      M_repeatedVolumeVector.
      \param _mesh - Mesh& - the unpartitioned mesh
     */
    void createRepeatedMap(Mesh &_mesh);
    //@}
    //! Private Data Members
    //@{
    mesh_type               M_mesh;
    std::vector<int>        M_vertexDist;
    std::vector<int>        M_iadj;
    std::vector<int>        M_jadj;
    std::vector<int>        M_localNodes;
    std::set<int>           M_localEdges;
    std::set<int>           M_localFaces;
    std::vector<int>        M_localVolumes;
    std::vector<int>        M_repeatedNodeVector;
    std::vector<int>        M_repeatedEdgeVector;
    std::vector<int>        M_repeatedFaceVector;
    std::vector<int>        M_repeatedVolumeVector;
    std::map<int, int>      M_globalToLocalNode;
    std::map<int, int>      M_globalToLocalEdge;
    std::map<int, int>      M_globalToLocalFace;
    std::map<int, int>      M_globalToLocalVolume;
    Epetra_Comm*            M_comm;
    UInt                    M_me;
    // The following are utility variables used throughout the partitioning
    // process
    UInt                                 M_elementNodes;
    UInt                                 M_elementFaces;
    UInt                                 M_elementEdges;
    UInt                                 M_faceNodes;
    boost::shared_ptr<std::vector<int> > M_repeatedFace;
    boost::shared_ptr<std::vector<int> > M_isOnProc;
    std::vector<int>                     M_part;
    std::vector<std::vector<int> >       M_locProc;
    UInt                                 M_nBoundaryPoints;
    UInt                                 M_nBoundaryEdges;
    UInt                                 M_nBoundaryFaces;
    //@}
}; // class partitionMesh



//
// IMPLEMENTATION
//

template<typename Mesh>
partitionMesh<Mesh>::partitionMesh(Mesh &_mesh, Epetra_Comm &_comm,
                                   Epetra_Map* interfaceMap,
                                   Epetra_Map* interfaceMapRep):
    M_mesh (new Mesh),
    M_comm (&_comm)
{
    execute(_mesh, _comm, interfaceMap, interfaceMapRep);
}

template<typename Mesh>
void partitionMesh<Mesh>::setElementParameters()
{
    switch (Mesh::ElementShape::Shape)
    {
    case HEXA:
        M_elementNodes = 8;
        M_elementFaces = 6;
        M_elementEdges = 12;
        M_faceNodes    = 4;
        break;
    case TETRA:
        M_elementNodes = 4;
        M_elementFaces = 4;
        M_elementEdges = 6;
        M_faceNodes    = 3;
        break;
    default:
        ERROR_MSG( "Face Shape not implemented in partitionMesh" );
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::distributeElements(UInt numElements)
{
    // ParMETIS is able to work in parallel: how many processors does it have at hand?
    int numProcessors = M_comm->NumProc();
    M_me              = M_comm->MyPID();

    // CAREFUL: ParMetis works on a graph abstraction.
    // A graph is built over the data structure to be split, each vertex being a mesh element
    // so hereby a "vertex" is actually a _graph_ vertex, i. e. a mesh element
    M_vertexDist.resize(numProcessors + 1);
    M_vertexDist[0] = 0;

    UInt k = numElements;

    // Evenly distributed graph vertices
    for (int i = 0; i < numProcessors; ++i)
    {
        UInt l = k / (numProcessors - i);
        M_vertexDist[i + 1] = M_vertexDist[i] + l;
        k -= l;
    }
    ASSERT(k == 0, "At this point we should have 0 volumes left") ;
}

template<typename Mesh>
void partitionMesh<Mesh>::findRepeatedFacesFSI(Mesh &_mesh, Epetra_Map *interfaceMap, 
                                               Epetra_Map *interfaceMapRep)
{
    std::vector<int>                     myRepeatedFace; // used for the solid partitioning
    boost::shared_ptr<std::vector<int> > myIsOnProc;     // used for the solid partitioning

    myIsOnProc.reset(new std::vector<int>(_mesh.numVolumes()));
    
    bool myFaceRep;
    bool myFace(false);
    short count;
    for(UInt h = 0; h < _mesh.numVolumes(); ++h)
    {
        (*myIsOnProc)[h] = -1;
    }

    // This loop is throughout the whole unpartitioned mesh,
    // it is expensive and not scalable.
    // Bad, this part should be done offline

    for (UInt ie = 1; ie <= _mesh.numVolumes(); ++ie)
    {
        for (UInt iface = 1; iface <= M_elementFaces; ++iface)
        {
            UInt face = _mesh.localFaceId(ie, iface);
            UInt vol  = _mesh.face(face).ad_first();
            if (vol == ie)
            {
                vol = _mesh.face(face).ad_second();
            }
            if (vol != 0)
            {
                myFace = false;
                myFaceRep = false;
                count = 0;                    
                for(int ipoint = 1; ipoint <= (int) M_faceNodes; ++ipoint) // vertex-based dofs
                {
                    myFaceRep = ((interfaceMap->LID(_mesh.face(face).point(ipoint).id())
                                  /* first is fluid */ == -1) && 
                                 (interfaceMapRep->LID(_mesh.face(face).point(ipoint).id())
                                  /* first is fluid */ != -1));
                    myFace = myFace || (interfaceMap->LID(_mesh.face(face).point(ipoint).id()) != -1);
                    if (myFaceRep)
                    {
                        ++count;
                    }
                }
                if (count > 1)
                {
                    myRepeatedFace.push_back(1);
                }
                else
                {
                    myRepeatedFace.push_back(0);
                }
            }
            if (myFace)
            {
                (*myIsOnProc)[ie-1] = M_me;
            }
        }
    }

    M_repeatedFace.reset(new std::vector<int> (myRepeatedFace.size()));
    M_isOnProc.reset(new std::vector<int> (*myIsOnProc));

    // Lot of communication here!!
    MPI_Allreduce(&myRepeatedFace[0], &(*M_repeatedFace)[0], myRepeatedFace.size(), 
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&(*myIsOnProc)[0], &(*M_isOnProc)[0], myIsOnProc->size(), 
                  MPI_INT, MPI_MAX, MPI_COMM_WORLD);
}

template<typename Mesh>
void partitionMesh<Mesh>::partitionConnectivityGraph(Mesh &_mesh,
                                                     Epetra_Map* interfaceMap)
{
    // This array's size is equal to the number of locally-stored vertices:
    // at the end of the partitioning process, "M_part" will contain the partitioning array:
    // M_part[m] = n; means that graph vertex m belongs to subdomain n
    M_part.resize(M_vertexDist[M_me + 1] - M_vertexDist[M_me]);

    // Now each processor will take care of its own graph vertices (i. e. mesh elements).
    // Nothing guarantees about the neighbor elements distribution across the processors,
    // since as of now we just split the set of volumes based on IDs.
    // Here we building up the neighbor arrays.

    UInt localStart = M_vertexDist[M_me] + 1;
    UInt localEnd   = M_vertexDist[M_me + 1] + 1;

    // this vector contains the weights for the edges of the graph,
    // it is set to null if it is not used.
    std::vector<int> adjwgt;

    M_iadj.resize(0);
    M_iadj.push_back(0);

    UInt sum = 0;

    for (UInt ie = localStart; ie < localEnd; ++ie)
    {
        for (UInt iface = 1; iface <= M_elementFaces; ++iface)
        {
            // global ID of the iface-th face in element ie
            UInt face = _mesh.localFaceId(ie, iface);
            // first adjacent element to face "face"
            UInt elem = _mesh.face(face).ad_first();
            if (elem == ie)
            {
                elem = _mesh.face(face).ad_second();
            }
            if (elem != 0)
            {
                // this is the list of adjacency
                // for each graph vertex, simply push back the ID of its neighbors
                M_jadj.push_back(elem - 1);
                ++sum;
                if (interfaceMap) // if I'm partitioning the solid in FSI
                {
                    if ((*M_repeatedFace)[sum])
                    {
                        adjwgt.push_back(0);
                    }
                    else
                    {
                        adjwgt.push_back(10);
                    }
                }
            }
        }
        // this is the list of "keys" to access M_jadj
        // graph element i has neighbors M_jadj[ k ],
        // with M_iadj[i] <= k < M_iadj[i+1]
        M_iadj.push_back(sum);
    }

    // **************
    // parMetis part

    // this array is to be used for weighted nodes on the graph:
    // usually we will set it to NULL

    int* vwgt = 0;

    int wgtflag;
    if (interfaceMap)
    {
        wgtflag = 1;
    }
    else
    {
        wgtflag = 0;
    }

    int ncon = 1;
    int numflag = 0;

    int edgecut; // here will be stored the number of edges cut in the partitioning process

    // additional options
    std::vector<int>  options(3,0);
    options[0] = 1; // means that additional options are actually passed
    options[1] = 3; // level of information to be returned during execution (see ParMETIS's defs.h file)
    options[2] = 1; // random number seed for the ParMETIS routine

    // number of desired subdomains: can be different from the number of procs
    int nparts = M_comm->NumProc();

    // fraction of vertex weight to be distributed to each subdomain.
    // here we want the subdomains to be of the same size
    std::vector<float> tpwgts(ncon * nparts, 1. / nparts);
    // imbalance tolerance for each vertex weight
    std::vector<float> ubvec (ncon, 1.05);

    Epetra_MpiComm* mpiComm = dynamic_cast <Epetra_MpiComm*> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();

    int nprocs;
    MPI_Comm_size(MPIcomm, &nprocs);

    /*
      (from ParMETIS v 3.1 manual)
      This routine is used to compute a k-way partitioning of a graph
      on p processors using the multilevel k-way multi-constraint
      partitioning algorithm.
    */

    int* adjwgtPtr(0);
    if (adjwgt.size() > 0)
    {
        adjwgtPtr = (int*) &adjwgt[0];
    }
    ParMETIS_V3_PartKway((int*) &M_vertexDist[0], (int*) &M_iadj[0], (int*) &M_jadj[0],
                         vwgt,  adjwgtPtr, &wgtflag, &numflag, &ncon, &nparts, &tpwgts[0], 
                         &ubvec[0], &options[0], &edgecut, &M_part[0], &MPIcomm);

    M_comm->Barrier();

    int nProc;
    nProc = M_comm->NumProc();

    // this is a vector of subdomains: each component is
    // the list of vertices belonging to the specific subdomain
    M_locProc.resize(nProc);

    // cycling on locally stored vertices
    for (UInt ii = 0; ii < M_part.size(); ++ii)
    {
        // here we are associating the vertex global ID to the subdomain ID
        M_locProc[M_part[ii]].push_back(ii + M_vertexDist[M_me]);
    }    
}

template<typename Mesh>
void partitionMesh<Mesh>::matchFluidPartitionsFSI()
{
    Epetra_MpiComm* mpiComm = dynamic_cast <Epetra_MpiComm*> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();       
    int nprocs;
    MPI_Comm_size(MPIcomm, &nprocs);

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

    for(UInt kk=0; kk<M_part.size(); ++kk)
    {
        if((*M_isOnProc)[kk+M_vertexDist[M_me]]!=-1)
        {
            ++mymatchesForProc[M_part[kk]][(*M_isOnProc)[kk+M_vertexDist[M_me]]];
        }
    }

    for(UInt j=0; (int)j<nprocs; ++j)
    {
        MPI_Allreduce(&mymatchesForProc[j][0], &matchesForProc[j][0], nprocs, 
                      MPI_INT, MPI_SUM, MPIcomm);
    }

    M_comm->Barrier();

    int suitableProcess = -1;
    UInt max = 0;

    for(int ii = 0; ii<nprocs; ++ii)
    {
        if(matchesForProc[M_me][ii] > max)
        {
            suitableProcess = ii;
            max = matchesForProc[M_me][ii];
        }
    }

    ASSERT(suitableProcess != -1, "one partition is without interface nodes!");
    procOrder[M_me] = suitableProcess;

    M_comm->Barrier();

    std::vector<UInt> maxs(nprocs);
    maxs[M_me] = max;
    for(int j = 0; j < nprocs ; ++j) // Allgather
    {
        MPI_Bcast(&maxs[j], 1, MPI_INT, j, MPIcomm); // perhaps generates errors
    }
        
    std::vector<pair<UInt, int> > procIndex(nprocs);
    for(int k = 0; k < nprocs; ++k)
    {
        procIndex[k] = std::make_pair( maxs[k], k);
    }

    std::sort(procIndex.begin(), procIndex.end() /*, &booleanCondition::reordering*/);

    for(int l=0;l<nprocs;++l)
    {
        for(int l=0;l<nprocs;++l)
        {
            for(int j=0; j<nprocs ; ++j) // Allgather
            {
                MPI_Bcast( &procOrder[j], 1, MPI_INT, j, MPIcomm); // perhaps generates errors
            }
        }
    }

    std::vector< std::vector<int> > locProc2(M_locProc);
    for(int j = nprocs; j > 0 ; --j)
    {
        if (orderingError[procOrder[procIndex[j - 1].second]] == false)
        {
            M_locProc[procOrder[procIndex[j - 1].second]] = locProc2[procIndex[j - 1].second];
        }
        else
        {
            std::cout << "Ordering error when assigning the processor"
                      << M_me << " to the partition," << std::endl 
                      << " parmetis did a bad job." << std::endl;
            for (int i = nprocs; i > 0; --i)
            {
                if(orderingError[procIndex[i - 1].second] == false) // means that i is the first proc not assigned
                {
                    procOrder[procIndex[j - 1].second] = procIndex[i - 1].second;
                    M_locProc[procIndex[i - 1].second] = locProc2[procIndex[j - 1].second];
                    break;
                }
            }
        }
        orderingError[procOrder[procIndex[j - 1].second]] = true;
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::redistributeElements()
{
    Epetra_MpiComm* mpiComm = dynamic_cast <Epetra_MpiComm*> (M_comm);
    MPI_Comm MPIcomm = mpiComm->Comm();       
    int nProc;
    MPI_Comm_size(MPIcomm, &nProc);

    int max_int (1000);
    int ssize[nProc];
    int rsize[nProc];
    // cycling on subdomains
    // TODO: Matteo please comment this part :)

    MPI_Status status;
    int size;

    MPI_Status  recv_status;
    MPI_Request send_request;

    for (int iproc = 0; iproc < nProc; ++iproc)
    {
        // all processes other than me are sending vertices
        // belonging to my subdomain
        if (int(M_me) != iproc)
        {                                                  
            size = M_locProc[iproc].size();
   
            // tell me how many vertices belonging to me you have to send me
            MPI_Isend(&size, 1, MPI_INT, iproc, 10, MPIcomm, &send_request);
            ssize[iproc] = size;
        }
        else
        {
            for (int jproc = 0; jproc < nProc; ++jproc)
            {
                if ((int)M_me != jproc)
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
                int size_part = size;
                
                // divide the whole data set into smaller packets
                while (size_part > max_int)
                {
                    incr += 1;
                    size_part = size / incr;
                }
                
                MPI_Send(&incr, 1, MPI_INT, iproc, 20, MPIcomm);
                MPI_Send(&size_part, 1, MPI_INT, iproc, 30, MPIcomm);
                
                for (int kk = 0; kk < incr; ++kk)
                {
                    MPI_Send(&pos, 1, MPI_INT, iproc, 100+kk, MPIcomm);
                    MPI_Send(&M_locProc[iproc][pos], size_part, MPI_INT, iproc, 5000000+kk, MPIcomm);
                    pos = pos + size_part;
                }

                int resto = size % incr;

                MPI_Send(&resto, 1, MPI_INT, iproc, 80, MPIcomm);

                if (resto != 0)
                {
                    MPI_Send(&pos, 1, MPI_INT, iproc, 40, MPIcomm);
                    MPI_Send(&M_locProc[iproc][pos], resto, MPI_INT, iproc, 50, MPIcomm);
                }
            }
            else
            {
                if (size != 0)
                {
                    MPI_Send(&M_locProc[iproc][0], size, MPI_INT, iproc, 60, MPIcomm);
                }
            }
        }
        else
        {
            for (int jproc = 0; jproc < nProc; ++jproc)
            {
                if (jproc != iproc)
                {
                    size = rsize[jproc];
                    std::vector<int> stack(size, 0);
                    
                    if (size > max_int)
                    {
                        int size_part, pos, incr;
                        
                        MPI_Recv(&incr, 1, MPI_INT, jproc, 20, MPIcomm, &status);
                        MPI_Recv(&size_part, 1, MPI_INT, jproc, 30, MPIcomm, &status);

                        for (int kk = 0; kk < incr; ++kk)
                        {
                            MPI_Recv(&pos, 1, MPI_INT, jproc, 100+kk, MPIcomm, &status);
                            MPI_Recv(&stack[pos], size_part , MPI_INT, jproc, 5000000+kk, MPIcomm, &status);
                        }
                        int resto = 0;
                        MPI_Recv(&resto, 1, MPI_INT, jproc, 80, MPIcomm, &status);

                        if (resto != 0)
                        {
                            MPI_Recv(&pos, 1, MPI_INT, jproc, 40, MPIcomm, &status);
                            MPI_Recv(&stack[pos],  resto, MPI_INT, jproc, 50, MPIcomm, &status);
                        }
                    }
                    else
                    {
                        if (size != 0)
                        {                            
                            MPI_Recv(&stack[0], size , MPI_INT, jproc, 60, MPIcomm, &status);
                        }
                    }
                    for (int jj = 0; jj < size; ++jj)
                    {
                        M_locProc[M_me].push_back(stack[jj]);
                    }
                }
            }
        }
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::constructLocalMesh(Mesh &_mesh)
{
    if (!M_me) 
    {    
        std::cout << "Building local mesh ..." << std::flush;
    }

    std::map<int, int>::iterator  im;
    std::set<int>::iterator       is;

    int count = 1;
    UInt ielem;
    UInt inode;

    // cycle on local element's ID
    for (UInt jj = 0; jj < M_locProc[M_me].size(); ++jj)
    {
        ielem = M_locProc[M_me][jj];
        M_localVolumes.push_back(ielem);

        // cycle on element's nodes
        for (UInt ii = 1; ii <= M_elementNodes; ++ii)
        {
            inode = _mesh.volume(ielem + 1).point(ii).id();
            im    = M_globalToLocalNode.find(inode);

            // if the node is not yet present in the list of local nodes, then add it
            // CAREFUL: also local numbering starts from 1 in RegionMesh
            if (im == M_globalToLocalNode.end())
            {
                M_globalToLocalNode.insert(std::make_pair(inode, count));
                ++count;
                // store here the global numbering of the node
                M_localNodes.push_back(_mesh.volume(ielem + 1).point(ii).id());
            }
        }

        // cycle on element's edges
        for (UInt ii = 1; ii <= M_elementEdges; ++ii)
        {
        	// store here the global numbering of the edge
          	M_localEdges.insert(_mesh.localEdgeId(ielem + 1, ii));
        }
        
        // cycle on element's faces
        for (UInt ii = 1; ii <= M_elementFaces; ++ii)
        {
            // store here the global numbering of the face
            M_localFaces.insert(_mesh.localFaceId(ielem + 1, ii));
        }
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::constructNodes(Mesh &_mesh)
{
    std::vector<int>::iterator it;

    M_nBoundaryPoints = 0;
    M_mesh->pointList.reserve(M_localNodes.size());
    // guessing how many boundary points on this processor.
    M_mesh->_bPoints.reserve(_mesh.numBPoints() * M_localNodes.size() / _mesh.numBPoints());

    UInt inode = 1;

    typename Mesh::PointType *pp = 0;

    // loop in the list of local nodes:
    // in this loop inode is the local numbering of the points
    for (it = M_localNodes.begin(); it != M_localNodes.end(); ++it, ++inode)
    {
        typename Mesh::PointType point = 0;

        // create a boundary point in the local mesh, if needed
        bool boundary = _mesh.isBoundaryPoint(*it);
        if (boundary)
        {
            ++M_nBoundaryPoints;
        }

        pp = &M_mesh->addPoint(boundary);
        pp->setMarker(_mesh.point(*it).marker());

        pp->x() = _mesh.point(*it).x();
        pp->y() = _mesh.point(*it).y();
        pp->z() = _mesh.point(*it).z();

        UInt id = _mesh.point(*it).id();

        pp->setId(id);
        pp->setLocalId(inode);

        M_mesh->localToGlobalNode().insert(std::make_pair(inode, id));
        M_mesh->globalToLocalNode().insert(std::make_pair(id, inode));
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::constructVolumes(Mesh &_mesh)
{
    std::map<int, int>::iterator im;
    std::vector<int>::iterator it;
    int count = 1;
    UInt inode;

    typename Mesh::VolumeType * pv = 0;

    M_mesh->volumeList.reserve(M_localVolumes.size());

    // loop in the list of local elements
    // CAREFUL! in this loop inode is the global numbering of the points
    // We insert the local numbering of the nodes in the local volume list
    for (it = M_localVolumes.begin(); it != M_localVolumes.end(); ++it, ++count)
    {
        pv = &M_mesh->addVolume();
        // CAREFUL! in ParMETIS data structures, numbering starts from 0
        pv->setId (_mesh.volume(*it + 1).id());
        pv->setLocalId(count);

        M_globalToLocalVolume.insert(make_pair(_mesh.volume(*it + 1).id(), count));

        for (ID id = 1; id <= M_elementNodes; ++id)
        {
            inode = _mesh.volume(*it + 1).point(id).id();
            // im is an iterator to a map element
            // im->first is the key (i. e. the global ID "inode")
            // im->second is the value (i. e. the local ID "count")
            im = M_globalToLocalNode.find(inode);
            pv->setPoint(id, M_mesh->pointList( (*im).second ));
        }

        int ibc = _mesh.volume(*it + 1).marker();

        pv->setMarker(EntityFlag( ibc ));
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::constructEdges(Mesh &_mesh)
{
    std::map<int, int>::iterator im;
    std::set<int>::iterator is;

    typename Mesh::EdgeType * pe;
    UInt inode;
    int count = 1;

    M_nBoundaryEdges = 0;
    M_mesh->edgeList.reserve(M_localEdges.size());

    // loop in the list of local edges
    for (is = M_localEdges.begin(); is != M_localEdges.end(); ++is, ++count)
    {
        // create a boundary edge in the local mesh, if needed
        bool boundary = (_mesh.isBoundaryEdge(*is));
        if (boundary)
        {
            // create a boundary edge in the local mesh, if needed
            ++M_nBoundaryEdges;
        }

        pe = &M_mesh->addEdge(boundary);

        pe->setId (_mesh.edge(*is).id());
        pe->setLocalId(count);

        for (ID id = 1; id <= 2; ++id)
        {
            inode = _mesh.edge(*is).point(id).id();
            // im is an iterator to a map element
            // im->first is the key (i. e. the global ID "inode")
            // im->second is the value (i. e. the local ID "count")
            im = M_globalToLocalNode.find(inode);
            pe->setPoint(id, M_mesh->pointList((*im).second));
        }
        pe->setMarker(_mesh.edge(*is).marker());
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::constructFaces(Mesh &_mesh)
{
    std::map<int, int>::iterator im;
    std::set<int>::iterator      is;

    typename Mesh::FaceType * pf = 0;
    
    UInt inode;
    int count = 1;

    M_nBoundaryFaces = 0;
    M_mesh->faceList.reserve(M_localFaces.size());

    // loop in the list of local faces
    for (is = M_localFaces.begin(); is != M_localFaces.end(); ++is, ++count)
    {
      // create a boundary face in the local mesh, if needed
        bool boundary = (_mesh.isBoundaryFace(*is));
        if (boundary)
        {
            ++M_nBoundaryFaces;
        }

        pf =  &M_mesh->addFace(boundary);

        pf->setId (_mesh.face(*is).id());
        pf->setLocalId(count);

        int elem1 = _mesh.face(*is).ad_first();
        int elem2 = _mesh.face(*is).ad_second();

        // find the mesh elements adjacent to the face
        im =  M_globalToLocalVolume.find(elem1);

        int localElem1;

        if (im == M_globalToLocalVolume.end())
        {
            localElem1 = 0;
        }
        else
        {
            localElem1 = (*im).second;
        }

        im =  M_globalToLocalVolume.find(elem2);

        int localElem2;
        if (im == M_globalToLocalVolume.end())
        {
            localElem2 = 0;
        }
        else
        {
            localElem2 = (*im).second;
        }

        // if this process does not own either of the adjacent elements
        // then the two adjacent elements and the respective face positions coincide in the local mesh
        //possible bug fixed: not only the two adjacent elements face, but also the face positions should coincide.
        //otherwise it could happen that a pair(element, position) is associated to different faces.
        //This can lead to a wrong treatment of the dofPerFace (in 2D of the dofPerEdge, as occurred with P2)
        //
        if ((localElem1 == 0) && !boundary)
        {
        	pf->ad_first() = localElem2;
        	pf->pos_first() = _mesh.face(*is).pos_second();
        }
        else
        {
        	pf->ad_first() = localElem1;
        	pf->pos_first() = _mesh.face(*is).pos_first();
        }

        if ((localElem2 == 0) && !boundary)
        {
        	pf->ad_second() = localElem1;
        	pf->pos_second() = _mesh.face(*is).pos_first();
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
        	pf->setPoint(id, M_mesh->pointList((*im).second));
        }

        pf->setMarker(_mesh.face(*is).marker());

        M_mesh->setLinkSwitch("HAS_ALL_FACES");
        M_mesh->setLinkSwitch("FACES_HAVE_ADIACENCY");
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::finalSetup(Mesh &_mesh)
{
    UInt nVolumes = M_localVolumes.size();
    UInt nNodes   = M_localNodes.size();
    UInt nEdges   = M_localEdges.size();
    UInt nFaces   = M_localFaces.size();

    M_mesh->setMaxNumPoints (nNodes, true);
    M_mesh->setMaxNumEdges  (nEdges, true);
    M_mesh->setMaxNumFaces  (nFaces, true);
    M_mesh->setMaxNumVolumes( nVolumes, true);

    M_mesh->setMaxNumGlobalPoints (_mesh.numPoints());
    M_mesh->setNumGlobalVertices  (_mesh.numPoints());
    M_mesh->setMaxNumGlobalEdges  (_mesh.numEdges());
    M_mesh->setMaxNumGlobalFaces  (_mesh.numFaces());

    M_mesh->setMaxNumGlobalVolumes(_mesh.numVolumes());
    M_mesh->setNumBFaces    (M_nBoundaryFaces);

    M_mesh->setNumBPoints   (M_nBoundaryPoints);
    M_mesh->setNumBEdges    (M_nBoundaryEdges);

    M_mesh->setNumVertices (nNodes );
    M_mesh->setNumBVertices(M_nBoundaryPoints);


    M_mesh->updateElementEdges();

    M_mesh->updateElementFaces();

    if (!M_me) 
    {
        std::cout << "done" << std::endl;
        std::cout << "Creating the map ... " << std::flush;
    }
}

template<typename Mesh>
void partitionMesh<Mesh>::createRepeatedMap(Mesh &_mesh)
{
    std::set<int>::iterator is;
    std::vector<int> elementList = M_locProc[M_me];

    UInt inode, ielem;

    // repeated element map creation

    // use sets to store each entity only once
    std::set<int>    repeatedNodeList;
    std::set<int>    repeatedEdgeList;
    std::set<int>    repeatedFaceList;

    for (UInt ii = 0; ii < elementList.size(); ++ii)
    {
        ielem = elementList[ii];
        M_repeatedVolumeVector.push_back(ielem + 1);
        for (UInt jj = 1; jj <= M_elementNodes; ++jj)
        {
        	inode = _mesh.volume(ielem + 1).point(jj).id();
        	repeatedNodeList.insert(inode);
        }
        for (UInt jj = 1; jj <= M_elementEdges; ++jj)
        {
             UInt iedge = _mesh.localEdgeId(ielem + 1, jj);
             repeatedEdgeList.insert((int) iedge);
        }
        for (UInt jj = 1; jj <= M_elementFaces; ++jj)
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

    if (!M_me) std::cout << "done" << std::endl;
}

template<typename Mesh>
void partitionMesh<Mesh>::execute(Mesh &_mesh, Epetra_Comm &_comm,
                                  Epetra_Map *interfaceMap,
                                  Epetra_Map *interfaceMapRep)
{
    // Set element parameters (number of nodes, faces, edges and number of nodes
    // on each face according to the type of mesh element used.
    setElementParameters();

    // Build graph vertex distribution vector. Graph vertex represents one element
    // in the mesh.
    distributeElements(_mesh.numElements());


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
    if(interfaceMap)
    {
        findRepeatedFacesFSI(_mesh, interfaceMap, interfaceMapRep);
    }
    //////////////////// END OF SOLID PARTITION PART ////////////////////////


    // Partition connectivity graph
    partitionConnectivityGraph(_mesh, interfaceMap);


    //////////////// BEGIN OF SOLID PARTITION PART ////////////////
    if(interfaceMap)
    {
        matchFluidPartitionsFSI();
    }
    ////////////////// END OF SOLID PARTITION PART /////////////////////


    // Redistribute elements to appropriate processors before building the
    // partitioned mesh.
    redistributeElements();

#ifdef DEBUG
    Debug(4000) << M_me << " has " << M_locProc[M_me].size() << " elements.\n";
#endif

    // ***********************
    // local mesh construction
    // ***********************
    constructLocalMesh(_mesh);

    // ******************
    // nodes construction
    // ******************
    constructNodes(_mesh);

    // ******************
    // volumes construction
    // ******************
    constructVolumes(_mesh);

    // ******************
    // edges construction
    // ******************
    constructEdges(_mesh);

    // ******************
    // faces construction
    // ******************
    constructFaces(_mesh);

    // ******************
    // final setup
    // ******************
    finalSetup(_mesh);

    // *********************
    // repeated map creation
    // *********************
    createRepeatedMap(_mesh);
}

template<typename Mesh>
partitionMesh<Mesh>::~partitionMesh()
{
}

}
#endif
