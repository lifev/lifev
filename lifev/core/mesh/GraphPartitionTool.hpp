//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010, 2011 EPFL, Politecnico di Milano, Emory University

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
    @brief Class that partitions the graph associated with a mesh

    @date 17-11-2011
    @author Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef GRAPH_PARTITION_TOOL_H
#define GRAPH_PARTITION_TOOL_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>
#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
#include <zoltan.h>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>

namespace LifeV {

//! This structure is used to pack objects in the Zoltan migration phase
/*!
    @author Radu Popescu <radu.popescu@epfl.ch>
 */
struct TransportBuffer
{
    // The gid of the object
    int gid;
    // The partition to which it belongs
    int part;
};

//! Class that partitions the graph associated with a mesh
/*!
    @author Radu Popescu <radu.popescu@epfl.ch>

    This class uses the Zoltan package to partition the graph associated
    with a mesh. This class builds the dual graph of the mesh, partitions
    it according to a set of parameters and the stores the partitioning
    in a table (vector of vectors).
    At the end of the partitioning process, each vector will contain the
    GID of the elements in a partition.

    While this class can be used stand-alone, it is contained inside the
    MeshPartitionTool class and used automatically during the mesh
    partitioning process.

    More on class functionality to follow. Stay tuned...
 */
template<typename MeshType>
class GraphPartitionTool
{
public:
    //! @name Public Types
    //@{
    typedef MeshType                                             mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                         meshPtr_Type;
    typedef std::map<Int, std::vector<Int> >                     internalTable_Type;
    typedef std::map<Int, boost::shared_ptr<std::vector<Int> > > exportTable_Type;
    //@}

    //! @name Constructor & Destructor
    //@{
    //! Empty Constructor
    GraphPartitionTool();

    //! Constructor taking the original mesh, the MPI comm and parameters
    /*!
        This constructor can be used to build the object and perform the
        graph partitioning in one shot.
 
        @param mesh The original mesh whose graph this object will partition
        @param comm The Epetra MPI comm object which contains the processes
                    which participate
        @param parameters The Teuchos parameter list which contains the
                          partitioning parameters
     */
    GraphPartitionTool(meshPtr_Type& mesh,
                       boost::shared_ptr<Epetra_Comm>& comm,
                       Teuchos::ParameterList& parameters);

    //! Destructor
    virtual ~GraphPartitionTool() {}
    //@}

    //! @name Methods
    //@{
    //! Configure the object (used with empty constructor)
    void setup(meshPtr_Type& mesh,
               boost::shared_ptr<Epetra_Comm>& comm,
               Teuchos::ParameterList& parameters);

    //! Set values for all the parameters, with default values where needed
    void setParameters(Teuchos::ParameterList& parameters);

    //! Performs the graph partitioning
    void run();
    //@}

    //! @name Get Methods
    //@{
    //! Get a pointer to one of the partitions
    const boost::shared_ptr<std::vector<Int> >& getPartition(const UInt i) const {return M_partitionTable.find(i)->second;}

    //! Get number of stored graph elements
    UInt numStoredElements() const {return static_cast<UInt>(M_elementList.size());}

    //! First global index that is initially assigned to process i
    Int firstIndex(const Int i) const {return M_indexBounds[i];}

    //! Last global index that is initially assigned to process i
    Int lastIndex(const Int i) const {return M_indexBounds[i + 1] - 1;}

    //! The internally stored dual graph of the mesh
    const internalTable_Type& graph() const {return M_graph;}

    //! The vector of stored element GIDs
    const std::vector<Int>& elementList() const {return M_elementList;}

    //! The vector of stored element GIDs (non-const)
    std::vector<Int>& elementList() {return M_elementList;}

    //! The vector of stored element partitions
    const std::vector<Int>& elementParts() const {return M_elementParts;}

    //! The vector of stored element partitions (non-const)
    std::vector<Int>& elementParts() {return M_elementParts;}

    //! The PID of the process
    Int myPID() const {return M_myPID;}

    //! The number of processes in the comm
    Int numProcessors() const {return M_numProcessors;}
    //@}

    //! @name Static methods
    //@{
    static int getNumElements(void *data, int *ierr);
    static void getElementList(void *data, int sizeGID, int sizeLID,
                               ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int wgt_dim, float *obj_wgts, int *ierr);
    static void getNumNeighboursList(void *data, int sizeGID, int sizeLID,
                                     int num_obj,
                                     ZOLTAN_ID_PTR globalID,
                                     ZOLTAN_ID_PTR localID,
                                     int *numEdges, int *ierr);
    static void getNeighbourList(void *data, int sizeGID, int sizeLID,
                                 int num_obj, ZOLTAN_ID_PTR globalID,
                                 ZOLTAN_ID_PTR localID, int *num_edges,
                                 ZOLTAN_ID_PTR nborGID, int *nborProc,
                                 int wgt_dim, float *ewgts, int *ierr);
    static void getTransferObjectSizes(void *data, int num_gid_entries,
                                      int num_lid_entries, int num_ids,
                                      ZOLTAN_ID_PTR global_ids,
                                      ZOLTAN_ID_PTR local_ids,
                                      int *sizes, int *ierr);
    static void packObjects(void *data, int num_gid_entries,
                            int num_lid_entries, int num_ids,
                            ZOLTAN_ID_PTR global_ids,
                            ZOLTAN_ID_PTR local_ids,
                            int *dest, int *sizes, int *idx,
                            char *buf, int *ierr);
    static void unpackObjects(void *data, int num_gid_entries,
                              int num_ids, ZOLTAN_ID_PTR global_ids,
                              int *sizes, int *idx, char *buf, int *ierr);

private:
    // Private copy constructor and assignment operator are disabled
    GraphPartitionTool(const GraphPartitionTool&);
    GraphPartitionTool& operator=(const GraphPartitionTool&);

    //! @name Private Methods
    //@{
    //! Distribute the partitions among available processors
    void distributePartitions();

    //! Build the dual graph of the unpartitioned mesh
    void buildGraph();

    //! Migrate elements between locally stored partitions
    void localMigrate(int numExport,
                      ZOLTAN_ID_PTR exportLocalGids,
                      int* exportProcs, int* exportToPart);

    //! Build the partition table that can be exported to the mesh builder
    void buildPartitionTable();

    //! Partition the graph
    void partitionGraph();
    //@}
 
    boost::shared_ptr<Epetra_Comm>             M_comm;
    Int                                        M_myPID;
    Int                                        M_numProcessors;
    Int                                        M_numPartitions;
    Int                                        M_numPartitionsPerProcessor;
    Int                                        M_myFirstPartition;
    Int                                        M_myLastPartition;
    std::string                                M_topology;
    std::vector<Int>                           M_indexBounds;
    Teuchos::ParameterList                     M_parameters;
    meshPtr_Type                               M_mesh;
    exportTable_Type                           M_partitionTable;
    internalTable_Type                         M_graph;
    // TODO: possible improvement (memory-wise) is to implement a bidirectional map
    // instead of using these two vectors and M_partitionTable
    // MeshPartitionTool expects a partition-to-element map, while in the graph
    // cutting routine an element-to-partition map is needed
    std::vector<Int>                           M_elementList;
    std::vector<Int>                           M_elementParts;

    struct Zoltan_Struct*                      M_zoltanStruct;
};

//
// IMPLEMENTATION
//

// =================================
// Constructors and Destructor
// =================================

template<typename MeshType>
GraphPartitionTool<MeshType>::GraphPartitionTool() :
    M_comm(),
    M_myPID(0),
    M_numProcessors(0),
    M_numPartitions(0),
    M_numPartitionsPerProcessor(0),
    M_myFirstPartition(0),
    M_myLastPartition(0),
    M_parameters(),
    M_mesh(),
    M_zoltanStruct(0)
{}

template<typename MeshType>
GraphPartitionTool<MeshType>::GraphPartitionTool(meshPtr_Type& mesh,
                                                 boost::shared_ptr<Epetra_Comm>& comm,
                                                 Teuchos::ParameterList& parameters) :
    M_comm(comm),
    M_myPID(M_comm->MyPID()),
    M_numProcessors(M_comm->NumProc()),
    M_numPartitions(0),
    M_numPartitionsPerProcessor(0),
    M_myFirstPartition(0),
    M_myLastPartition(0),
    M_parameters(),
    M_mesh(mesh),
    M_zoltanStruct(0)
{
    setParameters(parameters);
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::setup(meshPtr_Type& mesh,
                                         boost::shared_ptr<Epetra_Comm>& comm,
                                         Teuchos::ParameterList& parameters)
{
    M_comm = comm;
    M_myPID = M_comm->MyPID();
    M_numProcessors = M_comm->NumProc();
    M_mesh = mesh;

    setParameters(parameters);
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::setParameters(Teuchos::ParameterList& parameters)
{
    // Here put some default values for the parameters and then import
    // the user supplied list, overwriting the corresponding parameters
    M_parameters.set("num_partitions", static_cast<Int>(M_comm->NumProc()), "");
    M_parameters.set("topology", "1", "");

    M_parameters.setParameters(parameters);

    M_numPartitions = M_parameters.get<Int>("num_partitions");
    M_numPartitionsPerProcessor = M_numPartitions / M_numProcessors;
    M_topology = M_parameters.get<std::string>("topology");
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::run()
{ 
    distributePartitions();
    buildGraph();
    partitionGraph();
}

// =================
// Static methods
// =================

template<typename MeshType>
int GraphPartitionTool<MeshType>::getNumElements(void *data, int *ierr)
{
    GraphPartitionTool<MeshType>* object = (GraphPartitionTool<MeshType>*) data;

    *ierr = ZOLTAN_OK;
    return object->numStoredElements();
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::getElementList(void *data, int /*sizeGID*/, int /*sizeLID*/,
                                              ZOLTAN_ID_PTR globalID,
                                              ZOLTAN_ID_PTR localID,
                                              int /*wgt_dim*/, float* /*obj_wgts*/,
                                              int *ierr)
{
    GraphPartitionTool<MeshType>* object = (GraphPartitionTool<MeshType>*) data;

    UInt k = 0;
    for (UInt i = 0; i < object->numStoredElements(); ++i) {
        globalID[i] = object->elementList()[i];
        localID[i] = k;
        k++;
    }

    *ierr = ZOLTAN_OK;
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::getNumNeighboursList(void *data, int /*sizeGID*/, int /*sizeLID*/,
                                                    int num_obj,
                                                    ZOLTAN_ID_PTR globalID,
                                                    ZOLTAN_ID_PTR /*localID*/,
                                                    int *numEdges, int *ierr)
{
    GraphPartitionTool<MeshType>* object = (GraphPartitionTool<MeshType>*) data;

    for (int element = 0; element < num_obj; ++element) {
        numEdges[element] = object->graph().find(globalID[element])->second.size();
    }

    *ierr = ZOLTAN_OK;
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::getNeighbourList(void *data, int /*sizeGID*/, int /*sizeLID*/,
                                                int num_obj, ZOLTAN_ID_PTR globalID,
                                                ZOLTAN_ID_PTR /*localID*/, int *num_edges,
                                                ZOLTAN_ID_PTR nborGID, int *nborProc,
                                                int /*wgt_dim*/, float* /*ewgts*/, int *ierr)
{
    GraphPartitionTool<MeshType>* object = (GraphPartitionTool<MeshType>*) data;

    const std::vector<Int>& elementList = object->elementList();
    const Int numStoredElements = object->numStoredElements();

    std::vector<Int>::const_iterator iter;
    int pos = 0;
    for (int element = 0; element < num_obj; ++element) {
        iter = object->graph().find(globalID[element])->second.begin();
        for (int k = 0; k < num_edges[element]; ++k) {
            nborGID[pos] = *iter;
            int pid = 0;
            // TODO: do a binary search here for speed
            for (int i = 0; i < object->numProcessors(); ++i) {
                if ((*iter >= static_cast<Int>(elementList[0])) &&
                    (*iter <= static_cast<Int>(elementList[numStoredElements - 1]))) {
                    pid = i;
                    break;
                }
            }
            nborProc[pos] = pid;
            ++pos;
            ++iter;
        }
    }

    *ierr = ZOLTAN_OK;
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::getTransferObjectSizes(void* /*data*/, int /*num_gid_entries*/,
                                                          int /*num_lid_entries*/, int num_ids,
                                                          ZOLTAN_ID_PTR /*global_ids*/,
                                                          ZOLTAN_ID_PTR /*local_ids*/,
                                                          int* sizes, int *ierr)
{
    int sizeOfBuffer = sizeof(TransportBuffer);
    for (int i = 0; i < num_ids; ++i) {
        sizes[i] = sizeOfBuffer;
    }

    *ierr = ZOLTAN_OK;
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::packObjects(void *data, int /*num_gid_entries*/,
                                               int /*num_lid_entries*/, int num_ids,
                                               ZOLTAN_ID_PTR global_ids,
                                               ZOLTAN_ID_PTR local_ids,
                                               int *dest, int* /*sizes*/, int *idx,
                                               char *buf, int *ierr)
{
    GraphPartitionTool<MeshType>* object = (GraphPartitionTool<MeshType>*) data;

    // Pack gids and part numbers in the buffer
    for (int i = 0; i < num_ids; ++i) {
        TransportBuffer* buffer = (TransportBuffer*) &buf[idx[i]];
        buffer->gid = global_ids[i];
        buffer->part = dest[i];

        // Remove the objects from the local storage
        // TODO: find better solution. This will be slow and inefficient !!! 
        object->elementList()[local_ids[i]] = -1;
        object->elementParts()[local_ids[i]] = -1;
    }

    *ierr = ZOLTAN_OK;
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::unpackObjects(void *data, int /*num_gid_entries*/,
                                                 int num_ids, ZOLTAN_ID_PTR /*global_ids*/,
                                                 int* /*sizes*/, int *idx,
                                                 char *buf, int *ierr)
{
    GraphPartitionTool<MeshType>* object = (GraphPartitionTool<MeshType>*) data;

    // Unpack gids and part numbers from the buffer
    for (int i = 0; i < num_ids; ++i) {
        TransportBuffer* buffer = (TransportBuffer*) &buf[idx[i]];
        object->elementList().push_back(buffer->gid);
        object->elementParts().push_back(buffer->part);
    }

    *ierr = ZOLTAN_OK;
}

// =======================
// Private methods
// =======================

template<typename MeshType>
void GraphPartitionTool<MeshType>::distributePartitions()
{
    // The algorithm to distribute partitions isn't clever at all.
    // We assume the number of partitions is a multiple of the
    // number of processes.

    M_myFirstPartition = M_myPID * M_numPartitionsPerProcessor;
    M_myLastPartition = (M_myPID + 1) * M_numPartitionsPerProcessor - 1;

    for (Int i = M_myFirstPartition; i <= M_myLastPartition; ++i) {
        M_partitionTable[i].reset(new std::vector<Int>);
    }
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::buildGraph()
{
    // This next part computes the first and last element global index
    // that the processor has to handle and makes a local vector of
    // all the GIDs that the process owns

    Int k = M_mesh->numElements();

    M_indexBounds.resize(M_numProcessors + 1);
    M_indexBounds[0] = 0;

    for (Int i = 0; i < M_numProcessors; ++i) {
        UInt l = k / (M_numProcessors - i);
        M_indexBounds[i + 1] = M_indexBounds[i] + l;
        k -= l;
    }

    M_elementList.resize(M_indexBounds[M_myPID + 1] - M_indexBounds[M_myPID]);
    Int startIndex = firstIndex(M_myPID);
    for (UInt i = 0; i < numStoredElements(); ++i) {
        M_elementList[i] = startIndex + i;
    }

    M_elementParts.resize(numStoredElements());

    k = numStoredElements();
    std::vector<Int> partitionBounds(M_numPartitionsPerProcessor + 1);
    partitionBounds[0] = 0;
    for (Int i = 0; i < M_numPartitionsPerProcessor; ++i) {
        UInt l = k / (M_numPartitionsPerProcessor - i);
        partitionBounds[i + 1] = partitionBounds[i] + l;
        k -= l;
    }
    for (Int i = 0; i < M_numPartitionsPerProcessor; ++i) {
        for (Int lid = partitionBounds[i]; lid < partitionBounds[i + 1]; ++lid) {
            M_elementParts[lid] = M_myFirstPartition + i;
        }
    }

    // Build the graph of the local elements
    Int numDimensions = MeshType::elementShape_Type::S_nDimensions;
    int numNeighbours;
    switch (numDimensions) {
    case 2:
        numNeighbours = 3;
        break;
    case 3:
        numNeighbours = 4;
        break;
    default:
        numNeighbours = 0;
        break;
    }

    UInt numElementFaces = MeshType::elementShape_Type::S_numFaces;

    for (UInt i = 0; i < numStoredElements(); ++i) {
        UInt ie = M_elementList[i];
        M_graph.insert(std::pair<Int, std::vector<Int> >(ie, std::vector<Int>()));
        M_graph[ie].reserve(numNeighbours);
        for (UInt iface = 0; iface < numElementFaces; ++iface) {
            UInt face = M_mesh->localFaceId(ie, iface);
            UInt elem = M_mesh->face(face).firstAdjacentElementIdentity();
            if (elem == ie)
            {
                elem = M_mesh->face(face).secondAdjacentElementIdentity();
            }
            if (elem != NotAnId)
            {
                M_graph[ie].push_back(elem);
            }
        }
    }
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::localMigrate(int numExport,
                                                ZOLTAN_ID_PTR exportLocalGids,
                                                int* exportProcs,
                                                int* exportToPart)
{
    for (int i = 0; i < numExport; ++i) {
        if (exportProcs[i] == M_myPID) {
            // We shouldn't need to check this, still ...
            if (M_elementParts[exportLocalGids[i]] != exportToPart[i]) {
                M_elementParts[exportLocalGids[i]] = exportToPart[i];
            }
        }
    }
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::buildPartitionTable()
{
    for (exportTable_Type::iterator it = M_partitionTable.begin();
         it != M_partitionTable.end(); ++it) {
        it->second->reserve(numStoredElements() / M_numPartitionsPerProcessor);
    }

    for (UInt i = 0; i < numStoredElements(); ++i) {
        // We marked elements that were moved to a different processor with -1
        if (static_cast<Int>(M_elementList[i]) != -1) {
            M_partitionTable[M_elementParts[i]]->push_back(M_elementList[i]);
        }
    }
    for (exportTable_Type::iterator it = M_partitionTable.begin();
         it != M_partitionTable.end(); ++it) {
        std::sort(it->second->begin(), it->second->end());
    }
}

template<typename MeshType>
void GraphPartitionTool<MeshType>::partitionGraph()
{
    int argc = 1;
    char* argv;
    float ver;
    boost::shared_ptr<Epetra_MpiComm> mpiComm = boost::dynamic_pointer_cast<Epetra_MpiComm>(M_comm);

    Zoltan_Initialize(argc, &argv, &ver);
    M_zoltanStruct = Zoltan_Create(mpiComm->Comm());

    // Check these parameters. It may be wiser not to hard-code anything
    Zoltan_Set_Param(M_zoltanStruct, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(M_zoltanStruct, "HIER_DEBUG_LEVEL", "0");
    Zoltan_Set_Param(M_zoltanStruct, "LB_METHOD", "HIER");
    Zoltan_Set_Param(M_zoltanStruct, "HIER_ASSIST", "1");
    Zoltan_Set_Param(M_zoltanStruct, "NUM_GID_ENTRIES", "1"); 
    Zoltan_Set_Param(M_zoltanStruct, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(M_zoltanStruct, "RETURN_LISTS", "EXPORT");
    // We don't want remapping enabled since we need the best quality partitioning available
    Zoltan_Set_Param(M_zoltanStruct, "REMAP", "0");
    // Let the Zoltan_Migrate function only handle off processor transfers
    // We move the elements between same-processor parts manually, to avoid MPI calls (??)
    Zoltan_Set_Param(M_zoltanStruct, "MIGRATE_ONLY_PROC_CHANGES", "1");
    Zoltan_Set_Param(M_zoltanStruct, "TOPOLOGY",
                     M_parameters.get<std::string>("topology").c_str());
    Zoltan_Set_Param(M_zoltanStruct, "NUM_GLOBAL_PARTS",
                     (boost::lexical_cast<std::string>(M_numPartitions)).c_str());
    Zoltan_Set_Param(M_zoltanStruct, "NUM_LOCAL_PARTS",
                     (boost::lexical_cast<std::string>(M_numPartitionsPerProcessor)).c_str());

    Zoltan_Set_Num_Obj_Fn(M_zoltanStruct, GraphPartitionTool::getNumElements, this);
    Zoltan_Set_Obj_List_Fn(M_zoltanStruct, GraphPartitionTool::getElementList, this);
    Zoltan_Set_Num_Edges_Multi_Fn(M_zoltanStruct, GraphPartitionTool::getNumNeighboursList, this);
    Zoltan_Set_Edge_List_Multi_Fn(M_zoltanStruct, GraphPartitionTool::getNeighbourList, this);
    Zoltan_Set_Obj_Size_Multi_Fn(M_zoltanStruct, GraphPartitionTool::getTransferObjectSizes, this);
    Zoltan_Set_Pack_Obj_Multi_Fn(M_zoltanStruct, GraphPartitionTool::packObjects, this);
    Zoltan_Set_Unpack_Obj_Multi_Fn(M_zoltanStruct, GraphPartitionTool::unpackObjects, this);

    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;

    // Partition the graph
    Zoltan_LB_Partition(M_zoltanStruct,
                        &changes,
                        &numGidEntries,
                        &numLidEntries,
                        &numImport,
                        &importGlobalGids,
                        &importLocalGids,
                        &importProcs,
                        &importToPart,
                        &numExport,
                        &exportGlobalGids,
                        &exportLocalGids,
                        &exportProcs,
                        &exportToPart);

    // We first need to migrate elements between locally stored parts. Aftewards, Zoltan
    // can handle data movement between processors.
    localMigrate(numExport, exportLocalGids, exportProcs, exportToPart);

    M_comm->Barrier();

    // Migrate data after partitioning
    // WARNING! After Zoltan does the migration, the M_elementList and M_elementParts vectors
    // are NOT ordered in ascending order of LID and GID. Make no assumptions about the order
    // of local elements. Elements that have been migrated to other processors are marked with
    // -1 in M_elementList and M_elementParts
    Zoltan_Migrate(M_zoltanStruct,
                   -1,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   numExport,
                   exportGlobalGids,
                   exportLocalGids,
                   exportProcs,
                   exportToPart);

    // Clean up after partitioning and migration
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                        &exportProcs, &exportToPart);
    Zoltan_Destroy(&M_zoltanStruct);

    // Build the partition->element table that can be exported to the mesh partitioner
    buildPartitionTable();
}

} // Namespace LifeV

#endif // GRAPH_PARTITION_TOOL_H
