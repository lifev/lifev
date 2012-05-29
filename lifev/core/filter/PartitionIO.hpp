//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010, 2011, 2012 EPFL, Politecnico di Milano, Emory University

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
  @brief Class that handles I/O of mesh parts (for offline partitioning mode)

  @date 10-05-2012
  @author Radu Popescu <radu.popescu@epfl.ch>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef PARTITION_IO_H_
#define PARTITION_IO_H_

#include<lifev/core/LifeV.hpp>

#include <Epetra_config.h>

#ifdef HAVE_HDF5
#ifdef HAVE_MPI

#include <algorithm>
#include <string>
#include <vector>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <mpi.h>
#include <hdf5.h>

#include <boost/shared_ptr.hpp>
#include <Epetra_MpiComm.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV {

/*!
  @brief Class that handles I/O of mesh parts (for offline partitioning mode)
  @author Radu Popescu radu.popescu@epfl.ch

  This class is used to write the mesh parts produced by (offline) mesh
  partitioning into a single HDF5 container. This part is done with a single
  MPI process (usually on a workstation).

  The class is later used during an online simulation (multiple MPI processes,
  on a cluster, supercomputer etc.) to read the mesh parts.

  Usage:
  	  - call the constructor (a default one is also available)
  	  - if the default constructor was used, call setup method to initialize
  	  	the object
	  - call write method to store a set of mesh parts in the HDF5 file
	  - call read method to load the mesh part associate with the current
	    MPI process from the HDF5 file

  Creating, opening and closing the HDF5 file is done automatically by the
  object.

  // the description of the storage format could be useful here
*/
template<typename MeshType>
class PartitionIO {
public:
    //! @name Public Types
    //@{
    typedef MeshType mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef std::vector<meshPtr_Type> meshParts_Type;
    typedef boost::shared_ptr<meshParts_Type> meshPartsPtr_Type;
    // TODO: enable the following two typedefs if I decide to also
    // write the graph to the HDF5 container
    // typedef std::vector<std::vector<Int> > table_Type;
	// typedef boost::shared_ptr<table_Type> tablePtr_Type;
	typedef boost::shared_ptr<Epetra_Comm> commPtr_Type;
    //@}

    //! \name Constructors & Destructors
    //@{
    //! Default empty constructor
	PartitionIO() {}

    //! Constructor
	/*!
	 * Non-empty constructor
	 * \param fileName the name of the HDF5 file to be used
	 * \param comm pointer to Epetra_Comm
	 * \param transposeInFile (default false) write/read tables transposed
	 * 		  in the HDF5. This option is available because the transposed
	 * 		  format may be faster on certain machines. Leave as default
	 * 		  unless certain that it works for a given machine.
	 */
	PartitionIO(const std::string& fileName,
				const commPtr_Type& comm,
				const bool transposeInFile = false);

    //! Empty destructor
    virtual ~PartitionIO() {}
    //@}

    //! \name Public Methods
    //@{
    //! Initialization method
	/*!
	 * Initialization method that is used with the default constructor
	 * \param fileName the name of the HDF5 file to be used
	 * \param comm pointer to Epetra_Comm
	 * \param transposeInFile (default false) write/read tables transposed
	 * 		  in the HDF5. This option is available because the transposed
	 * 		  format may be faster on certain machines. Leave as default
	 * 		  unless certain that it works for a given machine.
	 */
    void setup(const std::string& fileName,
    		   const commPtr_Type& comm,
    		   const bool transposeInFile = false);
    //! Write method
    /*!
     * Call this method to write the mesh parts to disk
     * \param meshParts pointer to a vector containing pointers to the mesh
     * 		  parts (RegionMesh objects). This pointer is released after
     * 		  writing, so the mesh parts can be cleared from memory
     */
    void write(const meshPartsPtr_Type& meshParts);
    //! Read method
    /*!
     * Call this method to read from the HDF5 file the mesh part associated
     * with the rank of the current MPI process
     * \param meshPart pointer to a RegionMesh that will contain the mesh
     * 		  part. If the RegionMesh has been initialized and contains any
     * 		  state, it will be destroyed before reading
     */
    void read(meshPtr_Type& meshPart);
    //@}

private:
    // Private copy constructor and assignment operator are disabled
    PartitionIO(const PartitionIO&);
    PartitionIO& operator=(const PartitionIO&);

    //! Private Methods
    //@{
    // Methods for writing
    void createHDF5File();
    void writeStats();
    void writePoints();
    void writeEdges();
    void writeFaces();
    void writeElements();
    void closeHDF5File();
    // Methods for reading
    void openHDF5File();
    void readStats();
    void readPoints();
    void readEdges();
    void readFaces();
    void readElements();
	//@}

    //! Private Data Members
    //@{
    boost::shared_ptr<Epetra_Comm> M_comm;
    UInt M_myRank;
    UInt M_numProc;
    std::string M_fileName;
    bool M_transposeInFile;
    meshPartsPtr_Type M_meshPartsOut;
    meshPtr_Type M_meshPartIn;
    // Mesh geometry
    UInt M_elementNodes;
    UInt M_faceNodes;
    UInt M_maxNumPoints;
    UInt M_maxNumEdges;
    UInt M_maxNumFaces;
    UInt M_maxNumElements;
    UInt M_numParts;
    // Counters when loading the mesh part
    UInt M_numPoints;
    UInt M_numEdges;
    UInt M_numFaces;
    UInt M_numElements;
    // HDF5 file handle
    hid_t M_fileId;
    // Buffers for reading/writing
    std::vector<Int> M_intBuffer;
    std::vector<Real> M_realBuffer;
    //@}
}; // class PartitionIO

} /* namespace LifeV */

template<typename MeshType>
inline LifeV::PartitionIO<MeshType>::PartitionIO(const std::string& fileName,
												 const commPtr_Type& comm,
												 const bool transposeInFile) :
	M_comm(comm),
	M_fileName(fileName),
	M_transposeInFile(transposeInFile),
	M_maxNumPoints(0),
	M_maxNumEdges(0),
	M_maxNumFaces(0),
	M_maxNumElements(0)
{
	// Mesh geometry
	// we should directly use the data structures provided by the library
    // M_elementNodes = MeshType::elementShape_Type::S_numPoints;
    // M_faceNodes = MeshType::elementShape_Type::GeoBShape::S_numPoints;
    switch (MeshType::elementShape_Type::S_shape)
    {
    case HEXA:
        M_elementNodes = 8;
        M_faceNodes    = 4;
        break;
    case TETRA:
        M_elementNodes = 4;
        M_faceNodes    = 3;
        break;
    }
}

template<typename MeshType>
inline void LifeV::PartitionIO<MeshType>::setup(const std::string& fileName,
												const commPtr_Type& comm,
												const bool transposeInFile)
{
	M_comm = comm;
	M_fileName = fileName;
	M_transposeInFile = transposeInFile;
	M_maxNumPoints = 0;
	M_maxNumEdges = 0;
	M_maxNumFaces = 0;
	M_maxNumElements = 0;

    // Mesh geometry
	// we should directly use the data structures provided by the library
    // M_elementNodes = MeshType::elementShape_Type::S_numPoints;
    // M_faceNodes = MeshType::elementShape_Type::GeoBShape::S_numPoints;
    switch (MeshType::elementShape_Type::S_shape)
    {
    case HEXA:
        M_elementNodes = 8;
        M_faceNodes    = 4;
        break;
    case TETRA:
        M_elementNodes = 4;
        M_faceNodes    = 3;
        break;
    }
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::write(const meshPartsPtr_Type& meshParts)
{
	M_meshPartsOut = meshParts;
	M_numParts = M_meshPartsOut->size();

	createHDF5File();
	writeStats();
	writePoints();
	writeEdges();
	writeFaces();
	writeElements();
	closeHDF5File();

	M_meshPartsOut.reset();
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::read(meshPtr_Type& meshPart)
{
	meshPart.reset();
	M_meshPartIn.reset(new mesh_Type);

	openHDF5File();
	readStats();
	readPoints();
	readEdges();
	readFaces();
	readElements();
	closeHDF5File();

	meshPart = M_meshPartIn;
	M_meshPartIn.reset();
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::createHDF5File()
{
	hid_t plistId;
	MPI_Comm comm;
	// why is this cast needed?
	boost::shared_ptr<Epetra_MpiComm> tempComm =
			boost::dynamic_pointer_cast<Epetra_MpiComm>(M_comm);
	if (tempComm != 0) {
		comm = tempComm->Comm();
	}
	MPI_Info info = MPI_INFO_NULL;

    // Set up file access property list with parallel I/O access
    plistId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plistId, comm, info);

    // Create a new file collectively and release property list identifier.
    M_fileId = H5Fcreate(M_fileName.c_str(), H5F_ACC_TRUNC,
    				     H5P_DEFAULT, plistId);
    H5Pclose(plistId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::writeStats()
{
    // Write mesh partition stats (N = number of parts)
    // This is an N x 15 table of int
	hsize_t currentSpaceDims[2];
    // you could avoid the if statement:
    //    currentSpaceDims[!M_transposeInFile] = 15;
    //    currentSpaceDims[M_transposeInFile] = M_numParts;
	if (! M_transposeInFile) {
		currentSpaceDims[0] = M_numParts;
		currentSpaceDims[1] = 15;
	} else {
		currentSpaceDims[0] = 15;
		currentSpaceDims[1] = M_numParts;
	}
    hid_t filespace = H5Screate_simple(2, currentSpaceDims, currentSpaceDims);
    hid_t intDataset = H5Dcreate(M_fileId, "stats", H5T_STD_I32BE, filespace,
						  	  	 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t currentCount[2];
    // you could avoid the if statement:
    //    currentCount[M_transposeInFile] = 1;
    //    currentCount[!M_transposeInFile] = currentSpaceDims[!M_transposeInFile];
	if (! M_transposeInFile) {
		currentCount[0] = 1;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 1;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    // Fill buffer
	M_intBuffer.resize(15);
    for (UInt i = 0; i < M_numParts; ++i) {
    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
    	M_intBuffer[0] = M_numParts;
    	// Next one is unused. I'll put it to keep compatibility
    	// in case I decide to write the graphs
    	M_intBuffer[1] = 0;
    	M_intBuffer[2] = currentPart.numPoints();
    	if (M_intBuffer[2] > M_maxNumPoints) M_maxNumPoints = M_intBuffer[2];
        M_intBuffer[3] = currentPart.numBPoints();
        M_intBuffer[4] = currentPart.numVertices();
        M_intBuffer[5] = currentPart.numBVertices();
        M_intBuffer[6] = currentPart.numGlobalVertices();
        M_intBuffer[7] = currentPart.numEdges();
        if (M_intBuffer[7] > M_maxNumEdges) M_maxNumEdges = M_intBuffer[7];
        M_intBuffer[8] = currentPart.numBEdges();
        M_intBuffer[9] = currentPart.numGlobalEdges();
        M_intBuffer[10] = currentPart.numFaces();
        if (M_intBuffer[10] > M_maxNumFaces) M_maxNumFaces = M_intBuffer[10];
        M_intBuffer[11] = currentPart.numBFaces();
        M_intBuffer[12] = currentPart.numGlobalFaces();
        M_intBuffer[13] = currentPart.numVolumes();
        if (M_intBuffer[13] > M_maxNumElements)
        	M_maxNumElements = M_intBuffer[13];
		M_intBuffer[14] = currentPart.numGlobalVolumes();

		hsize_t currentOffset[2];
	    // you could avoid the if statement:
	    //    currentOffset[M_transposeInFile] = i;
	    //    currentOffset[!M_transposeInFile] = 0;
		if (! M_transposeInFile) {
			currentOffset[0] = i;
			currentOffset[1] = 0;
		} else {
			currentOffset[0] = 0;
			currentOffset[1] = i;
		}
		H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
          					currentCount, NULL);
        H5Dwrite(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
        		 &M_intBuffer[0]);
    }

	// HDF5 cleanup
	H5Dclose(intDataset);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plistId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::writePoints()
{
    // Write points (N = number of parts)
    // Two tables: one is (3 * N) x max_num_points of int
    //			   second is (3 * N) x max_num_points of real
	hsize_t currentSpaceDims[2];
	// you could avoid the if statement
	if (! M_transposeInFile) {
		currentSpaceDims[0] = 3 * M_numParts;
		currentSpaceDims[1] = M_maxNumPoints;
	} else {
		currentSpaceDims[0] = M_maxNumPoints;
		currentSpaceDims[1] = 3 * M_numParts;
	}
    hid_t filespace = H5Screate_simple(2, currentSpaceDims, currentSpaceDims);
    hid_t intDataset = H5Dcreate(M_fileId, "point_ids", H5T_STD_I32BE,
    							 filespace, H5P_DEFAULT, H5P_DEFAULT,
    							 H5P_DEFAULT);
    hid_t realDataset = H5Dcreate(M_fileId, "point_coords", H5T_IEEE_F64BE,
								  filespace, H5P_DEFAULT, H5P_DEFAULT,
								  H5P_DEFAULT);

	hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 3;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 3;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    // Fill buffer
	M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
    M_realBuffer.resize(currentCount[0] * currentCount[1], 0);
    UInt stride = currentCount[1];
	if (! M_transposeInFile) {
		for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numPoints(); ++j) {
	    		M_intBuffer[j] = currentPart.pointList[j].marker();
	    		M_intBuffer[stride + j] = currentPart.point(j).id();
	    		// we should save the flag (more informative) rather than a bool
	    		M_intBuffer[2 * stride + j] =
						static_cast<int>(currentPart.isBoundaryPoint(j));
	    		M_realBuffer[j] = currentPart.pointList[j].x();
	    		M_realBuffer[stride + j] = currentPart.pointList[j].y();
	    		M_realBuffer[2 * stride + j] = currentPart.pointList[j].z();
	    	}

	    	hsize_t currentOffset[2] = {i * currentCount[0], 0};
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
								currentCount, NULL);
			H5Dwrite(realDataset, H5T_NATIVE_DOUBLE, memspace, filespace,
					 plistId, &M_realBuffer[0]);
	    }
	} else {
	    for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numPoints(); ++j) {
	    		M_intBuffer[stride * j] = currentPart.pointList[j].marker();
	    		M_intBuffer[stride * j + 1] = currentPart.point(j).id();
                // we should save the flag (more informative) rather than a bool
	    		M_intBuffer[stride * j + 2] =
						static_cast<int>(currentPart.isBoundaryPoint(j));
	    		M_realBuffer[stride * j] = currentPart.pointList[j].x();
	    		M_realBuffer[stride * j + 1] = currentPart.pointList[j].y();
	    		M_realBuffer[stride * j + 2] = currentPart.pointList[j].z();
	    	}

	    	hsize_t currentOffset[2] = {0, i * currentCount[1]};
	        // common calls to H5 routines could be factored out of the if..else statement
	    	// to avoid code duplication
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
								currentCount, NULL);
			H5Dwrite(realDataset, H5T_NATIVE_DOUBLE, memspace, filespace,
					 plistId, &M_realBuffer[0]);
	    }
	}
    M_realBuffer.resize(0);

    // HDF5 cleanup
    H5Dclose(intDataset);
    H5Dclose(realDataset);
    H5Sclose(filespace);
	H5Sclose(memspace);
    H5Pclose(plistId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::writeEdges()
{
    // Write edges (N = number of parts)
    // Table is (5 * N) x max_num_edges of int
    hsize_t currentSpaceDims[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentSpaceDims[0] = 5 * M_numParts;
		currentSpaceDims[1] = M_maxNumEdges;
	} else {
		currentSpaceDims[0] = M_maxNumEdges;
		currentSpaceDims[1] = 5 * M_numParts;
	}
    hid_t filespace = H5Screate_simple(2, currentSpaceDims, currentSpaceDims);
    hid_t dataset = H5Dcreate(M_fileId, "edges", H5T_STD_I32BE, filespace,
							  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 5;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 5;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    // Fill buffer
	M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
	UInt stride = currentCount[1];
	if (! M_transposeInFile) {
	    for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numEdges(); ++j) {
	    		M_intBuffer[j] = currentPart.edgeList[j].point(0).localId();
	    		M_intBuffer[stride + j] =
						currentPart.edgeList[j].point(1).localId();
	    		M_intBuffer[2 * stride + j] = currentPart.edgeList[j].marker();
	    		M_intBuffer[3 * stride + j] = currentPart.edgeList[j].id();
                // we should save the flag (more informative) rather than a bool
	    		M_intBuffer[4 * stride + j] =
	    				static_cast<int>(currentPart.isBoundaryEdge(j));
			}

	    	hsize_t currentOffset[2] = {i * currentCount[0], 0};
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
		}
	} else {
	    for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numEdges(); ++j) {
	    		M_intBuffer[stride * j] =
	    				currentPart.edgeList[j].point(0).localId();
	    		M_intBuffer[stride * j + 1] =
						currentPart.edgeList[j].point(1).localId();
	    		M_intBuffer[stride * j + 2] = currentPart.edgeList[j].marker();
	    		M_intBuffer[stride * j + 3] = currentPart.edgeList[j].id();
                // we should save the flag (more informative) rather than a bool
	    		M_intBuffer[stride * j + 4] =
	    				static_cast<int>(currentPart.isBoundaryEdge(j));
			}

	    	hsize_t currentOffset[2] = {0, i * currentCount[1]};
            // common calls to H5 routines could be factored out of the if..else statement
            // to avoid code duplication
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
		}
	}

    // HDF5 cleanup
    H5Dclose(dataset);
    H5Sclose(filespace);
	H5Sclose(memspace);
    H5Pclose(plistId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::writeFaces()
{
    // Write faces (N = number of parts)
    // Table is ((7 + num_face_nodes * N) x max_num_faces of int
    hsize_t currentSpaceDims[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentSpaceDims[0] = (7 + M_faceNodes) * M_numParts;
		currentSpaceDims[1] = M_maxNumFaces;
	} else {
		currentSpaceDims[0] = M_maxNumFaces;
		currentSpaceDims[1] = (7 + M_faceNodes) * M_numParts;
	}
    hid_t filespace = H5Screate_simple(2, currentSpaceDims, currentSpaceDims);
    hid_t dataset = H5Dcreate(M_fileId, "faces", H5T_STD_I32BE, filespace,
							  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 7 + M_faceNodes;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 7 + M_faceNodes;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    // Fill buffer
	M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
	UInt stride = currentCount[1];
	if (! M_transposeInFile) {
	    for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numFaces(); ++j) {
	    		for (UInt k = 0; k < M_faceNodes; ++k) {
	    			M_intBuffer[k * stride + j] =
						currentPart.faceList[j].point(k).localId();
	    		}
	    		M_intBuffer[M_faceNodes * stride +j] =
						currentPart.faceList[j].marker();
	    		M_intBuffer[(M_faceNodes + 1) * stride + j] =
						currentPart.faceList[j].id();
	    		M_intBuffer[(M_faceNodes + 2) * stride + j] =
						currentPart.faceList[j].firstAdjacentElementIdentity();
	    		M_intBuffer[(M_faceNodes + 3) * stride + j] =
						currentPart.faceList[j].secondAdjacentElementIdentity();
	    		M_intBuffer[(M_faceNodes + 4) * stride + j] =
						currentPart.faceList[j].firstAdjacentElementPosition();
	    		M_intBuffer[(M_faceNodes + 5) * stride + j] =
						currentPart.faceList[j].secondAdjacentElementPosition();
                // we should save the flag (more informative) rather than a bool
	    		M_intBuffer[(M_faceNodes + 6) * stride + j] =
						static_cast<int>(currentPart.isBoundaryFace(j));
			}

	    	hsize_t currentOffset[2] = {i * currentCount[0], 0};
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
		}
	} else {
	    for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numFaces(); ++j) {
	    		for (UInt k = 0; k < M_faceNodes; ++k) {
	    			M_intBuffer[stride * j + k] =
						currentPart.faceList[j].point(k).localId();
	    		}
	    		M_intBuffer[stride * j + M_faceNodes] =
						currentPart.faceList[j].marker();
	    		M_intBuffer[stride * j + M_faceNodes + 1] =
						currentPart.faceList[j].id();
	    		M_intBuffer[stride * j + M_faceNodes + 2] =
						currentPart.faceList[j].firstAdjacentElementIdentity();
	    		M_intBuffer[stride * j + M_faceNodes + 3] =
						currentPart.faceList[j].secondAdjacentElementIdentity();
	    		M_intBuffer[stride * j + M_faceNodes + 4] =
						currentPart.faceList[j].firstAdjacentElementPosition();
	    		M_intBuffer[stride * j + M_faceNodes + 5] =
						currentPart.faceList[j].secondAdjacentElementPosition();
                // we should save the flag (more informative) rather than a bool
	    		M_intBuffer[stride * j + M_faceNodes + 6] =
						static_cast<int>(currentPart.isBoundaryFace(j));
			}

	    	hsize_t currentOffset[2] = {0, i * currentCount[1]};
            // common calls to H5 routines could be factored out of the if..else statement
            // to avoid code duplication
	    	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
		}
	}

    // HDF5 cleanup
    H5Dclose(dataset);
    H5Sclose(filespace);
	H5Sclose(memspace);
    H5Pclose(plistId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::writeElements()
{
    // Write faces (N = number of parts)
    // Table is ((2 + num_element_nodes * N) x max_num_elements of int
    hsize_t currentSpaceDims[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentSpaceDims[0] = (2 + M_elementNodes) * M_numParts;
		currentSpaceDims[1] = M_maxNumElements;
	} else {
		currentSpaceDims[0] = M_maxNumElements;
		currentSpaceDims[1] = (2 + M_elementNodes) * M_numParts;
	}
    hid_t filespace = H5Screate_simple(2, currentSpaceDims, currentSpaceDims);
    hid_t dataset = H5Dcreate(M_fileId, "elements", H5T_STD_I32BE, filespace,
    						  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 2 + M_elementNodes;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];;
		currentCount[1] = 2 + M_elementNodes;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    // Fill buffer
	M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
	UInt stride = currentCount[1];
	if (! M_transposeInFile) {
	    for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numVolumes(); ++j) {
	    		for (UInt k = 0; k < M_elementNodes; ++k) {
	    			M_intBuffer[k * stride + j] =
						currentPart.volumeList[j].point(k).localId();
	    		}
	    		M_intBuffer[M_elementNodes * stride +j] =
						currentPart.volumeList[j].marker();
	    		M_intBuffer[(M_elementNodes + 1) * stride + j] =
						currentPart.volumeList[j].id();
	    		}

	    	hsize_t currentOffset[2] = {i * currentCount[0], 0};
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
		}
	} else {
	    for (UInt i = 0; i < M_numParts; ++i) {
	    	mesh_Type& currentPart = (*(*M_meshPartsOut)[i]);
	    	for (UInt j = 0; j < currentPart.numVolumes(); ++j) {
	    		for (UInt k = 0; k < M_elementNodes; ++k) {
	    			M_intBuffer[stride * j + k] =
						currentPart.volumeList[j].point(k).localId();
	    		}
	    		M_intBuffer[stride * j + M_elementNodes] =
						currentPart.volumeList[j].marker();
	    		M_intBuffer[stride * j + M_elementNodes + 1] =
						currentPart.volumeList[j].id();
	    		}

	    	hsize_t currentOffset[2] = {0, i * currentCount[1]};
            // common calls to H5 routines could be factored out of the if..else statement
            // to avoid code duplication
	        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
	          					currentCount, NULL);
	        H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plistId,
	        		 &M_intBuffer[0]);
		}
	}

    // HDF5 cleanup
    H5Dclose(dataset);
    H5Sclose(filespace);
	H5Sclose(memspace);
    H5Pclose(plistId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::closeHDF5File()
{
	H5Fclose(M_fileId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::openHDF5File()
{
	hid_t plistId;
	MPI_Comm comm;
	boost::shared_ptr<Epetra_MpiComm> tempComm =
			boost::dynamic_pointer_cast<Epetra_MpiComm>(M_comm);
	if (tempComm != 0) {
		comm = tempComm->Comm();
		M_myRank = tempComm->MyPID();
		M_numProc = tempComm->NumProc();
	}
	MPI_Info info = MPI_INFO_NULL;

    // Set up file access property list with parallel I/O access
    plistId = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plistId, comm, info);

    // Create a new file collectively and release property list identifier.
    M_fileId = H5Fopen(M_fileName.c_str(), H5F_ACC_RDONLY, plistId);
    H5Pclose(plistId);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::readStats()
{
    // Write mesh partition stats (N = number of parts)
    // This is an N x 15 table of int
	hsize_t currentSpaceDims[2];
    hid_t intDataset = H5Dopen(M_fileId, "stats", H5P_DEFAULT);
    hid_t filespace = H5Dget_space(intDataset);
    H5Sget_simple_extent_dims(filespace, currentSpaceDims, NULL);

    hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 1;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 1;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    // Fill buffer
	M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
	hsize_t currentOffset[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentOffset[0] = M_myRank;
		currentOffset[1] = 0;
	} else {
		currentOffset[0] = 0;
		currentOffset[1] = M_myRank;
	}
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
						currentCount, NULL);
	H5Dread(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
			&M_intBuffer[0]);

	// HDF5 cleanup
	H5Dclose(intDataset);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plistId);

	// Insert stats into mesh partition object
    M_numPoints = M_intBuffer[2];
    M_meshPartIn->setMaxNumPoints(M_intBuffer[2], true);
    M_meshPartIn->setNumBPoints(M_intBuffer[3]);

    // Vertices
    M_meshPartIn->setNumVertices(M_intBuffer[4]);
    M_meshPartIn->setNumBVertices(M_intBuffer[5]);
    M_meshPartIn->setNumGlobalVertices(M_intBuffer[6]);

    // Edges
    M_numEdges = M_intBuffer[7];
    M_meshPartIn->setNumEdges(M_intBuffer[7]);
    M_meshPartIn->setMaxNumEdges(M_intBuffer[7], true);
    M_meshPartIn->setNumBEdges(M_intBuffer[8]);
    M_meshPartIn->setMaxNumGlobalEdges(M_intBuffer[9]);

    // Faces
    M_numFaces = M_intBuffer[10];
    M_meshPartIn->setNumFaces(M_intBuffer[10]);
    M_meshPartIn->setMaxNumFaces(M_intBuffer[10], true);
    M_meshPartIn->setNumBFaces(M_intBuffer[11]);
    M_meshPartIn->setMaxNumGlobalFaces(M_intBuffer[12]);

    // Volumes
    M_numElements = M_intBuffer[13];
    M_meshPartIn->setMaxNumVolumes(M_intBuffer[13], true);
    M_meshPartIn->setMaxNumGlobalVolumes(M_intBuffer[14]);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::readPoints()
{
    // Write mesh points (N = number of parts)
    // There are two tables: a (3 * N) x max_num_points table of int and
	// a (3 * N) x max_num_points table of real
	hsize_t currentSpaceDims[2];
    hid_t intDataset = H5Dopen(M_fileId, "point_ids", H5P_DEFAULT);
    hid_t realDataset = H5Dopen(M_fileId, "point_coords", H5P_DEFAULT);

    hid_t filespace = H5Dget_space(intDataset);
    H5Sget_simple_extent_dims(filespace, currentSpaceDims, NULL);

    hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 3;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 3;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
    M_realBuffer.resize(currentCount[0] * currentCount[1], 0);
    hsize_t currentOffset[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
	    currentOffset[0] = currentCount[0] * M_myRank;
	    currentOffset[1] = 0;
	} else {
	    currentOffset[0] = 0;
	    currentOffset[1] = currentCount[1] * M_myRank;
	}
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
    					currentCount, NULL);
	H5Dread(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
			&M_intBuffer[0]);
	H5Dread(realDataset, H5T_NATIVE_DOUBLE, memspace, filespace, plistId,
			&M_realBuffer[0]);

	// HDF5 cleanup
	H5Dclose(intDataset);
	H5Dclose(realDataset);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plistId);

	// Insert points into the mesh partition object
    M_meshPartIn->pointList.reserve(M_numPoints);
    M_meshPartIn->_bPoints.reserve(M_meshPartIn->numBPoints());

    typename MeshType::point_Type *pp = 0;

    UInt stride = currentCount[1];
	if (! M_transposeInFile) {
	    for (UInt j = 0; j < M_numPoints; ++j)
	    {
	    	pp = &(M_meshPartIn->addPoint(
					static_cast<bool>(M_intBuffer[2 * stride + j])));
	        pp->setMarker(M_intBuffer[j]);
	        pp->x() = M_realBuffer[j];
	        pp->y() = M_realBuffer[stride + j];
	        pp->z() = M_realBuffer[2 * stride + j];
	        pp->setLocalId(j);
	        pp->setId(M_intBuffer[stride + j]);
	    }
	} else {
	    for (UInt j = 0; j < M_numPoints; ++j)
	    {
	    	pp = &(M_meshPartIn->addPoint(
					static_cast<bool>(M_intBuffer[stride * j + 2])));
	        pp->setMarker(M_intBuffer[stride * j]);
	        pp->x() = M_realBuffer[stride * j];
	        pp->y() = M_realBuffer[stride * j + 1];
	        pp->z() = M_realBuffer[stride * j + 2];
	        pp->setLocalId(j);
	        pp->setId(M_intBuffer[stride * j + 1]);
	    }
	}
    M_realBuffer.resize(0);
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::readEdges()
{
    // Write mesh edges (N = number of parts)
    // Read a (5 * N) x max_num_edges table of int and
	hsize_t currentSpaceDims[2];
    hid_t intDataset = H5Dopen(M_fileId, "edges", H5P_DEFAULT);

    hid_t filespace = H5Dget_space(intDataset);
    H5Sget_simple_extent_dims(filespace, currentSpaceDims, NULL);

    hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 5;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 5;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
    hsize_t currentOffset[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
	    currentOffset[0] = currentCount[0] * M_myRank;
	    currentOffset[1] = 0;
	} else {
	    currentOffset[0] = 0;
	    currentOffset[1] = currentCount[1] * M_myRank;
	}
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
						currentCount, NULL);
	H5Dread(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
			&M_intBuffer[0]);

	// HDF5 cleanup
	H5Dclose(intDataset);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plistId);

    M_meshPartIn->edgeList.reserve(M_numEdges);

    typename MeshType::edge_Type *pe;

    UInt stride = currentCount[1];
	if (! M_transposeInFile) {
	    for (UInt j = 0; j < M_numEdges; ++j)
	    {
	        pe = &(M_meshPartIn->addEdge(
	        		static_cast<bool>(M_intBuffer[4 * stride + j])));
	        pe->setLocalId(j);
	        pe->setId(M_intBuffer[3 * stride + j]);
	        pe->setPoint(0, M_meshPartIn->point(M_intBuffer[j]));
	        pe->setPoint(1, M_meshPartIn->point(M_intBuffer[stride +j]));
	        pe->setMarker(M_intBuffer[2 * stride + j]);
	    }
	} else {
	    for (UInt j = 0; j < M_numEdges; ++j)
	    {
	        pe = &(M_meshPartIn->addEdge(
	        		static_cast<bool>(M_intBuffer[stride * j + 4])));
	        pe->setLocalId(j);
	        pe->setId(M_intBuffer[stride * j + 3]);
	        pe->setPoint(0, M_meshPartIn->point(M_intBuffer[stride * j]));
	        pe->setPoint(1, M_meshPartIn->point(M_intBuffer[stride * j + 1]));
	        pe->setMarker(M_intBuffer[stride * j + 2]);
	    }
	}
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::readFaces()
{
    // Write mesh faces (N = number of parts)
    // Read a ((7 + num_face_points) * N) x max_num_faces table of int
	hsize_t currentSpaceDims[2];
    hid_t intDataset = H5Dopen(M_fileId, "faces", H5P_DEFAULT);

    hid_t filespace = H5Dget_space(intDataset);
    H5Sget_simple_extent_dims(filespace, currentSpaceDims, NULL);

    hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 7 + M_faceNodes;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 7 + M_faceNodes;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
    hsize_t currentOffset[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
	    currentOffset[0] = currentCount[0] * M_myRank;
	    currentOffset[1] = 0;
	} else {
	    currentOffset[0] = 0;
	    currentOffset[1] = currentCount[1] * M_myRank;
	}
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
						currentCount, NULL);
	H5Dread(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
			&M_intBuffer[0]);

	// HDF5 cleanup
	H5Dclose(intDataset);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plistId);

    typename MeshType::face_Type *pf = 0;

    M_meshPartIn->faceList.reserve(M_numFaces);

    UInt stride = currentCount[1];
	if (! M_transposeInFile) {
	    for (UInt j = 0; j < M_numFaces; ++j)
	    {
	        pf = &(M_meshPartIn->addFace(
	        		static_cast<bool>(M_intBuffer[(6 + M_faceNodes)
	        		                              * stride + j])));
	        pf->setLocalId(j);
	        pf->setId(M_intBuffer[(M_faceNodes + 1) * stride + j]);
	        pf->firstAdjacentElementIdentity() =
	        		M_intBuffer[(M_faceNodes + 2) * stride + j];
	        pf->secondAdjacentElementIdentity() =
	        		M_intBuffer[(M_faceNodes + 3) * stride + j];
	        pf->firstAdjacentElementPosition() =
					M_intBuffer[(M_faceNodes + 4) * stride + j];
	        pf->secondAdjacentElementPosition() =
	        		M_intBuffer[(M_faceNodes + 5) * stride + j];
	        pf->setMarker(M_intBuffer[M_faceNodes * stride + j]);
	        for (UInt k = 0; k < M_faceNodes; ++k)
	        {
	            pf->setPoint(k, M_meshPartIn->point(
	            		M_intBuffer[k * stride + j]));
	        }
	    }
	} else {
	    for (UInt j = 0; j < M_numFaces; ++j)
	    {
	        pf = &(M_meshPartIn->addFace(
	        		static_cast<bool>(M_intBuffer[stride * j
	        		                              + M_faceNodes + 6])));
	        pf->setLocalId(j);
	        pf->setId(M_intBuffer[stride * j + M_faceNodes + 1]);
	        pf->firstAdjacentElementIdentity() =
	        		M_intBuffer[stride * j + M_faceNodes + 2];
	        pf->secondAdjacentElementIdentity() =
	        		M_intBuffer[stride * j + M_faceNodes + 3];
	        pf->firstAdjacentElementPosition() =
					M_intBuffer[stride * j + M_faceNodes + 4];
	        pf->secondAdjacentElementPosition() =
	        		M_intBuffer[stride * j + M_faceNodes + 5];
	        pf->setMarker(M_intBuffer[(7 + M_faceNodes) * j + M_faceNodes]);
	        for (UInt k = 0; k < M_faceNodes; ++k)
	        {
	            pf->setPoint(k, M_meshPartIn->point(
	            		M_intBuffer[stride * j + k]));
	        }
	    }
	}

    M_meshPartIn->setLinkSwitch("HAS_ALL_FACETS");
    M_meshPartIn->setLinkSwitch("FACETS_HAVE_ADIACENCY");
}

template<typename MeshType>
void LifeV::PartitionIO<MeshType>::readElements()
{
    // Write mesh faces (N = number of parts)
    // Read a ((2 + num_element_points) * N) x max_num_elements table of int
	hsize_t currentSpaceDims[2];
    hid_t intDataset = H5Dopen(M_fileId, "elements", H5P_DEFAULT);

    hid_t filespace = H5Dget_space(intDataset);
    H5Sget_simple_extent_dims(filespace, currentSpaceDims, NULL);

    hsize_t currentCount[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentCount[0] = 2 + M_elementNodes;
		currentCount[1] = currentSpaceDims[1];
	} else {
		currentCount[0] = currentSpaceDims[0];
		currentCount[1] = 2 + M_elementNodes;
	}
    hid_t memspace = H5Screate_simple(2, currentCount, currentCount);

    hid_t plistId = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

    M_intBuffer.resize(currentCount[0] * currentCount[1], 0);
    hsize_t currentOffset[2];
    // you could avoid the if statement
	if (! M_transposeInFile) {
		currentOffset[0] = currentCount[0] * M_myRank;
		currentOffset[1] = 0;
	} else {
		currentOffset[0] = 0;
		currentOffset[1] = currentCount[1] * M_myRank;
	}
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, currentOffset, NULL,
						currentCount, NULL);
	H5Dread(intDataset, H5T_NATIVE_INT, memspace, filespace, plistId,
			&M_intBuffer[0]);

	// HDF5 cleanup
	H5Dclose(intDataset);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plistId);

    M_meshPartIn->volumeList.reserve(M_numElements);

    typename MeshType::volume_Type *pv = 0;

    UInt stride = currentCount[1];
	if (! M_transposeInFile) {
	    for (UInt j = 0; j < M_numElements; ++j)
	    {
	        pv = &(M_meshPartIn->addVolume());
	        pv->setId(M_intBuffer[(M_elementNodes + 1) * stride + j]);
	        pv->setLocalId(j);
	        for (UInt k = 0; k < M_elementNodes; ++k)
	        {
	            pv->setPoint(k, M_meshPartIn->point(
	            		M_intBuffer[k * stride + j]));
	        }
	        pv->setMarker(M_intBuffer[M_elementNodes * stride + j]);
	    }
	} else {
	    for (UInt j = 0; j < M_numElements; ++j)
	    {
	        pv = &(M_meshPartIn->addVolume());
	        pv->setId(M_intBuffer[stride * j + M_elementNodes + 1]);
	        pv->setLocalId(j);
	        for (UInt k = 0; k < M_elementNodes; ++k)
	        {
	            pv->setPoint(k, M_meshPartIn->point(
	            		M_intBuffer[stride * j + k]));
	        }
	        pv->setMarker(M_intBuffer[stride * j + M_elementNodes]);
	    }
	}

    M_meshPartIn->updateElementEdges(false, false);
    M_meshPartIn->updateElementFaces(false, false);
}

#endif /* HAVE_MPI */
#endif /* HAVE_HDF5 */
#endif /* PARTITION_IO_H_ */
