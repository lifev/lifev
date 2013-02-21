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
  @brief Convenience wrapper for the C interface of the HDF5 library

  @date 8-06-2012
  @author Radu Popescu <radu.popescu@epfl.ch>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#include <lifev/core/filter/HDF5IO.hpp>

#ifdef HAVE_HDF5
#ifdef HAVE_MPI

// ===================================================
// Constructor
// ===================================================

LifeV::HDF5IO::HDF5IO (const std::string& fileName, const commPtr_Type& comm,
                       const bool& existing)
{
    openFile (fileName, comm, existing);
}

// ===================================================
// Public Methods
// ===================================================

void LifeV::HDF5IO::openFile (const std::string& fileName,
                              const commPtr_Type& comm,
                              const bool& existing)
{
    hid_t plistId;
    MPI_Comm mpiComm;
    boost::shared_ptr<Epetra_MpiComm> tempComm =
        boost::dynamic_pointer_cast<Epetra_MpiComm> (comm);
    if (tempComm != 0)
    {
        mpiComm = tempComm->Comm();
    }
    MPI_Info info = MPI_INFO_NULL;

    // Set up file access property list with parallel I/O access
    plistId = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio (plistId, mpiComm, info);

    // Create/open a file collectively and release property list identifier.
    if (existing)
    {
        M_fileId = H5Fopen (fileName.c_str(), H5F_ACC_RDONLY, plistId);
    }
    else
    {
        M_fileId = H5Fcreate (fileName.c_str(), H5F_ACC_TRUNC,
                              H5P_DEFAULT, plistId);
    }
    H5Pclose (plistId);
}

void LifeV::HDF5IO::createTable (const std::string& tableName,
                                 hid_t& fileDataType,
                                 hsize_t tableDimensions[])
{
    tableHandle& currentTable = M_tableList[tableName];

    currentTable.filespace = H5Screate_simple (2, tableDimensions,
                                               tableDimensions);
    currentTable.dataset = H5Dcreate (M_fileId, tableName.c_str(), fileDataType,
                                      currentTable.filespace, H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_COLLECTIVE);
}

void LifeV::HDF5IO::openTable (const std::string& tableName,
                               hsize_t tableDimensions[])
{
    tableHandle& currentTable = M_tableList[tableName];

    currentTable.dataset = H5Dopen (M_fileId, tableName.c_str(), H5P_DEFAULT);
    currentTable.filespace = H5Dget_space (currentTable.dataset);
    H5Sget_simple_extent_dims (currentTable.filespace, tableDimensions, NULL);
    currentTable.plist = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (currentTable.plist, H5FD_MPIO_COLLECTIVE);
}

void LifeV::HDF5IO::write (const std::string& tableName,
                           hid_t& memDataType, hsize_t currentCount[],
                           hsize_t currentOffset[], void* buffer)
{
    tableHandle& currentTable = M_tableList[tableName];

    hid_t memspace = H5Screate_simple (2, currentCount, currentCount);

    H5Sselect_hyperslab (currentTable.filespace, H5S_SELECT_SET, currentOffset,
                         NULL, currentCount, NULL);
    H5Dwrite (currentTable.dataset, memDataType, memspace,
              currentTable.filespace, currentTable.plist, buffer);

    H5Sclose (memspace);
}

void LifeV::HDF5IO::read (const std::string& tableName,
                          hid_t& memDataType, hsize_t currentCount[],
                          hsize_t currentOffset[], void* buffer)
{
    tableHandle& currentTable = M_tableList[tableName];

    hid_t memspace = H5Screate_simple (2, currentCount, currentCount);

    H5Sselect_hyperslab (currentTable.filespace, H5S_SELECT_SET, currentOffset,
                         NULL, currentCount, NULL);
    H5Dread (currentTable.dataset, memDataType, memspace,
             currentTable.filespace, currentTable.plist,
             buffer);

    H5Sclose (memspace);
}

void LifeV::HDF5IO::closeTable (const std::string& tableName)
{
    tableHandle& currentTable = M_tableList[tableName];
    H5Sclose (currentTable.filespace);
    H5Dclose (currentTable.dataset);
    H5Pclose (currentTable.plist);
    M_tableList.erase (tableName);
}

void LifeV::HDF5IO::closeFile()
{
    H5Fclose (M_fileId);
}

#endif /* HAVE_MPI */
#endif /* HAVE_HDF5 */
