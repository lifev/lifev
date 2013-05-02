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

#ifndef HDF5_IO_H_
#define HDF5_IO_H_

#include <lifev/core/LifeV.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_config.h>

#include <map>
#include <string>

#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

#include <hdf5.h>

#include <Epetra_MpiComm.h>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{

/*!
  @brief Convenience wrapper for the C interface of the HDF5 library
  @author Radu Popescu <radu.popescu@epfl.ch>

  This class provides an easy way to write and read data from an HDF5 container.
  It is designed to handle a single open file at one time, with multiple open
  datasets simultaneously.
  It does not make any assumptions about data types, or dataset dimensions.
  It is a very thin wrapper on top of the HDF5 library and it is designed for
  high performance parallel operations.

  Usage:
      - open (or create) a file with HDF5IO::openFile
      - create or open existing data table with HDF5IO::createTable or
        HDF5IO::openTable (you can have multiple open tables at one time)
      - read or write blocks of data with HDF5IO::read or HDF5IO::write
      - close the tables after you are finished, with HDF5::closeTable
      - close the file: HDF5::closeFile
*/
class HDF5IO
{
public:
    //! @name Public Types
    //@{
    typedef boost::shared_ptr<Epetra_MpiComm> commPtr_Type;
    //@}

    //! @name Constructors and Destructor
    //@{
    //! Default empty constructor
    HDF5IO() {}

    //! Constructor
    /*!
     * Constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param existing boolean flag indicating whether the file exists already
     *        or not. If it exists, data is appended
     */
    HDF5IO (const std::string& fileName, const commPtr_Type& comm,
            const bool& existing = false);

    //! Empty destructor
    virtual ~HDF5IO() {}
    //@}

    //! @name Public Methods
    //@{
    //! Open
    /*!
     * Create a file or open an existing file
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param existing boolean flag indicating whether the file exists already
     *        or not. If it exists, data is appended
     */
    void openFile (const std::string& fileName, const commPtr_Type& comm,
                   const bool& existing);
    //! Create a new table
    /*!
     * Create a new table in the open file
     * \param tableName a string containing the table name
     * \param fileDataType data type that is to be used in the HDF5 container
     *        should be a standard HDF5 type, not a machine native type;
     *        consult HDF5 documentation for more information
     * \param tableDimensions array of hsize_t of size 2 which holds the
     *        dimensions of the table
     */
    void createTable (const std::string& tableName, hid_t& fileDataType,
                      hsize_t tableDimensions[]);
    //! Open a new table
    /*!
     * Open a new table in the open file
     * \param tableName a string containing the table name
     * \param tableDimensions array of hsize_t of size 2 which will hold the
     *        dimensions of the table (output parameter)
     */
    void openTable (const std::string& tableName, hsize_t tableDimensions[]);
    //! Write
    /*!
     * \param tableName a string containing the table name
     * \param memDataType the type (described as an HDF5 machine native type)
     *        of the data in the buffer that is to be written
     * \param currentCount an array of hsize_t of size two describing the shape
     *        of the block to be written (see HDF5 documentation)
     * \param currentOffset an array of hsize_t of size two describing the
     *        stride of the block to be written (see HDF5 documentation)
     * \param buffer pointer to a memory region containing the data to be
     *        written
     */
    void write (const std::string& tableName,
                hid_t& memDataType, hsize_t currentCount[],
                hsize_t currentOffset[], void* buffer);
    //! Write
    /*!
     * \param tableName a string containing the table name
     * \param memDataType the type (described as an HDF5 machine native type)
     *        of the data in the destination buffer
     * \param currentCount an array of hsize_t of size two describing the shape
     *        of the block to be read (see HDF5 documentation)
     * \param currentOffset an array of hsize_t of size two describing the
     *        stride of the block to be read (see HDF5 documentation)
     * \param buffer pointer to a memory region that represents the destination
     *        of the read operation
     */
    void read (const std::string& tableName,
               hid_t& memDataType, hsize_t currentCount[],
               hsize_t currentOffset[], void* buffer);
    //! Write
    /*!
     * \param tableName a string containing the table name
     */
    void closeTable (const std::string& tableName);
    //! Close an open file
    /*!
     * Call this when finished operating with a file
     */
    void closeFile();
    //@}

private:
    // typedef for internal use
    typedef struct
    {
        hid_t filespace;
        hid_t dataset;
        hid_t plist;
    } tableHandle;

    // Copy constructor and assignment operator are disabled
    HDF5IO (const HDF5IO&);
    HDF5IO& operator= (const HDF5IO&);

    //! Private Data Members
    //@{
    // HDF5 handles
    std::map<std::string, tableHandle> M_tableList;
    hid_t M_fileId;
    //@}
}; // class HDF5IO

} /* namespace LifeV */

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

#endif /* HDF5_IO_H_ */
