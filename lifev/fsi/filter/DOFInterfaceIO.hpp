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
  @brief Class that handles I/O of fluid-solid DOF interfaces in FSI

  @date 2013-04-30
  @author Radu Popescu <radu.popescu@epfl.ch>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef DOF_INTERFACE_IO_
#define DOF_INTERFACE_IO_

#include <string>
#include <lifev/core/Core_config.h>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_config.h>

#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

#include <Epetra_MpiComm.h>

#include <boost/shared_ptr.hpp>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/HDF5IO.hpp>

#include <lifev/core/fem/DOFInterface3Dto3D.hpp>

namespace LifeV
{

/*!
  @brief Class that handles I/O of fluid-solid DOF interfaces in FSI
  @author Radu Popescu radu.popescu@epfl.ch

  TODO: write description (also storage format)
*/
class DOFInterfaceIO
{
public:
    //! @name Public Types
    //@{
    typedef boost::shared_ptr<Epetra_MpiComm> commPtr_Type;
    typedef DOFInterface3Dto3D interface_Type;
    typedef boost::shared_ptr<interface_Type> interfacePtr_Type;
    typedef std::vector<interfacePtr_Type> interfaceVector_Type;
    typedef boost::shared_ptr<interfaceVector_Type> interfaceVectorPtr_Type;
    typedef std::map<UInt, UInt>                     dofMap_Type;
    typedef boost::shared_ptr<dofMap_Type>           dofMapPtr_Type;

    //@}

    //! \name Constructors & Destructors
    //@{
    //! Constructor
    /*!
     * Non-empty constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     */
    DOFInterfaceIO (const std::string& fileName,
                    const commPtr_Type& comm);

    //! Empty destructor
    virtual ~DOFInterfaceIO() {}
    //@}

    //! \name Public Methods
    //@{
    //! Write method
    /*!
     * Call this method to write the DOF interfaces
     */
    void write (const interfaceVectorPtr_Type& interfaces);
    //! Read method
    /*!
     * Call this method to read from the HDF5 file the DOF interface part
     * with the rank of the current MPI process
     */
    void read (dofMapPtr_Type& interface);
    //@}

private:
    // Copy constructor and assignment operator are disabled
    DOFInterfaceIO (const DOFInterfaceIO&);
    DOFInterfaceIO& operator= (const DOFInterfaceIO&);

    //! Private Data Members
    //@{
    commPtr_Type M_comm;
    UInt M_myRank;
    std::string M_fileName;
    HDF5IO M_HDF5IO;
    // Buffer for reading/writing
    std::vector<UInt> M_keyBuffer;
    std::vector<UInt> M_valueBuffer;
    //@}
}; // class DOFInterfaceIO

} /* namespace LifeV */

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */
#endif /* DOF_INTERFACE_IO_ */
