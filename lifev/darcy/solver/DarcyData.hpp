//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the data for all the Darcy solver

     @date 05/2010
     @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

     @contributor M. Kern <michel.kern@inria.fr>
     @maintainer M. Kern <michel.kern@inria.fr>
 */

#ifndef _DATADARCY_H_
#define _DATADARCY_H_ 1

#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/TimeAdvanceData.hpp>

// LifeV namespace
namespace LifeV
{

//! @class DarcyData This class contain the basic data for the Darcy solver.
/*!
  In particolar it stores the data as GetPot object, the data for the mesh and the data for the time schemes.
*/
template < typename MeshType >
class DarcyData
{
public:

    // Policies.
    //! @name Public policies
    //@{

    //! Typedef for the mesh template.
    typedef MeshType mesh_Type;

    //! Typedef for the GetPot data.
    typedef GetPot data_Type;

    //! Shared pointer for the data.
    typedef boost::shared_ptr < data_Type > dataPtr_Type;

    //! Typedef for the time data.
    typedef TimeData timeData_Type;

    //! Shared pointer for the time data.
    typedef boost::shared_ptr < timeData_Type > timeDataPtr_Type;

    //! Typedef for the time advance data.
    typedef TimeAdvanceData timeAdvanceData_Type;

    //! Shared pointer for the time advance data.
    typedef boost::shared_ptr < timeAdvanceData_Type > timeAdvanceDataPtr_Type;

    //! Typedef for the mesh data.
    typedef MeshData meshData_Type;

    //! Shared pointer for the mesh data.
    typedef boost::shared_ptr < meshData_Type > meshDataPtr_Type;

    //! Self typedef.
    typedef DarcyData < mesh_Type > darcyData_Type;

    //@}

    // Constructors.
    //! @name Constructors
    //@{

    //! Empty Constructor.
    DarcyData ();

    //! Constructor using a data file.
    /*!
      @param dataFile GetPot data file for setup the problem.
      @param section the section for the Darcy data.
    */
    DarcyData ( const data_Type& dataFile,
                const std::string& section = "darcy" );

    //! Copy constructor.
    /*!
      @param darcyData object to take a copy.
    */
    DarcyData ( const darcyData_Type &darcyData );

    //@}

    // Methods.
    //! @name Methods
    //@{

    //! Overloading of the operator =.
    /*!
       @param darcyData The DarcyData to be copied.
    */
    darcyData_Type& operator= ( const darcyData_Type& darcyData );

    //! External setup.
    /*!
      @param dataFile The data file with all the data.
      @param section The global section.
    */
    void setup ( const data_Type& dataFile,
                 const std::string& section = "darcy"  );

    //@}

    // Set methods
    //! @name Set methods
    //@{

    //! Set data time container.
    /*!
      @param timeData Boost shared_ptr to timeData container
    */
    void setTimeData ( const timeDataPtr_Type& timeData )
    {
        M_time = timeData;
    }

    //! Set data time container.
    /*!
      @param timeAdvanceData Boost shared_ptr to timeAdvanceData container
    */
    void setTimeAdvanceData ( const timeDataPtr_Type& timeAdvanceData )
    {
        M_timeAdvance = timeAdvanceData;
    }

    //! Set mesh container.
    /*!
      @param meshData Boost shared_ptr to meshData container
    */
    void setMeshData ( const meshDataPtr_Type& meshData )
    {
        M_mesh = meshData;
    }

    // Get methods.
    //! @name Get methods
    //@{

    //! Get the level of verbosity of the problem.
    UInt verbose () const
    {
        return M_verbose;
    }

    //! Get the main section of the data file.
    std::string section () const
    {
        return M_section;
    }

    //! Get the data file of the problem.
    /*!
      @return shared_ptr to data container.
    */
    const dataPtr_Type& dataFilePtr () const
    {
        return M_data;
    }

    //! Get the data file of the problem.
    /*!
      @return shared_ptr to data container.
    */
    dataPtr_Type& dataFilePtr ()
    {
        return M_data;
    }

    //! Get data time container.
    /*!
      @return shared_ptr to TimeData container.
    */
    const timeDataPtr_Type& dataTimePtr () const
    {
        return M_time;
    }

    //! Get data time container.
    /*!
       @return shared_ptr to TimeData container.
    */
    timeDataPtr_Type& dataTimePtr ()
    {
        return M_time;
    }

    //! Get data time advance container.
    /*!
      @return shared_ptr to TimeAdvanceData container.
    */
    const timeAdvanceDataPtr_Type& dataTimeAdvancePtr () const
    {
        return M_timeAdvance;
    }

    //! Get data time advance container.
    /*!
       @return shared_ptr to TimeAdvanceData container.
    */
    timeAdvanceDataPtr_Type& dataTimeAdvancePtr ()
    {
        return M_timeAdvance;
    }

    //! Get mesh container
    /*!
      @return shared_ptr to meshData container.
    */
    const meshDataPtr_Type& meshDataPtr () const
    {
        return M_mesh;
    }

    //! Get mesh container
    /*!
      @return shared_ptr to meshData container.
    */
    meshDataPtr_Type& meshDataPtr ()
    {
        return M_mesh;
    }

    //@}


private:

    //! Data containers for time and mesh
    dataPtr_Type            M_data;
    timeDataPtr_Type        M_time;
    timeAdvanceDataPtr_Type M_timeAdvance;
    meshDataPtr_Type        M_mesh;

    //! Miscellaneous
    UInt                    M_verbose;
    std::string             M_section;

};

// ===================================================
// Constructors
// ===================================================

template < typename MeshType >
DarcyData < MeshType >::
DarcyData ():
        // Miscellaneous
        M_verbose       ( static_cast<UInt>(0) )
{}

// Copy constructor
template < typename MeshType >
DarcyData < MeshType >::
DarcyData ( const darcyData_Type &darcyData ):
        // Data containers
        M_data ( darcyData.M_data ),
        M_time ( darcyData.M_time ),
        M_timeAdvance ( darcyData.M_timeAdvance ),
        M_mesh ( darcyData.M_mesh ),
        // Miscellaneous
        M_verbose ( darcyData.M_verbose ),
        M_section ( darcyData.M_section )
{}

// Overloading of the operator =
template < typename MeshType >
DarcyData < MeshType >&
DarcyData < MeshType >::
operator= ( const darcyData_Type& darcyData )
{
    // Avoid auto-copy
    if ( this != &darcyData )
    {
        // Data containers
        M_data = darcyData.M_data;
        M_time = darcyData.M_time;
        M_timeAdvance = darcyData.M_timeAdvance;
        M_mesh = darcyData.M_mesh;
        // Mescellaneous
        M_verbose = darcyData.M_verbose;
    }

    return *this;
} // operator=


// External set up method
template < typename MeshType >
void
DarcyData < MeshType >::
setup ( const data_Type& dataFile, const std::string& section )
{
    M_section = section;

    // If data has not been set
    if ( !M_data.get() )
    {
        M_data.reset( new data_Type ( dataFile ) );
    }

    // If data time has not been set
    if ( !M_time.get() )
    {
        M_time.reset( new timeData_Type ( dataFile, M_section + "/time_discretization" ) );
    }

    // If data time has not been set
    if ( !M_timeAdvance.get() )
    {
        M_timeAdvance.reset( new timeAdvanceData_Type ( dataFile, M_section + "/time_discretization" ) );
    }

    // If data mesh has not been set
    if ( !M_mesh.get() )
    {
        M_mesh.reset( new meshData_Type ( dataFile, M_section + "/space_discretization" ) );
    }

    // Miscellaneous
    M_verbose = dataFile( ( M_section + "/miscellaneous/verbose" ).data(), 1 );
} // setup

} // Namespace LifeV

#endif // _DATADARCY_H_

// -*- mode: c++ -*-
