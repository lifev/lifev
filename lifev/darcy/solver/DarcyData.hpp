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

#ifndef _DATADARCY_HPP_
#define _DATADARCY_HPP_ 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/fem/TimeAdvanceData.hpp>

// LifeV namespace
namespace LifeV
{

//! @class DarcyData This class contain the basic data for the Darcy solver.
/*!
  In particular it stores the data as GetPot object, the data for the mesh and the data for the time schemes.
*/
template < typename MeshType >
class DarcyData
{
public:

    // Policies.
    //! @name Public Types
    //@{

    //! Typedef for the mesh template.
    typedef MeshType mesh_Type;

    //! Typedef for the GetPot data.
    typedef GetPot data_Type;

    //! Shared pointer for the data.
    typedef boost::shared_ptr < data_Type > dataPtr_Type;

    //! Teuchos parameter list.
    typedef Teuchos::ParameterList paramList_Type;

    //! Shared pointer for the Teuchos parameter list.
    typedef Teuchos::RCP< paramList_Type > paramListPtr_Type;

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
    DarcyData ():
        M_verbose(0)
    {}

    //! Constructor using a data file.
    /*!
      @param dataFile GetPot data file for set-up the problem.
      @param section the section for the Darcy data.
    */
    DarcyData ( const data_Type& dataFile,
                const std::string& section = "darcy" );

    // Methods.
    //! @name Methods
    //@{

    //! External set-up.
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

    //! Set Teuchos parameter list for linear algebra.
    /*!
      @param paramList Teuchos RCP with the parameter list.
      @param linearSolver Section of the linear solver in the paramList.
      @param precond Section of the preconditioner in the paramList.
    */
    void setLinearAlgebraList ( const paramListPtr_Type& linearAlgebraList,
                                const std::string& linearSolver = "Linear Solver",
                                const std::string& precond = "Preconditioner" )
    {
        M_linearAlgebraList = linearAlgebraList;
        M_linearSolverSection = linearSolver;
        M_precondSection = precond;
    } // setLinearAlgebraList

    //! Set data time container.
    /*!
      @param timeData Boost shared_ptr to timeData container
    */
    void setTimeData ( const timeDataPtr_Type& timeData )
    {
        M_time = timeData;
    } // setTimeData

    //! Set mesh container.
    /*!
      @param meshData Boost shared_ptr to meshData container
    */
    void setMeshData ( const meshDataPtr_Type& meshData )
    {
        M_mesh = meshData;
    } // setMeshData

    // Get methods.
    //! @name Get methods
    //@{

    //! Get the level of verbosity of the problem.
    UInt verbose () const
    {
        return M_verbose;
    } // verbose

    //! Get the main section of the data file.
    std::string section () const
    {
        return M_section;
    } // section

    //! Get the data file of the problem.
    /*!
      @return shared_ptr to data container.
    */
    const dataPtr_Type& dataFilePtr () const
    {
        return M_data;
    } // dataFilePtr

    //! Get the data file of the problem.
    /*!
      @return shared_ptr to data container.
    */
    dataPtr_Type& dataFilePtr ()
    {
        return M_data;
    } // dataFilePtr

    //! Get data time container.
    /*!
      @return shared_ptr to TimeData container.
    */
    const timeDataPtr_Type& dataTimePtr () const
    {
        return M_time;
    } // dataTimePtr

    //! Get data time container.
    /*!
       @return shared_ptr to TimeData container.
    */
    timeDataPtr_Type& dataTimePtr ()
    {
        return M_time;
    } // dataTimePtr

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
    } // meshDataPtr

    //! Get mesh container
    /*!
      @return shared_ptr to meshData container.
    */
    meshDataPtr_Type& meshDataPtr ()
    {
        return M_mesh;
    } // meshDataPtr

    //! Get Teuchos parameter list for the linear solver.
    /*!
      @return Teuchos RCP with the parameter list for the linear solver.
    */
    const paramList_Type& linearSolverList () const
    {
        ASSERT ( M_linearAlgebraList.get(), "Parameter list not set." );
        return M_linearAlgebraList->sublist( M_linearSolverSection );
    } // linearSolverList

    //! Get Teuchos parameter list for the preconditioner.
    /*!
      @return Teuchos RCP with the parameter list for the preconditioner.
    */
    const paramList_Type& preconditionerList () const
    {
        ASSERT ( M_linearAlgebraList.get(), "Parameter list not set." );
        return M_linearAlgebraList->sublist( M_precondSection );
    } // preconditionerList

    //@}

private:

    //! @name Private Constructors
    //@{

    //! Inhibited copy constructor.
    DarcyData ( const darcyData_Type & );

    //@}

    //! @name Private Operators
    //@{

    //! Inhibited assign operator.
    darcyData_Type& operator= ( const darcyData_Type& );

    //@}

    //! Inhibited assign operator.

    //! Data GetPot.
    dataPtr_Type M_data;

    //! Section in GetPot file.
    std::string M_section;

    //! Data container for time.
    timeDataPtr_Type M_time;

    //! Data container for time advance.
    timeAdvanceDataPtr_Type M_timeAdvance;

    //! Data container for mesh.
    meshDataPtr_Type M_mesh;

    //! Teuchos paramter list for linear algebra.
    paramListPtr_Type M_linearAlgebraList;

    //! Section in the parameter list for linear solver.
    std::string M_linearSolverSection;

    //! Section in the parameter list for preconditioner.
    std::string M_precondSection;

    //! Output verbose.
    UInt M_verbose;

};

// ===================================================
// Constructors
// ===================================================

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

#endif // _DATADARCY_HPP_

// -*- mode: c++ -*-
