/* -*- mode: c++ -*-
   This program is part of the LifeV library

   Autor(s): Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
       Date:

   Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2006-2010 EPFL, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/**
  @file dataDarcy.hpp
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @date 05/2010

  @brief This file contains the data for all the Darcy solver
*/
#ifndef _DATAHYPERBOLIC_H_
#define _DATAHYPERBOLIC_H_

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifefem/dataTime.hpp>

// LifeV namespace
namespace LifeV
{
/*!
  @class DataDarcy

  This class contain the basic data for the Darcy solver. In particoular it stores...
  @todo class not finished!
 */
template <typename Mesh>
class DataHyperbolic
{
public:

    // Policies
    //! @name Policies
    //@{

    typedef GetPot                          Data_Type;
    typedef boost::shared_ptr< Data_Type >  Data_ptrType;

    typedef DataTime                        Time_Type;
    typedef boost::shared_ptr< Time_Type >  Time_ptrType;

    typedef DataMesh                        Mesh_Type;
    typedef boost::shared_ptr< Mesh_Type >  Mesh_ptrType;

    //@}

    // Constructors.
	//! @name Constructors
	//@{

    //! Empty Constructor
    DataHyperbolic();

    /*!
    Constructor using a data file.
      @param dataFile GetPot data file for setup the problem
      @param section the section for the Darcy data
    */
    DataHyperbolic( const GetPot& dataFile, const std::string& section = "hyperbolic" );

    /*!
    Copy constructor.
      @param dataDarcy object to take a copy
    */
    DataHyperbolic( const DataHyperbolic &dataHyperbolic );

    //@}

    // Set methods
    //! @name Set methods
    //@{

    /*! Set data time container
        @param DataTime Boost shared_ptr to dataTime container
    */
    inline void setDataTime( const Time_ptrType DataTime )
    {
        M_time = DataTime;
    }

    /*! Set mesh container
        @param DataMesh Boost shared_ptr to dataMesh container
    */
    inline void setDataMesh( const Mesh_ptrType DataMesh )
    {
        M_mesh = DataMesh;
    }

    // Get methods.
    //! @name Get methods
    //@{

	//! Get the level of verbosity of the problem.
    inline const UInt verbose( void ) const
    {
        return M_verbose;
    }

    //! Get the main section of the data file.
    inline const std::string section( void ) const
    {
        return M_section;
    }

   	//! Get the data file of the problem.
   	inline Data_ptrType dataFile( void ) const
    {
        return M_data;
    }

    /*! Get data time container.
        @return shared_ptr to dataTime container
    */
    inline Time_ptrType dataTime( void ) const
    {
        return M_time;
    }

    /*! Get mesh container
       @return shared_ptr to dataMesh container
    */
    inline Mesh_ptrType dataMesh( void ) const
    {
        return M_mesh;
    }

    /*! Get the relaxation paramether for the CFL condition
      @return CFL relaxation paramether
    */
    inline const Real getCFLrelax( void ) const
    {
        return M_relaxCFL;
    }

    //@}

    // Methods.
	//! @name Methods
	//@{

    /*! Overloading of the operator =
        @param dataDarcy The DataDarcy to be copied.
    */
    DataHyperbolic& operator=( const DataHyperbolic& dataHyperbolic );

    /*! External setup
        @param dataFile The data file with all the data.
        @param section The global section.
    */
    void setup( const Data_Type& dataFile, const std::string& section = "hyperbolic"  );

	//@}



protected:

    //! Data containers for time and mesh
    Data_ptrType      M_data;
    Time_ptrType      M_time;
    Mesh_ptrType      M_mesh;

    //! Miscellaneous
    UInt              M_verbose;
    std::string       M_section;

    //! Relax paramether for the CFL condition
    Real              M_relaxCFL;

};

// ===================================================
// Constructors
// ===================================================

template < typename Mesh >
DataHyperbolic<Mesh>::DataHyperbolic( ):
    // Data containers
    M_data          ( ),
    M_time          ( ),
    M_mesh          ( ),
    // Miscellaneous
    M_verbose       ( static_cast<UInt>(0) ),
    M_section       ( ),
    // CFL
    M_relaxCFL      ( static_cast<Real>(0.) )
{
    CONSTRUCTOR( "DataHyperbolic" );
}

// Copy constructor
template < typename Mesh >
DataHyperbolic<Mesh>::DataHyperbolic( const DataHyperbolic &dataHyperbolic ):
    // Data containers
    M_data        ( dataHyperbolic.M_data ),
	M_time        ( dataHyperbolic.M_time ),
    M_mesh        ( dataHyperbolic.M_mesh ),
    // Miscellaneous
    M_verbose     ( dataHyperbolic.M_verbose ),
    M_section     ( dataHyperbolic.M_section ),
    // CFL
    M_relaxCFL    ( dataHyperbolic.M_relaxCFL )
{
    CONSTRUCTOR( "DataHyperbolic" );
}

// Overloading of the operator =
template < typename Mesh >
DataHyperbolic<Mesh>&
DataHyperbolic<Mesh>::operator=( const DataHyperbolic& dataHyperbolic )
{
    // Avoid auto-copy
    if ( this != &dataHyperbolic )
    {
        // Data containers
        M_data        = dataHyperbolic.M_data;
        M_time        = dataHyperbolic.M_time;
        M_mesh        = dataHyperbolic.M_mesh;
        // Mescellaneous
        M_verbose     = dataHyperbolic.M_verbose;
        // CFL
        M_relaxCFL    = dataHyperbolic.M_relaxCFL;
    }

    return *this;

}


// External set up method
template < typename Mesh >
void DataHyperbolic<Mesh>::setup( const Data_Type& dataFile, const std::string& section )
{
    M_section = section;

    // If data has not been set
    if( !M_data.get() )
        M_data.reset( new Data_Type( dataFile ) );

    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, M_section + "/time_discretization" ) );

    // If data mesh has not been set
    if ( !M_mesh.get() )
        M_mesh.reset( new Mesh_Type( dataFile, M_section + "/space_discretization" ) );

    // Miscellaneous
    M_verbose = dataFile( ( M_section + "/miscellaneous/verbose" ).data(), 1 );

    // CFL
    M_relaxCFL = dataFile( (M_section + "/numerical_flux/CFL/relax").data(), 0.9 );

}

}
#endif
