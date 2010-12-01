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
#ifndef _DATADARCY_H_
#define _DATADARCY_H_

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifefem/assemb.hpp>

// LifeV namespace
namespace LifeV
{
/*!
  @class DataDarcy

  This class contain the basic data for the Darcy solver. In particoular it stores...
  @todo class not finished!
 */
template <typename Mesh>
class DataDarcy
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
    DataDarcy();

    /*!
    Constructor using a data file.
      @param dataFile GetPot data file for setup the problem
      @param section the section for the Darcy data
    */
    DataDarcy( const GetPot& dataFile, const std::string& section = "darcy" );

    /*!
    Copy constructor.
      @param dataDarcy object to take a copy
    */
    DataDarcy( const DataDarcy &dataDarcy );

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


    //@}

    // Methods.
    //! @name Methods
    //@{

    /*! Overloading of the operator =
        @param dataDarcy The DataDarcy to be copied.
    */
    DataDarcy& operator=( const DataDarcy& dataDarcy );

    /*! External setup
        @param dataFile The data file with all the data.
        @param section The global section.
    */
    void setup( const Data_Type& dataFile, const std::string& section = "darcy"  );

    //@}



protected:

    //! Data containers for time and mesh
    Data_ptrType      M_data;
    Time_ptrType      M_time;
    Mesh_ptrType      M_mesh;

    //! Miscellaneous
    UInt              M_verbose;
    std::string       M_section;

};

// ===================================================
// Constructors
// ===================================================

template < typename Mesh >
DataDarcy<Mesh>::DataDarcy( ):
        // Data containers
        M_data          ( ),
        M_time          ( ),
        M_mesh          ( ),
        // Miscellaneous
        M_verbose       ( static_cast<UInt>(0) ),
        M_section       ( )
{
    CONSTRUCTOR( "DataDarcy" );
}

// Copy constructor
template < typename Mesh >
DataDarcy<Mesh>::DataDarcy( const DataDarcy &dataDarcy ):
        // Data containers
        M_data                ( dataDarcy.M_data ),
        M_time                ( dataDarcy.M_time ),
        M_mesh                ( dataDarcy.M_mesh ),
        // Miscellaneous
        M_verbose             ( dataDarcy.M_verbose ),
        M_section             ( dataDarcy.M_section )
{
    CONSTRUCTOR( "DataDarcy" );
}

// Overloading of the operator =
template < typename Mesh >
DataDarcy<Mesh>&
DataDarcy<Mesh>::operator=( const DataDarcy& dataDarcy )
{
    // Avoid auto-copy
    if ( this != &dataDarcy )
    {
        // Data containers
        M_data           = dataDarcy.M_data;
        M_time           = dataDarcy.M_time;
        M_mesh           = dataDarcy.M_mesh;
        // Mescellaneous
        M_verbose       = dataDarcy.M_verbose;
    }

    return *this;

}


// External set up method
template < typename Mesh >
void DataDarcy<Mesh>::setup( const Data_Type& dataFile, const std::string& section )
{
    M_section = section;

    // If data has not been set
    if ( !M_data.get() )
        M_data.reset( new Data_Type( dataFile ) );

    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, M_section + "/time_discretization" ) );

    // If data mesh has not been set
    if ( !M_mesh.get() )
        M_mesh.reset( new Mesh_Type( dataFile, M_section + "/space_discretization" ) );

    // Miscellaneous
    M_verbose      = dataFile( ( M_section + "/miscellaneous/verbose" ).data(), 1 );
}


/////////////////////////////////////////////////////////////////////////////////

template< typename Mesh,
typename SolverType = LifeV::SolverTrilinos >
class inversePermeability
{

public:

    // Policies.
    //! @name Policies
    //@{

    typedef boost::function<Matrix ( const Real&, const Real&,
                                     const Real&, const Real&,
                                     const std::vector<Real> & )>
    permeability_type;

    typedef typename SolverType::vector_type      vector_type;
    typedef boost::shared_ptr<vector_type>        vector_ptrtype;

    //@}

    // Constructors and destructor.
    //! @name Constructors and destructor
    //@{

    // Copy constructor
    inversePermeability ( const permeability_type& invPerm, FESpace<Mesh, EpetraMap>& fESpace ):
            M_inversePermeability ( invPerm ),
            M_fESpace             ( fESpace ),
            M_fields              ( std::vector< const vector_ptrtype* >(0) ) {};

    //! Virtual destructor.
    // virtual ~inversePermeability ();

    //@}

    // Set methods
    //! @name Set methods
    //@{

    // Add one field
    inline void setField ( const vector_ptrtype & field ) { M_fields.push_back( &field ); };

    inline void setFunction ( const permeability_type & invPerm ) { M_inversePermeability = invPerm; };

    //@}

    // Get methods
    //! @name Get methods
    //@{

    Matrix operator() ( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& iElem );

    //@}

private:

    // Vector of pointers for the dependences of the permeability to an external field.
    std::vector< const vector_ptrtype* > M_fields;

    // Inverse permeability function
    permeability_type             M_inversePermeability;

    // Finite element space
    FESpace<Mesh, EpetraMap>&     M_fESpace;

};

template < typename Mesh, typename SolverType >
Matrix
inversePermeability < Mesh, SolverType >::
operator() ( const Real& t, const Real& x, const Real& y, const Real& z, const UInt& iElem )
{
    std::vector<Real> values ( M_fields.size(), 0 );
    ElemVec value ( M_fESpace.refFE().nbDof(), 1 );

    // Update the value of the current element
    M_fESpace.fe().update( M_fESpace.mesh()->element( iElem ),
                           UPDATE_QUAD_NODES | UPDATE_WDET );

    for ( UInt i(static_cast<UInt>(0)); i < values.size(); ++i )
    {
        extract_vec ( *( *(M_fields)[i] ),
                      value,
                      M_fESpace.refFE(),
                      M_fESpace.dof(),
                      M_fESpace.fe().currentLocalId(), 0 );

        values[i] = value[0];

    }

    return M_inversePermeability ( t, x, y, z, values );
}

}
#endif
