//@HEADER
/*
*******************************************************************************
    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University
    This file is part of LifeV.
    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with LifeV. If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************
*/
//@HEADER
/*!
 * @file
 * @brief Data for hyperbolic scalar equations.
 *
 *
 * @date 30-09-2010
 *
 * @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
 * @author Michel Kern       <michel.kern@inria.fr>
 *
 * @contributor
 *
 * @mantainer Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
 *
 */

#ifndef _DATAHYPERBOLIC_H_
#define _DATAHYPERBOLIC_H_ 1

#include <life/lifearray/tab.hpp>

#include <life/lifemesh/dataMesh.hpp>

#include <life/lifefem/dataTime.hpp>

// LifeV namespace
namespace LifeV
{
/*!
  @class DataHyperbolic

  @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
  @author Michel Kern       <michel.kern@inria.fr>

  This class contain the basic data for the hyperbolic solver. In particoular it stores
  <ol>
  <li> the GetPot data file; </li>
  <li> the time class handler; </li>
  <li> the mesh class handler; </li>
  <li> the level of verbosity of the solver hyperbolic; </li>
  <li> the section for the GetPot data file; </li>
  <li> the relax parameter for the CFL condition. </li>
  </ol>
  @todo class not finished!
 */
template <typename Mesh>
class DataHyperbolic
{
public:

    //! @name Public Types
    //@{

    typedef GetPot                          Data_Type;
    typedef boost::shared_ptr< Data_Type >  DataPtr_Type;

    typedef DataTime                        Time_Type;
    typedef boost::shared_ptr< Time_Type >  TimePtr_Type;

    typedef DataMesh                        Mesh_Type;
    typedef boost::shared_ptr< Mesh_Type >  MeshPtr_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    DataHyperbolic();

    //! Constructor using a data file.
    /*!
      @param dataFile GetPot data file for setup the problem
      @param section the section for the Darcy data
    */
    DataHyperbolic( const GetPot& dataFile, const std::string& section = "hyperbolic" );

    //! Copy constructor.
    /*!
      @param dataDarcy object to take a copy
    */
    DataHyperbolic( const DataHyperbolic &dataHyperbolic );

    //! Virtual destructor
    virtual ~DataHyperbolic();

    //@}

    //! @name Operators
    //@{

    //! Assign operator overloading
    /*!
      @param dataDarcy The DataDarcy to be copied
    */
    DataHyperbolic& operator=( const DataHyperbolic& dataHyperbolic );

    //@}

    //! @name Methods
    //@{

    //! External setup
    /*!
      @param dataFile The data file with all the data.
      @param section The global section.
    */
    void setup( const Data_Type& dataFile, const std::string& section = "hyperbolic"  );

    //! Print attributes of the class
    /*!
      @param output Stream to put the output
    */
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{

    //! Set data time container
    /*!
      @param DataTime Boost shared_ptr to dataTime container
    */
    inline void setDataTime( const TimePtr_Type DataTime )
    {
        M_time = DataTime;
    }

    //! Set mesh container
    /*!
      @param DataMesh Boost shared_ptr to dataMesh container
    */
    inline void setDataMesh( const MeshPtr_Type DataMesh )
    {
        M_mesh = DataMesh;
    }

    //! @name Get Methods
    //@{

    //! Get the level of verbosity of the problem.
    inline UInt verbose () const
    {
        return M_verbose;
    }

    //! Get the main section of the data file.
    inline const std::string section () const
    {
        return M_section;
    }

    //! Get the data file of the problem.
    inline DataPtr_Type dataFile () const
    {
        return M_data;
    }

    //! Get data time container.
    /*!
      @return shared_ptr to dataTime container
    */
    inline TimePtr_Type dataTime () const
    {
        return M_time;
    }

    //! Get mesh container
    /*!
      @return shared_ptr to dataMesh container
    */
    inline MeshPtr_Type dataMesh () const
    {
        return M_mesh;
    }

    //! Get the relaxation parameter for the CFL condition
    /* Get the parameter to compute \f$ CFL_{\rm used} = \alpha CFL_{\rm computed} \f$
       @return CFL relaxation parameter
    */
    inline Real getCFLRelaxParameter () const
    {
        return M_relaxCFL;
    }

    //@}

    inline Real __attribute__ ((__deprecated__)) getCFLrelax () const
    {
        return getCFLRelaxParameter ();
    }

protected:

    //! Data containers for time and mesh
    DataPtr_Type      M_data;
    TimePtr_Type      M_time;
    MeshPtr_Type      M_mesh;

    //! Miscellaneous
    UInt              M_verbose;
    std::string       M_section;

    //! Relax parameter for the CFL condition
    Real              M_relaxCFL;

};

// ===================================================
// Constructors & Destructor
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

}

// Virtual destructor
template < typename Mesh >
DataHyperbolic<Mesh>::~DataHyperbolic()
{

}

// ===================================================
// Operators
// ===================================================

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

// ===================================================
// Methods
// ===================================================

// External set up method
template < typename Mesh >
void
DataHyperbolic<Mesh>::setup( const Data_Type& dataFile,
                             const std::string& section )
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
    M_verbose = dataFile( ( M_section + "/miscellaneous/verbose" ).data(), 1 );

    // CFL
    M_relaxCFL = dataFile( (M_section + "/numerical_flux/CFL/relax").data(), 0.9 );

}

// Print attiributes of the class
template < typename Mesh >
void
DataHyperbolic<Mesh>::showMe( std::ostream& output ) const
{
    output << "Class DataHyperbolic:" << std::endl;
    M_time->showMe( output );
    M_mesh->showMe( output );
    output << "Verbosity level     " << M_verbose << std::endl
           << "Section of GetPot   " << M_section << std::endl
           << "Relax CFL parameter " << M_relaxCFL << std::endl
           << std::flush;
}

} // Namespace LifeV

#endif /* _DATAHYPERBOLIC_H_ */
