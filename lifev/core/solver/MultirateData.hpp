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

#ifndef _MULTIRATEDATA_H_
#define _MULTIRATEDATA_H_ 1

#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/fem/TimeData.hpp>

// LifeV namespace
namespace LifeV
{
/*!
    @class HyperbolicData

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

template < typename MeshType >
class MultirateData: public HyperbolicData
{

	//! Self typedef.
	typedef MultirateData < mesh_Type > multirateData_Type;

    //! @name Constructor & Destructor.
    //@{

    //! Empty Constructor.
    MultirateData ();

    //! Constructor using a data file.
    /*!
      @param dataFile GetPot data file for setup the problem.
      @param section the section for the Hyperbolic data.
    */
   	MultirateData ( const data_Type& dataFile,
                     const std::string& section = "multirate" );

    //! Copy constructor.
    /*!
      @param hyperbolicData object to take a copy.
    */
    MultirateData ( const multirateData_Type &multirateData );

    //! Virtual destructor.
    virtual ~HyperbolicData () {};

    //@}

    //! @name Operators.
    //@{

    //! Assign operator overloading.
    /*!
      @param hyperbolicData The hyperbolicData to be copied.
    */
    hyperbolicData_Type& operator= ( const hyperbolicData_Type& hyperbolicData );

    //@}

    //! @name Methods
    //@{

    //! External setup.
    /*!
      @param dataFile The data file with all the data.
      @param section The global section.
    */
    void setup ( const data_Type& dataFile,
                 const std::string& section = "hyperbolic"  );

    //! Print attributes of the class.
    /*!
      @param output Stream to put the output.
    */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}

    //! @name Get Methods
    //@{

    //! Get the level multirate problem.
    UInt multirateLvel () const
    {
        return M_multirateLevel;
    }

    //@}

protected:

    

    //! @name Multirate
    //@{

    //! Multirate level.
    UInt M_multirateLevel;
    //@}

};

// ===================================================
// Constructors & Destructor
// ===================================================

template < typename MeshType >
HyperbolicData < MeshType >::
HyperbolicData ( ):
        // Miscellaneous
        M_verbose ( 0 ),
        // CFL
        M_relaxCFL ( 0. ),
        M_timeStepFrozen ( 1 ),
        M_computeCFL ( true ),
        // Inflow/outflow
        M_tolInflowOutflow ( 1e-10 )
{}

// Copy constructor
template < typename MeshType >
HyperbolicData < MeshType >::
HyperbolicData( const hyperbolicData_Type &hyperbolicData ):
        // Data containers
        M_data ( hyperbolicData.M_data ),
        M_time ( hyperbolicData.M_time ),
        M_mesh ( hyperbolicData.M_mesh ),
        // Miscellaneous
        M_verbose ( hyperbolicData.M_verbose ),
        M_section ( hyperbolicData.M_section ),
        // CFL
        M_relaxCFL ( hyperbolicData.M_relaxCFL ),
        M_timeStepFrozen ( hyperbolicData.M_timeStepFrozen ),
        M_computeCFL ( hyperbolicData.M_computeCFL ),
        M_tolInflowOutflow ( hyperbolicData.M_tolInflowOutflow )
{}

// ===================================================
// Operators
// ===================================================

// Overloading of the operator =
template < typename MeshType >
HyperbolicData < MeshType >&
HyperbolicData < MeshType >::
operator= ( const hyperbolicData_Type& hyperbolicData )
{
    // Avoid auto-copy
    if ( this != &hyperbolicData )
    {
        // Data containers
        M_data = hyperbolicData.M_data;
        M_time = hyperbolicData.M_time;
        M_mesh = hyperbolicData.M_mesh;
        // Mescellaneous
        M_verbose = hyperbolicData.M_verbose;
        // CFL
        M_relaxCFL = hyperbolicData.M_relaxCFL;
        M_timeStepFrozen = hyperbolicData.M_timeStepFroze;
        M_computeCFL = hyperbolicData.M_computeCFL;
        // Inflow/outflow
        M_tolInflowOutflow = hyperbolicData.M_tolInflowOutflow;
    }

    return *this;
} // operator=

// ===================================================
// Methods
// ===================================================

// External set up method
template < typename MeshType >
void
HyperbolicData < MeshType >::
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

    // If data mesh has not been set
    if ( !M_mesh.get() )
    {
        M_mesh.reset( new meshData_Type ( dataFile, M_section + "/space_discretization" ) );
    }

    // Miscellaneous
    M_verbose = dataFile ( ( M_section + "/miscellaneous/verbose" ).data(), 1 );

    // CFL
    M_relaxCFL = dataFile( ( M_section + "/flux/CFL_relax").data(), 0.9 );
    M_timeStepFrozen = dataFile( ( M_section + "/flux/CFL_frozen").data(), 1 );
    M_computeCFL = dataFile( ( M_section + "/flux/CFL_compute").data(), true );

    // Inflow/outflow
    M_tolInflowOutflow = dataFile( ( M_section + "/inflow/tol").data(), 1e-10 );

	// Multirate Level
    M_multirateLevel = dataFile( ( M_section + "/multirate_setting/multirate_level").data(), 4 );
} // setup

// Print attiributes of the class
template < typename MeshType >
void
HyperbolicData < MeshType >::
showMe ( std::ostream& output ) const
{
    output << "Class HyperbolicData:" << std::endl;
    M_time->showMe( output );
    M_mesh->showMe( output );
    output << "Verbosity level                      " << this->M_verbose << std::endl
           << "Section of GetPot                    " << this->M_section << std::endl
           << "Compute CFL condition?               " << this->M_computeCFL << std::endl
           << "Relax CFL parameter                  " << this->M_relaxCFL << std::endl
           << "Number of step with time step frozen " << this->M_timeStepFrozen << std::endl
           << "Inflow/outflow tolerance             " << this->M_tolInflowOutflow << std::endl
	<< "Multirate Level			" << this->M_multirateLevel << std::endl;
	
} // showMe

} // Namespace LifeV

#endif /* _HYPERBOLICDATA_H_ */
