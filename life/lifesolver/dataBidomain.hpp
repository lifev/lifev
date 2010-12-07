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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief File containing a class for handling Bidomain data with GetPot

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef _DATABIDOMAIN_H_
#define _DATABIDOMAIN_H_

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifesolver/heartFunctors.hpp>


namespace LifeV
{
/*!
  \class DataBidomain

  Base class which holds usual data for the Bidomain model solvers

*/

class DataBidomain:
    public DataMesh,
    public DataTime
{
public:


    //! @name Public Types
    //@{

    /*! @enum listOfTemplatesOptions
        Description of the purpose of the enumerator list.
    */

    typedef boost::function<Real ( Real const& x, Real const& y, Real const& z, Real const& t, ID const& id, Real const&)> region_Type;

    //@}



    //! @name Constructor & Destructor
    //@{

    //! Constructors
    DataBidomain();

    DataBidomain( boost::shared_ptr<HeartFunctors> heart);

    DataBidomain( const DataBidomain& dataBidomain );

    //@}


    //! @name Operators
    //@{

    DataBidomain& operator=( const DataBidomain& dataBidomain );

    //@}

    //! @name Methods
    //@{

    //! Output: show the data used for the simulation
    void showMe( std::ostream& output = std::cout );

    //@}



    //! @name Set Methods
    //@{

    //! external setup: set all the data for the simulation
    void setup( const GetPot& dataFile );

    //@}


    //! @name Get Methods
    //@{

    //! verbose
    Real        verbose()        const {return M_verbose;};

    //! FE space order
    std::string uOrder()         const {return M_uOrder;};

    //! Chi
    const Real&        Chi()            const {return M_volumeSurfaceRatio;}
    //! fiber File
    std::string fibers_file()    const {return M_fibersFile;}

    const Int&         heart_diff_fct() const {return M_heartDiffusionFactor;}

    const bool&        has_fibers()     const {return M_hasFibers;}

    //! format vct
    const bool&        fibers_format()  const {return M_fibersFormat;}

    //! sigma_l
    const Real&        sigmal_i()       const 	{return M_longitudinalInternalConductivity;}
    const Real&        sigmal_e()       const 	{return M_longitudinalExternalConductivity;}

    //! sigma_t
    const Real&        sigmat_i()       const 	{return M_transversalInternalConductivity;}
    const Real&        sigmat_e()       const 	{return M_transversalExternalConductivity;}

    //! Cm
    const Real&        Cm()             const 	{return M_membraneCapacitance;}
    //! D
    const Real&        D_i()            const 	{return M_internalDiffusivity;}
    //! Post_dir
    const Real&        D_e()            const 	{return M_externalDiffusivity;}
    //! Post_dir
    std::string Post_dir()       const {return M_postProcessingDirectory;}

    UInt        order_bdf()      const {return M_BDForder;}

    //@}


protected:
private:
    // format of fibers file
    bool        M_fibersFormat;
    bool        M_hasFibers;

    // order of time discretization BDF order
    UInt        M_BDForder;
    UInt        M_verbose;

    Int         M_heartDiffusionFactor;

    Real        M_externalDiffusivity;
    Real        M_internalDiffusivity;
    Real        M_longitudinalExternalConductivity;
    Real        M_longitudinalInternalConductivity;
    Real        M_membraneCapacitance;
    Real        M_transversalExternalConductivity;
    Real        M_transversalInternalConductivity;
    Real        M_volumeSurfaceRatio;

    std::string M_fibersDirectory;
    std::string M_fibersFile;
    std::string M_postProcessingDirectory; //! full name
    std::string M_uOrder;

    region_Type M_reducedConductivityBox;
    region_Type M_reducedConductivityCylinder;
    region_Type M_reducedConductivitySphere;






};

}
#endif
