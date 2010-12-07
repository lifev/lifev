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
    @brief File containing a class for handling Monodomain data with GetPot

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef _DATAMONODOMAIN_H_
#define _DATAMONODOMAIN_H_

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifesolver/heartFunctors.hpp>


namespace LifeV
{
/*!
  \class DataMonodomain

  Base class which holds usual data for the Monodomain model solvers

*/
class DataMonodomain:
    public DataMesh,
    public DataTime
{
public:

    //! @name Public Types
    //@{


    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const&, Real const&)> region_Type;

    //@}



    //! @name Constructor & Destructor
    //@{

    //! Constructors
    DataMonodomain();

    DataMonodomain( boost::shared_ptr<HeartFunctors> heart);

    DataMonodomain( const DataMonodomain& dataMonodomain );

    //@}

    //! @name Operators
    //@{

    DataMonodomain& operator=( const DataMonodomain& dataMonodomain );

    //@}

    //! @name Methods
    //@{

    //!output: show the data used for the simulation
    void showMe( std::ostream& output = std::cout );

    //@}



    //! @name Set Methods
    //@{


    //! external setup: set all the data for the simulation
    void        setup( const GetPot& dataFile );

    //@}



    //! @name Get Methods
    //@{

    //! verbose
    const Real&        verbose() const {return M_verbose;};

    //! FE space order
    std::string uOrder()         const {return M_uOrder;};

    //! Chi
    Real        Chi()            const {return M_volumeSurfaceRatio;}

    //! fiber File
    string      fibers_file()    const {return M_fibersFile;}

    //local change of the diffusivity
    Int         heart_diff_fct() const {return M_heartDiffusionFactor;}

    //would you consider the fibers?
    bool        has_fibers()     const {return M_hasFibers;}

    //! sigma_l
    Real        sigmal()         const 	{return M_longitudinalConductivity;}

    //! sigma_t
    Real        sigmat()         const 	{return M_transversalConductivity;}

    //! lambda
    //Real lambda() const 	{return M_lambda;}//not used?????

    //! Cm
    Real        Cm()             const 	{return M_membraneCapacitance;}

    //! D
    Real        D()              const 	{return M_diffusivity;}

    //! Post_dir
    std::string post_dir()       const {return M_postProcessingDirectory;}

    //@}

private:

    bool        M_hasFibers;

    UInt        M_verbose;

    Int         M_heartDiffusionFactor;

    Real        M_diffusivity;
    Real        M_longitudinalConductivity;
    Real        M_membraneCapacitance;
    Real        M_transversalConductivity;
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
