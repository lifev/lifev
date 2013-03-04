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

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/heart/solver/HeartFunctors.hpp>


namespace LifeV
{
/*!
  \class DataBidomain

  Base class which holds usual data for the Bidomain model solvers

*/

class HeartBidomainData:
    public MeshData,
    public TimeData
{
public:


    //! @name Public Types
    //@{

    /*! @enum listOfTemplatesOptions
        Description of the purpose of the enumerator list.
    */

    typedef boost::function < Real ( Real const& x,
                                     Real const& y,
                                     Real const& z,
                                     Real const& t,
                                     ID   const& id,
                                     Real const&) > region_Type;

    //@}



    //! @name Constructor & Destructor
    //@{

    //! Constructors
    HeartBidomainData();

    HeartBidomainData ( boost::shared_ptr<HeartFunctors> heart);

    HeartBidomainData ( const HeartBidomainData& dataBidomain );

    virtual ~HeartBidomainData() {}
    //@}


    //! @name Operators
    //@{

    HeartBidomainData& operator= ( const HeartBidomainData& dataBidomain );

    //@}

    //! @name Methods
    //@{

    //! Output: show the data used for the simulation
    void showMe ( std::ostream& output = std::cout );

    //@}



    //! @name Set Methods
    //@{

    //! external setup: set all the data for the simulation
    void setup ( const GetPot& dataFile );

    //@}


    //! @name Get Methods
    //@{

    const region_Type&      reducedConductivityBox()         const
    {
        return M_reducedConductivityBox;
    }

    const region_Type&      reducedConductivityCylinder()    const
    {
        return M_reducedConductivityCylinder;
    }

    const region_Type&      reducedConductivitySphere()      const
    {
        return M_reducedConductivitySphere;
    }

    //! verbose
    const bool&             verbose()                        const
    {
        return M_verbose;
    }

    //! FE space order
    std::string             uOrder()                         const
    {
        return M_uOrder;
    }

    //! Chi
    const Real&             volumeSurfaceRatio()             const
    {
        return M_volumeSurfaceRatio;
    }
    //! fiber File
    std::string             fibersFile()                     const
    {
        return M_fibersFile;
    }

    const Int&              heartDiffusionFactor()           const
    {
        return M_heartDiffusionFactor;
    }

    const bool&             hasFibers()                      const
    {
        return M_hasFibers;
    }

    //! format vct
    const bool&             fibersFormat()                   const
    {
        return M_fibersFormat;
    }

    //! sigma_l
    const Real&             longitudinalInternalConductivity() const
    {
        return M_longitudinalInternalConductivity;
    }
    const Real&             longitudinalExternalConductivity() const
    {
        return M_longitudinalExternalConductivity;
    }

    //! sigma_t
    const Real&             transversalInternalConductivity() const
    {
        return M_transversalInternalConductivity;
    }
    const Real&             transversalExternalConductivity() const
    {
        return M_transversalExternalConductivity;
    }

    //! Cm
    const Real&             membraneCapacitance()             const
    {
        return M_membraneCapacitance;
    }
    //! D
    const Real&             internalDiffusivity()             const
    {
        return M_internalDiffusivity;
    }
    //! Post_dir
    const Real&             externalDiffusivity()             const
    {
        return M_externalDiffusivity;
    }
    //! Post_dir
    std::string             postProcessingDirectory()         const
    {
        return M_postProcessingDirectory;
    }

    UInt                    BDForder()                        const
    {
        return M_BDForder;
    }

    //@}

private:

    region_Type M_reducedConductivityBox;
    region_Type M_reducedConductivityCylinder;
    region_Type M_reducedConductivitySphere;
    // format of fibers file
    bool        M_fibersFormat;
    bool        M_hasFibers;

    // order of time discretization BDF order
    UInt        M_BDForder;
    bool        M_verbose;

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

};

}
#endif
