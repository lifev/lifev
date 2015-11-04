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

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/heart/solver/HeartFunctors.hpp>


namespace LifeV
{
/*!
  \class DataMonodomain

  Base class which holds usual data for the Monodomain model solvers

*/
class HeartMonodomainData:
    public MeshData,
    public TimeData
{
public:

    //! @name Public Types
    //@{


    typedef std::function < Real ( Real const& x,
                                     Real const& y,
                                     Real const& z,
                                     Real const& t,
                                     ID const& id,
                                     Real const&) > region_Type;

    //@}



    //! @name Constructor & Destructor
    //@{

    //! Constructors
    HeartMonodomainData();

    HeartMonodomainData ( std::shared_ptr<HeartFunctors> heart);

    HeartMonodomainData ( const HeartMonodomainData& dataMonodomain );

    virtual           ~HeartMonodomainData() {}
    //@}

    //! @name Operators
    //@{

    HeartMonodomainData&    operator= ( const HeartMonodomainData& dataMonodomain );

    //@}

    //! @name Methods
    //@{

    //!output: show the data used for the simulation
    void               showMe ( std::ostream& output = std::cout );

    //@}

    //! @name Set Methods
    //@{

    //! external setup: set all the data for the simulation
    void               setup ( const GetPot& dataFile );

    //@}



    //! @name Get Methods
    //@{

    const region_Type& reducedConductivityBox()         const
    {
        return M_reducedConductivityBox;
    }

    const region_Type& reducedConductivityCylinder()    const
    {
        return M_reducedConductivityCylinder;
    }

    const region_Type& reducedConductivitySphere()      const
    {
        return M_reducedConductivitySphere;
    }

    //! verbose
    const UInt&        verbose() const
    {
        return M_verbose;
    };

    //! FE space order
    std::string        uOrder()                         const
    {
        return M_uOrder;
    }

    //! Chi
    const Real&        volumeSurfaceRatio()             const
    {
        return M_volumeSurfaceRatio;
    }

    //! lambda, key parameter in the derivation of Monodomain equations
    const Real&        conductivityRatio()              const
    {
        return M_conductivityRatio;
    }

    //! fiber File
    std::string             fibersFile()                     const
    {
        return M_fibersFile;
    }

    //local change of the diffusivity
    const Int&         heartDiffusionFactor()           const
    {
        return M_heartDiffusionFactor;
    }

    //would you consider the fibers?
    bool               hasFibers()                      const
    {
        return M_hasFibers;
    }

    //! sigma_l
    const Real&        longitudinalConductivity()       const
    {
        return M_longitudinalConductivity;
    }

    //! sigma_t
    const Real&        transversalConductivity()        const
    {
        return M_transversalConductivity;
    }

    //! Cm
    const Real&        membraneCapacitance()            const
    {
        return M_membraneCapacitance;
    }

    //! D
    const Real&        diffusivity()                    const
    {
        return M_diffusivity;
    }

    //! Post_dir
    std::string        postProcessingDirectory()        const
    {
        return M_postProcessingDirectory;
    }

    //@}







private:

    region_Type M_reducedConductivityBox;
    region_Type M_reducedConductivityCylinder;
    region_Type M_reducedConductivitySphere;

    bool        M_hasFibers;

    UInt        M_verbose;

    Int         M_heartDiffusionFactor;

    Real        M_diffusivity;
    Real        M_longitudinalConductivity;
    Real        M_membraneCapacitance;
    Real        M_transversalConductivity;
    Real        M_volumeSurfaceRatio;
    Real        M_conductivityRatio;

    std::string M_fibersDirectory;
    std::string M_fibersFile;
    std::string M_postProcessingDirectory; //! full name
    std::string M_uOrder;

};

}
#endif
