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
 *  @file
 *  @brief File containing the MultiScale Physical Data
 *
 *  @date 09-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleData_H
#define MultiscaleData_H 1

#include <lifemc/lifesolver/MultiscaleDefinitions.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleData - Global data container for the physical quantities of the problem
/*!
 *  @author Cristiano Malossi
 *
 *  Up to now, the container contains:
 *  <ul>
 *      <li> Fluid density
 *      <li> Fluid viscosity
 *      <li> Structure Young modulus
 *      <li> Structure Poisson coefficient
 *  </ul>
 */
class MultiscaleData
{
public:

    //! @name Type definitions
    //@{

    typedef DataTime                                     time_Type;
    typedef boost::shared_ptr< time_Type >               timePtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleData();

    //! Copy constructor
    /*!
     * @param PhysicalData MultiscaleData
     */
    explicit MultiscaleData( const MultiscaleData& data );

    //! Destructor
    virtual ~MultiscaleData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param PhysicalData MultiscaleData
     * @return reference to a copy of the class
     */
    MultiscaleData& operator=( const MultiscaleData& data );

    //@}


    //! @name Methods
    //@{

    //! Read the physical quantities from a GetPot file
    /*!
     * @param dataFile GetPot file
     */
    void readData( const GetPot& dataFile );

    //! Display some information about the physical quantities
    void showMe();

    //@}


    //! @name Get Methods
    //@{

    //! Get the time container.
    /*!
     * @return time container
     */
    timePtr_Type dataTime() const { return M_dataTime; }

    //! Get the global fluid density.
    /*!
     * @return density of the fluid.
     */
    const Real& fluidDensity() const { return M_fluidDensity; }

    //! Get the global fluid viscosity.
    /*!
     * @return viscosity of the fluid.
     */
    const Real& fluidViscosity() const { return M_fluidViscosity; }

    //! Get the global fluid reference pressure (used by 1D model).
    /*!
     * @return reference pressure of the fluid.
     */
    const Real& fluidReferencePressure() const { return M_fluidReferencePressure; }

    //! Get the global structural Poisson coefficient.
    /*!
     * @return Poisson coefficient of the structure.
     */
    const Real& structureDensity() const { return M_structureDensity; }

    //! Get the global structural density.
    /*!
     * @return density of the structure.
     */
    const Real& structurePoissonCoefficient() const { return M_structurePoissonCoefficient; }

    // //! Get the global structural thickness.
    // /*!
    // * @return thickness of the structure.
    // */
    // const Real& GetStructureThickness() const { return M_structureThickness; }

    //! Get the global structural Young modulus.
    /*!
     * @return Young modulus of the structure.
     */
    const Real& structureYoungModulus() const { return M_structureYoungModulus; }

    //@}

private:

    timePtr_Type                        M_dataTime;

    Real                                M_fluidDensity;
    Real                                M_fluidViscosity;
    Real                                M_fluidReferencePressure;

    Real                                M_structureDensity;
    Real                                M_structurePoissonCoefficient;
//    Real                                M_structureThickness;
    Real                                M_structureYoungModulus;

};

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleData_H */
