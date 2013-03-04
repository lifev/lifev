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
 *  @brief File containing the Multiscale Physical Data
 *
 *  @date 09-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleGlobalData_H
#define MultiscaleGlobalData_H 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// Boost classes
#include <boost/shared_ptr.hpp>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/TimeData.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleGlobalData - Global data container for the physical quantities of the problem
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  At the present time the global container has the following parameters:
 *  <ul>
 *      <li> Fluid density
 *      <li> Fluid viscosity
 *      <li> Fluid venous pressure
 *      <li> Solid external pressure
 *      <li> Solid density
 *      <li> Solid Young's modulus
 *      <li> Solid Poisson's ratio
 *      <li> ScalingFactor Resistance
 *      <li> ScalingFactor Compliance
 *  </ul>
 */
class MultiscaleGlobalData
{
public:

    //! @name Type definitions
    //@{

    typedef TimeData                                     time_Type;
    typedef boost::shared_ptr< time_Type >               timePtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleGlobalData();

    //! Copy constructor
    /*!
     * @param PhysicalData MultiscaleGlobalData
     */
    explicit MultiscaleGlobalData ( const MultiscaleGlobalData& data );

    //! Destructor
    virtual ~MultiscaleGlobalData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param PhysicalData MultiscaleGlobalData
     * @return reference to a copy of the class
     */
    MultiscaleGlobalData& operator= ( const MultiscaleGlobalData& data );

    //@}


    //! @name Methods
    //@{

    //! Read the physical quantities from a GetPot file
    /*!
     * @param dataFile GetPot file
     */
    void readData ( const GetPot& dataFile );

    //! Display some information about the physical quantities
    void showMe();

    //@}


    //! @name Get Methods
    //@{

    //! Get the time container.
    /*!
     * @return time container
     */
    timePtr_Type dataTime() const
    {
        return M_timeData;
    }

    //! Get the global fluid density.
    /*!
     * @return density of the fluid.
     */
    const Real& fluidDensity() const
    {
        return M_fluidDensity;
    }

    //! Get the global fluid viscosity.
    /*!
     * @return viscosity of the fluid.
     */
    const Real& fluidViscosity() const
    {
        return M_fluidViscosity;
    }

    //! Get the global fluid venous pressure.
    /*!
     * @return venous pressure of the fluid.
     */
    const Real& fluidVenousPressure() const
    {
        return M_fluidVenousPressure;
    }

    //! Get the global fluid reference pressure (used by 1D model).
    /*!
     * @return reference pressure of the fluid.
     */
    const Real& solidExternalPressure() const
    {
        return M_solidExternalPressure;
    }

    //! Get the global structural Poisson coefficient.
    /*!
     * @return Poisson coefficient of the solid.
     */
    const Real& solidDensity() const
    {
        return M_solidDensity;
    }

    //! Get the global structural density.
    /*!
     * @return density of the solid.
     */
    const Real& solidPoissonCoefficient() const
    {
        return M_solidPoissonCoefficient;
    }

    // //! Get the global structural thickness.
    // /*!
    // * @return thickness of the solid.
    // */
    // const Real& solidThickness() const { return M_structureThickness; }

    //! Get the global structural Young modulus.
    /*!
     * @return Young modulus of the solid.
     */
    const Real& solidYoungModulus() const
    {
        return M_solidYoungModulus;
    }

    //! Get the global resistance scaling factor.
    /*!
     * @return resistance scaling factor.
     */
    const Real& scalingFactorResistance() const
    {
        return M_scalingFactorResistance;
    }

    //! Get the global compliance scaling factor.
    /*!
     * @return compliance scaling factor.
     */
    const Real& scalingFactorCompliance() const
    {
        return M_scalingFactorCompliance;
    }

    //@}

private:

    timePtr_Type                        M_timeData;

    Real                                M_fluidDensity;
    Real                                M_fluidViscosity;
    Real                                M_fluidVenousPressure;

    Real                                M_solidExternalPressure;
    Real                                M_solidDensity;
    Real                                M_solidPoissonCoefficient;
    //    Real                                M_structureThickness;
    Real                                M_solidYoungModulus;

    Real                                M_scalingFactorResistance;
    Real                                M_scalingFactorCompliance;

};

} // Namespace Multiscale
} // Namespace LifeV

#endif /* MultiscaleGlobalData_H */
