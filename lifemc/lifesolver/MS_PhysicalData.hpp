//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Physical Data
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 09-09-2009
 */

#ifndef MS_PhysicalData_H
#define MS_PhysicalData_H 1

#include <lifemc/lifesolver/MS_Definitions.hpp>

namespace LifeV {

//! MS_PhysicalData - Global data container for the physical quantities of the problem
/*!
 *  @author Cristiano Malossi
 *
 *  Up to now, the container contains:
 *  <ul>
 *      <li> Fluid density
 *      <li> Fluid viscosity
 *  </ul>
 */
class MS_PhysicalData
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_PhysicalData();

    //! Copy constructor
    /*!
     * @param PhysicalData MS_PhysicalData
     */
    MS_PhysicalData( const MS_PhysicalData& PhysicalData );

    //! Destructor
    ~MS_PhysicalData() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param PhysicalData MS_PhysicalData
     * @return reference to a copy of the class
     */
    MS_PhysicalData& operator=( const MS_PhysicalData& PhysicalData );

    //@}


    //! @name Methods
    //@{

    //! Read the physical quantities from a GetPot file
    /*!
     * @param dataFile GetPot file
     */
    void ReadData( const GetPot& dataFile );

    //! Display some information about the physical quantities
    void ShowMe();

    //@}


    //! @name Get Methods
    //@{

    //! Get the global fluid density
    /*!
     * @return density of the fluid
     */
    const Real& GetFluidDensity() const;

    //! Get the global fluid viscosity
    /*!
     * @return viscosity of the fluid
     */
    const Real& GetFluidViscosity() const;

    //@}

private:

    Real M_fluidDensity;
    Real M_fluidViscosity;

};

} // Namespace LifeV

#endif /* MS_PhysicalData_H */
