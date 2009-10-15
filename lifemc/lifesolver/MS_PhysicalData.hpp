/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-09-09

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MultiScale_DataModel.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-09-09
 */

#ifndef __MS_PhysicalData_H
#define __MS_PhysicalData_H 1

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>

#include <string>

namespace LifeV {

//todo Allow specific value for the models

//! MS_PhysicalData - Global data container for the physical quantities of the problem
/*!
 *  @author Cristiano Malossi
 */
class MS_PhysicalData
{
public:

    //! @name Constructors, Destructor
    //@{

    //! Constructor
    MS_PhysicalData();

    //! Copy constructor
    /*!
     * \param PhysicalData - MS_PhysicalData
     */
    MS_PhysicalData( const MS_PhysicalData& PhysicalData );

    //! Destructor
    virtual ~MS_PhysicalData() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param PhysicalData - MS_PhysicalData
     */
    MS_PhysicalData& operator=( const MS_PhysicalData& PhysicalData );

    //! Read the physical quantities from a GetPot file
    /*!
     * \param dataFile - GetPot file
     */
    void ReadData( const GetPot& dataFile );

    //! Display some information about the physical quantities
    void ShowMe( void );

    //@}


    //! @name Get functions
    //@{

    //! Get the global fluid density
    const Real& GetFluidDensity() const
    {
        return M_fluidDensity;
    }

    //! Get the global fluid viscosity
    const Real& GetFluidViscosity() const
    {
        return M_fluidViscosity;
    }

    //@}

private:

    Real M_fluidDensity;
    Real M_fluidViscosity;

};

} // Namespace LifeV

#endif /* __MS_PhysicalData_H */
