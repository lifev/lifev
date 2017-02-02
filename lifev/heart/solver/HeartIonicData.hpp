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
    @brief File containing a class for handling Ionic model data with GetPot

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>, J.Castelneau (INRIA)
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef _DATAIONIC_H_
#define _DATAIONIC_H_

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/fem/TimeData.hpp>

namespace LifeV
{
/*!
  \class DataIonic

  Base class which holds usual data for the ionic model solvers

*/
class HeartIonicData:
    public MeshData,
    public TimeData
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Constructors
    HeartIonicData();


    HeartIonicData ( const GetPot& dataFile );

    HeartIonicData ( const HeartIonicData& dataIonic );

    virtual ~HeartIonicData() {}
    //@}


    //! @name Operators
    //@{



    HeartIonicData& operator= ( const HeartIonicData& dataIonic );

    //@}

    //! @name Methods
    //@{

    //! output: show the data used for the simulation
    void showMe ( std::ostream& output = std::cout );

    //@}


    //! @name Set Methods
    //@{

    //!external setup: set all the data for the simulation
    void setup ( const GetPot& dataFile );

    const Real& RMCParameterA() const
    {
        return M_RMCParameterA;
    }

    const Real& RMCParameterB() const
    {
        return M_RMCParameterB;
    }

    const Real& RMCParameterC1() const
    {
        return M_RMCParameterC1;
    }

    const Real& RMCParameterC2() const
    {
        return M_RMCParameterC2;
    }

    const Real& RMCParameterD() const
    {
        return M_RMCParameterD;
    }

    const Real& RMCTimeUnit() const
    {
        return M_RMCTimeUnit;
    }

    const Real& RMCPotentialAmplitude() const
    {
        return M_RMCPotentialAmplitude;
    }

    const Real& RMCRestPotential() const
    {
        return M_RMCRestPotential;
    }

    const Real& RMCInitialRepolarization() const
    {
        return M_RMCInitialRepolarization;
    }

    const Real& MSTauIn() const
    {
        return M_MSTauIn;
    }

    const Real& MSTauOut() const
    {
        return M_MSTauOut;
    }

    const Real& MSTauOpen() const
    {
        return M_MSTauOpen;
    }

    const Real& MSTauClose() const
    {
        return M_MSTauClose;
    }

    const Real& MSCriticalPotential() const
    {
        return M_MSCriticalPotential;
    }

    const Real& MSPotentialMinimum() const
    {
        return M_MSPotentialMinimum;
    }

    const Real& MSPotentialMaximum() const
    {
        return M_MSPotentialMaximum;
    }

    const Real& MSReactionAmplitude() const
    {
        return M_MSReactionAmplitude;
    }

    const Real& MSInitialTime() const
    {
        return M_MSInitialTime;
    }

    const Real& MSTend() const
    {
        return M_MSTend;
    }

    const Real& MSBDForder() const
    {
        return M_MSBDForder;
    }

    const bool& MSHasHeterogeneousTauClose() const
    {
        return M_MSHasHeterogeneousTauClose;
    }
    //@}

    /*//! End time
    Real endtime() const;

    //! FE space order
    std::string wOrder() const;
    */


private:



    std::string M_meshFile;
    std::string M_meshDirectory;

    UInt        M_verbose;
    //! RogersMcCulloch (RMC) 1994 Ionic Model parameters
    Real        M_RMCParameterA;
    Real        M_RMCParameterB;
    Real        M_RMCParameterC1;
    Real        M_RMCParameterC2;
    Real        M_RMCParameterD;
    Real        M_RMCTimeUnit;
    Real        M_RMCPotentialAmplitude;
    Real        M_RMCRestPotential;
    Real        M_RMCInitialRepolarization;

    //!Mitchell & Schaeffer (MS)
    Real        M_MSTauIn;   // = 0.8
    Real        M_MSTauOut;  // = 18.0
    Real        M_MSTauOpen; // = 300.0
    Real        M_MSTauClose;// = 100.0
    Real        M_MSCriticalPotential;    // =  -67.0
    Real        M_MSPotentialMinimum;
    Real        M_MSPotentialMaximum;
    Real        M_MSReactionAmplitude;
    Real        M_MSInitialTime;
    Real        M_MSTend;
    Real        M_MSBDForder;       //= 1

    bool        M_MSHasHeterogeneousTauClose;


};

}
#endif
