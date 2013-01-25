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
    @brief Contains the data container for the LevelSetSolver

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 10-11-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */

#include <lifev/level_set/solver/LevelSetData.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
DataLevelSet::
DataLevelSet():
        M_time          ( ),
        M_timeAdvance   ( ),
        M_stabilization ( ),
        M_IPTreatment   ( ),
        M_IPCoef        ( )
{}

// ===================================================
// Methods
// ===================================================
void
DataLevelSet::
setup( const GetPot& dataFile, const std::string& section)
{
    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new time_Type( dataFile, section + "/time_discretization" ) );

    if ( !M_timeAdvance.get() )
        M_timeAdvance.reset( new timeAdvance_Type( dataFile, section + "/time_discretization" ) );

    std::string stabName = dataFile((section+"/stabilization").data(),"none");
    setStabilization(stabName);
    std::string ipName = dataFile((section+"/ip/treatment").data(),"implicit");
    setIPTreatment(ipName);
    M_IPCoef = dataFile((section+"/ip/coefficient").data(),0.0);
}

void
DataLevelSet::
showMe(std::ostream& output) const
{
    output << " Time data : " << std::endl;
    M_time->showMe( output );
    M_timeAdvance->showMe( output );

    output << " Stabilization : ";
    if (M_stabilization == NONE) output << "none" << std::endl;
    if (M_stabilization == IP) output << "ip" << std::endl;

    output << " IP Treatment  : ";
    if (M_IPTreatment == IMPLICIT) output << "implicit" << std::endl;
    if (M_IPTreatment == SEMI_IMPLICIT) output << "semi-implicit" << std::endl;
    if (M_IPTreatment == EXPLICIT) output << "explicit" << std::endl;

    output << " IP coefficient : " << M_IPCoef << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
DataLevelSet::
setStabilization(const std::string& stab)
{
    if (stab.compare("ip") ==0)
    {
        M_stabilization = IP;
    }
    else
    {
        ASSERT( stab.compare("none") ==0, " Unknown stabilization! ");
        M_stabilization = NONE;
    }
}

void
DataLevelSet::
setIPTreatment(const std::string& treat)
{
    if (treat.compare("implicit") ==0)
    {
        M_IPTreatment = IMPLICIT;
    }
    else if (treat.compare("semi-implicit") ==0)
    {
        M_IPTreatment = SEMI_IMPLICIT;
    }
    else
    {
        ASSERT( treat.compare("explicit") ==0, " Unknown IP treatment! ");
        M_IPTreatment = EXPLICIT;
    }
}


} // Namespace LifeV
