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
 *  @brief File containing the MultiScale Explicit Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MS_Algorithm_Explicit.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Algorithm_Explicit::MS_Algorithm_Explicit() :
        MS_Algorithm_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MS_Algorithm_Explicit::MS_Algorithm_Explicit() \n";
#endif

    M_type = Explicit;
}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MS_Algorithm_Explicit::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8012 ) << "MS_Algorithm_Explicit::SetupData( fileName ) \n";
#endif

    MS_Algorithm_Type::setupData( fileName );
}

void
MS_Algorithm_Explicit::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8012 ) << "MS_Algorithm_Explicit::SubIterate() \n";
#endif

    toleranceSatisfied();
}

void
MS_Algorithm_Explicit::updateCouplingVariables()
{
    // We use the initialize method for updating the coupling
    M_multiscale->InitializeCouplingVariables();
}

void
MS_Algorithm_Explicit::showMe()
{
    if ( M_displayer->isLeader() )
    {
        MS_Algorithm_Type::showMe();
    }
}

} // Namespace LifeV
