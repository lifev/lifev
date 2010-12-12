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

#include <lifemc/lifesolver/MultiscaleAlgorithmExplicit.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithmExplicit::MultiscaleAlgorithmExplicit() :
        MS_Algorithm_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MultiscaleAlgorithmExplicit::MultiscaleAlgorithmExplicit() \n";
#endif

    M_type = Explicit;
}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmExplicit::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8012 ) << "MultiscaleAlgorithmExplicit::SetupData( fileName ) \n";
#endif

    MS_Algorithm_Type::setupData( fileName );
}

void
MultiscaleAlgorithmExplicit::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8012 ) << "MultiscaleAlgorithmExplicit::SubIterate() \n";
#endif

    toleranceSatisfied();
}

void
MultiscaleAlgorithmExplicit::updateCouplingVariables()
{
    // We use the initialize method for updating the coupling
    M_multiscale->initializeCouplingVariables();
}

void
MultiscaleAlgorithmExplicit::showMe()
{
    if ( M_displayer->isLeader() )
    {
        MS_Algorithm_Type::showMe();
    }
}

} // Namespace LifeV
