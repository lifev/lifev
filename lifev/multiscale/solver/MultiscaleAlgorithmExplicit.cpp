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
 *  @brief File containing the Multiscale Explicit Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleAlgorithmExplicit.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithmExplicit::MultiscaleAlgorithmExplicit() :
        multiscaleAlgorithm_Type    ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8011 ) << "MultiscaleAlgorithmExplicit::MultiscaleAlgorithmExplicit() \n";
#endif

    M_type = Explicit;
}

// ===================================================
// Multiscale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmExplicit::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8013 ) << "MultiscaleAlgorithmNewton::setupData( fileName ) \n";
#endif

    // Read parameters
    multiscaleParameterListPtr_Type solverParametersList = Teuchos::rcp( new Teuchos::ParameterList );
    solverParametersList = Teuchos::getParametersFromXmlFile( fileName );

    setAlgorithmName( solverParametersList->sublist( "Multiscale", true, "" ) );
    setAlgorithmParameters( solverParametersList->sublist( "Multiscale Algorithm", true, "" ) );
}

void
MultiscaleAlgorithmExplicit::subIterate()
{
    checkResidual( 0 );
}

} // Namespace Multiscale
} // Namespace LifeV
