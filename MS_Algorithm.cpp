/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-23

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
 \file MS_Algorithm.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-23
 */

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV {

std::map< std::string, algorithmsTypes > algorithmMap;

// ===================================================
//! Constructors
// ===================================================
MS_Algorithm::MS_Algorithm() :
    M_type                       (),
    M_multiscale                 (),
    M_couplingVariables          (),
    M_couplingResiduals          (),
    M_comm                       (),
    M_displayer                  (),
    M_SubiterationsMaximumNumber (),
    M_Tolerance                  ()
{

#ifdef DEBUG
    Debug( 8010 ) << "MS_Algorithm::MS_Algorithm() \n";
#endif

}

MS_Algorithm::MS_Algorithm( const MS_Algorithm& algorithm ) :
    M_type                       ( algorithm.M_type ),
    M_multiscale                 ( algorithm.M_multiscale ),
    M_couplingVariables          ( algorithm.M_couplingVariables ),
    M_couplingResiduals          ( algorithm.M_couplingResiduals ),
    M_comm                       ( algorithm.M_comm ),
    M_displayer                  ( algorithm.M_displayer ),
    M_SubiterationsMaximumNumber ( algorithm.M_SubiterationsMaximumNumber ),
    M_Tolerance                  ( algorithm.M_Tolerance )
{

#ifdef DEBUG
    Debug( 8010 ) << "MS_Algorithm::MS_Algorithm( algorithm ) \n";
#endif

}

// ===================================================
//! Methods
// ===================================================
MS_Algorithm&
MS_Algorithm::operator=( const MS_Algorithm& algorithm )
{
    if ( this != &algorithm )
    {
        M_type                       = algorithm.M_type;
        M_multiscale                 = algorithm.M_multiscale;
        M_couplingVariables          = algorithm.M_couplingVariables;
        M_couplingResiduals          = algorithm.M_couplingResiduals;
        M_comm                       = algorithm.M_comm;
        M_displayer                  = algorithm.M_displayer;
        M_SubiterationsMaximumNumber = algorithm.M_SubiterationsMaximumNumber;
        M_Tolerance                  = algorithm.M_Tolerance;
    }
    return *this;
}

void
MS_Algorithm::SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm )
{

#ifdef DEBUG
    Debug( 8100 ) << "MS_Algorithm::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm.get() ) );
}

void
MS_Algorithm::SetMultiScaleProblem( const boost::shared_ptr< MS_Model_MultiScale > multiscale )
{

#ifdef DEBUG
    Debug( 8010 ) << "MS_Algorithm::SetMultiScaleProblem( multiscale ) \n";
#endif

    M_multiscale = multiscale;

    // Build coupling variables and residuals vectors
    std::vector<int> MyGlobalElements; //MyGlobalElements.resize( M_multiscale.GetCouplingVariablesNumber() );
    UInt couplingVariablesGlobalNumber = 0;
    M_multiscale->BuildCouplingVectorsMap( couplingVariablesGlobalNumber, MyGlobalElements );

    EpetraMap map( static_cast<int> ( MyGlobalElements.size() ), static_cast<int> ( MyGlobalElements.size() ), &MyGlobalElements[0], 0, *M_comm );
    M_couplingVariables.reset( new EpetraVector( map, Unique ) );
    M_couplingResiduals.reset( new EpetraVector( map, Unique ) );
}

// ===================================================
//! Virtual Methods
// ===================================================
void
MS_Algorithm::SetupData( const GetPot& DataFile )
{

#ifdef DEBUG
    Debug( 8010 ) << "MS_Algorithm::SetupData( DataFile ) \n";
#endif

    M_SubiterationsMaximumNumber = DataFile( "Solver/Algorithm/subITMax", 100 );
    M_Tolerance                  = DataFile( "Solver/Algorithm/tolerance", 1.e-10 );
}

void
MS_Algorithm::ShowMe( void )
{
    std::cout << "=================== Algorithm Information ===================" << std::endl << std::endl;

    std::cout << "Algorithm type      = " << Enum2String( M_type, algorithmMap ) << std::endl
              << "Max Sub-iterations = " << M_SubiterationsMaximumNumber << std::endl
              << "Tolerance          = " << M_Tolerance << std::endl << std::endl;
    std::cout << std::endl << std::endl;
}

} // Namespace LifeV
