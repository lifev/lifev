/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-09-02

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
   \file MS_Coupling_BoundaryCondition.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-09-02
 */

#include <lifemc/lifesolver/MS_Coupling_BoundaryCondition.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition( ) :
	super					( ),
	M_list					( ),
	M_listSize				( )
{

#ifdef DEBUG
	Debug( 8210 ) << "MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition() \n";
#endif

	M_type = BoundaryCondition;

}

MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition( const MS_Coupling_BoundaryCondition& BoundaryCondition ) :
	super					( BoundaryCondition ),
	M_list					( BoundaryCondition.M_list ),
	M_listSize				( BoundaryCondition.M_listSize )
{

#ifdef DEBUG
	Debug( 8210 ) << "MS_Coupling_BoundaryCondition::MS_Coupling_BoundaryCondition( BoundaryCondition ) \n";
#endif

}



// ===================================================
//! Methods
// ===================================================
MS_Coupling_BoundaryCondition&
MS_Coupling_BoundaryCondition::operator=( const MS_Coupling_BoundaryCondition& boundaryCondition )
{
    if ( this != &boundaryCondition )
    {
    	super::operator=( boundaryCondition );
    	M_list				= boundaryCondition.M_list;
    	M_listSize			= boundaryCondition.M_listSize;
    }
	return *this;
}



// ===================================================
//! MultiScale Physical Coupling
// ===================================================
void
MS_Coupling_BoundaryCondition::SetupData( void )
{

#ifdef DEBUG
	Debug( 8210 ) << "MS_Coupling_BoundaryCondition::SetupData() \n";
#endif

	//Load the list of boundary conditions
    M_listSize = M_dataFile.vector_variable_size( "boundary_conditions/list" );

    M_list.reserve( M_listSize );
    for ( UInt i(0) ; i < M_listSize ; ++i )
    	M_list.push_back( M_dataFile( "boundary_conditions/list", " ", i ) );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_BoundaryCondition::SetupCoupling( void )
{

#ifdef DEBUG
	Debug( 8210 ) << "MS_Coupling_BoundaryCondition::SetupCoupling() \n";
#endif

	for ( UInt i(0) ; i < modelsNumber() ; ++i )
		switch ( M_models[i]->GetType() )
		{
			case Fluid3D :

				ApplyBoundaryConditions<MS_Model_Fluid3D>( M_models[i] );

				break;

			default :

				if ( M_displayer->isLeader() )
				    switchErrorMessage( M_models[i] );
		}

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_BoundaryCondition::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
    	super::ShowMe();

    	std::cout 	<< "List size         = " << M_listSize	<< std::endl;
    	std::cout 	<< "List              = ";
    	for ( UInt i(0) ; i < M_listSize ; ++i )
			std::cout << M_list[i] << " ";
    	std::cout	<< std::endl << std::endl << std::endl << std::endl;
    }

    //MPI Barrier
    M_comm->Barrier();
}

} // Namespace LifeV
