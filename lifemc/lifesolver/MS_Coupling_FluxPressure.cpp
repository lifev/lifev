/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-08-24

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
   \file MS_Coupling_FluxPressure.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-08-24
 */

#include <lifemc/lifesolver/MS_Coupling_FluxPressure.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
MS_Coupling_FluxPressure::MS_Coupling_FluxPressure( ) :
	super					( ),
	M_baseFlux				( ),
	M_basePressure			( ),
	M_flux					( 0. ),
	M_pressure				( 0. )
{

#ifdef DEBUG
	Debug( 8220 ) << "MS_Coupling_FluxPressure::MS_Coupling_FluxPressure() \n";
#endif

	M_type = FluxPressure;
}

MS_Coupling_FluxPressure::MS_Coupling_FluxPressure( const MS_Coupling_FluxPressure& FluxPressure ) :
	super					( FluxPressure ),
	M_baseFlux				( FluxPressure.M_baseFlux ),
	M_basePressure			( FluxPressure.M_basePressure ),
	M_flux					( FluxPressure.M_flux ),
	M_pressure				( FluxPressure.M_pressure )
{

#ifdef DEBUG
	Debug( 8220 ) << "MS_Coupling_FluxPressure::MS_Coupling_FluxPressure( FluxPressure ) \n";
#endif

}



// ===================================================
//! Methods
// ===================================================
MS_Coupling_FluxPressure&
MS_Coupling_FluxPressure::operator=( const MS_Coupling_FluxPressure& fluxPressure )
{
    if ( this != &fluxPressure )
    {
    	super::operator=( fluxPressure );
    	M_baseFlux				= fluxPressure.M_baseFlux;
    	M_basePressure			= fluxPressure.M_basePressure;
    	M_flux					= fluxPressure.M_flux;
    	M_pressure				= fluxPressure.M_pressure;
    }
	return *this;
}



// ===================================================
//! MultiScale Physical Coupling
// ===================================================
void
MS_Coupling_FluxPressure::SetupData( void )
{

#ifdef DEBUG
	Debug( 8220 ) << "MS_Coupling_FluxPressure::SetupData() \n";
#endif

	//Set Functions
	M_baseFlux.setFunction		( boost::bind( &MS_Coupling_FluxPressure::functionFlux,		this, _1, _2, _3, _4, _5 ) );
	M_basePressure.setFunction	( boost::bind( &MS_Coupling_FluxPressure::functionPressure,	this, _1, _2, _3, _4, _5 ) );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_FluxPressure::SetupCoupling( void )
{

#ifdef DEBUG
	Debug( 8220 ) << "MS_Coupling_FluxPressure::SetupCoupling() \n";
#endif

	// Impose flux
	switch ( M_models[0]->GetType() )
	{
		case Fluid3D :

			imposeFlux<MS_Model_Fluid3D>( M_models[0] );

			break;

		default :

			if ( M_displayer->isLeader() )
			    switchErrorMessage( M_models[0] );
	}

	// Impose pressure
	for ( UInt i(1) ; i < modelsNumber() ; ++i )
		switch ( M_models[i]->GetType() )
		{
			case Fluid3D :

				imposePressure<MS_Model_Fluid3D>( M_models[i], i );

				break;

			default :

				if ( M_displayer->isLeader() )
				    switchErrorMessage( M_models[i] );
		}

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_FluxPressure::UpdateCoupling( void )
{

#ifdef DEBUG
	Debug( 8220 ) << "MS_Coupling_FluxPressure::UpdateCoupling() \n";
#endif

	// Compute the flux
	M_flux = 0;

	for ( UInt i(1) ; i < modelsNumber() ; ++i )
		switch ( M_models[i]->GetType() )
		{
			case Fluid3D :
			{
				MS_Model_Fluid3D *Model = dynamic_cast<MS_Model_Fluid3D *>( &( *M_models[i] ) );

				M_flux -= Model->GetFlux( M_flags[i] );

				break;
			}

			default :

				if ( M_displayer->isLeader() )
				    switchErrorMessage( M_models[i] );
		}


	// Compute the pressure
	M_pressure = 0;

	switch ( M_models[0]->GetType() )
	{
		case Fluid3D :
		{
			MS_Model_Fluid3D *Model = dynamic_cast<MS_Model_Fluid3D *>( &( *M_models[0] ) );

			M_pressure = -( Model->GetPressure( M_flags[0] ) + 0.5 * Model->GetDensity( M_flags[0] ) * Model->GetFlux( M_flags[0] ) * Model->GetFlux( M_flags[0] ) / ( Model->GetArea( M_flags[0] ) * Model->GetArea( M_flags[0] ) ) );

			break;
		}

		default :

			if ( M_displayer->isLeader() )
			    switchErrorMessage( M_models[0] );
	}

#ifdef DEBUG
	Debug( 8220 ) << "                                            Coupling Flux     = " << M_flux		<< "\n";
	Debug( 8220 ) << "                                            Coupling Pressure = " << M_pressure	<< "\n";
#endif

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_FluxPressure::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
    	super::ShowMe();

    	std::cout 	<< "Coupling Flux     = " << M_flux		<< std::endl
					<< "Coupling Pressure = " << M_pressure	<< std::endl << std::endl;

    	std::cout	<< std::endl << std::endl;
    }

    //MPI Barrier
    M_comm->Barrier();
}



// ===================================================
//! Private Methods
// ===================================================
Real
MS_Coupling_FluxPressure::functionFlux( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ )
{
	return M_flux;
}

Real
MS_Coupling_FluxPressure::functionPressure( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ )
{
	return M_pressure;
}

} // Namespace LifeV
