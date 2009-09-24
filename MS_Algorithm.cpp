/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-03-12

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
   \date 2009-03-12
 */

#include <lifemc/lifesolver/MS_Algorithm.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
MS_Algorithm::MS_Algorithm() :
	M_comm						( ),
	M_displayer					( ),
	M_dataFile					( "Undefined" ),
	M_dataTime					( ),
	M_dataPhysics				( new MS_PhysicalData( ) ),
	M_models					( ),
	M_couplings					( ),
	M_modelsNumber				( ),
	M_couplingsNumber			( ),
	M_chrono					( )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::MS_Algorithm() \n";
#endif

}

MS_Algorithm::MS_Algorithm( const MS_Algorithm& algorithm ) :
	M_comm						( algorithm.M_comm ),
	M_displayer					( algorithm.M_displayer ),
	M_dataFile					( algorithm.M_dataFile ),
	M_dataTime					( algorithm.M_dataTime ),
	M_dataPhysics				( algorithm.M_dataPhysics ),
	M_models					( algorithm.M_models ),
	M_couplings					( algorithm.M_couplings ),
	M_modelsNumber				( algorithm.M_modelsNumber ),
	M_couplingsNumber			( algorithm.M_couplingsNumber ),
	M_chrono					( algorithm.M_chrono )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::MS_Algorithm( algorithm ) \n";
#endif

}



// ===================================================
//! Methods
// ===================================================
MS_Algorithm&
MS_Algorithm::operator=( const MS_Algorithm& algorithm )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::operator=( algorithm ) \n";
#endif

    if ( this != &algorithm )
    {
    	M_comm						= algorithm.M_comm;
    	M_displayer					= algorithm.M_displayer;
    	M_dataFile					= algorithm.M_dataFile;
    	M_dataTime					= algorithm.M_dataTime;
    	M_dataPhysics				= algorithm.M_dataPhysics;
    	M_models					= algorithm.M_models;
    	M_couplings					= algorithm.M_couplings;
    	M_modelsNumber				= algorithm.M_modelsNumber;
    	M_couplingsNumber			= algorithm.M_couplingsNumber;
    	M_chrono					= algorithm.M_chrono;
    }
	return *this;
}

void
MS_Algorithm::SetDataFile( const std::string& dataFile )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::SetDataFile( dataFile ) \n";
#endif

	M_dataFile	= GetPot( dataFile );
}

void
MS_Algorithm::SetCommunicator( const boost::shared_ptr<Epetra_Comm>& comm )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::SetCommunicator( comm ) \n";
#endif

	M_comm = comm;
	M_displayer.reset( new Displayer( M_comm.get() ) );
}

void
MS_Algorithm::SetupData( void )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::SetupData() \n";
#endif

	//Models & Couplings mapCreation
	modelsMap["Fluid3D"] 				= Fluid3D;
	couplingsMap["BoundaryCondition"]	= BoundaryCondition;
	couplingsMap["FluxPressure"]		= FluxPressure;

	//Models & Couplings Factory registration
	FactoryModels::instance().registerProduct	(	Fluid3D,			&createFluid3D );
	FactoryCouplings::instance().registerProduct(	FluxPressure,		&createFluxPressure );
	FactoryCouplings::instance().registerProduct(	BoundaryCondition,	&createBoundaryCondition );

	// Load
	loadModels( M_dataFile );
	loadCouplings( M_dataFile );
	loadGeometry( M_dataFile );

	// Time & Physics containers
	M_dataTime.reset( new DataTime ( M_dataFile, "Algorithm/time_discretization" ) );
	M_dataPhysics->ReadData( M_dataFile );
}

void
MS_Algorithm::SetupProblem( void )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::SetupProblem() \n";
#endif

	// SetData
	for ( UInt i(0) ; i < M_modelsNumber ; ++i )
		M_models[i]->SetData( M_dataPhysics, M_dataTime );



	// SetCommunicator
	for ( UInt i(0) ; i < M_modelsNumber ; ++i )
		M_models[i]->SetCommunicator( M_comm );

	for ( UInt i(0) ; i < M_couplingsNumber ; ++i )
		M_couplings[i]->SetCommunicator( M_comm );



	// SetupData
	for ( UInt i(0) ; i < M_modelsNumber ; ++i )
		M_models[i]->SetupData();

	for ( UInt i(0) ; i < M_couplingsNumber ; ++i )
		M_couplings[i]->SetupData();



	// SetupCoupling
	for ( UInt i(0) ; i < M_couplingsNumber ; ++i )
		M_couplings[i]->SetupCoupling();

	// SetupModel
	for ( UInt i(0) ; i < M_modelsNumber ; ++i )
		M_models[i]->SetupModel();
}

void
MS_Algorithm::SolveProblem( void )
{

#ifdef DEBUG
	Debug( 8000 ) << "MS_Algorithm::SolveProblem() \n";
#endif

	// Time loop
	for ( ; M_dataTime->canAdvance() ; M_dataTime->updateTime() )
	{
		M_chrono.start();

		if ( M_displayer->isLeader() )
		{
			std::cout << std::endl;
			std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
			std::cout << "            MULTISCALE SIMULATION TIME: " << M_dataTime->getTime() << " s" << std::endl;
			std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
		}

		//MPI Barrier
	    M_comm->Barrier();

		// Build or Update System
		for ( UInt i(0) ; i < M_modelsNumber ; ++i )
			if ( M_dataTime->isFirstTimeStep() )
				M_models[i]->BuildSystem();
			else
				M_models[i]->UpdateSystem();


		// UpdateCouplings
		for ( UInt i(0) ; i < M_couplingsNumber ; ++i )
			M_couplings[i]->UpdateCoupling();



		// SolveSystem
		for ( UInt i(0) ; i < M_modelsNumber ; ++i )
			M_models[i]->SolveSystem();

		M_chrono.stop();

		if ( M_displayer->isLeader() )
			std::cout << "      Total iteration time:                    " << M_chrono.diff() << " s" << std::endl;
	}
}

void
MS_Algorithm::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
    	std::cout 	<< std::endl << std::endl
					<< "=================== MultiScale Information =================="
					<< std::endl << std::endl;

    	std::cout 	<< "Models number    = " << M_modelsNumber		<< std::endl
					<< "Couplings number = " << M_couplingsNumber	<< std::endl << std::endl;

        std::cout 	<< "Initial time     = " << M_dataTime->getInitialTime()	<< std::endl
					<< "End time         = " << M_dataTime->getEndTime()		<< std::endl
					<< "TimeStep         = " << M_dataTime->getTimeStep()		<< std::endl << std::endl; //MOVE THESE!!!
        std::cout 	<< std::endl << std::endl;

        std::cout 	<< "==================== Models Information ====================="
    				<< std::endl << std::endl;

		for ( UInt i(0) ; i < M_modelsNumber ; ++i )
			M_models[i]->ShowMe();

		std::cout 	<< "=================== Couplings Information ==================="
					<< std::endl << std::endl;

		for ( UInt i(0) ; i < M_couplingsNumber ; ++i )
			M_couplings[i]->ShowMe();

		std::cout 	<< "============================================================="
					<< std::endl << std::endl;
    }
}



// ===================================================
//! Private Methods
// ===================================================
inline void
MS_Algorithm::loadModels( const GetPot& dataFile )
{
	UInt				id;
	modelsTypes			model;

	UInt columnNumber	= 3.0;
	UInt linesNumber	= dataFile.vector_variable_size( "Problem/models" ) / columnNumber;

	std::string path = dataFile( "Problem/modelsPath", "./");
	for ( UInt i(0) ; i < linesNumber ; ++i )
	{
		id		= 				dataFile( "Problem/models", 0,				i*columnNumber );
		model	= modelsMap[	dataFile( "Problem/models", "undefined",	i*columnNumber + 1 ) ];

		M_models[ id ] = PhysicalModel_ptr( FactoryModels::instance().createObject( model ) );
		M_models[ id ]->SetID( id );
		M_models[ id ]->SetDataFile( path + dataFile( "Problem/models", "undefined", i*columnNumber + 2 ) + ".dat" );
	}

	M_modelsNumber = static_cast<UInt> ( M_models.size() );
}

inline void
MS_Algorithm::loadCouplings( const GetPot& dataFile )
{
	UInt				id;
	couplingsTypes		coupling;

	std::vector<UInt>	modelsIDVector;
	std::vector<UInt>	flagsIDVector;

	UInt columnNumber	= 5.0;
	UInt linesNumber	= dataFile.vector_variable_size( "Problem/couplings" ) / columnNumber;

	std::string path = dataFile( "Problem/couplingsPath", "./");
	for ( UInt i(0) ; i < linesNumber ; ++i )
	{
		id			= 				dataFile( "Problem/couplings", 0,			i*columnNumber );
		coupling	= couplingsMap[	dataFile( "Problem/couplings", "undefined",	i*columnNumber + 1 ) ];

		M_couplings[ id ] = PhysicalCoupling_ptr( FactoryCouplings::instance().createObject( coupling ) );
		M_couplings[ id ]->SetID( id );
		M_couplings[ id ]->SetDataFile( path + dataFile( "Problem/couplings", "undefined", i*columnNumber + 2 ) + ".dat" );

		modelsIDVector = string2numVect<UInt>	( dataFile( "Problem/couplings", "undefined", i*columnNumber + 3 ) );
		flagsIDVector  = string2numVect<UInt>	( dataFile( "Problem/couplings", "undefined", i*columnNumber + 4 ) );
		for ( UInt j(0) ; j < static_cast<UInt> ( modelsIDVector.size() ) ; ++j )
		{
			M_couplings[ id ]->AddModel( M_models[ modelsIDVector[j] ] );
			M_couplings[ id ]->AddFlagID( flagsIDVector[j] );
		}
	}

	M_couplingsNumber = static_cast<UInt> ( M_couplings.size() );
}

inline void
MS_Algorithm::loadGeometry( const GetPot& dataFile )
{
	UInt					id;

	boost::array<Real,3>	geometryScale;
	boost::array<Real,3>	geometryRotate;
	boost::array<Real,3>	geometryTranslate;

	UInt columnNumber	= 10.0;
	UInt linesNumber	= dataFile.vector_variable_size( "Geometry/offset" ) / columnNumber;

	for ( UInt i(0) ; i < linesNumber ; ++i )
	{
		id						= dataFile( "Geometry/offset", 0,	i*columnNumber );

		geometryScale[0]		= dataFile( "Geometry/offset", 1.,	i*columnNumber + 1 );
		geometryScale[1]		= dataFile( "Geometry/offset", 1.,	i*columnNumber + 2 );
		geometryScale[2]		= dataFile( "Geometry/offset", 1.,	i*columnNumber + 3 );

		geometryRotate[0]		= dataFile( "Geometry/offset", 0.,	i*columnNumber + 4 ) * Pi / 180;
		geometryRotate[1]		= dataFile( "Geometry/offset", 0.,	i*columnNumber + 5 ) * Pi / 180;
		geometryRotate[2]		= dataFile( "Geometry/offset", 0.,	i*columnNumber + 6 ) * Pi / 180;

		geometryTranslate[0]	= dataFile( "Geometry/offset", 0.,	i*columnNumber + 7 );
		geometryTranslate[1]	= dataFile( "Geometry/offset", 0.,	i*columnNumber + 8 );
		geometryTranslate[2]	= dataFile( "Geometry/offset", 0.,	i*columnNumber + 9 );

		M_models[ id ]->SetGeometry( geometryScale, geometryRotate, geometryTranslate );
	}
}

template <typename number>
inline std::vector<number>
MS_Algorithm::string2numVect( const std::string& string )
{
	//Split the string
	std::vector<std::string> stringVector;
	boost::split( stringVector, string, boost::is_any_of(",") );

	//Convert to the right type
	std::vector< number > numberVector;
	for ( UInt i(0) ; i < static_cast<UInt>( stringVector.size() ) ; ++i )
		numberVector.push_back( static_cast<number>( std::atoi( stringVector[i].c_str() ) ) );

	return numberVector;
}

} // Namespace LifeV
