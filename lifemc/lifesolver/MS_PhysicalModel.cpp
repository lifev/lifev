/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-09-03

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
   \file MS_PhysicalModel.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-09-03
 */

#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

namespace LifeV {

std::map<std::string, modelsTypes> modelsMap;

std::string getModelTypeString( const modelsTypes& type )
{
	for ( std::map<std::string, modelsTypes>::iterator j = modelsMap.begin() ; j != modelsMap.end() ; ++j )
		if ( j->second == type )
			return j->first;

	return "NOTYPE";
}

// ===================================================
//! Constructors
// ===================================================
MS_PhysicalModel::MS_PhysicalModel( ) :
	M_ID					( ),
	M_type					( ),
	M_dataFile				( ),
	M_modelName				( ),
	M_flags					( ),
	M_geometryScale			( ),
	M_geometryRotate		( ),
	M_geometryTranslate		( ),
	M_dataPhysics			( ),
	M_dataTime				( ),
	M_comm					( ),
	M_displayer				( )
{

#ifdef DEBUG
	Debug( 8100 ) << "MS_PhysicalModel::MS_PhysicalModel() \n";
#endif

	//Initialization of geometry arrays
	for ( UInt i(0); i < nDimensions; ++i )
	{
		M_geometryScale[i]		= 1.;
		M_geometryRotate[i]		= 0.;
		M_geometryTranslate[i]	= 0.;
	}
}

MS_PhysicalModel::MS_PhysicalModel( const MS_PhysicalModel& model ) :
	M_ID					( model.M_ID ),
	M_type					( model.M_type ),
	M_dataFile				( model.M_dataFile ),
	M_modelName				( model.M_modelName ),
	M_flags					( model.M_flags ),
	M_geometryScale			( model.M_geometryScale ),
	M_geometryRotate		( model.M_geometryRotate ),
	M_geometryTranslate		( model.M_geometryTranslate ),
	M_dataPhysics			( model.M_dataPhysics ),
	M_dataTime				( model.M_dataTime ),
	M_comm					( model.M_comm ),
	M_displayer				( model.M_displayer )
{

#ifdef DEBUG
	Debug( 8100 ) << "MS_PhysicalModel::MS_PhysicalModel( model ) \n";
#endif

}



// ===================================================
//! Methods
// ===================================================
MS_PhysicalModel&
MS_PhysicalModel::operator=( const MS_PhysicalModel& model )
{
    if ( this != &model )
    {
    	M_ID				= model.M_ID;
    	M_type				= model.M_type;
    	M_dataFile			= model.M_dataFile;
    	M_modelName			= model.M_modelName;
    	M_flags				= model.M_flags;
    	M_geometryScale		= model.M_geometryScale;
    	M_geometryRotate	= model.M_geometryRotate;
    	M_geometryTranslate	= model.M_geometryTranslate;
    	M_dataPhysics		= model.M_dataPhysics;
    	M_dataTime			= model.M_dataTime;
    	M_comm				= model.M_comm;
    	M_displayer			= model.M_displayer;
    }
	return *this;
}



// ===================================================
//! Set Methods
// ===================================================
void
MS_PhysicalModel::SetDataFile( const std::string& dataFile )
{

#ifdef DEBUG
	Debug( 8100 ) << "MS_PhysicalModel::SetDataFile( dataFile ) \n";
#endif

	M_dataFile	= GetPot( dataFile );

	// Read modelName
	M_modelName = M_dataFile( "multiscale/modelName", "modelName" );

	// Read flags
    UInt componentSize = M_dataFile.vector_variable_size( "multiscale/couplingFlags" );
    M_flags.reserve( componentSize );
    for ( UInt j(0) ; j < componentSize ; ++j )
    	M_flags.push_back( M_dataFile( "multiscale/couplingFlags", 0, j) );
}

void
MS_PhysicalModel::SetGeometry( 	const boost::array<Real,nDimensions>& scale,
								const boost::array<Real,nDimensions>& rotate,
								const boost::array<Real,nDimensions>& translate )
{

#ifdef DEBUG
	Debug( 8100 ) << "MS_PhysicalModel::SetGeometry( scale, rotate, translate ) \n";
#endif

	M_geometryScale		= scale;
	M_geometryRotate	= rotate;
	M_geometryTranslate	= translate;
}


void
MS_PhysicalModel::SetData( 	const boost::shared_ptr<MS_PhysicalData>& dataPhysics,
							const boost::shared_ptr<DataTime>& dataTime )
{

#ifdef DEBUG
	Debug( 8100 ) << "MS_PhysicalModel::SetData( dataPhysics, dataTime ) \n";
#endif

	M_dataPhysics	= dataPhysics;
	M_dataTime		= dataTime;
}

void
MS_PhysicalModel::SetCommunicator( const boost::shared_ptr<Epetra_Comm>& comm )
{

#ifdef DEBUG
	Debug( 8100 ) << "MS_PhysicalModel::SetCommunicator( comm ) \n";
#endif

	M_comm = comm;
	M_displayer.reset( new Displayer( M_comm.get() ) );
}



// ===================================================
//! Virtual Functions
// ===================================================
void
MS_PhysicalModel::ShowMe( void )
{
	std::cout 	<< "Model id           = " << M_ID << std::endl
				<< "Model name         = " << M_modelName	<< std::endl
				<< "Model type         = " << getModelTypeString( M_type )	<< std::endl;

	std::cout	<< "Flags list         = ";
	for ( UInt i(0) ; i < flagsNumber() ; ++i )
		std::cout << M_flags[i] << " ";
	std::cout	<< std::endl << std::endl;

	std::cout	<< "Geometry scale     = ";
	for ( UInt i(0) ; i < nDimensions ; ++i )
		std::cout << M_geometryScale[i] << " ";
	std::cout	<< std::endl;
	std::cout	<< "Geometry rotate    = ";
	for ( UInt i(0) ; i < nDimensions ; ++i )
		std::cout << M_geometryRotate[i] << " ";
	std::cout	<< std::endl;
	std::cout	<< "Geometry translate = ";
	for ( UInt i(0) ; i < nDimensions ; ++i )
		std::cout << M_geometryTranslate[i] << " ";
	std::cout	<< std::endl << std::endl;
}

} // Namespace LifeV
