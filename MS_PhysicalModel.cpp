//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Physical Model
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-03-2009
 */

#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

namespace LifeV {

std::map< std::string, modelsTypes > modelsMap;

UInt MS_PhysicalModel::M_modelsNumber = 0;

// ===================================================
// Constructors & Destructor
// ===================================================
MS_PhysicalModel::MS_PhysicalModel() :
    M_ID                (),
    M_type              (),
    M_dataFile          (),
    M_couplings         (),
    M_modelName         (),
    M_flags             (),
    M_geometryScale     (),
    M_geometryRotate    (),
    M_geometryTranslate (),
    M_dataPhysics       (),
    M_dataTime          (),
    M_comm              (),
    M_displayer         ()
{

#ifdef DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::MS_PhysicalModel() \n";
#endif

    M_ID = M_modelsNumber++;

    //Initialization of geometry arrays
    for ( UInt i( 0 ); i < nDimensions; ++i )
    {
        M_geometryScale[i]     = 1.;
        M_geometryRotate[i]    = 0.;
        M_geometryTranslate[i] = 0.;
    }
}

MS_PhysicalModel::MS_PhysicalModel( const MS_PhysicalModel& model ) :
    M_ID                ( model.M_ID ),
    M_type              ( model.M_type ),
    M_dataFile          ( model.M_dataFile ),
    M_couplings         ( model.M_couplings ),
    M_modelName         ( model.M_modelName ),
    M_flags             ( model.M_flags ),
    M_geometryScale     ( model.M_geometryScale ),
    M_geometryRotate    ( model.M_geometryRotate ),
    M_geometryTranslate ( model.M_geometryTranslate ),
    M_dataPhysics       ( model.M_dataPhysics ),
    M_dataTime          ( model.M_dataTime ),
    M_comm              ( model.M_comm ),
    M_displayer         ( model.M_displayer )
{

#ifdef DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::MS_PhysicalModel( model ) \n";
#endif

    M_ID = M_modelsNumber++;
}

// ===================================================
// Operators
// ===================================================
MS_PhysicalModel&
MS_PhysicalModel::operator=( const MS_PhysicalModel& model )
{
    if ( this != &model )
    {
        M_ID                = model.M_ID;
        M_type              = model.M_type;
        M_dataFile          = model.M_dataFile;
        M_couplings         = model.M_couplings;
        M_modelName         = model.M_modelName;
        M_flags             = model.M_flags;
        M_geometryScale     = model.M_geometryScale;
        M_geometryRotate    = model.M_geometryRotate;
        M_geometryTranslate = model.M_geometryTranslate;
        M_dataPhysics       = model.M_dataPhysics;
        M_dataTime          = model.M_dataTime;
        M_comm              = model.M_comm;
        M_displayer         = model.M_displayer;
    }
    return *this;
}

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_PhysicalModel::ShowMe()
{
    std::cout << "Model id            = " << M_ID << std::endl
              << "Model name          = " << M_modelName << std::endl
              << "Model type          = " << Enum2String( M_type, modelsMap ) << std::endl;

    std::cout << "Couplings number    = " << GetCouplingsNumber() << std::endl;
    std::cout << "Couplings type(s)   = ";
    for ( UInt i( 0 ); i < GetCouplingsNumber(); ++i )
        std::cout << Enum2String( M_couplings[i]->GetType(), couplingsMap ) << " ";
    std::cout << std::endl;
    std::cout << "Flags list          = ";
    for ( UInt i( 0 ); i < GetCouplingsNumber(); ++i )
        std::cout << M_flags[i] << " ";
    std::cout << std::endl << std::endl;

    std::cout << "Geometry scale      = ";
    for ( UInt i( 0 ); i < nDimensions; ++i )
        std::cout << M_geometryScale[i] << " ";
    std::cout << std::endl;
    std::cout << "Geometry rotate     = ";
    for ( UInt i( 0 ); i < nDimensions; ++i )
        std::cout << M_geometryRotate[i] << " ";
    std::cout << std::endl;
    std::cout << "Geometry translate  = ";
    for ( UInt i( 0 ); i < nDimensions; ++i )
        std::cout << M_geometryTranslate[i] << " ";
    std::cout << std::endl << std::endl;
}

// ===================================================
// Methods
// ===================================================
void
MS_PhysicalModel::ClearCouplingsList()
{
    M_couplings.clear();
}

// ===================================================
// Set Methods
// ===================================================
void
MS_PhysicalModel::SetID( const UInt& id )
{
    M_ID = id;
}

void
MS_PhysicalModel::SetDataFile( const std::string& dataFile )
{

#ifdef DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::SetDataFile( dataFile ) \n";
#endif

    M_dataFile = GetPot( dataFile );

    // Read modelName
    M_modelName = M_dataFile( "MultiScale/modelName", "modelName" );

    // Read flags
    UInt componentSize = M_dataFile.vector_variable_size( "MultiScale/couplingFlags" );
    M_flags.reserve( componentSize );
    for ( UInt j( 0 ); j < componentSize; ++j )
        M_flags.push_back( M_dataFile( "MultiScale/couplingFlags", 0, j ) );
}

void
MS_PhysicalModel::AddCoupling( const Coupling_ptrType& coupling )
{
    M_couplings.push_back( coupling );
}

void
MS_PhysicalModel::SetGeometry( const boost::array< Real, NDIM >& scale, const boost::array< Real, NDIM >& rotate, const boost::array< Real, NDIM >& translate )
{

#ifdef DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::SetGeometry( scale, rotate, translate ) \n";
#endif

    M_geometryScale     = scale;
    M_geometryRotate    = rotate;
    M_geometryTranslate = translate;
}

void
MS_PhysicalModel::SetData( const boost::shared_ptr< MS_PhysicalData >& dataPhysics,
                           const boost::shared_ptr< DataTime >& dataTime )
{

#ifdef DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::SetData( dataPhysics, dataTime ) \n";
#endif

    M_dataPhysics = dataPhysics;
    M_dataTime    = dataTime;
}

void
MS_PhysicalModel::SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm )
{

#ifdef DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm.get() ) );
}

// ===================================================
// Get Methods
// ===================================================
const UInt&
MS_PhysicalModel::GetID() const
{
    return M_ID;
}

const modelsTypes&
MS_PhysicalModel::GetType() const
{
    return M_type;
}

const BCFlag&
MS_PhysicalModel::GetFlag( const UInt& id ) const
{
    return M_flags[id];
}

const std::vector< BCFlag >&
MS_PhysicalModel::GetFlags() const
{
    return M_flags;
}

const std::string&
MS_PhysicalModel::GetModelName() const
{
    return M_modelName;
}

UInt
MS_PhysicalModel::GetCouplingsNumber() const
{
    return static_cast< UInt > ( M_couplings.size() );
}

UInt
MS_PhysicalModel::GetCouplingLocalID( const UInt& ID ) const
{
    for ( UInt localID( 0 ); localID < GetCouplingsNumber(); ++localID )
        if ( M_couplings[localID]->GetID() == ID )
            return localID;

    return -1;
}

Coupling_ptrType
MS_PhysicalModel::GetCoupling( const UInt& LocalID ) const
{
    return M_couplings[LocalID];
}

} // Namespace LifeV
