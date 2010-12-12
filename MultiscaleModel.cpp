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
 *  @brief File containing the MultiScale Physical Model
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleModel.hpp>

namespace LifeV
{

std::map< std::string, models_Type > MS_modelsMap;

UInt MS_PhysicalModel::M_modelsNumber = 0;

// ===================================================
// Constructors & Destructor
// ===================================================
MS_PhysicalModel::MS_PhysicalModel() :
        M_ID                (),
        M_type              (),
        M_couplings         (),
        M_modelName         (),
        M_flags             (),
        M_globalData        (),
        M_geometryScale     (),
        M_geometryRotate    (),
        M_geometryTranslate (),
        M_comm              (),
        M_displayer         ()
{

#ifdef HAVE_LIFEV_DEBUG
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

// ===================================================
// MultiScale PhysicalModel Virtual Methods
// ===================================================
void
MS_PhysicalModel::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::SetupData( fileName ) \n";
#endif

    GetPot dataFile( fileName );

    // Read modelName
    M_modelName = dataFile( "MultiScale/modelName", "modelName" );

    // Read flags
    UInt componentSize = dataFile.vector_variable_size( "MultiScale/couplingFlags" );
    M_flags.reserve( componentSize );
    for ( UInt j( 0 ); j < componentSize; ++j )
        M_flags.push_back( dataFile( "MultiScale/couplingFlags", 0, j ) );
}

void
MS_PhysicalModel::showMe()
{
    std::cout << "Model id            = " << M_ID << std::endl
              << "Model name          = " << M_modelName << std::endl
              << "Model type          = " << Enum2String( M_type, MS_modelsMap ) << std::endl;

    std::cout << "Couplings number    = " << couplingsNumber() << std::endl;
    std::cout << "Couplings ID(s)     = ";
    for ( UInt i( 0 ); i < couplingsNumber(); ++i )
        std::cout << M_couplings[i]->ID() << " ";
    std::cout << std::endl;
    std::cout << "Couplings type(s)   = ";
    for ( UInt i( 0 ); i < couplingsNumber(); ++i )
        std::cout << Enum2String( M_couplings[i]->type(), MS_couplingsMap ) << " ";
    std::cout << std::endl;
    std::cout << "Flags list          = ";
    for ( UInt i( 0 ); i < couplingsNumber(); ++i )
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
MS_PhysicalModel::clearCouplingsList()
{
    M_couplings.clear();
}

// ===================================================
// Set Methods
// ===================================================
void
MS_PhysicalModel::setID( const UInt& id )
{
    M_ID = id;
}

void
MS_PhysicalModel::addCoupling( const MS_Coupling_PtrType& coupling )
{
    M_couplings.push_back( coupling );
}

void
MS_PhysicalModel::setGlobalData( const MS_GlobalDataContainer_PtrType& globalData )
{
    M_globalData = globalData;
}

void
MS_PhysicalModel::setGeometry( const boost::array< Real, NDIM >& scale,
                               const boost::array< Real, NDIM >& rotate,
                               const boost::array< Real, NDIM >& translate )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::SetGeometry( scale, rotate, translate ) \n";
#endif

    M_geometryScale     = scale;
    M_geometryRotate    = rotate;
    M_geometryTranslate = translate;
}

void
MS_PhysicalModel::setCommunicator( const MS_Comm_PtrType& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MS_PhysicalModel::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

// ===================================================
// Get Methods
// ===================================================
const UInt&
MS_PhysicalModel::ID() const
{
    return M_ID;
}

const models_Type&
MS_PhysicalModel::type() const
{
    return M_type;
}

const BCFlag&
MS_PhysicalModel::flag( const UInt& id ) const
{
    return M_flags[id];
}

const std::vector< BCFlag >&
MS_PhysicalModel::flags() const
{
    return M_flags;
}

const std::string&
MS_PhysicalModel::modelName() const
{
    return M_modelName;
}

UInt
MS_PhysicalModel::couplingsNumber() const
{
    return static_cast< UInt > ( M_couplings.size() );
}

UInt
MS_PhysicalModel::couplingLocalID( const UInt& ID ) const
{
    for ( UInt localID( 0 ); localID < couplingsNumber(); ++localID )
        if ( M_couplings[localID]->ID() == ID )
            return localID;

    return -1;
}

MS_Coupling_PtrType
MS_PhysicalModel::coupling( const UInt& localID ) const
{
    return M_couplings[localID];
}

const MS_GlobalDataContainer_PtrType&
MS_PhysicalModel::globalData() const
{
    return M_globalData;
}

} // Namespace LifeV
