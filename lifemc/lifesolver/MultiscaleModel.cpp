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
namespace multiscale
{

std::map< std::string, models_Type > multiscaleModelsMap;

UInt MultiscaleModel::M_modelsNumber = 0;

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModel::MultiscaleModel() :
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
    Debug( 8100 ) << "MultiscaleModel::MultiscaleModel() \n";
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
MultiscaleModel::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MultiscaleModel::SetupData( fileName ) \n";
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
MultiscaleModel::showMe()
{
    std::cout << "Model id            = " << M_ID << std::endl
              << "Model name          = " << M_modelName << std::endl
              << "Model type          = " << Enum2String( M_type, multiscaleModelsMap ) << std::endl;

    std::cout << "Couplings number    = " << couplingsNumber() << std::endl;
    std::cout << "Couplings ID(s)     = ";
    for ( UInt i( 0 ); i < couplingsNumber(); ++i )
        std::cout << M_couplings[i]->ID() << " ";
    std::cout << std::endl;
    std::cout << "Couplings type(s)   = ";
    for ( UInt i( 0 ); i < couplingsNumber(); ++i )
        std::cout << Enum2String( M_couplings[i]->type(), multiscaleCouplingsMap ) << " ";
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
// Set Methods
// ===================================================
void
MultiscaleModel::setGeometry( const boost::array< Real, NDIM >& scale,
                               const boost::array< Real, NDIM >& rotate,
                               const boost::array< Real, NDIM >& translate )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MultiscaleModel::SetGeometry( scale, rotate, translate ) \n";
#endif

    M_geometryScale     = scale;
    M_geometryRotate    = rotate;
    M_geometryTranslate = translate;
}

void
MultiscaleModel::setCommunicator( const multiscaleCommPtr_Type& comm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8100 ) << "MultiscaleModel::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_displayer.reset( new Displayer( M_comm ) );
}

// ===================================================
// Get Methods
// ===================================================
UInt
MultiscaleModel::couplingLocalID( const UInt& ID ) const
{
    for ( UInt localID( 0 ); localID < couplingsNumber(); ++localID )
        if ( M_couplings[localID]->ID() == ID )
            return localID;

    return -1;
}

} // Namespace multiscale
} // Namespace LifeV
